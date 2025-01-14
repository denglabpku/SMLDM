# import debugpy

# # 5678 is the default attach port in the VS Code debug configurations. Unless a host and port are specified, host defaults to 127.0.0.1
# debugpy.listen(5678)
# print("Waiting for debugger attach")
# debugpy.wait_for_client()
# debugpy.breakpoint()
# print('break on this line')

import torch.multiprocessing
torch.multiprocessing.set_sharing_strategy('file_system') # these two lines are important to avoid file descriptor errors

import argparse
import logging
import os
# import random
# import sys
import torch
import torch.nn as nn
import torch.nn.functional as F
import torchvision.transforms as v2
# import torchvision.transforms.functional as TF
from pathlib import Path
from torch import optim
from torch.utils.data import DataLoader, random_split
from tqdm import tqdm
import torch.distributed as dist
from scipy.ndimage import distance_transform_edt
from skimage.morphology import binary_erosion, label

# from torchvision.transforms import Compose, ToTensor
# from torchvision.transforms import RandomHorizontalFlip, RandomVerticalFlip

import wandb
from evaluate import para_evaluate
from unet import UNet
from utils.data_loading import ParaBasicDataset
from utils.dice_score import dice_loss
from utils.dice_score import multiclass_dice_coeff, dice_coeff

def find_objects(binary_mask):
    labeled_array = label(binary_mask.cpu().numpy(),connectivity=1)
    num_objects = labeled_array.max()
    return labeled_array, num_objects

def find_nearest_border(binary_mask):
    eroded = binary_erosion(binary_mask)
    border = binary_mask ^ eroded  # XOR operation to find border pixels
    return border

def create_weight_matrix(y, w0=10.0, sigma=5.0):
    weight_matrices = []
    fb_ratios = []

    device = y.device
    for i in range(y.size(0)):
        binary_mask = y[i]
        labeled_array, num_objects = find_objects(binary_mask)

        # if binary_mask.sum() == 0:
        #     weight_matrix = torch.ones(binary_mask.shape, device=device)
        #     fb_ratio = torch.tensor(1.0, device=device)
        # elif num_objects <= 4:
        #     weight_matrix = torch.ones(binary_mask.shape, device=device)
        #     fb_ratio = torch.tensor(1.0, device=device)
        # else:
        dist_stack = torch.zeros((32, 32, num_objects), device=device)  # Tensor to store distances for all objects

        for obj_label in range(1, num_objects + 1):
            object_mask = (labeled_array == obj_label)

            # Find the nearest border 1 pixel for the current object
            border_pixels = find_nearest_border(object_mask)

            # Calculate distance of 0 pixels to the nearest border 1 pixel for each object
            dist_transform = torch.tensor(distance_transform_edt(1 - border_pixels), device=device)

            dist_transform[binary_mask == 1] = 0

            dist_stack[..., obj_label-1] = dist_transform

        # Sort the distances along the last axis to get A1 and A2
        sorted_distances, _ = torch.sort(dist_stack, dim=-1)

        A1 = sorted_distances[..., 0]  # Smallest pixel values
        A2 = sorted_distances[..., 1]  # Next smallest pixel values

        foreground = torch.sum(binary_mask)
        background = binary_mask.numel() - foreground
        fb_ratio = background / foreground

        distance_component = w0 * torch.exp(-((A1 ** 2 + A2 ** 2) / (2 * sigma**2)))
        distance_component[binary_mask == 1] = 0

        num_balance_component = torch.ones(binary_mask.shape, device=device)
        num_balance_component[binary_mask == 1] = fb_ratio

        weight_matrix = distance_component + num_balance_component

        weight_matrices.append(weight_matrix)
        fb_ratios.append(fb_ratio)

    return torch.stack(weight_matrices), torch.tensor(fb_ratios, device=y.device)

def train_model(
        model,
        device,
        dir_img,
        dir_mask,
        dir_checkpoint,
        run,
        epochs: int = 5,
        batch_size: int = 1,
        learning_rate: float = 1e-5,
        val_percent: float = 0.1,
        save_checkpoint: bool = True,
        img_scale: float = 0.5,
        amp: bool = False,
        weight_decay: float = 1e-8,
        momentum: float = 0.999,
        gradient_clipping: float = 1.0,

        
):   
    # Check to see if local_rank is 0
    is_master = device == 0
    do_log = run is not None

    # initialize PyTorch distributed using environment variables, world_size equals GPU cards
    dist.init_process_group(backend="nccl", init_method="env://", rank=device, world_size=2)
    torch.cuda.set_device(device)

    model = nn.parallel.DistributedDataParallel(model, device_ids=[device],output_device=device)

    # watch gradients only for rank 0
    if is_master and do_log:
        run.watch(model)
        run.config.update(
        dict(epochs=epochs, batch_size=batch_size, learning_rate=learning_rate,
             val_percent=val_percent, save_checkpoint=save_checkpoint, img_scale=img_scale, amp=amp)
    )

    # 1. Create dataset
    transform = v2.Compose([
        v2.RandomHorizontalFlip(),
        v2.RandomVerticalFlip(),
        v2.RandomRotation(90), # random_rotation_range [-90, 90]
        v2.RandomAffine(degrees=0, translate=(0.1, 0.3)),
    ])

    dataset = ParaBasicDataset(dir_img, dir_mask, img_scale, transform=transform, device=device)

    # 2. Split into train / validation partitions
    n_val = int(len(dataset) * val_percent)
    n_train = len(dataset) - n_val
    train_set, val_set = random_split(dataset, [n_train, n_val], generator=torch.Generator().manual_seed(0))

    # 3. Create data loaders
    num_workers = min(os.cpu_count(), 35)  # Limit to a maximum of 8 workers
    loader_args = dict(batch_size=batch_size, num_workers=num_workers, pin_memory=True)
    train_loader = DataLoader(train_set, shuffle=True, **loader_args)
    val_loader = DataLoader(val_set, shuffle=False, drop_last=False, **loader_args)

    # Inform user training begun
    print('Training model...')

    # 4. Set up the optimizer, the loss, the learning rate scheduler
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, 'min', patience=5)  # goal: minimize loss score
    grad_scaler = torch.cuda.amp.GradScaler(enabled=amp)    

    global_train_step = 0
    global_val_step = 0

    # 5. Begin training
    for epoch in range(1, epochs + 1):
        
        # Training round
        model.train()
        epoch_loss = 0
        t = len(train_loader)

        with tqdm(total=n_train, desc=f'Epoch {epoch}/{epochs}', unit='img') as pbar:
            for batch in train_loader:
                images, true_masks = batch['image'], batch['mask']

                assert images.shape[1] == model.module.n_channels, \
                    f'Network has been defined with {model.module.n_channels} input channels, ' \
                    f'but loaded images have {images.shape[1]} channels. Please check that ' \
                    'the images are loaded correctly.'

                images = images.to(device=device, dtype=torch.float32, memory_format=torch.channels_last)
                true_masks = true_masks.to(device=device, dtype=torch.long)

                masks_pred = model(images)
                if model.module.n_classes == 1:
                    y = true_masks.float()
                    weight_matrix, _ = create_weight_matrix(y,w0=10.0, sigma=2.0)                    
                    criterion = nn.BCEWithLogitsLoss(reduction="none")
                    raw_loss = criterion(masks_pred.squeeze(1), true_masks.float())
                    # Apply the weights from weight_matrix_tensor to the raw loss
                    loss = (raw_loss * weight_matrix).mean()                    
                else:
                    loss = criterion(masks_pred, true_masks)
                    loss += dice_loss(
                        F.softmax(masks_pred, dim=1).float(),
                        F.one_hot(true_masks, model.module.n_classes).permute(0, 3, 1, 2).float(),
                        multiclass=True
                    )

                optimizer.zero_grad(set_to_none=True)
                grad_scaler.scale(loss).backward()
                torch.nn.utils.clip_grad_norm_(model.parameters(), gradient_clipping)
                grad_scaler.step(optimizer)
                grad_scaler.update()                

                pbar.update(images.shape[0])
                global_train_step += 1
                t -= 1
                epoch_loss += loss.item()
                if do_log:
                    run.log({
                        'train loss': loss.item(),
                        'train step': global_train_step,
                        'epoch': epoch
                    })
                pbar.set_postfix(**{'loss (batch)': loss.item()})

        # Evaluation round       
        model.eval()
        dice_score = 0
        t = len(val_loader)

        with tqdm(total=len(val_loader), desc=f'Validation', unit='batch') as pbar:
            for batch in val_loader:
                image, mask_true = batch['image'], batch['mask']

                # move images and labels to correct device and type
                image = image.to(device=device, dtype=torch.float32, memory_format=torch.channels_last)
                mask_true = mask_true.to(device=device, dtype=torch.long)

                # predict the mask
                mask_pred = model(image)

                if model.module.n_classes == 1:
                    assert mask_true.min() >= 0 and mask_true.max() <= 1, 'True mask indices should be in [0, 1]' 
                    y = mask_true.float()
                    weight_matrix, _ = create_weight_matrix(y,w0=10.0, sigma=2.0)                    
                    criterion = nn.BCEWithLogitsLoss(reduction="none")
                    raw_loss = criterion(mask_pred.squeeze(1), mask_true.float())
                    # Apply the weights from weight_matrix_tensor to the raw loss
                    loss = (raw_loss * weight_matrix).mean()                    
                    
                    mask_pred = mask_pred.squeeze()
                    mask_pred = (torch.sigmoid(mask_pred) > 0.5).float()
                    
                else:
                    assert mask_true.min() >= 0 and mask_true.max() < model.module.n_classes, 'True mask indices should be in [0, n_classes]'
                    # convert to one-hot format
                    mask_true = F.one_hot(mask_true, model.module.n_classes).permute(0, 3, 1, 2).float()
                    mask_pred = F.one_hot(mask_pred.argmax(dim=1), model.module.n_classes).permute(0, 3, 1, 2).float()
                    # compute the Dice score, ignoring background
                    batch_dice_score = multiclass_dice_coeff(mask_pred[:, 1:], mask_true[:, 1:], reduce_batch_first=False)
                    dice_score += batch_dice_score

                # Update the learning rate using the scheduler based on the current dice score                
                scheduler.step(loss)

                pbar.update()

                global_val_step += 1

                t -= 1

                # logging.info('Validation Dice score: {}'.format(val_score))
                if t <= 3 and do_log:
                    run.log({
                        'learning rate': optimizer.param_groups[0]['lr'],
                        'validation loss': loss,
                        'images': wandb.Image(image[0].cpu()),
                        'masks': {
                            'true': wandb.Image(mask_true[0].float().cpu()),
                            'pred': wandb.Image(mask_pred[0].float().cpu()),
                        },
                        'val step': global_val_step,
                        'epoch': epoch
                    })
                    
        if is_master and save_checkpoint:
            Path(dir_checkpoint).mkdir(parents=True, exist_ok=True)
            state_dict = model.module.state_dict()
            state_dict['mask_values'] = dataset.mask_values
            torch.save(state_dict, str(dir_checkpoint / 'checkpoint_epoch{}.pth'.format(epoch)))            
            print(f'Checkpoint {epoch} saved!')
        
        dist.barrier()

    cleanup()

    # Finish the run
    wandb.finish()
        
    # Inform user training ended
    print('Training Completed!')

    return



def get_args():
    parser = argparse.ArgumentParser(description='Train the UNet on images and target masks')
    parser.add_argument('--epochs', '-e', metavar='E', type=int, default=5, help='Number of epochs')
    parser.add_argument('--batch-size', '-b', dest='batch_size', metavar='B', type=int, default=1, help='Batch size')
    parser.add_argument('--learning-rate', '-l', metavar='LR', type=float, default=1e-3,
                        help='Learning rate', dest='lr')
    parser.add_argument('--load', '-f', type=str, default=False, help='Load model from a .pth file')
    parser.add_argument('--scale', '-s', type=float, default=1, help='Downscaling factor of the images')
    parser.add_argument('--validation', '-v', dest='val', type=float, default=10.0,
                        help='Percent of the data that is used as validation (0-100)')
    parser.add_argument('--amp', action='store_true', default=False, help='Use mixed precision')
    parser.add_argument('--bilinear', action='store_true', default=False, help='Use bilinear upsampling')
    parser.add_argument('--classes', '-c', type=int, default=1, help='Number of classes')
    parser.add_argument('--channels','-C',type=int, default=1, help='Number of channels')
    parser.add_argument('--dir_imgs', '-I', type=str, default=False, help='directory of ground truth images',required=True)
    parser.add_argument('--dir_masks', '-M', type=str, default=False, help='directory of ground truth masks',required=True)
    parser.add_argument('--dir_checkpoints', '-check', type=str, default='./checkpoints/', help='directory of saving checkpoint model',required=True)

    return parser.parse_args()

def remove_module_key(state_dict):
    new_state_dict = {}
    for key, value in state_dict.items():
        if key.startswith('module.'):
            new_key = key[7:] # remove 'module.' from the key
        else:
            new_key = key
        new_state_dict[new_key] = value
    return new_state_dict

def setup_run(local_rank):
    if local_rank == 0:
        run = wandb.init(
            project='Parallel U-Net',
            anonymous='must',
        )
    else:
        run = None
    
    return run

def cleanup():
    dist.destroy_process_group()


if __name__ == '__main__':

    local_rank = int(os.environ["LOCAL_RANK"])
    # local_rank = 0

    args = get_args()
    dir_img = Path(args.dir_imgs)
    dir_mask = Path(args.dir_masks)
    dir_checkpoint = Path(args.dir_checkpoints)

    if local_rank == 0:
        logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
        # device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        # logging.info(f'Using device {device}')

    # Change here to adapt to your data
    # n_channels=3 for RGB images
    # n_classes is the number of probabilities you want to get per pixel
    model = UNet(n_channels=args.channels, n_classes=args.classes, bilinear=args.bilinear)
    model = model.to(memory_format=torch.channels_last)

    # logging.info(f'Network:\n'
    #              f'\t{model.n_channels} input channels\n'
    #              f'\t{model.n_classes} output channels (classes)\n'
    #              f'\t{"Bilinear" if model.bilinear else "Transposed conv"} upscaling')

    if args.load:                
        state_dict = torch.load(args.load, map_location=torch.device('cuda', local_rank))
        state_dict = remove_module_key(state_dict) # Remove module. prefix from keys
        del state_dict['mask_values']
        model.load_state_dict(state_dict)
        logging.info(f'Model loaded from {args.load}')

    model.to(local_rank)
    model.n_channels = args.channels
    model.n_classes = args.classes    

    # wandb.init a run if logging, otherwise return None
    run = setup_run(local_rank)
    # run = None

    try:
        train_model(
            model=model,
            dir_img = dir_img,
            dir_mask = dir_mask,
            dir_checkpoint = dir_checkpoint,
            run = run,
            epochs=args.epochs,
            batch_size=args.batch_size,
            learning_rate=args.lr,
            device=local_rank,
            img_scale=args.scale,
            val_percent=args.val / 100,
            amp=args.amp
        )
    except torch.cuda.OutOfMemoryError:
        logging.error('Detected OutOfMemoryError!')
        torch.cuda.empty_cache()
        model.use_checkpointing()
        train_model(
            model=model,
            dir_img = dir_img,
            dir_mask = dir_mask,
            dir_checkpoint = dir_checkpoint,
            run = run,
            epochs=args.epochs,
            batch_size=args.batch_size,
            learning_rate=args.lr,
            device=local_rank,
            img_scale=args.scale,
            val_percent=args.val / 100,
            amp=args.amp
        )
