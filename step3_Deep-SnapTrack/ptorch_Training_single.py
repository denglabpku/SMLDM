# import debugpy

# # 5678 is the default attach port in the VS Code debug configurations. Unless a host and port are specified, host defaults to 127.0.0.1
# debugpy.listen(5678)
# print("Waiting for debugger attach")
# debugpy.wait_for_client()
# debugpy.breakpoint()
# print('break on this line')

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
# from torchvision import transforms
# from torchvision.transforms.functional import InterpolationMode
from torch.utils.data import DataLoader, Dataset,random_split
import argparse
import os
from os import listdir
from os.path import splitext, isfile, join
from ptorch_CNN_model import L1L2loss, Model as buildModel
import wandb
from pathlib import Path
from PIL import Image
from tqdm import tqdm

class IndividualNormalize(object):

    def __call__(self, sample):
        img = np.array(sample)

        # Define a function that projects an image to the range [0,1]
        img = np.squeeze(img)
        min_val = img.min()
        max_val = img.max()
        img = (img - min_val) / (max_val - min_val)

        # Normalize image given mean and std
        mean_img = img.mean()
        std_img = img.std()
        img = (img - mean_img) / std_img    

        img = img.reshape(1, img.shape[0], img.shape[1])  # Reshape to (1, height, width)    
        
        return img

class BasicDataset(Dataset):
    def __init__(self, images_dir: str, mask_dir: str, mask_suffix: str = '_trackHeatmap', transform=None):
        self.images_dir = Path(images_dir)
        self.mask_dir = Path(mask_dir)
        self.mask_suffix = mask_suffix
        self.transform = transform

        self.ids = [splitext(file)[0] for file in listdir(images_dir) if isfile(join(images_dir, file)) and not file.startswith('.')]
        if not self.ids:
            raise RuntimeError(f'No input file found in {images_dir}, make sure you put your images there')
        
    def __len__(self):
        return len(self.ids)

    def __getitem__(self, idx):
        name = self.ids[idx]
        mask_file = list(self.mask_dir.glob(name + self.mask_suffix + '.*'))
        img_file = list(self.images_dir.glob(name + '.*'))

        assert len(img_file) == 1, f'Either no image or multiple images found for the ID {name}: {img_file}'
        assert len(mask_file) == 1, f'Either no mask or multiple masks found for the ID {name}: {mask_file}'
        mask = Image.open(mask_file[0]).convert('I')  # Convert to grayscale
        img = Image.open(img_file[0]).convert('I')  # Convert to grayscale

        # Upsample the image to match the size of the mask using 'nearest' interpolation
        img = img.resize(mask.size, resample=Image.NEAREST)
        
        if self.transform is not None:
            img = self.transform(img)
            mask = self.transform(mask)
       
        return {
            'image': torch.as_tensor(img.copy()).float().contiguous(),
            'mask': torch.as_tensor(mask.copy()).float().contiguous(),
        }

# define a function that trains a model for a given data SNR and density
def train_model(
        dir_img,
        dir_mask,
        dir_checkpoint,
        weights_name, 
        batch_size,
        epochs,
        lr,
        weight_file=None):
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f'Using device {device}')

    transform_norm = IndividualNormalize()
 
    # 1. Create dataset
    dataset = BasicDataset(dir_img, dir_mask, transform=transform_norm)

    # 2. Split into train / validation partitions
    val_percent = 0.3
    n_val = int(len(dataset) * val_percent)
    n_train = len(dataset) - n_val
    train_set, val_set = random_split(dataset, [n_train, n_val], generator=torch.Generator().manual_seed(0))
    # train_set = dataset
    # val_set = dataset

    # 3. Create data loaders
    loader_args = dict(batch_size=batch_size, num_workers=os.cpu_count(), pin_memory=True)
    train_loader = DataLoader(train_set, shuffle=True, **loader_args)
    val_loader = DataLoader(val_set, shuffle=False, drop_last=True, **loader_args)

    # Example input dimensions
    input_dim = (batch_size, 1, 320, 320)

    # Instantiate the model and set up criterion, optimizer, and scheduler
    model = buildModel()
    criterion = L1L2loss(input_shape=input_dim, device=device)
    optimizer = optim.Adam(model.parameters(), lr=lr)
    factor=0.5
    patience=5
    min_lr=0.00005
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=factor, patience=patience, min_lr=min_lr)
    # change_lr = ReduceLROnPlateau(monitor='val_loss', factor=0.1, patience=5, min_lr=0.00005) # from Colab Deep Storm

    if weight_file is not None:
        # Ensure that the model and state_dict share similar device placement
        state_dict = torch.load(weight_file, map_location=device)
        
        # If the model is wrapped with DistributedDataParallel, use model.module
        if isinstance(model, nn.DataParallel) or isinstance(model, nn.parallel.DistributedDataParallel):
            model.module.load_state_dict(state_dict)
        else:
            model.load_state_dict(state_dict)
            
        print(f'Model loaded from {weight_file}')


    # Move the model to the selected device
    model.to(device)

    # Inform user training begun
    print('Training model...')

    # (Initialize logging)
    experiment = wandb.init(project='Deep-STORM', resume='allow', anonymous='must')
    experiment.config.update(
        dict(epochs=epochs, batch_size=batch_size, learning_rate=lr, val_percent=val_percent,
             factor = factor, patience = patience, min_learning_rate = min_lr)
    )

    global_train_step = 0
    global_val_step = 0

    # New: Save the model weights after each epoch if the validation loss decreased
    best_val_loss = float('inf')

    for epoch in range(1, epochs + 1):
    
        # Training mode
        model.train()

        batch_train_loss = []
        t = len(train_loader)

        with tqdm(total=len(train_loader), desc=f'Epoch {epoch}/{epochs}', unit='batch') as pbar:        
            for batch in train_loader:
                # X_train_norm, Y_train = batch['image'], batch['mask']
                X_train_norm, Y_train = batch['image'].to(device), batch['mask'].to(device)

                optimizer.zero_grad()
                output = model(X_train_norm)
                # loss = criterion(output, Y_train)
                loss = criterion(Y_train, output)
                loss.backward()
                optimizer.step()

                batch_train_loss.append(loss.item())

                global_train_step += 1
                pbar.update()

                t -= 1
                if t == 1:
                    experiment.log({
                    "step train loss": loss.item(),
                    "train step": global_train_step,
                    'learning rate': optimizer.param_groups[0]['lr'],  # Log the current learning rate
                    })
        
        experiment.log({"epoch": epoch, 
                        "epoch train_loss": np.mean(batch_train_loss)})

        # Evaluation mode
        model.eval()
        val_losses = []  # List to store validation losses
        t = len(val_loader)
        
        with tqdm(total=len(val_loader), desc=f'Validation', unit='batch') as pbar:
            for batch in val_loader:
                # X_test_norm, Y_test = batch['image'], batch['mask']
                X_test_norm, Y_test = batch['image'].to(device), batch['mask'].to(device)
                val_output = model(X_test_norm)
                # val_loss = criterion(val_output, Y_test)
                val_loss = criterion(Y_test, val_output)
                val_losses.append(val_loss.item())  # Store the validation loss for logging

                scheduler.step(val_loss)  # Step the scheduler after each batch to update the learning rate

                pbar.update()

                global_val_step += 1

                t -= 1

                if t == 1:
                    experiment.log({
                        "step val loss": val_loss.item(),
                        "val step": global_val_step,
                        'images': [wandb.Image(X_test_norm[i].cpu()) for i in range(min(5, X_test_norm.shape[0]))],
                        'masks': {
                            'true': [wandb.Image(Y_test[i].float().cpu()) for i in range(min(5, Y_test.shape[0]))],
                            'pred': [wandb.Image(val_output[i, 0].cpu()) for i in range(min(5, val_output.shape[0]))],# Assuming val_output is a single-channel image
                            },
                        })
                    
        mean_val_loss = sum(val_losses) / len(val_losses)
        
        # Save model weights if validation loss improved
        if mean_val_loss < best_val_loss:
            Path(dir_checkpoint).mkdir(parents=True, exist_ok=True)
            best_val_loss = mean_val_loss
            torch.save(model.state_dict(), f'{dir_checkpoint}/{weights_name}_epoch{epoch}.pth')
            experiment.log({"epoch": epoch, 
                            "epoch val loss": mean_val_loss})
        
    # Finish the run
    wandb.finish()
        
    # Inform user training ended
    print('Training Completed!')

    return

if __name__ == '__main__':

    # start a parser
    parser = argparse.ArgumentParser()

    # path of the training data: patches and heatmaps, created in MATLAB using
    # the function "GenerateTrainingExamples.m"
    parser.add_argument('--dir_img', '-I', type=str, default=False, help='directory of ground truth images',required=True)
    parser.add_argument('--dir_mask', '-M', type=str, default=False, help='directory of ground truth masks',required=True)
    parser.add_argument('--dir_checkpoint', '-check', type=str, default='./checkpoints/', help='directory of saving checkpoint model',required=True)

    # path for saving the optimal model weights and normalization factors after
    # training with the function "train_model.py" is completed.
    parser.add_argument('--weights_name', type=str, help="path to save model weights as hdf5-file")
    parser.add_argument('--epochs', '-e', metavar='E', type=int, default=5, help='Number of epochs')
    parser.add_argument('--batch-size', '-b', dest='batch_size', metavar='B', type=int, default=1, help='Batch size')
    parser.add_argument('--learning-rate', '-l', metavar='LR', type=float, default=1e-3, help='Learning rate', dest='lr')
    parser.add_argument('--weight_file', type=str, help="path to load previous model weights")

    # parse the input arguments
    args = parser.parse_args()
    dir_img = Path(args.dir_img)
    dir_mask = Path(args.dir_mask)
    dir_checkpoint = Path(args.dir_checkpoint)
    weights_name = Path(args.weights_name)
    
    # run the training process
    train_model(dir_img = dir_img,
                dir_mask = dir_mask,
                dir_checkpoint = dir_checkpoint,
                weights_name = weights_name,
                batch_size = args.batch_size,
                epochs = args.epochs,
                lr = args.lr,
                weight_file = args.weight_file)
    

