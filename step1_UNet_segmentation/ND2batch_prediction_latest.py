# import debugpy

# # 5678 is the default attach port in the VS Code debug configurations. Unless a host and port are specified, host defaults to 127.0.0.1
# debugpy.listen(5678)
# print("Waiting for debugger attach")
# debugpy.wait_for_client()
# debugpy.breakpoint()
# print('break on this line')



import argparse
import logging
import os

from tqdm import tqdm
import nd2
import tifffile as tiff
import numpy as np
import torch
import torch.nn.functional as F
from PIL import Image
from torchvision import transforms

from utils.data_loading import ParaBasicDataset
from unet import UNet
from utils.utils import plot_img_and_mask

def predict_img(net,
                full_img,
                device,
                scale_factor=1,
                out_threshold=0.5):
    net.eval()
    img = torch.from_numpy(ParaBasicDataset.preprocess(None, full_img, scale_factor, is_mask=False))
    img = img.unsqueeze(0)
    img = img.to(device=device, dtype=torch.float32)

    with torch.no_grad():
        output = net(img).cpu()
        output = F.interpolate(output, (full_img.size[1], full_img.size[0]), mode='bilinear')
        if net.n_classes > 1:
            mask = output.argmax(dim=1)
        else:
            mask = torch.sigmoid(output) > out_threshold

    return mask[0].long().squeeze().numpy()


def get_args():
    parser = argparse.ArgumentParser(description='Predict masks from input images')
    parser.add_argument('--model', '-m', default='MODEL.pth', metavar='FILE',
                        help='Specify the file in which the model is stored')
    parser.add_argument('--input', '-i', metavar='INPUT', nargs='+', help='Filenames of input images', required=True)
    parser.add_argument('--output', '-o', metavar='OUTPUT', nargs='+', help='Filenames of output images', required=True)
    parser.add_argument('--viz', '-v', action='store_true',
                        help='Visualize the images as they are processed')
    parser.add_argument('--no-save', '-n', action='store_true', help='Do not save the output masks')
    parser.add_argument('--mask-threshold', '-t', type=float, default=0.5,
                        help='Minimum probability value to consider a mask pixel white')
    parser.add_argument('--scale', '-s', type=float, default=1,
                        help='Scale factor for the input images')
    parser.add_argument('--bilinear', action='store_true', default=False, help='Use bilinear upsampling')
    parser.add_argument('--classes', '-c', type=int, default=1, help='Number of classes')
    parser.add_argument('--channels','-C',type=int, default=1, help='Number of channels')
    parser.add_argument('--model_version','-mv',help='Model versions that use as prefix of generated mask file',required=True)
    parser.add_argument('--mpalm_channel','-MC',type=int, default=2, help='Specify which channel for UNet processing if read a multi channel ND2')
    
    return parser.parse_args()

def get_output_filenames(args):
    output_folder = []  
    for input_path in args.input:
        input_ND2name = os.path.split(input_path)
        # output_path = os.path.join(args.output[0],f'{os.path.splitext(input_ND2name[1])[0]}',f"UNet_mask_{args.model_version}.tif")
        # output_path = os.path.join(args.output,f'{os.path.splitext(input_ND2name[1])[0]}'+f"_UNet_mask_{args.model_version}_Ch{args.mpalm_channel}.tif")
        output_path = os.path.join(args.output[0],f'{os.path.splitext(input_ND2name[1])[0]}',f"UNet_mask_{args.model_version}_Ch{args.mpalm_channel}.tif")
        output_folder.append(output_path)
    return output_folder


def mask_to_image(mask: np.ndarray, mask_values):
    if isinstance(mask_values[0], list):
        out = np.zeros((mask.shape[-2], mask.shape[-1], len(mask_values[0])), dtype=np.uint8)
    elif mask_values == [0, 1]:
        out = np.zeros((mask.shape[-2], mask.shape[-1]), dtype=bool)
    else:
        out = np.zeros((mask.shape[-2], mask.shape[-1]), dtype=np.uint8)

    if mask.ndim == 3:
        mask = np.argmax(mask, axis=0)

    for i, v in enumerate(mask_values):
        out[mask == i] = v

    return Image.fromarray(out)

def remove_module_key(state_dict):
    new_state_dict = {}
    for key, value in state_dict.items():
        if key.startswith('module.'):
            new_key = key[7:] # remove 'module.' from the key
        else:
            new_key = key
        new_state_dict[new_key] = value
    return new_state_dict


if __name__ == '__main__':
    os.environ["CUDA_VISIBLE_DEVICES"] = "0"  # Set to use the second GPU (indexing starts from 0)
    args = get_args()
    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

    in_files = args.input
    out_files = get_output_filenames(args)

    net = UNet(n_channels=args.channels, n_classes=args.classes, bilinear=args.bilinear)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    logging.info(f'Loading model {args.model}')
    logging.info(f'Using device {device}')

    net.to(device=device)
    state_dict = torch.load(args.model, map_location=device)
    state_dict = remove_module_key(state_dict) # Remove module. prefix from keys
    mask_values = state_dict.pop('mask_values', [0, 1])
    net.load_state_dict(state_dict)

    logging.info('Model loaded!')

    

    for i, filename in enumerate(in_files):
        logging.info(f'Predicting image {filename} ...')
        # ND2_img = tiff.imread(filename)
        ND2_img = nd2.imread(filename) # numpy array
        if ND2_img.ndim == 4:
            ND2_img = ND2_img[:, args.mpalm_channel-1, :, :]

        final_result = []
        for iFrame in tqdm(range(ND2_img.shape[0])): # loop through each frame
            # Check if frame is valid, sometimes frame may all zero, a bug in Nikon software
            # return a empty mask if true
            if np.all(ND2_img[iFrame, :] == 0):
                # Assuming ND2_img is a NumPy array
                mask_shape = (ND2_img.shape[1], ND2_img.shape[2])  # Shape matching the second and third dimensions

                # Create a mask with all zeros using numpy.zeros
                mask = np.zeros(mask_shape, dtype=np.int64)
            else:
                # Convert rescaled array into Pillow image object
                img = Image.fromarray(ND2_img[iFrame,...])
                mask = predict_img(net=net,
                                full_img=img,
                                scale_factor=args.scale,
                                out_threshold=args.mask_threshold,
                                device=device)

            if not args.no_save:
                result = mask_to_image(mask, mask_values)
                final_result.append(result)

            if args.viz:
                logging.info(f'Visualizing results for image {filename}, close to continue...')
                plot_img_and_mask(img, mask)
        
        # save GIF file from list of frames
        if not args.no_save:
            out_filename = out_files[i]
            out_dirname = os.path.dirname(out_filename)
            if not os.path.exists(out_dirname):
                os.makedirs(out_dirname)
                print(f"Mask output directory '{out_dirname}' created successfully.")
            else:
                print(f"Mask output directory '{out_dirname}' already exists.")
            final_result[0].save(out_filename, format="tiff", append_images=final_result[1:],save_all=True,compression="tiff_lzw") # tiff_lzw is a lossless compression, safe to use
            logging.info(f'Mask saved to {out_filename}')

