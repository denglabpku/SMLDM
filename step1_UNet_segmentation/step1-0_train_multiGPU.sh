# Script for training U-Net with multiple GPU
# used for training U-Net with two GPU in parallel, modify your GPU number using nproc_per_node
OMP_NUM_THREADS=1 torchrun \
    --standalone \
    --nnodes 1 \
    --nproc_per_node 2 \
train_multiGPU_parallel_latest.py \
    --batch-size 32 \
    --epochs 20 \
    --do_log \
    --dir_imgs /path/to/save/simulated_dataset_forUNet/images/ \
    --dir_masks /path/to/save/simulated_dataset_forUNet/masks/ \
    --dir_checkpoints ./checkpoints/    

# -- uncomment below -- #
# Script for training U-Net with one GPU
# python train_multiGPU_single_latest.py \
#     --batch-size 32 \
#     --epochs 20 \
#     --do_log \
#     --dir_imgs /path/to/save/simulated_dataset_forUNet/images/ \
#     --dir_masks /path/to/save/simulated_dataset_forUNet/masks/ \
#     --dir_checkpoints ./checkpoints/