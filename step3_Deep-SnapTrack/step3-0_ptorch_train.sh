# Training Deep-SnapTrack with multiple GPUs (e.g. 2 GPUs)
OMP_NUM_THREADS=1 torchrun \
    --standalone \
    --nnodes 1 \
    --nproc_per_node 2 \
ptorch_Training_parallel.py \
    --project Deep-SnapTrack \
    --dir_img /path/to/save/simulated_dataset_forDeepSnapTrack/imgs \
    --dir_mask /path/to/save/simulated_dataset_forDeepSnapTrack/trackHeatmap \
    --dir_checkpoint ./checkpoints/ \
    --weights_name MBX_20231220_110nmPix_rep2 \
    --epochs 30 \
    --batch-size 12 \
    --gpus 2 \
    --log_all \
    --RunName MBX_20231220_110nmPix_rep2

# # Training Deep-SnapTrack with single GPUs
# python ptorch_Training_single.py \
#     --project Deep-SnapTrack \
#     --dir_img /path/to/save/simulated_dataset_forDeepSnapTrack/imgs \
#     --dir_mask /path/to/save/simulated_dataset_forDeepSnapTrack/trackHeatmap \
#     --dir_checkpoint ./checkpoints/ \
#     --weights_name MBX_20231220_110nmPix_rep2 \
#     --epochs 30 \
#     --batch-size 12 \
#     --log_all \
#     --RunName MBX_20231220_110nmPix_rep2

