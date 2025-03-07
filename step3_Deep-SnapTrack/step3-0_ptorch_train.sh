# Training Deep-SnapTrack
OMP_NUM_THREADS=1 torchrun \
    --standalone \
    --nnodes 1 \
    --nproc_per_node 2 \
ptorch_Training_parallel.py \
    --project Deep-STORM \
    --dir_img /path/to/save/simulated_dataset_forDeepSnapTrack/imgs \
    --dir_mask /path/to/save/simulated_dataset_forDeepSnapTrack/trackHeatmap \
    --dir_checkpoint ./checkpoints/ \
    --weights_name MBX_20231220_110nmPix_rep2 \
    --epochs 30 \
    --batch-size 12 \
    --gpus 2 \
    --RunName MBX_20231220_110nmPix_rep2

