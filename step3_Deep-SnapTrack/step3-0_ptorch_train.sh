# python ptorch_Training_single.py \
#     --dir_img /mnt/f58069a5-1cf3-43b8-bb9b-ea74327327c9/WZH-DataCenter/PROCESS-SPT/2023/simPSF_results/DeepSTORMDataset/20231219_SNR19-35_NoLocError_110nmPixSize_deepStorm/imgs \
#     --dir_mask /mnt/f58069a5-1cf3-43b8-bb9b-ea74327327c9/WZH-DataCenter/PROCESS-SPT/2023/simPSF_results/DeepSTORMDataset/20231219_SNR19-35_NoLocError_110nmPixSize_deepStorm/locHeatmap \
#     --dir_checkpoint ./checkpoints/ \
#     --weights_name MBX_20231219_110nmPix_loc \
#     --epochs 100 \
#     --batch-size 12

# python ptorch_Training_single.py \
#     --dir_img /mnt/f58069a5-1cf3-43b8-bb9b-ea74327327c9/WZH-DataCenter/PROCESS-SPT/2023/simPSF_results/DeepSTORMDataset/20231219_SNR19-35_NoLocError_110nmPixSize_deepStorm/imgs \
#     --dir_mask /mnt/f58069a5-1cf3-43b8-bb9b-ea74327327c9/WZH-DataCenter/PROCESS-SPT/2023/simPSF_results/DeepSTORMDataset/20231219_SNR19-35_NoLocError_110nmPixSize_deepStorm/trackHeatmap \
#     --dir_checkpoint ./checkpoints/ \
#     --weights_name MBX_20231219_110nmPix_track \
#     --epochs 100 \
#     --batch-size 12

# OMP_NUM_THREADS=1 torchrun \
#     --standalone \
#     --nnodes 1 \
#     --nproc_per_node 2 \
# ptorch_Training_parallel.py \
#     --project Deep-STORM \
#     --dir_img /mnt/f58069a5-1cf3-43b8-bb9b-ea74327327c9/WZH-DataCenter/PROCESS-SPT/2023/simPSF_results/DeepSTORMDataset/20231220_SNR19-35_NoLocError_110nmPixSize_deepStorm/imgs \
#     --dir_mask /mnt/f58069a5-1cf3-43b8-bb9b-ea74327327c9/WZH-DataCenter/PROCESS-SPT/2023/simPSF_results/DeepSTORMDataset/20231220_SNR19-35_NoLocError_110nmPixSize_deepStorm/trackHeatmap \
#     --dir_checkpoint ./checkpoints/ \
#     --weights_name MBX_20231220_110nmPix_rep3 \
#     --epochs 30 \
#     --batch-size 12

# # No transfer learning
# OMP_NUM_THREADS=1 torchrun \
#     --standalone \
#     --nnodes 1 \
#     --nproc_per_node 2 \
# ptorch_Training_parallel.py \
#     --project Deep-STORM \
#     --dir_img /mnt/f58069a5-1cf3-43b8-bb9b-ea74327327c9/WZH-DataCenter/PROCESS-SPT/2023/simPSF_results/DeepSTORMDataset/20240319_SNR19-35_NoLocError_110nmPixSize_log10D/imgs \
#     --dir_mask /mnt/f58069a5-1cf3-43b8-bb9b-ea74327327c9/WZH-DataCenter/PROCESS-SPT/2023/simPSF_results/DeepSTORMDataset/20240319_SNR19-35_NoLocError_110nmPixSize_log10D/trackHeatmap \
#     --dir_checkpoint ./checkpoints/ \
#     --weights_name MBX_20240319_110nmPix_track \
#     --epochs 30 \
#     --batch-size 12

# # With transfer learning
# OMP_NUM_THREADS=1 torchrun \
#     --standalone \
#     --nnodes 1 \
#     --nproc_per_node 2 \
# ptorch_Training_parallel.py \
#     --project Deep-STORM \
#     --dir_img /mnt/disk1/WZH-DataCenter/PROCESS-SPT/2023/simPSF_results/DeepSTORMDataset/20240323_SNR19-35_NoLocError_110nmPixSize_log10D/imgs \
#     --dir_mask /mnt/disk1/WZH-DataCenter/PROCESS-SPT/2023/simPSF_results/DeepSTORMDataset/20240323_SNR19-35_NoLocError_110nmPixSize_log10D/trackHeatmap \
#     --dir_checkpoint ./checkpoints/ \
#     --weights_name MBX_20240323_110nmPix_track \
#     --epochs 30 \
#     --batch-size 12

# # With transfer learning
# OMP_NUM_THREADS=1 torchrun \
#     --standalone \
#     --nnodes 1 \
#     --nproc_per_node 2 \
# ptorch_Training_parallel.py \
#     --project Deep-STORM \
#     --dir_img /mnt/disk1/WZH-DataCenter/PROCESS-SPT/2023/simPSF_results/DeepSTORMDataset/20240521_SNR19-35_NoLocError_160nmPixSize_log10D/imgs \
#     --dir_mask /mnt/disk1/WZH-DataCenter/PROCESS-SPT/2023/simPSF_results/DeepSTORMDataset/20240521_SNR19-35_NoLocError_160nmPixSize_log10D/trackHeatmap \
#     --dir_checkpoint ./checkpoints/ \
#     --weights_name MBX_20240521_160nmPix_track \
#     --epochs 30 \
#     --batch-size 12

# With transfer learning
OMP_NUM_THREADS=1 torchrun \
    --standalone \
    --nnodes 1 \
    --nproc_per_node 2 \
ptorch_Training_parallel.py \
    --project Deep-STORM \
    --dir_img /mnt/disk1/WZH-DataCenter/PROCESS-SPT/2023/simPSF_results/DeepSTORMDataset/20240808_SNR19-35_NoLocError_110nmPixSize_log10D_10ms/imgs \
    --dir_mask /mnt/disk1/WZH-DataCenter/PROCESS-SPT/2023/simPSF_results/DeepSTORMDataset/20240808_SNR19-35_NoLocError_110nmPixSize_log10D_10ms/trackHeatmap \
    --dir_checkpoint ./checkpoints/ \
    --weights_name MBX_20240808_110nmPix_trackHeatmap \
    --epochs 30 \
    --batch-size 12 \
    --gpus 2 \
    --RunName MBX_20240808_110nmPix_trackHeatmap

