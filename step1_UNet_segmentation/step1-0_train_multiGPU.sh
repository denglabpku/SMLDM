# used for denoised image training
OMP_NUM_THREADS=1 torchrun \
    --standalone \
    --nnodes 1 \
    --nproc_per_node 2 \
train_multiGPU_parallel_latest.py \
    --batch-size 32 \
    --epochs 20 \
    --dir_imgs /dataE/WZH-DataCenter/PROCESS-SPT/2023/simPSF_results/UNetDetectionDataset/20240620_SNR19-35_diffBackground_NoLocError_110nmPixSize_3nd/imgs/ \
    --dir_masks /dataE/WZH-DataCenter/PROCESS-SPT/2023/simPSF_results/UNetDetectionDataset/20240620_SNR19-35_diffBackground_NoLocError_110nmPixSize_3nd/masks/ \
    --dir_checkpoints ./checkpoints/    
