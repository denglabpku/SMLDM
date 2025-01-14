
# ===== U-Net segmentation (110 nm pixel size)  ====== 
model_version=MBX_20240620_2035_epoch20
python ND2batch_prediction_latest.py \
    --model /mnt/disk1/WZH-DataCenter/PROCESS-SPT/2023/simPSF_results/UNetDetectionDataset/pretrain_models/MBX_20240620_2035/checkpoint_epoch20.pth \
    --model_version=${model_version} \
    --mpalm_channel=1 \
    --input \
    '/dataB/SPT_raw_backup/ZuhuiWang/2024/20240602_U2OS_FOXA2_MPALM/20240602_Clust01_U2OS_FOXA2_25uMPA646_30p5ms_10kframe_001.nd2'\
    '/dataB/SPT_raw_backup/ZuhuiWang/2024/20240602_U2OS_FOXA2_MPALM/20240602_Clust02_U2OS_FOXA2_25uMPA646_30p5ms_10kframe_002.nd2'\
    '/dataB/SPT_raw_backup/ZuhuiWang/2024/20240602_U2OS_FOXA2_MPALM/20240602_Clust03_U2OS_FOXA2_25uMPA646_30p5ms_5kframe_003.nd2'\
    '/dataB/SPT_raw_backup/ZuhuiWang/2024/20240602_U2OS_FOXA2_MPALM/20240602_Clust04_U2OS_FOXA2_25uMPA646_30p5ms_5kframe_004.nd2'\
    '/dataB/SPT_raw_backup/ZuhuiWang/2024/20240602_U2OS_FOXA2_MPALM/20240602_Clust05_U2OS_FOXA2_25uMPA646_30p5ms_5kframe_005.nd2'\
    '/dataB/SPT_raw_backup/ZuhuiWang/2024/20240602_U2OS_FOXA2_MPALM/20240602_Clust06_U2OS_FOXA2_25uMPA646_30p5ms_5kframe_006.nd2'\
    '/dataB/SPT_raw_backup/ZuhuiWang/2024/20240602_U2OS_FOXA2_MPALM/20240602_Clust07_U2OS_FOXA2_25uMPA646_30p5ms_5kframe_007.nd2'\
    '/dataB/SPT_raw_backup/ZuhuiWang/2024/20240602_U2OS_FOXA2_MPALM/20240602_Clust08_U2OS_FOXA2_25uMPA646_30p5ms_5kframe_008.nd2'\
    '/dataB/SPT_raw_backup/ZuhuiWang/2024/20240602_U2OS_FOXA2_MPALM/20240602_Clust09_U2OS_FOXA2_25uMPA646_30p5ms_5kframe_009.nd2'\
    '/dataB/SPT_raw_backup/ZuhuiWang/2024/20240602_U2OS_FOXA2_MPALM/20240602_Clust10_U2OS_FOXA2_25uMPA646_30p5ms_5kframe_010.nd2'\
    '/dataB/SPT_raw_backup/ZuhuiWang/2024/20240602_U2OS_FOXA2_MPALM/20240602_Clust11_U2OS_FOXA2_25uMPA646_30p5ms_5kframe_011.nd2'\
    '/dataB/SPT_raw_backup/ZuhuiWang/2024/20240602_U2OS_FOXA2_MPALM/20240602_Clust12_U2OS_FOXA2_25uMPA646_30p5ms_5kframe_012.nd2'\
    '/dataB/SPT_raw_backup/ZuhuiWang/2024/20240602_U2OS_FOXA2_MPALM/20240602_Clust13_U2OS_FOXA2_25uMPA646_30p5ms_2kframe_001.nd2'\
    '/dataB/SPT_raw_backup/ZuhuiWang/2024/20240602_U2OS_FOXA2_MPALM/20240602_Clust13_U2OS_FOXA2_25uMPA646_30p5ms_2kframe_002.nd2'\
    '/dataB/SPT_raw_backup/ZuhuiWang/2024/20240602_U2OS_FOXA2_MPALM/20240602_Clust13_U2OS_FOXA2_25uMPA646_30p5ms_2kframe_003.nd2'\
    '/dataB/SPT_raw_backup/ZuhuiWang/2024/20240602_U2OS_FOXA2_MPALM/20240602_Clust13_U2OS_FOXA2_25uMPA646_30p5ms_2kframe_004.nd2'\
    '/dataB/SPT_raw_backup/ZuhuiWang/2024/20240602_U2OS_FOXA2_MPALM/20240602_Clust13_U2OS_FOXA2_25uMPA646_30p5ms_2kframe_005.nd2'\
    '/dataB/SPT_raw_backup/ZuhuiWang/2024/20240602_U2OS_FOXA2_MPALM/20240602_Clust13_U2OS_FOXA2_25uMPA646_30p5ms_2kframe_006.nd2'\
    '/dataB/SPT_raw_backup/ZuhuiWang/2024/20240602_U2OS_FOXA2_MPALM/20240602_Clust13_U2OS_FOXA2_25uMPA646_30p5ms_2kframe_007.nd2'\
    '/dataB/SPT_raw_backup/ZuhuiWang/2024/20240602_U2OS_FOXA2_MPALM/20240602_Clust13_U2OS_FOXA2_25uMPA646_30p5ms_2kframe_008.nd2'\
    '/dataB/SPT_raw_backup/ZuhuiWang/2024/20240602_U2OS_FOXA2_MPALM/20240602_Clust13_U2OS_FOXA2_25uMPA646_30p5ms_2kframe_009.nd2'\
    '/dataB/SPT_raw_backup/ZuhuiWang/2024/20240602_U2OS_FOXA2_MPALM/20240602_Clust13_U2OS_FOXA2_25uMPA646_30p5ms_2kframe_010.nd2'\
    --output /mnt/disk1/WZH-DataCenter/PROCESS-SPT/2024/20240602_U2OS_FOXA2_MPALM/











