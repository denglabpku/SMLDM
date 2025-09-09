# ===== U-Net segmentation with pretrained model on the 110 nm pixel size EMCCD camera  ====== 
model_version=MBX_20240620_epoch20
python ND2batch_prediction_latest.py \
    --model /path/to/your/model/checkpoint_UNet_epoch20.pth \
    --model_version=${model_version} \
    --mpalm_channel=1 \
    --input \
    '/path/to/your/Paxillin_raw_data/20240712_Clust01_U2OS_Paxillin_30p5ms_2kframe_001.nd2'\
    --output /path/to/save/Paxillin_results/











