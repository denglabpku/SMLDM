# ===== U-Net segmentation with pretrained model on the 110 nm pixel size EMCCD camera  ====== 
model_version=MBX_20240620_epoch13
python ND2batch_prediction_latest.py \
    --model /path/to/your/model/checkpoint_epoch13.pth \
    --model_version=${model_version} \
    --mpalm_channel=1 \
    --input \
    '/path/to/your/data.nd2'\
    --output /path/to/save/masks/











