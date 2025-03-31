// Below is the ImageJ Macro to perform Thunder-STORM localization,
// which the results are used in step2-2 to keep segmentation containing only one single-molecule localization.

// version number: batchThunderStorm_ZW_latest_20231221
// Below params is designed for 110 nm pixel size EMCCD camera.
// You may need to adjust params if using camera with other pixel sizes.

dir_processing = getDirectory("Choose a Directory to process");
dir_saving = getDirectory("Choose a Directory to save");
list = getFileList(dir_processing);
channel_number = 1;

for (i = 0; i < list.length; i++) { //list.length
	if (endsWith(list[i], ".nd2") && list[i].indexOf("30p5ms") != -1) {
	// if (endsWith(list[i], ".tif")) {
		img_dir = dir_processing + list[i];
		sav_img_dir = dir_saving + list[i];
		run("Bio-Formats Importer", "open="+img_dir+" color_mode=Default rois_import=[ROI manager] split_channels view=Hyperstack stack_order=XYCZT");
		
		if (channel_number == 1) {
			// For single channel (Ch1: Halo-PAJF646), Ch1 use tighter cutoff for UNet segmentation QC. Also pixel size use 110 nm.
			// Do not allow multi fit, because from 1xNLS data, some single-molecule blurry signal splited into two locs.
			// Also set bigger sigma and fitradius to match the histogram distribution from trial run (mean value).
			selectWindow(list[i]+" - C=0");
			run("Camera setup", "offset=500.0 isemgain=true photons2adu=2 gainem=300.0 pixelsize=110.0");
			run("Run analysis", "filter=[Wavelet filter (B-Spline)] scale=2.0 order=3 detector=[Centroid of connected components] watershed=false threshold=2*std(Wave.F1) estimator=[PSF: Integrated Gaussian] sigma=3 fitradius=9 method=[Least squares] full_image_fitting=false mfaenabled=false renderer=[No Renderer]");
			run("Export results", "filepath="+sav_img_dir.substring(0, sav_img_dir.length-4)+"_C1.csv fileformat=[CSV (comma separated)] sigma=true intensity=true chi2=true offset=true saveprotocol=true x=true y=true bkgstd=true id=true uncertainty=true frame=true");
		}
				
		close("*");
	}
}
run("Quit");


