// Mitopits quantification 
// Written by Meghane Sittewelle 3 dec 2021
// WARNING - need to adapt channel depending of acquisition!


close("*");
roiManager("reset");

#@ File    (label = "File to analyse", style = "file") namefile
#@ File    (label = "Saving directory", style = "directory") dirsave
dir=dirsave+File.separator

open(namefile);
namefilewoext=File.getNameWithoutExtension(namefile);

print("namefilewoext=",namefilewoext); 

//Duplicate, chose the cell, and split channels, close nucleus, rename beta2 and mito
run("Duplicate...", "duplicate")
namefile=getTitle();
Stack.setChannel(2);
run("Enhance Contrast", "saturated=0.5");
Stack.setChannel(1);
run("Enhance Contrast", "saturated=0.35");

waitForUser("Waiting for user to do draw selection around the cell to measure. Press Okay to continue....");	
roiManager("Add");
roiManager("select",0); 
roiManager("rename","entireCell"); 

roiManager("Save", dir+File.separator+namefilewoext+"_ROIentirecell.roi");

run("Split Channels");

//NEED TO ADAPT CHANNEL DEPENDING OF EXPERIMENTS
//Nucleus channel
selectWindow("C3-"+namefile); 		
close();

//Beta2 channel
selectWindow("C1-"+namefile); 
rename("beta2_"+namefilewoext);
beta2="beta2_"+namefilewoext

//Mitochondria channel
selectWindow("C2-"+namefile); 
rename("mito_"+namefilewoext);
mito="mito_"+namefilewoext

//If 4th non used channel
//selectWindow("C4-"+namefile); 
//close();

//Mitochondria masks
selectWindow(mito); 
run("Threshold...");
waitForUser("Waiting for user to do threshold on mitochondria signal. Press Okay to continue....");
run("Convert to Mask");
rename("mask_"+mito);
mask_mito="mask_"+mito

roiManager("select","entireCell"); 
run("Make Inverse");
run("Fill", "slice");

saveAs("Tiff", dir+File.separator+mask_mito+".tiff");
mask_mito="mask_"+mito+".tiff"

//Beta2 masks
selectWindow(beta2); 

//Comdet- need to be adapted for each acquisition
run("Detect Particles");

run("Create Mask");
rename("mask_"+beta2);
mask_beta2="mask_"+beta2
saveAs("Tiff", dir+File.separator+mask_beta2+".tiff");
mask_beta2="mask_"+beta2+".tiff"

roiManager("select","entireCell"); 
run("Analyze Particles...", "size=0-4 clear include summarize add"); 	//Count particles only into selection of the entire cell


//multiply masks
imageCalculator("Multiply create", mask_beta2, mask_mito);
rename("mitoSpots_"+namefilewoext);
mitospots="mitoSpots_"+namefilewoext
saveAs("Tiff", dir+File.separator+mitospots+".tiff");
mitospots="mitoSpots_"+namefilewoext+".tiff"

selectWindow(mitospots);
run("Analyze Particles...", "size=0-3 clear include summarize add");

//Save summary
selectWindow("Summary");
saveAs("Results", dir+File.separator+"Summary_"+namefilewoext+".csv");
