
run("Make Substack...", "slices=20-32");
run("StackReg ", "transformation=[Rigid Body]");
//saveAs("Tiff", "F:/8-20230810/1-Mouse9/1-M9_N2-FOV2_001/8-branch4-stack.tif");

run("Gamma...", "value=1.50 stack");
run("Gaussian Blur...", "sigma=1 stack");
run("Z Project...", "projection=[Max Intensity]");

//saveAs("Tiff", "F:/5-20230807/4-M7/1-M7-N1-FOV1/5-branch7-max.tif");


