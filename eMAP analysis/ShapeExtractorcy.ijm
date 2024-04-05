for (i=0;i<10;i=i+1) 
{
open("/Users/cyaeger/Desktop/DVharnett/09092021/Series033/Series033.tif");
roiManager("Open", "/Users/cyaeger/Desktop/DVharnett/09092021/Series033/RoiSet.zip");
roiManager("Select", i);
run("Median...", "radius=2 slice");
setAutoThreshold("Percentile dark");
setThreshold(350, 65535);
run("Analyze Particles...", "size=80-Infinity pixel show=Outlines display clear add slice");
roiManager("Select", 0);
roiManager("Rename", i);

roiManager("Save","/Users/cyaeger/Desktop/DVharnett/09092021/Series033/"+i+".roi");
selectWindow("Series033.tif");
close();
}

{while (nImages>0) { 
          selectImage(nImages); 
          close(); 
}
