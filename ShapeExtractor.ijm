for (i=0;i<50;i=i+1)
{
open("/Users/dimvardal/Desktop/MAPV1/22122020/Series016/Series016.tif");
roiManager("Open", "/Users/dimvardal/Desktop/MAPV1/22122020/Series016/RoiSet.zip");
roiManager("Select", i);
run("Median...", "radius=2 slice");
setAutoThreshold("Percentile dark");
setThreshold(755, 65535);
run("Analyze Particles...", "size=80-Infinity pixel show=Outlines display clear add slice");
roiManager("Select", 0);
roiManager("Rename", i);

roiManager("Save","/Users/dimvardal/Desktop/MAPV1/22122020/Series016/b/"+i+".roi");
selectWindow("Series016.tif");
close();
}

