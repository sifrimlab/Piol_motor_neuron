title = getTitle();
run("*To ROI Manager [F4]");

output_dir = "/path/to/dir/"
width = getWidth();
height = getHeight();

close();

newImage(title, "16-bit black", width, height, 1);

for (index = 0; index < roiManager("count"); index++) {
	roiManager("select", index);
	setColor(index+1);
	fill();
}

resetMinAndMax();
saveAs("Tiff", output_dir +  title + "_labeled.tif");

