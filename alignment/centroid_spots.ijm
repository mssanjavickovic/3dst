input = getDirectory("Input directory");
output = getDirectory("Output directory");

suffix = "spots.jpg";  //you only want to apply to Mask.jpg images, no need to ask

processFolder(input);

function processFolder(input) {
    list = getFileList(input);
    for (i = 0; i < list.length; i++) {
        if(File.isDirectory(input + list[i]))   //if it's a directory, go to subfolder
            processFolder("" + input + list[i]);   
        if(endsWith(list[i], suffix))   //if it's a jpg image, process it
            processFile(input, output, list[i]);
            close();  //close image
        //if it's neither a tiff nor a directory, do nothing
    }
}

function processFile(input, output, file) {
    //here, you put the real processing
    print("Processing: " + input + file);
    open(input + file);  //open image
	setAutoThreshold("Default");
	//run("Threshold...");
	setAutoThreshold("Default");
	//setThreshold(0, 143);
	setOption("BlackBackground", false);
	run("Convert to Mask");
	setThreshold(255, 255);
	run("Set Measurements...", "centroid redirect=None decimal=3");
	run("Analyze Particles...", "size=3000-Infinity show=Overlay display");
	saveAs("results", output + file + "_spots.tsv");
	run("Clear Results");
}




