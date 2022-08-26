Beginner tutorial
=====================================
This tutorial will introduce some of the basic functionalities of HAT and generate some basic plots. 

.. toctree::
 :maxdepth: 1

Downloading and extracting example dataset
````````````````

An example dataset can be downloaded from the following addesss:

https://figshare.com/articles/dataset/Au_111_in_NaOH_High-energy_surface_x-ray_diffraction_dataset/20160632

The files are compressed with the 7zip format, in windows you will need software that can extract such files. Such a tool can be downloaded from here:

https://www.7-zip.org/download.html

In linux the following terminal command will download the dataset.

``wget https://figshare.com/ndownloader/files/36048185``

 in linux you can usually use the 7za command to uncompress the file. Firstly check it is available::

  whereis 7za
  7za: /usr/bin/7za /usr/share/man/man1/7za.1.gz
  
You should then extract the data using the following command (replacing the number by the actualy filename):
 
``7za e ./36048185``

If it isn't found then you need to install it.

Install p7zip to unzip *.7z files on Debian and Ubuntu
......

``$ sudo apt-get install p7zip``

Install p7zip to unzip *.7z files on Fedora
......

``$ yum install p7zip``


Loading the data
````````````````
Start HAT, if you installed hat with pip (and you're in the correct enviroment) then at the terminal you can simply type:

``xrayhat``

Otherwise from the source directory you can run (you mayneed to replace python with python3):

``python main.py``

Next you should change the beamline preset in the parameter tree (right side of the window) to P07 2019. It is generally recommend you use the manual beamline preset as beamlines seem to be constantly chaning their file formats and HAT is unlikely to keep up. However, the chosen preset is known to work with this data. 

Next click on |File -> Select Data Images| and then find and select the .TIF images that you should have extracted above.

They should have names such as: "au_111_naoh_2_00023-00001.tif", one of the files is a dark image and this can be ignored. The compressed archive also contains .metadata files for each image, the P07 2019 preset will automatically handel this metadata. In particular it will read the angle values and an intensity monitor value. 

Next several parameters should be chosen under Experiment in the parameter tree.
The energy = 73.7 keV
Sample-Detector Dist. = 1.6m

Center Pixel X = 1003
Center Pixel Y = 2024

The number of pixels and size the pixel should be chosen correctly since we alreay selected the correct preset. The angle offset is to be determined later. But we no from the experimental geomtery that the direction is negative. 

Under crystall we can choose the Au(111) preset since the data was collected from a Au(111) sample. These will set the lattice parameters to the correct values for Au(111) in surface coordinates. 

In Data Proccessing we should select 180 degrees for Image rotation as the images are rotated. The binning value will depend on how much system RAM  is available, but since we are only quickly exploring the data - setting it to value of 4 or higher is a good idea.

Unless you know you have a GPU that is correctly setup acceleration should be set to "Numba (cpu)"

Now we can load the data, although it could be a good idea to save the file selection and parameters first so you can restore them if needed (File -> Quick Save). 

To load the data click on (File -> Load)


Detector View
````````````````
While the images are loading the view will update after every 100 images. After loading you should have a view that looks like the screenshot below. This is called the detecor view mode and can be accessed at anytime by clicking (View -> Detector View) or pressing CTRL+1.

.. image:: ./images/screenshot1.png


In this view the image is construsted from maximum intensity of each pixel accross the whole rotation. If instead you'd like the average intensity across the whole rotation you can do this by selecting "Mean Images instead of Max" under Data Processing in the parameter tree. 

From the toolbar, masks can be added. Mask are attached the view they are created in, for detector view the masks are exclusionary and can be used to reject pixels from later binning. This could be useful for excluding gaps between detector pixels or areas of the detector that are covered. 

If the "Toggle ROI" button is pressed on the toolbar a box profile is created and the profile shown in a panel below the detector view. The profile runs between the two draggable circles. "Save ROI" will save the shape and position of the ROI vox so the same profile can be used for multiple datasets. "Extract ROI" will save the profile to a csv file.

Transformed Detector View
````````````````



