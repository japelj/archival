# Archival

A script to cross-check whether a transient source coincides with a galaxy, star, or a previously known transient source. The final product is an html document with a detailed summary.

### Prerequisites and installation

The script is written in Python 3. The following modules are required:
- astropy
- astroquery
- lxml
- matplotlib
- numpy
- scipy
- pandas
- pillow
- requests
- yattag

The easiest way to prepare the proper environment is with [anaconda](https://www.anaconda.com/download). After installing anaconda, create Python 3 environment
```
conda create -n py36 python=3.6
```
and then install the modules
```
conda install astropy, astroquery, matplotlib, numpy, pandas, pillow, requests
```
The final package has to be installed via pip:
```
pip install yattag
```
Alternatively, if you already have Python 3 installed, you can install each of the necessary packages with pip.

### How to use it 

The script is run from the folder with the input list as
```
python /path_to_script/archival.py -i input_list.txt [-p] [-d]
```
Input list is a list of lines (Name RA DEC [DATE]), such as  
GWT1	03:55:50.1	-45:06:20.0	2019-01-06.70  
GWT2	62.89175	-9.122556

DATE corresponds to the time of the detection of the source and is necessary for the search of coincident Solar system objects. An example of the input file is included.

Two optional inputs are a parameter file (eg. -p par_file.txt), to change the default searching parameters (such as the radius around the source in which objects are looked for; an example parameter file is included) and luminosity distance and error in Mpc (eg. -d 80 20).

The information (eg. plots, images, numbers) is collected and a summary is provided in an html document (default: output/summary.html). 

