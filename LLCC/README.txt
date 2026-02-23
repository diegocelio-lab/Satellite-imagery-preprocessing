LOGLOCAL COLOR CORRECTION ALGORITHM (LLCC)

Version 2.0 - April 22, 2020
by Jose-Luis Lisani <joseluis.lisani@uib.es>


ASSOCIATED PUBLICATION

This code is associated to the following IPOL publication:

J.L. Lisani, 
Local Contrast Enhancement based on Adaptive Logarithmic Mappings
Image Processing On Line (IPOL), www.ipol.im

An online demo that uses the code is available at IPOL


ORGANIZATION OF THE CODE

The code is organized as follows:

- a folder 'library' containing the following: 
	- io_RGB.h, io_RGB.cpp: functions for the management of
	the RGB channels of an image and the computation of the 
	intensity component
	- io_png folder: containing code for reading/writing
	PNG images
	- gaussianGetreuer: code that implements the Gaussian filter,
	created by Pascal Getreuer and obtained from 
	http://www.ipol.im/pub/art/2013/87/
	- MCMMondelli: code that implements the MCM filter,
	created by Marco Mondelli and obtained from 
	http://www.ipol.im/pub/art/2013/53/
	- TRUNCATED_KERNEL_BF: code that implements the fast bilateral
	filter, created by Sylvain Paris and Frédo Durand and
	obtained from http://people.csail.mit.edu/sparis/bf/
- llcc.cpp: the main function
- makefile: for compilation
- README.txt: this file
- agpl-3.0.txt: GNU Affero General Public License 
- iris.png: an image to test the algorithm

The algorithm described in the associated paper is implemented 
in llcc.cpp. This code have been submitted for review at IPOL.


COMPILATION

1) Decompress code and access to folder:
tar xvzf llcc.tgz
cd llcc

2) Compilation:
make 

The executable file obtained after compilation is 'llcc'.



USAGE

Usage: llcc input.png output.png [option=1] [sigmaS=5] [sigmaI=70/gradth=10]
       option == 1 --> Use Gaussian filter
       option == 2 --> Use MCM filter
       option == 3 --> Use Bilinear filter
       Parameter of Gaussian filter: sigmaS (spatial scale)
       Parameter of MCM filter: sigmaS (spatial scale)
                              : gradth (gradient threshold)
       Parameter of Bilinear filter: sigmaI (intensity scale)
                                     sigmaS (spatial scale)


IMPORTANT NOTES

- The program only admits PNG input images.
The output images are also in this format.
In Linux or MacOS you can use the ImageMagicks'convert' 
function to convert any graphic file to PNG format.
- The code has been compiled and tested using a Linux environment
 

EXAMPLES OF USE

The following commands apply MLHE to the provided image 
using the three possible processing options for the weight map:

1) Gaussian weight map, with scale parameter 5

llcc iris.png output.png 1 5


2) MCM weight map, with normalized scale parameter 5 and 
gradient threshold 10

llcc iris.png output.png 2 5 10


3) Bilateral weight map, with scale parameter 5 and 
range parameter 70

llcc iris.png output.png 3 5 70



COPYRIGHT AND LICENSE

Copyright (c) 2020 J.L. Lisani

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.


THANKS

We would be grateful to receive any comment, especially about errors, bugs,
or strange results. Address them to joseluis.lisani@uib.es





