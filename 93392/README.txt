INTRO
======================================================================

Infer full compartmental model given only access to the voltage in the
compartments. This code is released in conjunction with the paper 

      Huys QJM, Ahrens M and Paninski L (2006): Efficient estimation of
      detailed single-neurone models

and can be downloaded from 

      http://www.gatsby.ucl.ac.uk/~qhuys/code.html

The paper can be downloaded from

      http://www.gatsby.ucl.ac.uk/~qhuys/pub/hap06.pdf

Copyright Quentin Huys 2006
qhuys@gatsby.ucl.ac.uk


HOW TO
======================================================================

To use the code, unzip it, eg on a linux machine type
	
	gunzip hap06_code.gz

Which will create a directory with all the files. 

From within Matlab, change to that directory by typing eg
	
	cd hap06_code

Edit the file PARAM.M to change any parameters you want, like the number of
compartments in the cell, the size of each compartment, the amount of noise etc.
You should not have to edit any of the other files. 

To run the inference, simply type
	
	main

and hit ENTER. Enjoy. 


BUGS
======================================================================

I don't know of any bugs at the moment, but please do let me know if you find
any (qhuys@gatsby.ucl.ac.uk)
