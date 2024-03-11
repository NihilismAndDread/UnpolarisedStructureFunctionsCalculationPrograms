The code used to create the unpolarised iso-singlet nucleon structure functions within the NJL soliton model.

This was the code I used during my masters project. It is quite lengthy and complicated. However, I will give a brief description in this readme.

I have also added my thesis, this will give a good description of the theory behind the programs.

****************************************
cutoff.py
****************************************
This is a basic python script used to attain the physically accurate lambda cut-off when given the pion decay constant fpi and the pion mass mpi and the number of colours Nc.

****************************************
constants.h
****************************************
This file defines the numerical and physical parameters used during the experiments, these can be changed as you see fit.

****************************************
SCS_builder.c
****************************************
Next we build the NJL soliton self-consistently. This can be quite an involved process but what we end up with is a file for the chiral angle (Data/expnum**/chiralangle_**.b where ** is the relevant experiment number) as an unformatted file. In addtion to this we get the eigen values and eigen vectors for the quark wave functions, the number of momenta  per grand spin and the energy values for wave functions with a non-soliton background. These are just spherical bessel functions.

*****************************************
routines.h
*****************************************
There is a variety of functions that need to be used throughout the program, this file stores some of them.

*****************************************
SCS_fourier_old.c
*****************************************
While this piece of code went through multiple revisions this was the final product, so to speak. This codes simply takes the quark wave functions and performs a fourier transformation on them. The information is stored in Data/expnum**/Fdata

*****************************************
SCS_UNPOLARISED_STRUCTURE_FUNCTION_OLD.C
*****************************************
This file is used to create the unpolarised structure functions as well as the sum rules
 at every quark energy level. These are then saved in the directories Data/expnum**/USFdata and Data/expnum**/SRdata respectively.

