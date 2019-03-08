The program "damfits" is a simple python program for doing simple
operations of fits files as they are handled in the damic project.

Please see INSTALL for prerequisites and installation instructions.

Examples:

Simple invocation:

$ damfits

usage: damfits[options] afile.fits# use -h for help

Getting a little help:

$ damfits -h

damfits [-h|--help]

damfits file.fits #displays fits file info

damfits -p [extraOptions] file.fits #plots

extraOptions:

	-i:	Needs one integer as argument
 		so it can select the corresponding image

	-b:	Needs one integer for defining the bining

	--side:	Needs one argument specifying the left or right
		side of the image

	--noLog:	A simple flag for turning of the ylog scale

	--color:	Needs one argument defining the color

For seeing some info about the fits file:

$ damfits d44_snolab_Int-200_Exp-100000_4283.fits

Filename: d44_snolab_Int-200_Exp-100000_4283.fits

No.    Name      Ver    Type      Cards   Dimensions   Format

  0  PRIMARY       1 PrimaryHDU     105   ()

  1                1 ImageHDU       110   (8544, 4298)   int16 (rescales to uint16)

  2                1 ImageHDU       110   (8544, 4298)   int16 (rescales to uint16)

  3                1 ImageHDU       110   (8544, 4298)   int16 (rescales to uint16)

  4                1 ImageHDU       110   (8544, 4298)   int16 (rescales to uint16)

  5                1 ImageHDU       110   (8544, 4298)   int16 (rescales to uint16)

  6                1 ImageHDU       110   (8544, 4298)   int16 (rescales to uint16)

  7                1 ImageHDU       110   (8544, 4298)   int16 (rescales to uint16)

  8                1 ImageHDU       110   (8544, 4298)   int16 (rescales to uint16)

  9                1 ImageHDU       110   (8544, 4298)   int16 (rescales to uint16)

 10                1 ImageHDU       110   (8544, 4298)   int16 (rescales to uint16)

 11                1 ImageHDU       110   (8544, 4298)   int16 (rescales to uint16)

 12                1 ImageHDU       110   (8544, 4298)   int16 (rescales to uint16)


Try it out and post suggestions and comments.
