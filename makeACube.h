/*
File:			makeACube.h

Description:	Header file for makeACube.cxx. 

				Function effectively sets the square from index value (0,0,0) to (50,50,50) of the input image to 0 and sets the square from index value (10,10,10) to (40,40,40) to 255.

Inputs:			OutputImage - (void *) Void pointer that points to the ITK::image type pointer that holds the image that will be modified by this function. Due to how this function was programmed, the pointer has to point to a "itk::Image< float, 3 >"  image.

Author:			Frankie (Hoi-Ki) Tong, Sep 6 2014

Changes:		Frankie (Hoi-Ki) Tong, Sep 6 2014
				Initial creation of the file.
*/

#include "itkImage.h"

#ifndef __MAKE_A_CUBE_H_
#define __MAKE_A_CUBE_H_

int makeACube(void * OutputImage);

#endif