/*
File:			mainTransform.h

Description:	Header file for mainTransform.cxx.  Input parameters based on input processing from calling mainTransform as a command line function, ignoring the first 3 input parameters. Refer to mainTransform.cxx for more detailed input information.

Inputs:			argc - (int) Number of input arguments into the function
				argv - (char *) Character pointer containing the inputs of this function

Author:			Frankie (Hoi-Ki) Tong, Jan 12 2015

Changes:		Frankie (Hoi-Ki) Tong, Jan 12 2015
				Initial creation of the file.
*/

#include "ImageStats.h"
#include "itkImageIORaw.h"
#include "ImageIOData.h"

#include "itkVectorLinearInterpolateImageFunction.h"
#include "itkVectorResampleImageFilter.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkWarpImageFilter.h"

//#include "itkVectorLinearInterpolateImageFunction.h"
//#include "itkDisplacementFieldTransform.h"

#include "itkResampleImageFilter.h"
#include "itkConstantBoundaryCondition.h"
#include "itkWindowedSincInterpolateImageFunction.h"

#include <iostream>
#include <string>

#ifndef __MAINTRANSFORM__H__
#define __MAINTRANSFORM__H__

bool mainTransform (int argc, char *argv[]);

#endif