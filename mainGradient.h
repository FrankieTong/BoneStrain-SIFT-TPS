/*
File:			mainGradient.h

Description:	Header file for mainGradient.cxx.  Input parameters based on input processing from calling mainGradient as a command line function, ignoring the first 3 input parameters. Refer to mainGradient.cxx for more detailed input information.

Inputs:			argc - (int) Number of input arguments into the function
				argv - (char *) Character pointer containing the inputs of this function

Author:			Frankie (Hoi-Ki) Tong, Apr 11 2015

Changes:		Frankie (Hoi-Ki) Tong, Apr 11 2015
				Initial creation of the file.
*/

#include "ImageStats.h"
#include "itkImageIORaw.h"
#include "ImageIOData.h"


#include "itkBSplineInterpolateImageFunction.h"

#include "itkAffineTransform.h"
#include "itkIdentityTransform.h"
#include "itkResampleImageFilter.h"

#include "itkVectorGradientMagnitudeImageFilter.h"
#include "itkGradientImageFilter.h"
#include "itkVectorIndexSelectionCastImageFilter.h"
#include "itkCropImageFilter.h"

#include "itkConstantBoundaryCondition.h"
#include "itkWindowedSincInterpolateImageFunction.h"

#include "itkImageRegionIterator.h"

#include <iostream>
#include <string>
#include <cmath>

#ifndef __MAINGRADIENT__H__
#define __MAINGRADIENT__H__

bool mainGradient (int argc, char *argv[]);

#endif