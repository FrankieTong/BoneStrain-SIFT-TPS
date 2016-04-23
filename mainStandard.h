/*
File:			mainStandard.h

Description:	Header file for mainStandard.cxx.  Input parameters based on input processing from calling mainStandard as a command line function, ignoring the first 3 input parameters. Refer to mainStandard.cxx for more detailed input information.

Inputs:			argc - (int) Number of input arguments into the function
				argv - (char *) Character pointer containing the inputs of this function

Author:			Frankie (Hoi-Ki) Tong, Nov 9 2013

Changes:		Frankie (Hoi-Ki) Tong, Nov 9 2013
				Initial creation of the file.
*/

#include "ImageStats.h"
#include "itkImageIORaw.h"
#include "ImageIOData.h"
#include "KeypointsIO.h"
#include "MatchedPoints.h"

#include "itkIdentityTransform.h"
#include "itkThinPlateSplineKernelTransform.h"

#include "itkResampleImageFilter.h"
#include "itkConstantBoundaryCondition.h"
#include "itkWindowedSincInterpolateImageFunction.h"

#include "itkTransformToDisplacementFieldFilter.h"

#include "itknSift/itkScaleInvariantFeatureImageFilter.h"

#include <iostream>
#include <string>

#ifndef __MAINSTANDARD__H__
#define __MAINSTANDARD__H__

bool mainStandard (int argc, char *argv[]);


#define TIMER

#ifdef TIMER
#include <ctime>
#endif

#endif