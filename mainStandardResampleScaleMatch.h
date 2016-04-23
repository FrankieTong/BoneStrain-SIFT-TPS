/*
File:			mainStandardResampleScaleMatch.h

Description:	Header file for mainStandardResampleScaleMatch.cxx.  Input parameters based on input processing from calling mainStandardResampleScaleMatch as a command line function, ignoring the first 3 input parameters. Refer to mainStandardResampleScaleMatch.cxx for more detailed input information.

Inputs:			argc - (int) Number of input arguments into the function
				argv - (char *) Character pointer containing the inputs of this function

Author:			Frankie (Hoi-Ki) Tong, May 17 2015

Changes:		Frankie (Hoi-Ki) Tong, May 17 2015
				Initial creation of the file.
*/

#include "ImageStats.h"
#include "itkImageIORaw.h"
#include "ImageIOData.h"
#include "KeypointsIO.h"
#include "MatchedPoints.h"

#include "itkIdentityTransform.h"
#include "itkThinPlateSplineKernelTransform.h"

#include "itkAffineTransform.h"
#include "itkNearestNeighborExtrapolateImageFunction.h"

#include "itkResampleImageFilter.h"
#include "itkConstantBoundaryCondition.h"
#include "itkWindowedSincInterpolateImageFunction.h"

#include "itkBSplineInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

#include "itkIntensityWindowingImageFilter.h"

#include "itkTransformToDisplacementFieldFilter.h"

#include "itknSift/itkScaleInvariantFeatureImageFilter.h"

#include <iostream>
#include <string>
#include <cmath>




#ifndef __MAINSTANDARDRESAMPLESCALEMATCH__H__
#define __MAINSTANDARDRESAMPLESCALEMATCH__H__

bool mainStandardResampleScaleMatch (int argc, char *argv[]);


#define TIMER

#ifdef TIMER
#include <ctime>
#endif

#endif