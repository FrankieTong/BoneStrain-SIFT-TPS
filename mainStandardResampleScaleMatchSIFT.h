/*
File:			mainStandardResampleScaleMatchSIFT.h

Description:	Header file for mainStandardResampleScaleMatchSIFT.cxx.  Input parameters based on input processing from calling mainStandardResampleScaleMatchSIFT as a command line function, ignoring the first 3 input parameters. Refer to mainStandardResampleScaleMatchSIFT.cxx for more detailed input information.
				NOTE: This function operates exactly the same as mainStandardResampleScaleMatch() except it does not attempt to calculate and output a displacement field.

Inputs:			argc - (int) Number of input arguments into the function
				argv - (char *) Character pointer containing the inputs of this function

Author:			Frankie (Hoi-Ki) Tong, Aug 20 2015

Changes:		Frankie (Hoi-Ki) Tong, Aug 20 2015
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




#ifndef __MAINSTANDARDRESAMPLESCALEMATCHSIFT__H__
#define __MAINSTANDARDRESAMPLESCALEMATCHSIFT__H__

bool mainStandardResampleScaleMatchSIFT (int argc, char *argv[]);


#define TIMER

#ifdef TIMER
#include <ctime>
#endif

#endif