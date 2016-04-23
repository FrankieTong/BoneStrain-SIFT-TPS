/*
File:			mainTPSPerfectMatch.h

Description:	Header file for mainTPSPerfectMatch.cxx.  Input parameters based on input processing from calling mainTPSPerfectMatch as a command line function, ignoring the first 3 input parameters. Refer to mainTPSPerfectMatch.cxx for more detailed input information.

Inputs:			argc - (int) Number of input arguments into the function
				argv - (char *) Character pointer containing the inputs of this function

Author:			Frankie (Hoi-Ki) Tong, Apr 4 2015

Changes:		Frankie (Hoi-Ki) Tong, Apr 4 2015
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

#include "itkTransformToDisplacementFieldFilter.h"
#include "itkDisplacementFieldTransform.h"

#include "itknSift/itkScaleInvariantFeatureImageFilter.h"

#include <iostream>
#include <string>
#include <cmath>




#ifndef __MAINTPSPERFECTMATCH__H__
#define __MAINTPSPERFECTMATCH__H__

bool mainTPSPerfectMatch (int argc, char *argv[]);

#endif