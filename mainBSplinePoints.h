/*
File:			mainBSplinePoints.h

Description:	Header file for mainBSplinePoints.cxx.  Input parameters based on input processing from calling mainBSplinePoints as a command line function, ignoring the first 3 input parameters. Refer to mainBSplinePoints.cxx for more detailed input information.

Inputs:			argc - (int) Number of input arguments into the function
				argv - (char *) Character pointer containing the inputs of this function

Author:			Frankie (Hoi-Ki) Tong, Feb 24 2016

Changes:		Frankie (Hoi-Ki) Tong, Feb 24 2016
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
#include "itkBSplineScatteredDataPointSetToImageFilter.h"

#include "itkTransformToDisplacementFieldFilter.h"
#include "itkDisplacementFieldTransform.h"

#include "itknSift/itkScaleInvariantFeatureImageFilter.h"

#include "itkMesh.h"
#include "itkTransformMeshFilter.h"

#include <iostream>
#include <string>
#include <cmath>




#ifndef __MAINBSPLINEPOINTS__H__
#define __MAINBSPLINEPOINTS__H__

bool mainBSplinePoints (int argc, char *argv[]);

#endif