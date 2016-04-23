/*
File:			mainTPS.h

Description:	Header file for mainTPS.cxx.  Input parameters based on input processing from calling mainTPS as a command line function, ignoring the first 3 input parameters. Refer to mainTPS.cxx for more detailed input information.

Inputs:			argc - (int) Number of input arguments into the function
				argv - (char *) Character pointer containing the inputs of this function

Author:			Frankie (Hoi-Ki) Tong, Aug 13 2015

Changes:		Frankie (Hoi-Ki) Tong, Aug 13 2015
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

#include "itkMesh.h"
#include "itkTransformMeshFilter.h"

#include <iostream>
#include <string>
#include <cmath>




#ifndef __MAINTPS__H__
#define __MAINTPS__H__

bool mainTPS (int argc, char *argv[]);

#endif