/*
File:			mainWarp.h

Description:	Header file for mainWarp.cxx.  Input parameters based on input processing from calling mainWarp as a command line function, ignoring the first 3 input parameters. Refer to mainWarp.cxx for more detailed input information.

Inputs:			argc - (int) Number of input arguments into the function
				argv - (char *) Character pointer containing the inputs of this function

Author:			Frankie (Hoi-Ki) Tong, Nov 9 2014

Changes:		Frankie (Hoi-Ki) Tong, Nov 9 2014
				Initial creation of the file.
*/

#include "itkImageDuplicator.h"

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
//#include "itkNearestNeighborExtrapolateImageFunction.h"

#include "itkTransformToDisplacementFieldFilter.h"

#include <iostream>
#include <string>

#ifndef __MAINWARP__H__
#define __MAINWARP__H__

bool mainWarp(int argc, char *argv[]);

#define TIMER

#ifdef TIMER
#include <ctime>
#endif

#endif