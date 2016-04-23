/*
File:			mainResample.h

Description:	Header file for mainResample.cxx.  Input parameters based on input processing from calling mainResample as a command line function, ignoring the first 3 input parameters. Refer to mainResample.cxx for more detailed input information.

Inputs:			argc - (int) Number of input arguments into the function
				argv - (char *) Character pointer containing the inputs of this function

Author:			Frankie (Hoi-Ki) Tong, Feb 16 2015

Changes:		Frankie (Hoi-Ki) Tong, Feb 16 2015
				Initial creation of the file.
*/

#include "ImageStats.h"
#include "itkImageIORaw.h"
#include "ImageIOData.h"


#include "itkBSplineInterpolateImageFunction.h"

#include "itkAffineTransform.h"
#include "itkIdentityTransform.h"
#include "itkResampleImageFilter.h"

#include "itkNearestNeighborExtrapolateImageFunction.h"

#include "itkConstantBoundaryCondition.h"
#include "itkWindowedSincInterpolateImageFunction.h"

#include <iostream>
#include <string>
#include <sstream>
#include <cmath>

#ifndef __MAINRESAMPLE__H__
#define __MAINRESAMPLE__H__

bool mainResample (int argc, char *argv[]);

#endif