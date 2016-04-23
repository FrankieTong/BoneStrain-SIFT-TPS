/*
File:			mainRotate.h

Description:	Header file for mainRotate.cxx.  Input parameters based on input processing from calling mainRotate as a command line function, ignoring the first 3 input parameters. Refer to mainRotate.cxx for more detailed input information.

Inputs:			argc - (int) Number of input arguments into the function
				argv - (char *) Character pointer containing the inputs of this function

Author:			Frankie (Hoi-Ki) Tong, Nov 1 2014

Changes:		Frankie (Hoi-Ki) Tong, Nov 1 2014
				Initial creation of the file.
*/

#include "itkRawImageIO.h"
#include "itkNumericSeriesFileNames.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSimpleFilterWatcher.h"
#include "itkShrinkImageFilter.h"
#include "itkComposeImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkMeanImageFilter.h"
#include "itkStatisticsImageFilter.h"

#include "itkImageDuplicator.h"

#include "itkAffineTransform.h"
#include "itkIdentityTransform.h"
#include "itkResampleImageFilter.h"

#include "itkNearestNeighborExtrapolateImageFunction.h"

#include "itkConstantBoundaryCondition.h"
#include "itkWindowedSincInterpolateImageFunction.h"

#include "itknSift/itkScaleInvariantFeatureImageFilter.h"

#include "ImageStats.h"
#include "itkImageIORaw.h"
#include "ImageIOData.h"
#include "KeypointsIO.h"
#include "MatchedPoints.h"

#include "mainRotateSupportClass.h"

#include "makeACube.h"

#include <iostream>
#include <string>

#ifndef __MAINROTATE__H__
#define __MAINROTATE__H__

bool mainRotate(int argc, char *argv[]);

//#define RESAMPLE

#ifdef RESAMPLE
#include "itkLinearInterpolateImageFunction.h"
#define RESAMPLE_FACTOR 2
#endif

#define TIMER

#ifdef TIMER
#include <ctime>
#endif

#endif