/*
File:			ImageStats.h

Description:	Header file for ImageStats.txx. Function prints out some image statistics onto the command line that is useful as a sanity to see if the image has been modified properly.

Inputs:			image - (ImageType *) Pointer that holds the ITK::Image that is to be analyzed

Author:			Frankie (Hoi-Ki) Tong, Sep 6 2014

Changes:		Frankie (Hoi-Ki) Tong, Sep 6 2014
				Initial creation of the file.
*/


#include "itkMeanImageFilter.h"
#include "itkStatisticsImageFilter.h"

#include <iostream>

#pragma once
#ifndef __IMAGEISTATS_H__
#define __IMAGEISTATS_H__

template <typename ImageType>
void ImageStats(ImageType * image);

#include "ImageStats.txx"

#endif