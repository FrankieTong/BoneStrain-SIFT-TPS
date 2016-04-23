/*
File:			ImageImport.h

Description:	Header file for ImageImport.cxx. Import and export of ITK image data in a raw data format

Author:			Frankie (Hoi-Ki) Tong, Sep 6 2014

Changes:		Frankie (Hoi-Ki) Tong, Sep 6 2014
				Initial creation of the file.
				
				Frankei (Hoi-Ki) Tong, Oct 4 2014
				Depreciated in exchange for ImageIOConfig structure due to awkard method to store configuration parameters
*/

#include "itkRawImageIO.h"
#include "itkNumericSeriesFileNames.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"

#include <string>
#include <list>
#include <sstream>
#include <iostream>
#include <new>
#include <limits>
#include <iterator>

#pragma once
#ifndef __IMAGEIMPORT_H__
#define __IMAGEIMPORT_H__

//Abandonded function prototype. Using ImageImportReadImage as a temporary fix...
template <class OutputImageType, class OutputPixelType> 
void ImageImportScalor(std::list<std::string> * ConfigFileParam, 
							std::list<std::string> * ConfigFileContent, 
							std::string InputFile,
							void * OutputImage);

//Namespace defined for identifying if PixelType is a proper data type
enum EnumImageType {eChar, eUnsignedChar, eInt, eUnsignedInt,
					eShort, eUnsignedShort, eFloat, eDouble, eUnknownType};
EnumImageType hashitType (std::string inString);

//Read in ITK image from raw data file using configuration parameters
template <class PixelType, unsigned int ImageDimensionality> 
void ImageImportReadImage (std::list<std::string> * ConfigFileParam,  std::list<std::string> * ConfigFileContent, std::string InputFile, void * Image) ;

//Write out ITK iamge as raw data file using configuration parameters
template <class PixelType, unsigned int ImageDimensionality> 
void ImageImportWriteImage (std::list<std::string> * ConfigFileParam,  std::list<std::string> * ConfigFileContent, std::string OutputFile, void * Image) ;

//Support function to assist in identifing cofiguration parameter contents
int FindParamIndex (std::string Param, std::list<std::string> * ConfigFileParam );
std::string ReturnConfigIndex (int Index, std::list<std::string> * ConfigFileContent);

#include "ImageImport.txx"

#endif