/*
File:			ImageConfigFileCheck.h

Description:	Header file for ImageConfigFileCheck.cxx. Runs a checking protocol to see if all the image configuration parameters are properly set.
				
Author:			Frankie (Hoi-Ki) Tong, Sep 6 2014

Changes:		Frankie (Hoi-Ki) Tong, Sep 6 2014
				Initial creation of the file.
				
				Frankie (Hoi-Ki) Tong, Oct 4 2014
				Depreciated in exchange for ImageIOConfig structure due to awkard method to store configuration parameters
*/

#include <string>
#include <list>
#include <sstream>
#include <iostream>
#include <new>
#include <limits>

#pragma once
#ifndef __IMAGECONFIGFILECHECK_H__
#define __IMAGECONFIGFILECHECK_H__

//Namespace defined for identifying if PixelType is a proper data type
enum EnumConfigFile {ePixelType, eImageDimensionality, eDimensionSize, eSpacingSize,
					eHeaderSize, eByteOrder, eOrigin, ePixelDimensionality, eUnknownParam};
EnumConfigFile hashitParam (std::string inString);

//Check to see if the image parameters needed to import an image are all present in ConfigFileParam and ConfigFileContent
bool ImageConfigFileCheck (std::list<std::string> * ConfigFileParam, std::list<std::string> * ConfigFileContent);
#endif