/*
File:			ImageConfigFileImport.h

Description:	Header file for ImageConfigFileImport.cxx. Reading and writing image configuration parameters.

Author:			Frankie (Hoi-Ki) Tong, Sep 6 2014

Changes:		Frankie (Hoi-Ki) Tong, Sep 6 2014
				Initial creation of the file.
				
				Frankie (Hoi-Ki) Tong, Oct 4 2014
				Depreciated in exchange for ImageIOConfig structure due to awkard method to store configuration parameters
*/

#pragma once
#ifndef __IMAGECONFIGFILEIMPORT_H__
#define __IMAGECONFIGFILEIMPORT_H__

#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <limits>

//Pass out image io pointer
int ImageConfigFileImport ( std::string ConfigFileName, std::list<std::string> * ConfigFileParam, std::list<std::string> * ConfigFileContent);

//endif statement here
#endif