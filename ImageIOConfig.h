/*
File:			ImageIOConfig.h

Description:	Header file for ImageIOConfig.cxx. Adds aditional reading and writing functionality to/from text files for class ImageIOData
				
				Performs some input data verification for data stored in ImageIOData classes when reading or writing to/from text files.
				
				The public functions defined here are not designed well generally and should be rewritten if the functions need to be modified.

Author:			Frankie (Hoi-Ki) Tong, Oct 4 2014

Changes:		Frankie (Hoi-Ki) Tong, Oct 4 2014
				Initial creation of the file.
*/

#include "ImageIOData.h"
#include <fstream>
#include <sstream>
#include <string>
#include <list>
#include <limits>
#include <iterator>
#include <algorithm> 

#ifndef __ImageIOConfig_H__
#define __ImageIOConfig_H__

//Namespace defined for identifying if PixelType is a proper data type
namespace returnType {
	enum e {ePixelType, eImageDimensionality, eDimensionSize, eSpacingSize,
							eHeaderSize, eByteOrder, eOrigin, ePixelDimensionality, eUnknownParam};
};

class ImageIOConfig{
public:

	//Standard constructors and descructors
	ImageIOConfig::ImageIOConfig();															// constructor; initialize everything to empty values
	ImageIOConfig::ImageIOConfig(std::string ConfigFileName, ImageIOData inputImageIOData);	// constructor; initialize with preset values

	ImageIOConfig::~ImageIOConfig(); //Deconstructor

	ImageIOConfig::ImageIOConfig(const ImageIOConfig &obj);	//Copy constructor
	ImageIOConfig & operator=(const ImageIOConfig & obj);	//Assignment constructor

	
	//Standard get and set commands for data values in the class
	std::string getConfigFileName() const;
	bool setConfigFileName(std::string inputConfigFileName);
  
	ImageIOData getImageIOData() const;
	bool setImageIOData(ImageIOData inputImageIOData);
  
  
	//Primary functions called for reading and writing configurations files to/from an ImageIOData class to/from a text file
	bool readConfigFile();
	bool writeConfigFile();
	
protected:

	//Protected Data Classes
	std::string ConfigFileName;	//Path and name of text file where image configuration parameters are stored
	ImageIOData ImageMetaData;	//ImageMetaData class where image configuration parameters are stored for reading and writing to/from text files

private:

	//Subfinctions that will only be used in this class (they are ugly so I am hiding them...)
	bool ImageConfigFileImport ( std::string ConfigFileName, std::list<std::string> * ConfigFileParam, std::list<std::string> * ConfigFileContent);
	bool ImageConfigFileCheck (std::list<std::string> * ConfigFileParam, std::list<std::string> * ConfigFileContent);
	bool setConfigFileParam(std::list<std::string> * ConfigFileParam, std::list<std::string> * ConfigFileContent);
	returnType::e hashitParam (std::string inString);
	int FindParamIndex (std::string Param, std::list<std::string> * ConfigFileParam );
	std::string ReturnConfigIndex (int Index, std::list<std::string> * ConfigFileContent);
  
};

#endif