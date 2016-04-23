/*
File:			ImageIOConfig.cxx

Description:	Adds aditional reading and writing functionality to/from text files for class ImageIOData
				
				Performs some input data verification for data stored in ImageIOData classes when reading or writing to/from text files.

Author:			Frankie (Hoi-Ki) Tong, Oct 4 2014

Changes:		Frankie (Hoi-Ki) Tong, Oct 4 2014
				Initial creation of the file.
*/

#include "ImageIOConfig.h"

//Description:	Basic constructor with default values set to empty 
ImageIOConfig::ImageIOConfig(): 
	ConfigFileName(""){
}

//Description:	Basic constructor with default values set to parameters defined on construction by the user 
ImageIOConfig::ImageIOConfig(std::string inputConfigFileName, ImageIOData inputImageIOData): 
	ConfigFileName(inputConfigFileName), ImageMetaData(inputImageIOData) {
}

//Description: Basic deconstructor
ImageIOConfig::~ImageIOConfig(){
}

//Description: Copy constructor
ImageIOConfig::ImageIOConfig(const ImageIOConfig &obj) {
	ConfigFileName = obj.getConfigFileName();
	ImageMetaData = obj.getImageIOData();
}

//Description: Assignment constructor
ImageIOConfig & ImageIOConfig::operator=(const ImageIOConfig & obj) {
	ConfigFileName = obj.getConfigFileName();
	ImageMetaData = obj.getImageIOData();

	return *this;
}


//Description: 	Get command for the ConfigFileName protected data set	
//Return:		(std::string) - String with configuration file and path name
std::string ImageIOConfig::getConfigFileName() const { 
	return ConfigFileName; 
}

//Description: 	Set command for the ConfigFileName protected data set
//Input:		inputConfigFileName - (std::string) String with configuration file and path name to set data member as
//Return:		(bool) 0
bool ImageIOConfig::setConfigFileName(std::string inputConfigFileName){
	ConfigFileName = inputConfigFileName;
	return 0;
}


//Description: 	Get command for the ImageMetaData protected data set	
//Return:		(ImageIOData) - Copy of ImageMetaData 
ImageIOData ImageIOConfig::getImageIOData() const { 
	return  ImageMetaData;
}

//Description: 	Set command for the inputImageIOData protected data set
//Input:		inputImageIOData - (ImageIOData) ImageIOData structure containing image configuration data to set data member as
//Return:		(bool) 0
bool ImageIOConfig::setImageIOData(ImageIOData inputImageIOData){
	ImageMetaData = inputImageIOData;
	return 0;
}

//Description: 	Read in the text file defined in ConfigFileName, check the input parameters are correct, and store them in ConfigFileParam
//Return:		(bool) 0 if successful
//				(bool) 1 otherwise
bool ImageIOConfig::readConfigFile() {

	//Open config file and read in the parameters and their content
	std::list<std::string> * ConfigFileParam = new std::list<std::string>;
	std::list<std::string> * ConfigFileContent = new std::list<std::string>;
	if (this->ImageConfigFileImport(ConfigFileName, ConfigFileParam, ConfigFileContent)) {
		return 1;
	}

	//Check to make sure input parameters are correct
	if (this->ImageConfigFileCheck(ConfigFileParam,ConfigFileContent)) {
		return 1;
	}
	
	//Convert config parmeters to their proper values in ImageIOData
	if (this->setConfigFileParam(ConfigFileParam, ConfigFileContent)) {
		return 1;
	}
	
	return 0;

}

//Description: 	Write to the text file defined in ConfigFileName the configuration file using the parameters in ImageMetaData
//Return:		(bool) 0 if successful
//				(bool) 1 otherwise
bool ImageIOConfig::writeConfigFile() {

	//Check if config file exists, if it does delete it and get ready to make new config file
	std::ofstream outfile;
	outfile.open(ConfigFileName.c_str(),std::ofstream::trunc);

	if (!outfile.is_open()) {
		return 1;
	}

	//Start writing config file
	outfile << "//Contains image configuration parameters" <<std::endl;
	outfile << std::endl;
	outfile << "PixelType " << ImageMetaData.getPixelType() << std::endl;
	outfile << "PixelDimensionality " << ImageMetaData.getPixelDimensionality() << std::endl;
	outfile << "ImageDimensionality " << ImageMetaData.getImageDimensionality() << std::endl;

	outfile << "DimensionSize";
	int * DimensionSize = ImageMetaData.getDimensionSize();
	for (int i = 0; i < ImageMetaData.getImageDimensionality() ; i++) {
		outfile << " " <<DimensionSize[i];
	} outfile << std::endl;

	outfile << "SpacingSize";
	float * SpacingSize = ImageMetaData.getSpacingSize();
	for (int i = 0; i < ImageMetaData.getImageDimensionality() ; i++) {
		outfile << " " <<SpacingSize[i];
	} outfile << std::endl;

	outfile << "Origin";
	float * Origin = ImageMetaData.getOrigin();
	for (int i = 0; i < ImageMetaData.getImageDimensionality() ; i++) {
		outfile << " " <<Origin[i];
	} outfile << std::endl;

	outfile << "HeaderSize " << ImageMetaData.getHeaderSize() << std::endl;
	outfile << "ByteOrder " << ImageMetaData.getByteOrder() <<std::endl;

	outfile.close();
	return 0;
}



//Subfunctions that will be used only in this class (they are ugly which is why I am hiding them...)

//Description: 	Read in a text file where each new line in the text file consists of "(parameter name)\s(parameter contents)\n". C/C++ style content lines are ignored and each new line is a new parameter.
//Input:		ConfigFileName - (std::string) Path and file name of the text file to be read
//				ConfigFileParam - (std::list<std::string> *) String list pointer containing the first word of each new line in the text file which is the parameter name
//				ConfigFileContent - (std::list<std::string> *) String list pointer containing each new line in the text file excluding the first word which then together consist of the parameter contents
//Return:		(bool) 0 if successful
//				(bool) 1 otherwise
bool ImageIOConfig::ImageConfigFileImport (std::string ConfigFileName, std::list<std::string> * ConfigFileParam, std::list<std::string> * ConfigFileContent) {

	std::ifstream infile;
	infile.open(ConfigFileName.c_str());

	//Check if file is open, if not return 1
	if (!infile.is_open()) {
		return 1;
	}

	std::string line;

	//Split the file into param name param content lists
	while (!infile.eof())
	{
		std::getline(infile, line);
		std::istringstream iss(line);
		std::string ConfigFileParamName;
		std::string ConfigFileParamContent = "";

		//Grab first parameter
		iss >> ConfigFileParamName;

		//Check if line is a comment (denoted by "//" in front)  or blank line and ignore line
		std::size_t comment = ConfigFileParamName.find("//");
		if (comment==std::string::npos && ConfigFileParamName.length() > 0) {

			//Store parameter name and its contents in seperate lists
			ConfigFileParam->push_back(ConfigFileParamName);

			while ( iss >> ConfigFileParamName ) {
				ConfigFileParamContent += (ConfigFileParamName + " ");
			}
			
			//Remove last space
			ConfigFileParamContent = ConfigFileParamContent.substr(0,ConfigFileParamContent.size()-1);
			ConfigFileContent->push_back(ConfigFileParamContent);
		}

		//Flush istringstream
		iss.str("");
		iss.clear();
	}

	return 0;
}

//Description: 	Hash parameter to be used to determine if the input data type used in PixelType is correct
//Input:		inString - (std::string) input string to be used for comparison with all allowable pixel types
//Return:		returnType::(valid PixelType) if inString matches one of the defined PixelTypes
//				returnType::eUnknownParam if inString does not match any defined PixelTypes	
returnType::e ImageIOConfig::hashitParam (std::string inString) {
	if (inString == "PixelType") {return returnType::ePixelType;}
	else if (inString == "ImageDimensionality") {return returnType::eImageDimensionality;}
	else if (inString == "DimensionSize") {return returnType::eDimensionSize;}
	else if (inString == "SpacingSize") {return returnType::eSpacingSize;}
	else if (inString == "HeaderSize") {return returnType::eHeaderSize;}
	else if (inString == "ByteOrder") {return returnType::eByteOrder;}
	else if (inString == "Origin") {return returnType::eOrigin;}
	else if (inString == "PixelDimensionality") {return returnType::ePixelDimensionality;}
	else {return returnType::eUnknownParam;}
}

//Description: 	Verifies the values inside ImageIOData is correct and complete
//Input:		ConfigFileParam - (std::list<std:;string> *) string list pointer storing the configuration parameters for ITK images
//				ConfigFileContent - (std::list<std:;string> *) string list pointer storing the configuration parameters content/values for ITK images
//Return:		(bool) 0 if all data in ImageIOData is correct and complete
//				(bool) 1 otherwise
bool ImageIOConfig::ImageConfigFileCheck (std::list<std::string> * ConfigFileParam, std::list<std::string> * ConfigFileContent) {

	bool finishedChecking = false;

	bool checkPixelType = false;
	bool checkImageDimensionality = false;
	int ImageDimensionality = 0;
	bool checkDimensionSize = false;
	bool checkSpacingSize = false;
	bool checkHeaderSize = false;
	bool checkByteOrder = false;
	bool checkOrigin = false;
	bool checkPixelDimensionality = false;

	while (!finishedChecking) {
		
		finishedChecking = true;

		//Go through list of parameter names and pick out the ones you want
		std::list<std::string>::iterator itParam;
		std::list<std::string>::iterator itContent;
		std::istringstream testParam;
		int testInt = 0;
		float testFloat = 0;

		for (itParam = ConfigFileParam->begin(), itContent = ConfigFileContent->begin(); itParam != ConfigFileParam->end() && itContent != ConfigFileContent->end(); ++itParam, ++itContent) {
			switch(this->hashitParam(*itParam)) {
				case returnType::ePixelType:
					//check if already checked
					if (checkPixelType) {
						break;
					}
					//check for valid type
					if (!itContent->compare("char") || !itContent->compare("unsigned char") || 
						!itContent->compare("int") || !itContent->compare("unsigned int") ||
						!itContent->compare("short") || !itContent->compare("unsigned short") ||
						!itContent->compare("float") || !itContent->compare("double")) {
							checkPixelType = true;
							finishedChecking = false;
					} else {
						//Invalid data type
						return 1;
					}
					break;

				case returnType::eImageDimensionality:
					//check if already checked
					if (checkImageDimensionality) {
						break;
					}
					//check for valid dimension number
					testParam.str(*itContent);
					testParam >> ImageDimensionality;
					if (!testParam.fail() && ImageDimensionality > 0) {
						checkImageDimensionality = true;
						finishedChecking = false;
					} else {
						//Invalid ImageDimensionality
						return 1;
					}
					testParam.str("");
					testParam.clear();
					break;

				case returnType::eDimensionSize:
					//check if already checked
					if (checkDimensionSize) {
						break;
					}
					//Check if dimensionailty has been set. If not, check later
					if (!checkImageDimensionality) {
						finishedChecking = false;
					} else {
						//ImageDimensionality is correct, find if dimension size is correct
						int numDimensions = 0;

						size_t start = 0;
						size_t pos = 0;
						std::string token;

						//Split into dimension sizes
						while (start <= itContent->length()) {

							pos = itContent->find(" ",start);

							if (itContent->find(" ",start) == std::string::npos) {
								token = itContent->substr(start, itContent->length() - start);
								pos = itContent->length();
							} else {
								token = itContent->substr(start, pos - start);
							}

							//Does the dimension size make sense?
							testParam.str(token);
							testParam >> testInt;

							if (testParam.fail() || testParam.peek() != EOF || testInt <= 0) {
								//Invalid dimension size
								return 1;
							}
							//Count number of dimension sizes
							numDimensions++;

							start = pos + 1;

							testParam.str("");
							testParam.clear();
						}
						
						//Does number of dimension sizes match ImageDimensionality?
						if (numDimensions == ImageDimensionality) {
							checkDimensionSize = true;
							finishedChecking = false;
						} else {
							return 1;
						}
					}
					testParam.str("");
					testParam.clear();
					break;

				case returnType::eSpacingSize:
					//check if already checked
					if (checkSpacingSize) {
						break;
					}
					//Check if dimensionailty has been set. If not, check later
					if (!checkImageDimensionality) {
						finishedChecking = false;
					} else {
						//ImageDimensionality is correct, find if spacing size is correct
						int numSpacings = 0;

						size_t start = 0;
						size_t pos = 0;
						std::string token;

						//Split into spacing sizes
						while (start <= itContent->length()) {

							pos = itContent->find(" ",start);

							if (itContent->find(" ",start) == std::string::npos) {
								token = itContent->substr(start, itContent->length() - start);
								pos = itContent->length();
							} else {
								token = itContent->substr(start, pos - start);
							}

							//Does the spacing size make sense?
							testParam.str(token);
							testParam >> testFloat;

							if (testParam.fail() || testParam.peek() != EOF || testFloat <= 0) {
								//Invalid dimension size
								return 1;
							}

							//Count number of spacing sizes
							numSpacings++;

							start = pos + 1;
							testParam.str("");
							testParam.clear();
						}
						
						//Does number of spacing sizes match ImageDimensionality?
						if (numSpacings == ImageDimensionality) {
							checkSpacingSize = true;
							finishedChecking = false;
						} else {
							return 1;
						}
					}
					testParam.str("");
					testParam.clear();
					break;

				case returnType::eHeaderSize:
					//check if already checked
					if (checkHeaderSize) {
						break;
					}
					//check for valid header size
					testParam.str(*itContent);
					testParam >> testInt;
					if (!testParam.fail() && testParam.peek() == EOF && testInt >= 0) {
						checkHeaderSize = true;
						finishedChecking = false;
					} else {
						//Invalid header size
						return 1;
					}
					testParam.str("");
					testParam.clear();
					break;

				case returnType::eByteOrder:
					//check if already checked
					if (checkByteOrder) {
						break;
					}
					//check for byte order
					if (!itContent->compare("little") || !itContent->compare("big")) {
							checkByteOrder = true;
							finishedChecking = false;
					} else {
						//Invalid byte order
						return 1;
					}
					break;

				case returnType::eOrigin:
					if (checkOrigin) {
						break;
					}
					//Check if dimensionailty has been set. If not, check later
					if (!checkImageDimensionality) {
						finishedChecking = false;
					} else {
						//ImageDimensionality is correct, find if origin is correct
						int numOrigin = 0;

						size_t start = 0;
						size_t pos = 0;
						std::string token;

						//Split into origin numbers
						while (start <= itContent->length()) {

							pos = itContent->find(" ",start);

							if (itContent->find(" ",start) == std::string::npos) {
								token = itContent->substr(start, itContent->length() - start);
								pos = itContent->length();
							} else {
								token = itContent->substr(start, pos - start);
							}

							//Does the origin numbers make sense?
							testParam.str(token);
							testParam >> testFloat;

							if (testParam.fail() || testParam.peek() != EOF) {
								//Invalid origin numbers
								return 1;
							}
							//Count number of origin numbers
							numOrigin++;

							start = pos + 1;

							testParam.str("");
							testParam.clear();
						}
						
						//Does number of origin numbers match ImageDimensionality?
						if (numOrigin == ImageDimensionality) {
							checkOrigin = true;
							finishedChecking = false;
						} else {
							return 1;
						}
					}
					testParam.str("");
					testParam.clear();
					break;

				case returnType::ePixelDimensionality:
					//check if already checked
					if (checkPixelDimensionality) {
						break;
					}
					//check for valid dimension number
					testParam.str(*itContent);
					testParam >> testInt;
					if (!testParam.fail() && testInt > 0) {
						checkPixelDimensionality = true;
						finishedChecking = false;
					} else {
						//Invalid dimensionily
						return 1;
					}
					testParam.str("");
					testParam.clear();
					break;

				default:
					break;
			}
		}
	}

	//Finished checking. Find if we meet all the parameters
	if (checkPixelType == true && checkImageDimensionality == true && 
		checkDimensionSize == true && checkSpacingSize == true &&  
		checkHeaderSize == true && checkByteOrder == true &&  
		checkOrigin == true) {
			return 0;
	} else {
		return 1;
	}
}

//Description: 	Places the parameter contents in ConfigFileContent to thier respective lcoations in ImageMetaData
//Input:		ConfigFileParam - (std::list<std:;string> *) string list pointer storing the configuration parameters for ITK images
//				ConfigFileContent - (std::list<std:;string> *) string list pointer storing the configuration parameters content/values for ITK images
//Return:		(bool) 0 if all data in ImageIOData is correct and complete
//				(bool) 1 otherwise
bool ImageIOConfig::setConfigFileParam(std::list<std::string> * ConfigFileParam, std::list<std::string> * ConfigFileContent) {

	std::string  ConfigContentString;
	std::istringstream stringConvert;

	//Find and set pixel type
	ConfigContentString = this->ReturnConfigIndex (this->FindParamIndex("PixelType",ConfigFileParam),ConfigFileContent);
	if (ImageMetaData.setPixelType(ConfigContentString)) {return 1;};

	//Find and set pixel dimensionsionality
	ConfigContentString = this->ReturnConfigIndex (this->FindParamIndex("PixelDimensionality",ConfigFileParam),ConfigFileContent);
	stringConvert.str(ConfigContentString);
	unsigned int PixelDimensionality;
	stringConvert >> PixelDimensionality;
	if (ImageMetaData.setPixelDimensionality(PixelDimensionality)) {return 1;};
	stringConvert.str("");
	stringConvert.clear();
	
	//Find and set image dimensionsionality
	ConfigContentString = this->ReturnConfigIndex (this->FindParamIndex("ImageDimensionality",ConfigFileParam),ConfigFileContent);
	stringConvert.str(ConfigContentString);
	int ImageDimensionality;
	stringConvert >> ImageDimensionality;
	if (ImageMetaData.setImageDimensionality(ImageDimensionality)) {return 1;};
	stringConvert.str("");
	stringConvert.clear();

	//Find and set dimension size
	ConfigContentString = this->ReturnConfigIndex (this->FindParamIndex("DimensionSize",ConfigFileParam),ConfigFileContent);
	stringConvert.str(ConfigContentString);
	int* DimensionSize = new int[ImageDimensionality];
	for (int i = 0; i < ImageDimensionality; i++ ) {
		stringConvert >> DimensionSize[i];
	}
	if (ImageMetaData.setDimensionSize(DimensionSize)) {return 1;};
	stringConvert.str("");
	stringConvert.clear();

	//Find and set spacing size
	ConfigContentString = this->ReturnConfigIndex (this->FindParamIndex("SpacingSize",ConfigFileParam),ConfigFileContent);
	stringConvert.str(ConfigContentString);
	float* SpacingSize = new float[ImageDimensionality];
	for (int i = 0; i < ImageDimensionality; i++ ) {
		stringConvert >> SpacingSize[i];
	}
	if (ImageMetaData.setSpacingSize(SpacingSize)) {return 1;};
	stringConvert.str("");
	stringConvert.clear();

	//Find and set origin
	ConfigContentString = this->ReturnConfigIndex (this->FindParamIndex("Origin",ConfigFileParam),ConfigFileContent);
	stringConvert.str(ConfigContentString);
	float* Origin = new float[ImageDimensionality];
	for (int i = 0; i < ImageDimensionality; i++ ) {
		stringConvert >> Origin[i];
	}
	if (ImageMetaData.setOrigin(Origin)) {return 1;};
	stringConvert.str("");
	stringConvert.clear();

	//Find and set header size
	ConfigContentString = this->ReturnConfigIndex (this->FindParamIndex("HeaderSize",ConfigFileParam),ConfigFileContent);
	stringConvert.str(ConfigContentString);
	int HeaderSize;
	stringConvert >> HeaderSize;
	if (ImageMetaData.setHeaderSize(HeaderSize)) {return 1;};
	stringConvert.str("");
	stringConvert.clear();

	//Find and set byte order
	ConfigContentString = this->ReturnConfigIndex (this->FindParamIndex("ByteOrder",ConfigFileParam),ConfigFileContent);
	if (ImageMetaData.setByteOrder(ConfigContentString)) {return 1;};
	
	return 0;
}

//Description: 	Finds and returns the index value of the parameter name found in the parameter list
//Input:		Param - (std::string) String with the parameter name we are searching for
//				ConfigFileParam - (std::list<std:;string> *) string list pointer storing the configuration parameters for ITK images
//Return:		(int) index of the string list where the parameter name was found
int ImageIOConfig::FindParamIndex (std::string Param, std::list<std::string> * ConfigFileParam ) {
	std::list<std::string>::iterator itParam = std::find(ConfigFileParam->begin(),ConfigFileParam->end(), Param);
	int index = (int) std::distance(ConfigFileParam->begin(),itParam);
	return index;
}

//Description: 	When given an index calue, finds and returns the contents of the string list it is given
//Input:		Index - (int) Index value that you want to return the contents of ConfigFileContent from
//				ConfigFileContent - (std::list<std:;string> *) string list pointer storing the configuration parameters content/values for ITK images
//Return:		(std::sting) Content of the string list at the given index valuw
std::string ImageIOConfig::ReturnConfigIndex (int Index, std::list<std::string> * ConfigFileContent) {
	std::list<std::string>::iterator itContent = ConfigFileContent->begin();
	std::advance(itContent,Index);
	return *itContent;
}
