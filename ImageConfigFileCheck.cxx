/*
File:			ImageConfigFileCheck.h

Description:	Runs a checking protocol to see if all the image configuration parameters are properly set.
				
Author:			Frankie (Hoi-Ki) Tong, Sep 6 2014

Changes:		Frankie (Hoi-Ki) Tong, Sep 6 2014
				Initial creation of the file.
				
				Frankie (Hoi-Ki) Tong, Oct 4 2014
				Depreciated in exchange for ImageIOConfig structure due to awkard method to store configuration parameters
*/


#include "imageConfigFileCheck.h"

//Description: 	Hash parameter to be used to determine if the input data type used in PixelType is correct
//Input:		inString - (std::string) input string to be used for comparison with all allowable pixel types
//Return:		returnType::(valid PixelType) if inString matches one of the defined PixelTypes
//				returnType::eUnknownParam if inString does not match any defined PixelTypes	
EnumConfigFile hashitParam (std::string inString) {
	if (inString == "PixelType") {return ePixelType;}
	else if (inString == "ImageDimensionality") {return eImageDimensionality;}
	else if (inString == "DimensionSize") {return eDimensionSize;}
	else if (inString == "SpacingSize") {return eSpacingSize;}
	else if (inString == "HeaderSize") {return eHeaderSize;}
	else if (inString == "ByteOrder") {return eByteOrder;}
	else if (inString == "Origin") {return eOrigin;}
	else if (inString == "PixelDimensionality") {return ePixelDimensionality;}
	else {return eUnknownParam;}
}

//Description: 	Verifies the values inside ImageIOData is correct and complete
//Input:		ConfigFileParam - (std::list<std:;string> *) string list pointer storing the configuration parameters for ITK images
//				ConfigFileContent - (std::list<std:;string> *) string list pointer storing the configuration parameters content/values for ITK images
//Return:		(bool) 0 if all data in ImageIOData is correct and complete
//				(bool) 1 otherwise
bool ImageConfigFileCheck (std::list<std::string> * ConfigFileParam, std::list<std::string> * ConfigFileContent) {

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
			switch(hashitParam(*itParam)) {
				case ePixelType:
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
						return false;
					}
					break;

				case eImageDimensionality:
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
						return false;
					}
					testParam.str("");
					testParam.clear();
					break;

				case eDimensionSize:
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
								return false;
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
							return false;
						}
					}
					testParam.str("");
					testParam.clear();
					break;

				case eSpacingSize:
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
								return false;
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
							return false;
						}
					}
					testParam.str("");
					testParam.clear();
					break;

				case eHeaderSize:
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
						return false;
					}
					testParam.str("");
					testParam.clear();
					break;

				case eByteOrder:
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
						return false;
					}
					break;

				case eOrigin:
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
								return false;
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
							return false;
						}
					}
					testParam.str("");
					testParam.clear();
					break;

				case ePixelDimensionality:
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
						return false;
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
			return true;
	} else {
		return false;
	}
}

/*
	EnumConfigFile ParamIndex;
	
	//Do we want to pass this back directly back up to the previous function to make the ITK image containers?

	//Properly format input lines into proper containers
	std::string rPixelType;
	std::string rDimensionality;
	std::string rDimensionSize;
	std::string rSpacingSize;
	std::string rHeaderSize;
	std::string rByteOrder;

	while (std::getline(infile, line))
	{
		std::istringstream iss(line);
		std::string rParamName;

		if (!(iss >> rParamName)) { return 0;}

		switch(hashit(rParamName)) {
			case ePixelType:
				rPixelType = iss.str();
				break;
			case eDimensionality:
				rDimensionality = iss.str();
				break;
			case eDimensionSize:
				rDimensionSize = iss.str();
				break;
			case eSpacingSize:
				rSpacingSize = iss.str();
				break;
			case eHeaderSize:
				rHeaderSize = iss.str();
				break;
			case eByteOrder:
				rByteOrder = iss.str();
				break;
			default:
				break;
		}

		//Will typedef disappear out of while loop?

		//int a, b;
		//if (!(iss >> a >> b)) { break; } // error

		// process pair (a,b)
	}
	
	std::ifstream infile(ConfigFileName);

	std::string line;

	EnumConfigFile ParamIndex;

	std::string rPixelType;
	int Dimensionality;

	while (std::getline(infile, line))
	{
		std::istringstream iss(line);
		std::string rParamName;

		if (!(iss >> rParamName)) { return 0;}

		switch(hashit(rParamName)) {
			case ePixelType:
				if (!(iss >> rPixelType)) { return 0; }
				break;
			case eDimensionality:
				if (!(iss >> Dimensionality)) { return 0; }
				break;
			default:
				break;
		}

		//Will typedef disappear out of while loop?

		//int a, b;
		//if (!(iss >> a >> b)) { break; } // error

		// process pair (a,b)
	}

	//Try to make itk image structure

	if (rPixelType == "unsigned int") { typedef unsigned int PixelType; }
	if (rPixelType == "int") { typedef int PixelType; }
	if (rPixelType == "char") { typedef char PixelType; }
	if (rPixelType == "unsigned char") { typedef unsigned char PixelType; }
	if (rPixelType == "float") { typedef float PixelType; }
	if (rPixelType == "double") { typedef double PixelType; }

	try
	{
		typedef itk::Image<PixelType, Dimensionality>		ImageType;
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << "Exception thrown while making ImageType container" << std::endl;
		std::cerr << excp << std::endl;
	}

}*/