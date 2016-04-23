/*
File:			ImageConfigFileImport.h

Description:	Reading and writing image configuration parameters.

Author:			Frankie (Hoi-Ki) Tong, Sep 6 2014

Changes:		Frankie (Hoi-Ki) Tong, Sep 6 2014
				Initial creation of the file.
				
				Frankei (Hoi-Ki) Tong, Oct 4 2014
				Depreciated in exchange for ImageIOConfig structure due to awkard method to store configuration parameters
*/


#include "ImageConfigFileImport.h"

//Description: 	Import configuration file given the path name of the image configuration file and 2 std::list<std::string> for configuration parameter names and contents
//Input:		ConfigFileName - (std::string) String with configuration file and path name 
//				ConfigFileParam - (std::list<std::string>) String list to store configuration parameter names
//				ConfigFileParam - (std::list<std::string>) String list to store configuration parameter contents
//Return:		(bool) 0 if successful
//				(bool) 1 if there was a problem importing the configuration parameters
int ImageConfigFileImport ( std::string ConfigFileName, std::list<std::string> * ConfigFileParam, std::list<std::string> * ConfigFileContent) {

	//Changing approach, making a list of variable names and stuff and sending it back up for configuration purposes
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