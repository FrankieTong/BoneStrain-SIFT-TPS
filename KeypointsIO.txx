/*
File:			KeypointsIO.txx

Description:	Class that defines methods for reading and writing key point sets from ITK to and from text files. 
				Also provides functionality for exporting point data to MatLab and Amira.
				
Author:			Frankie (Hoi-Ki) Tong, Jun 22 2015

Changes:		Frankie (Hoi-Ki) Tong, Jun 22 2015
				Initial creation of the file.
*/


#include "KeypointsIO.h"

//Description:	Basic constructor with default values set to empty 
template <class FeatureType, unsigned int ImageDimensionality>
KeypointsIO<FeatureType, ImageDimensionality>::KeypointsIO(): 
	keypointSetPtr(NULL), fileName(""){
}

//Description:	Basic constructor with default values set to parameters defined on construction by the user 
template <class FeatureType, unsigned int ImageDimensionality>
KeypointsIO<FeatureType, ImageDimensionality>::KeypointsIO(PointSetTypePointer * inputkeypointSetPtr, std::string inputfileName): 
	fileName(inputfileName){
	
	this->setkeypointSetPtr(inputkeypointSetPtr);
	
}

//Description: Basic deconstructor
template <class FeatureType, unsigned int ImageDimensionality>
KeypointsIO<FeatureType, ImageDimensionality>::~KeypointsIO(){
}

//Description: 	Copy constructor
//				Note that keypointSetPtr is a pointer to an ITK key point set and is not a clean copy
template <class FeatureType, unsigned int ImageDimensionality>
KeypointsIO<FeatureType, ImageDimensionality>::KeypointsIO(const KeypointsIO &obj) {
	fileName = obj.getFileName();
	keypointSetPtr = obj.getKeypointSetPtr();
}

//Description: 	Assignment constructor
//				Note that keypointSetPtr is a pointer to an ITK key point set and is not a clean copy
template <class FeatureType, unsigned int ImageDimensionality>
KeypointsIO<FeatureType, ImageDimensionality>& KeypointsIO<FeatureType, ImageDimensionality>::operator=(const KeypointsIO & obj){
	fileName = obj.getFileName();
	keypointSetPtr = obj.getKeypointSetPtr();
	return *this;
}

//Description: 	Set command for the keypointSetPtr private data set	
//Input:		inputkeypointSetPtr -(PointSetType::Pointer) Pointer to the key point set to be manipulated by this class
//Return:		(bool) 0
template <class FeatureType, unsigned int ImageDimensionality>
bool KeypointsIO<FeatureType, ImageDimensionality>::setKeypointSetPtr(PointSetTypePointer * inputkeypointSetPtr){
	keypointSetPtr = inputkeypointSetPtr;
	return 0;
}

//Description: 	Get command for the keypointSetPtr private data set	
//Return:		(PointSetType::Pointer) - Pointer to the key point set stored in this class
template <class FeatureType, unsigned int ImageDimensionality>
typename KeypointsIO<FeatureType,ImageDimensionality>::PointSetTypePointerPtr 
KeypointsIO<FeatureType, ImageDimensionality>::getKeypointSetPtr() const{
	return keypointSetPtr;
}

//Description: 	Get command for the fileName private data set	
//Input:		inputfileName - (std::string) Path of the file name that the key point set data is stored in
//Return:		(bool 0)
template <class FeatureType, unsigned int ImageDimensionality>
bool KeypointsIO<FeatureType, ImageDimensionality>::setFileName(std::string inputfileName){
	fileName = inputfileName;
	return 0;
}

//Description: 	Get command for the fileName private data set	
//Return:		(std::string) - String of the file name parameter stored in this class
template <class FeatureType, unsigned int ImageDimensionality>
std::string KeypointsIO<FeatureType, ImageDimensionality>::getFileName() const{
	return fileName;
}

//Description: 	Print out current points in the key point set on screen
template <class FeatureType, unsigned int ImageDimensionality>
void KeypointsIO<FeatureType, ImageDimensionality>::printKeypointSet(){

	//Check if keypointSetPtr is set
	if (keypointSetPtr == NULL) {
		std::cout<<"Keypoint set has not been set. Nothing to print."<<std::endl;
		return;
	}

	//Iterate through the ITK key point set
	PointSetTypePointer keypointSet = *keypointSetPtr;

	unsigned int numberOfPoints = keypointSet->GetNumberOfPoints();
	PointType tmpp;
	
	typedef itk::Vector<double, ImageDimensionality> VectorType;
	VectorType tmpv;

	for (int i = 0; i < numberOfPoints; i++) {
		
		//For each key point in the set, print out each index in the point seperated by a comma
		keypointSet->GetPoint(i,&tmpp);
		tmpv = tmpp.GetVectorFromOrigin();

		for (int j = 0; j < ImageDimensionality; j++) {
			if(j != 0) {
				std::cout<<", ";
			}
			std::cout<<tmpv[j];
		}

		//Each point is seperated by a new line on screen
		std::cout<<std::endl;

	}
	return;
}

//Description:	Read in key point set data from a text file. Keypoint file is a comma seperated file with each point on a new line and each index of the point seperated by commas
//Return:		(bool) 0 if successful
//				(bool) 1 otherwise
template <class FeatureType, unsigned int ImageDimensionality>
bool KeypointsIO<FeatureType, ImageDimensionality>::readKeypointSet(){

	//Check if we have a file name
	if (fileName == "") {
		std::cout<<"No file name to read keypoints into."<<std::endl;
		return 1;
	}

	//Open file
	std::ifstream infile(fileName.c_str());
	if (!infile.is_open()) {
		std::cout<<"Could not find file "<<fileName<<"."<<std::endl;
		return 1;
	}

	//Check if we keypointSetPtr is empty
	if (keypointSetPtr == NULL) {
		std::cout<<"No keypoints to read."<<std::endl;
		return 1;
	}

	//Make new keypoint list and store number of points in keypoint list
	//(*keypointSetPtr) = PointSetType::New();
	PointSetType::PointsContainerPointer keypoints = (*keypointSetPtr)->GetPoints();
	int keypointsSize = 0;

	//Read in each line in the file until we reach end of file
	std::string line;
	while (getline(infile,line)) {

		//Determine if the line is a comment or blank line and skip line
		if (line.find("//") == std::string::npos && line != "") {
			
			//String is not comment, time to do some parsing (delimiter is ","). We only read until all index points are filled
			typedef itk::Vector<double, ImageDimensionality> VectorType;
			VectorType pointIndex;

			size_t currentPosition = 0;
			std::string token;

			for (int i = 0; i < ImageDimensionality; i++) {
				
				//grab characters up to delimeter ","
				currentPosition = line.find(",");

				if ((currentPosition == std::string::npos) && (i < ImageDimensionality - 1)) {
					std::cout<<"Keypoint number "<< i <<" in file "<< fileName <<" is not formatted properly. Stopping read in operation." << std::endl;
					infile.close();
					return 1;
				}

				//If we are reading in last index number that has no delimetor at the end. Read in rest of line as token. Else read up to delimeter
				if (currentPosition == std::string::npos) {
					token = line;
				} else {
					token = line.substr(0, currentPosition);
					line.erase(0,currentPosition + 1);
				}

				//Convert token to double using stringstream and store in keypointIndex
				std::stringstream convert(token);

				if ( !(convert >> pointIndex[i]) ) {
					std::cout<<"Keypoint number "<<i<<" in file "<<fileName<<" is not formatted properly. Stopping read in operation."<<std::endl;
					infile.close();
					return 1;
				} 
			}

			//Insert point into keypointSetPtr
			keypoints->InsertElement(keypointsSize, pointIndex);
			keypointsSize++;
		}

		//clear line to get ready to read in next line
		line = "";

	}

	//Finished reading file. Close file.
	infile.close();
	return 0;
}

//Description:	Write out key point set data to a text file. Keypoint file is a comma seperated file with each point on a new line and each index of the point seperated by commas
//Return:		(bool) 0 if successful
//				(bool) 1 otherwise
template <class FeatureType, unsigned int ImageDimensionality>
bool KeypointsIO<FeatureType, ImageDimensionality>::writeKeypointSet(){

	//Check if we have a file name
	if (fileName == "") {
		std::cout<<"No file name to write keypoints into."<<std::endl;
		return 1;
	}

	//Check if we have keypoints to write
	if (keypointSetPtr == NULL) {
		std::cout<<"No keypoints to write."<<std::endl;
		return 1;
	}

	//Check if config file exists, if it does delete it and get ready to make new config file
	std::ofstream outfile;
	outfile.open(fileName.c_str(),std::ofstream::trunc);

	if (!outfile.is_open()) {
		return 1;
	}

	//Start writing keypoints file
	outfile << "//Contains keypoint locations with each point's index comma deliminated and each line being a new point:" <<std::endl;
	outfile << std::endl;

	if (keypointSetPtr == NULL) {
		std::cout<<"Keypoint set has not been set. Nothing to write."<<std::endl;
		return 1;
	}

	//Iterate through every point in the keypoint set
	PointSetTypePointer keypointSet = *keypointSetPtr;

	unsigned int numberOfPoints = keypointSet->GetNumberOfPoints();
	PointType tmpp;
	
	typedef itk::Vector<double, ImageDimensionality> VectorType;
	VectorType tmpv;

	for (int i = 0; i < numberOfPoints; i++) {
		keypointSet->GetPoint(i,&tmpp);
		tmpv = tmpp.GetVectorFromOrigin();

		//Each index in the point is seperated by a comma and a space
		for (int j = 0; j < ImageDimensionality; j++) {
			outfile<<tmpv[j]<<", ";
		}
		outfile<<std::endl;
	}

	outfile.close();
	return 0;
	
}

//Description:	Write out key point set data to an Amira Mesh file format
//Return:		(bool) 0 if successful
//				(bool) 1 otherwise
template <class FeatureType, unsigned int ImageDimensionality>
bool KeypointsIO<FeatureType, ImageDimensionality>::writeAmiraKeypointSet(){

	//Check if we have a file name
	if (fileName == "") {
		std::cout<<"No file name to write keypoints into."<<std::endl;
		return 1;
	}

	//Check if we have keypoints to write
	if (keypointSetPtr == NULL) {
		std::cout<<"No keypoints to write."<<std::endl;
		return 1;
	}

	//Check if config file exists, if it does delete it and get ready to make new config file
	std::ofstream outfile;
	outfile.open(fileName.c_str(),std::ofstream::trunc);

	if (!outfile.is_open()) {
		return 1;
	}

	if (keypointSetPtr == NULL) {
		std::cout<<"Keypoint set has not been set. Nothing to write."<<std::endl;
		return 1;
	}

	//Just some formatting code to make 2D images easier to render in Amira
	unsigned int dimensions = ImageDimensionality;
	if (ImageDimensionality == 2) {
		dimensions = 3;
	}

	//Start writing keypoints file

	//Amira file header format
	outfile << "# AmiraMesh 3D ASCII 2.0" <<std::endl;
	outfile << std::endl;
	outfile << std::endl;

	PointSetTypePointer keypointSet = *keypointSetPtr;

	unsigned int numberOfPoints = keypointSet->GetNumberOfPoints();

	outfile << "define Markers "<< numberOfPoints <<std::endl;
	outfile << std::endl;
	outfile << "Parameters {" << std::endl;
	outfile << "    ContentType \"LandmarkSet\"," << std::endl;
	outfile << "    NumSets 1" << std::endl;
	outfile << "}" << std::endl;
	outfile << std::endl;
	outfile << "Markers { float[" << dimensions <<"] Coordinates } @1" << std::endl;
	outfile << std::endl;
	outfile << "# Data section follows" << std::endl;
	outfile << "@1" << std::endl;

	outfile.precision(15);
	outfile << std::scientific;

	//Iterate through each point in the key point set and write data to file
	PointType tmpp;
	
	typedef itk::Vector<double, ImageDimensionality> VectorType;
	VectorType tmpv;

	for (int i = 0; i < numberOfPoints; i++) {
		keypointSet->GetPoint(i,&tmpp);
		tmpv = tmpp.GetVectorFromOrigin();

		//Iterate through each point index and write out coordinate positions
		for (int j = 0; j < ImageDimensionality; j++) {
			outfile<<tmpv[j]<<" ";
		}

		//Just some formatting code to make 2D images easier to render in Amira
		if (ImageDimensionality == 2) {
			outfile<<"0.000000000000000"<<" ";
		}

		outfile<<std::endl;
	}

	outfile.close();
	return 0;
	
}

//Description:	Write out key point set data to an MatLab script file format
//Input:		variableName - (std::string) Variable name of the key point list/array in MatLab
//Return:		(bool) 0 if successful
//				(bool) 1 otherwise
template <class FeatureType, unsigned int ImageDimensionality>
bool KeypointsIO<FeatureType, ImageDimensionality>::writeMatlabKeypointSet(std::string variableName){

	//Check if we have a file name
	if (fileName == "") {
		std::cout<<"No file name to write keypoints into."<<std::endl;
		return 1;
	}

	//Check if we have keypoints to write
	if (keypointSetPtr == NULL) {
		std::cout<<"No keypoints to write."<<std::endl;
		return 1;
	}

	//Check if config file exists, if it does delete it and get ready to make new config file
	std::ofstream outfile;
	outfile.open(fileName.c_str(),std::ofstream::trunc);

	if (!outfile.is_open()) {
		return 1;
	}

	if (keypointSetPtr == NULL) {
		std::cout<<"Keypoint set has not been set. Nothing to write."<<std::endl;
		return 1;
	}

	//Start writing keypoints file

	outfile<<variableName<<" = [" <<std::endl;

	PointSetTypePointer keypointSet = *keypointSetPtr;

	unsigned int numberOfPoints = keypointSet->GetNumberOfPoints();

	outfile.precision(15);
	outfile << std::scientific;

	//Iterate through each point in the key point set and write data to file
	PointType tmpp;
	
	typedef itk::Vector<double, ImageDimensionality> VectorType;
	VectorType tmpv;

	for (int i = 0; i < numberOfPoints; i++) {
		keypointSet->GetPoint(i,&tmpp);
		tmpv = tmpp.GetVectorFromOrigin();

		//Iterate through each point index and write out coordinate positions
		for (int j = 0; j < ImageDimensionality; j++) {
			outfile<<tmpv[j]<<" ";
		}

		//Just some formatting code to make 2D images easier to render in Amira
		if (ImageDimensionality == 2) {
			outfile<<"0.000000000000000"<<" ";
		}

		outfile<<";"<<std::endl;
	}

	outfile << "];" << std::endl;
	outfile << std::endl;

	outfile.close();
	return 0;
	
}
