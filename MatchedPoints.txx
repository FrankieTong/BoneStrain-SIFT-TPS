/*
File:			KeypointsIO.txx

Description:	Class that defines methods for reading and writing key point sets from ITK to and from text files. 
				Also provides functionality for exporting point data to MatLab and Amira.
				
				Matched key point sets are basically 2 sets of key points. Each point in the first key point set is linked to a point in the second key point set by index value.
				Thus the number of key points in both key point sets are the same. This formulation is useful when describing matched feature points where points from the first
				image are matched to a corresponding point in the second image in order to describe how the feature points moved between the 2 key point sets.
				
Author:			Frankie (Hoi-Ki) Tong, Jun 22 2015

Changes:		Frankie (Hoi-Ki) Tong, Jun 22 2015
				Initial creation of the file.
*/

#include "MatchedPoints.h"

//Description:	Basic constructor with default values set to empty 
template <class FeatureType, unsigned int ImageDimensionality>
MatchedPoints<FeatureType, ImageDimensionality>::MatchedPoints(){
}

//Description: Basic deconstructor
template <class FeatureType, unsigned int ImageDimensionality>
MatchedPoints<FeatureType, ImageDimensionality>::~MatchedPoints(){
}

//Description: 	Copy constructor
template <class FeatureType, unsigned int ImageDimensionality>
MatchedPoints<FeatureType, ImageDimensionality>::MatchedPoints(const MatchedPoints &obj) {
	obj.getMatchedPoints(&matchedpoints1,&matchedpoints2);
}

//Description: 	Assignment constructor
template <class FeatureType, unsigned int ImageDimensionality>
MatchedPoints<FeatureType, ImageDimensionality>& MatchedPoints<FeatureType, ImageDimensionality>::operator=(const MatchedPoints & obj){
	obj.getMatchedPoints(&matchedpoints1,&matchedpoints2);
	return *this;
}

//Description: 	Set command for the matchedpoints1 and matchedpoints2 private data set	
//Input:		inputmatchedpoints1 -(PointSetTypePointer *) Pointer to store the first matched keypointIO set to be manipulated by this class
//				inputmatchedpoints2 -(PointSetTypePointer *) Pointer to store the second matched keypointIO set to be manipulated by this class
//Return:		(bool) 0
template <class FeatureType, unsigned int ImageDimensionality>
bool MatchedPoints<FeatureType, ImageDimensionality>::setMatchedPoints(PointSetTypePointer * inputmatchedpoints1, PointSetTypePointer * inputmatchedpoints2) {

	matchedpoints1.setKeypointSetPtr(inputmatchedpoints1);
	matchedpoints2.setKeypointSetPtr(inputmatchedpoints2);

	return 0;
}

//Description: 	Set command for the matchedpoints1 and matchedpoints2 private data set	
//Input:		inputmatchedpoints1 -(KeypointsIOType *) Pointer to store the first matched keypointIO set to be manipulated by this class
//				inputmatchedpoints2 -(KeypointsIOType *) Pointer to store the second matched keypointIO set to be manipulated by this class
//Return:		(bool) 0
template <class FeatureType, unsigned int ImageDimensionality>
bool MatchedPoints<FeatureType, ImageDimensionality>::setMatchedPoints(KeypointsIOType * inputmatchedpoints1, KeypointsIOType * inputmatchedpoints2) {

	matchedpoints1 = (*inputmatchedpoints1);
	matchedpoints2 = (*inputmatchedpoints2);

	return 0;
}

//Description: 	Get command for the matchedpoints1 and matchedpoints2 private data set	
//Input:		inputmatchedpoints1 -(PointSetTypePointer *) Pointer to the first matched keypointIO set of this class
//				inputmatchedpoints2 -(PointSetTypePointer *) Pointer to the second matched keypointIO set of this class
//Return:		(bool) 0
template <class FeatureType, unsigned int ImageDimensionality>
bool MatchedPoints<FeatureType, ImageDimensionality>::getMatchedPoints(PointSetTypePointer * inputmatchedpoints1, PointSetTypePointer * inputmatchedpoints2) const{

	inputmatchedpoints1 = matchedpoints1.getKeypointSetPtr();
	inputmatchedpoints2 = matchedpoints2.getKeypointSetPtr();

	return 0;
}

//Description: 	Get command for the matchedpoints1 and matchedpoints2 private data set	
//Input:		inputmatchedpoints1 -(KeypointsIOType *) Pointer to the first matched keypointIO set of this class
//				inputmatchedpoints2 -(KeypointsIOType *) Pointer to the second matched keypointIO set of this class
//Return:		(bool) 0
template <class FeatureType, unsigned int ImageDimensionality>
bool MatchedPoints<FeatureType, ImageDimensionality>::getMatchedPoints(KeypointsIOType * inputmatchedpoints1, KeypointsIOType * inputmatchedpoints2)  const{

	(*inputmatchedpoints1) = matchedpoints1;
	(*inputmatchedpoints2) = matchedpoints2;

	return 0;
}

//Description:	Read in matched key point sets from 2 text files. Keypoint files are comma seperated files with each point on a new line and each index of the point seperated by commas
//Input:		matchedpoints1filename - (std::string) Path name where first set of matched key point data is stored
//				matchedpoints1filename - (std::string) Path name where second set of matched key point data is stored
//Return:		(bool) 0 if successful
//				(bool) 1 otherwise
template <class FeatureType, unsigned int ImageDimensionality>
bool MatchedPoints<FeatureType, ImageDimensionality>::readMatchedPoints(std::string matchedpoints1filename, std::string matchedpoints2filename) {

	matchedpoints1.setFileName(matchedpoints1filename);
	matchedpoints2.setFileName(matchedpoints1filename);

	if (matchedpoints1.readKeypointSet() || matchedpoints2.readKeypointSet()) {
		std::cout<<"Error encountered when reading matched points from files: "<<matchedpoints1filename<<" "<<matchedpoints2filename<<std::endl;
		return 1;
	}

	return 0;
}

//Description:	Write out matched key point sets to 2 text files. Keypoint files are comma seperated files with each point on a new line and each index of the point seperated by commas
//Return:		(bool) 0 if successful
//				(bool) 1 otherwise
template <class FeatureType, unsigned int ImageDimensionality>
bool MatchedPoints<FeatureType, ImageDimensionality>::writeMatchedPoints(std::string matchedpoints1filename, std::string matchedpoints2filename) {

	matchedpoints1.setFileName(matchedpoints1filename);
	matchedpoints2.setFileName(matchedpoints1filename);

	if (matchedpoints1.writeKeypointSet() || matchedpoints2.writeKeypointSet()) {
		std::cout<<"Error encountered when writing matched points to files: "<<matchedpoints1filename<<" "<<matchedpoints2filename<<std::endl;
		return 1;
	}

	return 0;
}

//Description:	Write out matched key point set data to an Amira Mesh format file
//Input:		amiraLandmarksASCIIFileName - (std::string) Path and file name of the Amira mesh file that will store the matched key points
//Return:		(bool) 0 if successful
//				(bool) 1 otherwise
template <class FeatureType, unsigned int ImageDimensionality>
bool MatchedPoints<FeatureType, ImageDimensionality>::writeAmiraMatchedPoints(std::string amiraLandmarksASCIIFileName) {

	//Check if we have a file name
	if (amiraLandmarksASCIIFileName == "") {
		std::cout<<"No file name to write keypoints into."<<std::endl;
		return 1;
	}

	PointSetTypePointer * matchedKeypointSetPtr1 = matchedpoints1.getKeypointSetPtr();
	PointSetTypePointer * matchedKeypointSetPtr2 = matchedpoints2.getKeypointSetPtr();

	//Check if we have keypoints to write
	if (matchedKeypointSetPtr1 == NULL || matchedKeypointSetPtr2 == NULL) {
		std::cout<<"No keypoints to write."<<std::endl;
		return 1;
	}

	//Check if amiraMatchedPointsFileName file exists, if it does delete it and get ready to make new config file
	std::ofstream outfile;
	outfile.open(amiraLandmarksASCIIFileName.c_str(),std::ofstream::trunc);

	if (!outfile.is_open()) {
		return 1;
	}

	//Just some formatting code to make 2D images easier to render in Amira
	unsigned int dimensions = ImageDimensionality;
	if (ImageDimensionality == 2) {
		dimensions = 3;
	}

	//Start writing Amira LandmarksASCII file
	
	//Amira file header format
	outfile << "# AmiraMesh 3D ASCII 2.0" <<std::endl;
	outfile << std::endl;
	outfile << std::endl;

	PointSetTypePointer matchedKeypointSet1 = *matchedKeypointSetPtr1;
	PointSetTypePointer matchedKeypointSet2 = *matchedKeypointSetPtr2;

	unsigned int numberOfPoints = matchedKeypointSet1->GetNumberOfPoints();

	outfile << "define Markers "<< numberOfPoints <<std::endl;
	outfile << std::endl;
	outfile << "Parameters {" << std::endl;
	outfile << "    ContentType \"LandmarkSet\"," << std::endl;
	outfile << "    NumSets 2" << std::endl;
	outfile << "}" << std::endl;
	outfile << std::endl;
	outfile << "Markers { float[" << dimensions <<"] Coordinates } @1" << std::endl;
	outfile << "Markers { float[" << dimensions <<"] Coordinates2 } @2" << std::endl;
	outfile << std::endl;
	outfile << "# Data section follows" << std::endl;
	outfile << "@1" << std::endl;

	outfile.precision(15);
	outfile << std::scientific;

	PointType tmpp;
	
	typedef itk::Vector<double, ImageDimensionality> VectorType;
	VectorType tmpv;

	//Iterate through each point in the first key point set and write data to file
	for (int i = 0; i < numberOfPoints; i++) {
		matchedKeypointSet1->GetPoint(i,&tmpp);
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

	outfile << std::endl;
	outfile << "@2" << std::endl;

	//Iterate through each point in the second key point set and write data to file
	for (int i = 0; i < numberOfPoints; i++) {
		matchedKeypointSet2->GetPoint(i,&tmpp);
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

	outfile << std::endl;
	outfile << std::endl;

	outfile.close();
	return 0;
}

//Description: 	Print out matched key point set points on screen
template <class FeatureType, unsigned int ImageDimensionality>
void MatchedPoints<FeatureType, ImageDimensionality>::printMatchedPoints() {

	PointSetTypePointer * matchedPointSetPtr1 = NULL;
	PointSetTypePointer * matchedPointSetPtr2 = NULL;
	matchedPointSetPtr1 = matchedpoints1.getKeypointSetPtr();
	matchedPointSetPtr2 = matchedpoints2.getKeypointSetPtr();

	//Check if keypointSetPtr is set
	if (matchedPointSetPtr1 == NULL || matchedPointSetPtr2 == NULL) {
		std::cout<<"Keypoint sets has not been set."<<std::endl;
		return;
	}

	//Iterate through the ITK key point sets
	PointSetTypePointer matchedPointSet1 = (*matchedPointSetPtr1);
	PointSetTypePointer matchedPointSet2 = (*matchedPointSetPtr2);

	int numberOfPoints = matchedPointSet1->GetNumberOfPoints();
	PointType tmpp1, tmpp2;
	
	typedef itk::Vector<double, ImageDimensionality> VectorType;
	VectorType tmpv1, tmv2;

	//Each matched key point is printed in the following manner on screen:
	//
	//(x1 ,y1 ,z1) -> (x2, y2, z2)
	//
	for (int i = 0; i < numberOfPoints; i++) {
		matchedPointSet1->GetPoint(i,&tmpp1);
		matchedPointSet2->GetPoint(i,&tmpp2);
		tmpv1 = tmpp1.GetVectorFromOrigin();
		tmpv2 = tmpp2.GetVectorFromOrigin();

		//Print coordinatees of first matched point
		std::cout<<"(";
		for (int j = 0; j < ImageDimensionality; j++) {
			if(j != 0) {
				std::cout<<", ";
			}
			std::cout<<tmpv1[j];
		}

		std::cout<<") -> (";

		//Print coordinatees of second matched point
		std::cout<<"(";
		for (int j = 0; j < ImageDimensionality; j++) {
			if(j != 0) {
				std::cout<<", ";
			}
			std::cout<<tmpv2[j];
		}

		std::cout<<")"<<std::endl;

	}

	return;
}

//Description: 	Print out matched key point set points as well as the euclidean distance between the matched keypoints on screen
template <class FeatureType, unsigned int ImageDimensionality>
void MatchedPoints<FeatureType, ImageDimensionality>::printMatchedPointsWithDistance() {

	PointSetTypePointer * matchedPointSetPtr1 = NULL;
	PointSetTypePointer * matchedPointSetPtr2 = NULL;
	matchedPointSetPtr1 = matchedpoints1.getKeypointSetPtr();
	matchedPointSetPtr2 = matchedpoints2.getKeypointSetPtr();

	//Check if keypointSetPtr is set
	if (matchedPointSetPtr1 == NULL || matchedPointSetPtr2 == NULL) {
		std::cout<<"Keypoint sets has not been set."<<std::endl;
		return;
	}

	//Iterate through the ITK key point sets
	PointSetTypePointer matchedPointSet1 = (*matchedPointSetPtr1);
	PointSetTypePointer matchedPointSet2 = (*matchedPointSetPtr2);

	int numberOfPoints = matchedPointSet1->GetNumberOfPoints();
	PointType tmpp1, tmpp2;
	
	typedef itk::Vector<double, ImageDimensionality> VectorType;
	VectorType tmpv1, tmv2;

	//Calculate the euclidean distance between all matched key points in the class object
	float * distance = this->getDistance();

	if (distance == NULL) {
		std::cout<<"Could not calculate distance when executing printMatchedPointsWithDistance()."<<std::endl;
		return;
	}

	//Each matched key point is printed in the following manner on screen:
	//
	//(x1 ,y1 ,z1) -> (x2, y2, z2) Distance: (distance between the 2 points)
	//
	for (int i = 0; i < numberOfPoints; i++) {
		matchedPointSet1->GetPoint(i,&tmpp1);
		matchedPointSet2->GetPoint(i,&tmpp2);
		tmpv1 = tmpp1.GetVectorFromOrigin();
		tmpv2 = tmpp2.GetVectorFromOrigin();

		//Print coordinatees of first matched point
		std::cout<<"(";
		for (int j = 0; j < ImageDimensionality; j++) {
			if(j != 0) {
				std::cout<<", ";
			}
			std::cout<<tmpv1[j];
		}

		std::cout<<") -> (";

		//Print coordinatees of second matched point
		std::cout<<"(";
		for (int j = 0; j < ImageDimensionality; j++) {
			if(j != 0) {
				std::cout<<", ";
			}
			std::cout<<tmpv2[j];
		}

		std::cout<<") Distance: ";

		//Print distance between the 2 points
		std::cout<<distance[i]<<std::endl;

	}

	//Clean up float pointer used to store distance values
	delete[] distance;

	return;

}

//Description: 	Calculate the euclidean distance between the matched key points in the class object
//Return:		(float *) - A float pointer containing a 1 x n array with the distances between the key points in the matched key point sets defined by array/set index value. 
//				(float *) - NULL if distances could not be calculated.
template <class FeatureType, unsigned int ImageDimensionality>
float * MatchedPoints<FeatureType, ImageDimensionality>::getDistance(){

	PointSetTypePointer * matchedPointSetPtr1 = NULL;
	PointSetTypePointer * matchedPointSetPtr2 = NULL;
	matchedPointSetPtr1 = matchedpoints1.getKeypointSetPtr();
	matchedPointSetPtr2 = matchedpoints2.getKeypointSetPtr();

	//Check to make sure the matched key point sets are defined
	if (matchedPointSetPtr1 == NULL || matchedPointSetPtr2 == NULL) {
		std::cout<<"Keypoint sets has not been set."<<std::endl;
		return NULL;
	}

	PointSetTypePointer matchedPointSet1 = (*matchedPointSetPtr1);
	PointSetTypePointer matchedPointSet2 = (*matchedPointSetPtr2);

	int numberOfPoints = matchedPointSet1->GetNumberOfPoints();
	PointType tmpp1, tmpp2;

	float * distance = new float[numberOfPoints]; 

	//For each point in the matched key point sets, calculate the euclidean distance between the 2 points.
	for (int i = 0; i < numberOfPoints; i++) {
		matchedPointSet1->GetPoint(i,&tmpp1);
		matchedPointSet2->GetPoint(i,&tmpp2);
	
		distance[i] = tmpp1.EuclideanDistanceTo(tmpp2);
	}

	return distance;

}

//Description: 	Returns the total number of matched points in the matched key point set
//Return:		(int) -  Total number of matched points in the matched key point set
template <class FeatureType, unsigned int ImageDimensionality>
int MatchedPoints<FeatureType, ImageDimensionality>::getDistanceIndex(){

	PointSetTypePointer * matchedPointSetPtr1 = NULL;
	PointSetTypePointer * matchedPointSetPtr2 = NULL;
	matchedPointSetPtr1 = matchedpoints1.getKeypointSetPtr();
	matchedPointSetPtr2 = matchedpoints2.getKeypointSetPtr();

	if (matchedPointSetPtr1 == NULL || matchedPointSetPtr2 == NULL) {
		std::cout<<"Keypoint sets has not been set."<<std::endl;
		return 0;
	}

	PointSetTypePointer matchedPointSet1 = (*matchedPointSetPtr1);
	PointSetTypePointer matchedPointSet2 = (*matchedPointSetPtr2);

	int numberOfPoints = matchedPointSet1->GetNumberOfPoints();

	return numberOfPoints;

}

//Description: 	Print out the euclidean distance between the matched keypoints on screen
template <class FeatureType, unsigned int ImageDimensionality>
void MatchedPoints<FeatureType, ImageDimensionality>::printDistance() {

	//Calculate the euclidean distance between matched key point sets
	float * distance = this->getDistance();

	if (distance == NULL) {
		std::cout<<"Could not calculate distance when executing printDistance()."<<std::endl;
		return;
	}

	//Display distance values
	std::cout<<"Distance Values:"<<std::endl;

	for (int i = 0; i < this->getDistanceIndex(); i++) {
		std::cout<<distance[i]<<std::endl;
	}

	std::cout<<std::endl;

	//Clean up float pointer used to store distance values
	delete[] distance;

	return;

}

//Description:	Write out the euclidean distance values to a text file. Each distance value is given on a new line in the file
//Input:		distancefilename - (std::string) Path and file name of the text file that will store euclidean distance values
//Return:		(bool) 0 if successful
//				(bool) 1 otherwise
template <class FeatureType, unsigned int ImageDimensionality>
bool MatchedPoints<FeatureType, ImageDimensionality>::writeDistance(std::string distancefilename) {

	//Calculate euclidean distance
	float * distance = this->getDistance();

	if (distance == NULL) {
		std::cout<<"Could not calculate distance when executing writeDistance()."<<std::endl;
		return 1;
	}

	//Check if we have a file name
	if (distancefilename == "") {
		std::cout<<"No file name to write keypoints into."<<std::endl;
		return 1;
	}

	//Check if output file exists, if it does delete it and get ready to make new output file
	std::ofstream outfile;
	outfile.open(distancefilename.c_str(),std::ofstream::trunc);

	if (!outfile.is_open()) {
		return 1;
	}

	//Start writing to file
	outfile<<"//Contains distance values of matched points:"<<std::endl;
	outfile<<std::endl;

	for (int i = 0; i < this->getDistanceIndex(); i++) {
		outfile<<distance[i]<<std::endl;
	}

	outfile<<std::endl;
	outfile.close();

	delete[] distance;

	return 0;

}

//Description:	Write out matched key point set data to a MatLab script file. The script file will store the matched key points in "fixed_points" and "moving_points" array variables respectively.
//Input:		matlabScriptFileName - (std::string) Path and file name of the MatLab script file that will store the matched key points
//Return:		(bool) 0 if successful
//				(bool) 1 otherwise
template <class FeatureType, unsigned int ImageDimensionality>
bool MatchedPoints<FeatureType, ImageDimensionality>::writeMatlabMatchedPoints(std::string matlabScriptFileName) {

	//Check if we have a file name
	if (matlabScriptFileName == "") {
		std::cout<<"No file name to write keypoints into."<<std::endl;
		return 1;
	}

	PointSetTypePointer * matchedKeypointSetPtr1 = matchedpoints1.getKeypointSetPtr();
	PointSetTypePointer * matchedKeypointSetPtr2 = matchedpoints2.getKeypointSetPtr();

	//Check if we have keypoints to write
	if (matchedKeypointSetPtr1 == NULL || matchedKeypointSetPtr2 == NULL) {
		std::cout<<"No keypoints to write."<<std::endl;
		return 1;
	}

	//Check if matlabScriptFileName file exists, if it does delete it and get ready to make new config file
	std::ofstream outfile;
	outfile.open(matlabScriptFileName.c_str(),std::ofstream::trunc);

	if (!outfile.is_open()) {
		return 1;
	}

	//Just some formatting code to make 2D images easier to render in Matlab
	unsigned int dimensions = ImageDimensionality;
	if (ImageDimensionality == 2) {
		dimensions = 3;
	}

	//Start writing Matlab script file

	//Define "fixed_points" variable
	outfile<<"fixed_points = [" <<std::endl;

	PointSetTypePointer matchedKeypointSet1 = *matchedKeypointSetPtr1;
	PointSetTypePointer matchedKeypointSet2 = *matchedKeypointSetPtr2;

	unsigned int numberOfPoints = matchedKeypointSet1->GetNumberOfPoints();

	outfile.precision(15);
	outfile << std::scientific;

	PointType tmpp;
	
	typedef itk::Vector<double, ImageDimensionality> VectorType;
	VectorType tmpv;

	//Iterate through each point in the first list of matched key points
	for (int i = 0; i < numberOfPoints; i++) {
		matchedKeypointSet1->GetPoint(i,&tmpp);
		tmpv = tmpp.GetVectorFromOrigin();

		//Write the coordinate values of each point
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
	
	//Define "moving_points" variable
	outfile << "moving_points =[" <<std::endl;

	//Iterate through each point in the second list of matched key points
	for (int i = 0; i < numberOfPoints; i++) {
		matchedKeypointSet2->GetPoint(i,&tmpp);
		tmpv = tmpp.GetVectorFromOrigin();

		//Write the coordinate values of each point
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