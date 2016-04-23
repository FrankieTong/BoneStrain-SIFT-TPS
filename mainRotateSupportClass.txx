/*
File:			mainRotateSupportClass.txx

Description:	Provides support for analyzing feature point registration accuracy.

				Given a set of matched feature points and the ideal transform used to map the first point set to the second point set,
				this class applies the inverse of the transform to map the second point set back to its original position where it can
				then be used to meaasure "error" using the first point set as reference.
				
				The class also includes the ability to match 2 sets of points based on minimum euclidean distance.
				
				Overall, this class itself is not ideal as concepts used in this class would usually be deemed flawed. Consider using another
				class when trying to determine positional errors in matching feature points.

Inputs:			argc - (int) Number of input arguments into the function
				argv - (char *) Character pointer containing the inputs of this function

Author:			Frankie (Hoi-Ki) Tong, Nov 1 2014

Changes:		Frankie (Hoi-Ki) Tong, Nov 1 2014
				Initial creation of the file.
				
				Frankie (Hoi-Ki) Tong, Apr 20 2016
				File and class renamed to mainRotateSupportClass
*/

#include "mainRotateSupportClass.h"

//Description:	Basic constructor with default values set to empty 
template <class PixelType, unsigned int ImageDimensionality>
mainRotateSupportClass<PixelType, ImageDimensionality>::mainRotateSupportClass() 
{
	distanceBetweenPoints = NULL;
	distanceBetweenPointsIndex = 0;
}

//Description: Basic deconstructor
template <class PixelType, unsigned int ImageDimensionality>
mainRotateSupportClass<PixelType, ImageDimensionality>::~mainRotateSupportClass(){
	if (distanceBetweenPoints != NULL) {delete[] distanceBetweenPoints;}
}

//Description: 	Copy constructor
template <class PixelType, unsigned int ImageDimensionality>
mainRotateSupportClass<PixelType, ImageDimensionality>::mainRotateSupportClass(const mainRotateSupportClass &obj) {
	distanceBetweenPoints = obj.getDistanceBetweenPoints();
}

//Description: 	Assignment constructor
template <class PixelType, unsigned int ImageDimensionality>
mainRotateSupportClass<PixelType, ImageDimensionality>& mainRotateSupportClass<PixelType, ImageDimensionality>::operator=(const mainRotateSupportClass & obj){
	distanceBetweenPoints = obj.getDistanceBetweenPoints();
	return *this;
}

//Description: 	Set command for the inputKeypoints1 and inputKeypoints2 private data set	
//Input:		inputKeypoints1 -(PointSetTypePointer *) Pointer to store the first key point set to be used for analysis
//				inputKeypoints2 -(PointSetTypePointer *) Pointer to store the second key point set to be used for analysis
//Return:		(bool) 0
template <class PixelType, unsigned int ImageDimensionality>	
bool mainRotateSupportClass<PixelType, ImageDimensionality>::setImageKeyPoints(PointSetTypePointer * inputKeypoints1, PointSetTypePointer * inputKeypoints2) {
	keypointsptr1 = inputKeypoints1;
	keypointsptr2 = inputKeypoints2;

	return 0;
}

//Description: 	Set command for the inputmatchedkeypoints1 and inputmatchedkeypoints2 private data set	
//Input:		inputmatchedkeypoints1 -(PointSetTypePointer *) Pointer to store the first key point set to be used for retrival of matched points
//				inputmatchedkeypoints2 -(PointSetTypePointer *) Pointer to store the second key point set to be used for retrival of matched points
//Return:		(bool) 0
template <class PixelType, unsigned int ImageDimensionality>	
bool mainRotateSupportClass<PixelType, ImageDimensionality>::setMatchedKeyPoints(PointSetTypePointer * inputmatchedkeypoints1, PointSetTypePointer * inputmatchedkeypoints2) {
	
	matchedkeypointsptr1 = inputmatchedkeypoints1;
	matchedkeypointsptr2 = inputmatchedkeypoints2;
	
	return 0;
}

//Description: 	Get command for distanceBetweenPoints private data set	
//Return:		(float *) Float array containing the distance error between each index of matchedkeypointsptr1 and inputmatchedkeypoints2
//				(float *) NULL if updateDistance() or updateMatched() have not been used yet
template <class PixelType, unsigned int ImageDimensionality>	
float* mainRotateSupportClass<PixelType, ImageDimensionality>::getDistanceBetweenPoints() const {
	if(distanceBetweenPoints == NULL) { return NULL; }
	float * tmp = new float[this->getDistanceBetweenPointsIndex()];
	for (int i = 0; i<this->getDistanceBetweenPointsIndex(); i++) {
		tmp[i] = distanceBetweenPoints[i];
	}
	return tmp;
}

//Description: 	Set command for the inverse_transform private data set	
//Input:		input_inverse_transform -(TransformType::Pointer) ITK image pointer storing the Affine transform mapping keypointsptr2 back to keypointsptr1
//Return:		(bool) 0
template <class PixelType, unsigned int ImageDimensionality>	
bool mainRotateSupportClass<PixelType, ImageDimensionality>::setInverseTransform(typename TransformType::Pointer input_inverse_transform) {
	inverse_transform = input_inverse_transform;
	return 0;
}

//Description: 	Get command for distanceBetweenPointsIndex private data set	
//Return:		(int) - int with the size of the distanceBetweenPoints array
template <class PixelType, unsigned int ImageDimensionality>
int mainRotateSupportClass<PixelType, ImageDimensionality>::getDistanceBetweenPointsIndex() const {
	return distanceBetweenPointsIndex;
}
//Description: 	Set command for distanceBetweenPointsIndex private data set	
//Input:		(int) - int with the size of the distanceBetweenPoints array
//Return:		(bool) 0
template <class PixelType, unsigned int ImageDimensionality>
bool mainRotateSupportClass<PixelType, ImageDimensionality>::setDistanceBetweenPointsIndex(int inputdistanceBetweenPointsIndex) {
	distanceBetweenPointsIndex = inputdistanceBetweenPoints;
	return 0;
}

//Description: 	Apply the inverse transform stored in inverse_transform onto keypointsptr2 and compare agaisnt the original values in keypointsptr1.
//				Then generate a set of matched points between keypointsptr1 and keypointsptr2 using a greedy approach with the euclidean distance as
//				the minimization factor and store these matched points in matchedkeypointsptr1 and matchedkeypointsptr2. Calculate the positional error between
//				matchedkeypointsptr1 and the transformed matchedkeypointsptr2 and store them in distanceBetweenPoints.
//Return:		(bool) 0
template <class PixelType, unsigned int ImageDimensionality>
bool mainRotateSupportClass<PixelType, ImageDimensionality>::updateDistance() {

	//Make sure all inputs are in place
	if(distanceBetweenPoints != NULL) { delete[] distanceBetweenPoints; }
	
	//Compare keypoints after inverse rotation and find closest matched points.
	//List out keypoint matches (based on distance) and distance between the keypoints
	//No real way to verify other than through visual checking really...
	unsigned long numpoints1, numpoints2;
    numpoints1 = (*keypointsptr1)->GetNumberOfPoints();
    std::cout << "Keypoints1 Found: " << numpoints1 << std::endl;
    numpoints2 = (*keypointsptr2)->GetNumberOfPoints();
    std::cout << "Keypoints2 Found: " << numpoints2 << std::endl;

	//Make spare pointer for keypoints1 so we don't overwrite original keypoints1
	PointSetType::Pointer copykeypoints1 = PointSetType::New();
	PointsContainerPointer copypoints1 = copykeypoints1->GetPoints();

	for (unsigned int i = 0;  i < numpoints1; ++i) {
		//Reverse the transform for keypoints2
		PointType pp;
		(*keypointsptr1)->GetPoint(i, &pp);
		
		copypoints1->InsertElement(i, pp);
	}

	//Transform the second set of keypoints with the inverse transform to prepare for matching
	PointSetType::Pointer transformedkeypoints2 = PointSetType::New();
	PointsContainerPointer transformedpoints2 = transformedkeypoints2->GetPoints();

	for (unsigned int i = 0;  i < numpoints2; ++i) {
		//Reverse the transform for keypoints2
		PointType pp, tmp;
		(*keypointsptr2)->GetPoint(i, &tmp);
		pp = inverse_transform->TransformPoint(tmp);
		
		transformedpoints2->InsertElement(i, pp);
	}
	
	//Make 2 new point conatianer pointers to store matched points in their proper containers (by matching indicies)
	PointSetType::Pointer * matchedkeypoints1 = (PointSetTypePointer *) matchedkeypointsptr1;
	PointSetType::Pointer * matchedkeypoints2 = (PointSetTypePointer *) matchedkeypointsptr2;
	
	*matchedkeypoints1 = PointSetType::New();
	*matchedkeypoints2 = PointSetType::New();
	PointsContainerPointer matchedpoints1 = (*matchedkeypoints1)->GetPoints();
	PointsContainerPointer matchedpoints2 = (*matchedkeypoints2)->GetPoints();

	//Want to match the smaller keypoint set to the larger keypoint set (for simplicity sake)
	bool flipKeypoints = 0;
	if (numpoints1 > numpoints2) {
	
		PointSetTypePointer tempKeypoints = copykeypoints1;
		copykeypoints1 = transformedkeypoints2;
		transformedkeypoints2 = tempKeypoints;
		
		int tempNum = numpoints1;
		numpoints1 = numpoints2;
		numpoints2 = tempNum;

		flipKeypoints = 1;
		
	}
	
	//Make new float array to store distance values
	distanceBetweenPoints = new float[numpoints1];
	distanceBetweenPointsIndex = numpoints1;
	
	//Make list of pointers that have already been matched for list 2
	bool* matchedPoint = new bool[numpoints2];
	for (unsigned int i = 0; i < numpoints2; ++i) { matchedPoint[i] = false; }
	
	for (unsigned int i = 0; i < numpoints1; ++i) {
		
		PointType pp1;
		PointType tmpp, pp2;
		copykeypoints1->GetPoint(i, &pp1);	
	
		float minDistance = std::numeric_limits<float>::max();
		unsigned int minDistanceIndex = -1;
		
		for(unsigned int j = 0; j < numpoints2; ++j) {
		
			//Reverse the transform for keypoints2
			transformedkeypoints2->GetPoint(j, &pp2);
			
			//Find the point in keypooints2 that is closet to current point
			if (minDistance > pp1.EuclideanDistanceTo(pp2) && matchedPoint[j] == false ) {
				minDistance = pp1.EuclideanDistanceTo(pp2);
				minDistanceIndex = j;
			}
			
		}
		
		//Once closest keypoint has been found, record the match
		matchedpoints1->InsertElement(i,pp1);
		transformedkeypoints2->GetPoint(minDistanceIndex, &pp2);
		matchedpoints2->InsertElement(i,pp2);
		distanceBetweenPoints[i] = minDistance;
		
		//Remove matched keypoint from search pool
		matchedPoint[minDistanceIndex] = true;
	}

	delete[] matchedPoint;
	
	//Flip back the matched points if flip for keypoints occured
	if (flipKeypoints) {
	
		PointSetTypePointer tempKeypoints = *matchedkeypoints1;
		*matchedkeypoints1 = *matchedkeypoints2;
		*matchedkeypoints2 = tempKeypoints;
		
	}
	
	return 0;
	
}

//Description: 	Apply the inverse transform stored in inverse_transform onto keypointsptr2 and compare agaisnt the original values in keypointsptr1.
//				Store the values of keypointsptr1 and the transformed keypointsptr2 in matchedkeypointsptr1 and matchedkeypointsptr2 respectively.
//				Calculate the positional error between matchedkeypointsptr1 and the transformed matchedkeypointsptr2 and store them in distanceBetweenPoints.
//Return:		(bool) 0
template <class PixelType, unsigned int ImageDimensionality>
bool mainRotateSupportClass<PixelType, ImageDimensionality>::updateMatched() {
	//Make sure all inputs are in place
	if(distanceBetweenPoints != NULL) { delete[] distanceBetweenPoints; }
	
	//Compare keypoints after inverse rotation and find closest matched points.
	unsigned long numpoints1, numpoints2;
    numpoints1 = (*keypointsptr1)->GetNumberOfPoints();
    std::cout << "Keypoints1 Found: " << numpoints1 << std::endl;
    numpoints2 = (*keypointsptr2)->GetNumberOfPoints();
    std::cout << "Keypoints2 Found: " << numpoints2 << std::endl;

	//Make spare pointer for keypoints1 so we don't overwrite original keypoints1
	PointSetType::Pointer copykeypoints1 = PointSetType::New();
	PointsContainerPointer copypoints1 = copykeypoints1->GetPoints();

	for (unsigned int i = 0;  i < numpoints1; ++i) {
		//Reverse the transform for keypoints2
		PointType pp;
		(*keypointsptr1)->GetPoint(i, &pp);
		
		copypoints1->InsertElement(i, pp);
	}

	//Transform the second set of keypoints with the inverse transform to prepare for matching
	PointSetType::Pointer transformedkeypoints2 = PointSetType::New();
	PointsContainerPointer transformedpoints2 = transformedkeypoints2->GetPoints();

	for (unsigned int i = 0;  i < numpoints2; ++i) {
		//Reverse the transform for keypoints2
		PointType pp, tmp;
		(*keypointsptr2)->GetPoint(i, &tmp);
		pp = inverse_transform->TransformPoint(tmp);
		
		transformedpoints2->InsertElement(i, pp);
	}
	
	//Make 2 new point conatianer pointers to store matched points in their proper containers (by matching indicies)
	PointSetType::Pointer * matchedkeypoints1 = (PointSetTypePointer *) matchedkeypointsptr1;
	PointSetType::Pointer * matchedkeypoints2 = (PointSetTypePointer *) matchedkeypointsptr2;
	
	*matchedkeypoints1 = PointSetType::New();
	*matchedkeypoints2 = PointSetType::New();
	PointsContainerPointer matchedpoints1 = (*matchedkeypoints1)->GetPoints();
	PointsContainerPointer matchedpoints2 = (*matchedkeypoints2)->GetPoints();

	//Make new float array to store distance values
	distanceBetweenPoints = new float[numpoints1];
	distanceBetweenPointsIndex = numpoints1;

	//calculate distance between the 2 points and store each point into the matchedpoints container
	for (unsigned int i = 0; i < distanceBetweenPointsIndex; ++i) {
		
		PointType pp1;
		PointType pp2;
		copykeypoints1->GetPoint(i, &pp1);
		transformedkeypoints2->GetPoint(i, &pp2);
	
		//Calculate distance and store match
		matchedpoints1->InsertElement(i,pp1);
		matchedpoints2->InsertElement(i,pp2);
		distanceBetweenPoints[i] = pp1.EuclideanDistanceTo(pp2);
	}

	return 0;

}