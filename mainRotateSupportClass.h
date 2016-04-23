/*
File:			mainRotateSupportClass.h

Description:	Header file for mainRotateSupportClass.txx.  Provides support for analyzing feature point registration accuracy.

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

#include "itkAffineTransform.h"
#include "itkIdentityTransform.h"
#include "itkResampleImageFilter.h"

#include "itkPointSet.h"
#include "itkVector.h"
#include <limits>

//This is going to be a pipeline with all inputs feed into and out of this managed by the user using this pipeline
//Be careful of memory leaks!!!
#ifndef __MAINROTATESUPPORTCLASS_H__
#define __MAINROTATESUPPORTCLASS_H__

template <class PixelType, unsigned int ImageDimensionality> 
class mainRotateSupportClass {
	public:
		
		//Typedefs to make ITK point sets easier to manipulate
		typedef itk::Array<float> FeatureType;

		typedef itk::PointSet< FeatureType, ImageDimensionality,
			itk::DefaultStaticMeshTraits< FeatureType, ImageDimensionality, ImageDimensionality, double > > PointSetType;
					  
		typedef typename PointSetType::PointType PointType;
		typedef typename PointSetType::Pointer PointSetTypePointer;
		typedef typename PointSetType::PointsContainerPointer PointsContainerPointer;
		typedef typename itk::ScalableAffineTransform< double, ImageDimensionality > TransformType;
	  
		//Standard constructors and descructors
		mainRotateSupportClass();                         // constructor; initialize the list to be empty
		~mainRotateSupportClass();						// destructor;

		mainRotateSupportClass(const mainRotateSupportClass &obj);	//copy constuctor
		mainRotateSupportClass& operator=(const mainRotateSupportClass & obj); //assignment constructor
		
		//Set matched key point sets to be analyzed for error
		bool setImageKeyPoints(PointSetTypePointer * inputKeypoints1, PointSetTypePointer * inputKeypoints2);
		
		//Set key point set pointers to retrieve analyzed key point sets after the transform has been applied to the second set
		bool setMatchedKeyPoints(PointSetTypePointer * inputmatchedkeypoints1, PointSetTypePointer * inputmatchedkeypoints2);
		
		//Retireve distance between matched key points from private data set
		float * getDistanceBetweenPoints() const;

		//Set ITK inverse transform to be applied to the second feature point set
		bool setInverseTransform(typename TransformType::Pointer input_inverse_transform);
		
		//Get and set method to determine how many matched points are in the list
		int getDistanceBetweenPointsIndex() const;
		bool setDistanceBetweenPointsIndex(int inputdistanceBetweenPointsIndex);

		//Match keypoints in keypointsptr1 and keypointsptr2 based on Ecludian distance after inverse_transform has been applied to keypointsptr2
		bool updateDistance();

		//Calculates error between keypointsptr1 and keypointsptr2 after inverse_transform has been applied to keypointsptr2
		bool updateMatched();
		
	private:

		//Private data classes
		PointSetTypePointer * keypointsptr1;				//pointer to store the first key point set
		PointSetTypePointer * keypointsptr2;				//pointer to store the second key point set
		PointSetTypePointer * matchedkeypointsptr1;			//pointer to store the first key point set for retrieval by the main function
		PointSetTypePointer * matchedkeypointsptr2;			//pointer to store the second key point set for retrieval by the main function
		float * distanceBetweenPoints;						//float pointer to store the distance values between matchedkeypointsptr1 and matchedkeypointsptr2
		typename TransformType::Pointer inverse_transform;	//pointer to store the inverse image transform used to map keypointsptr2 back to keypointsptr1
		int distanceBetweenPointsIndex;						//int to store how many matched points were found
  
};

#include "mainRotateSupportClass.txx"

#endif