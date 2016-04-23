/*
File:			MatchedPoints.h

Description:	Header file for MatchedPoints.hxx. Class that defines methods for reading and writing matched key point sets from ITK to and from text files. 
				Also provides functionality for exporting matched point data to MatLab and Amira.
				
				Matched key point sets are basically 2 sets of key points. Each point in the first key point set is linked to a point in the second key point set by index value.
				Thus the number of key points in both key point sets are the same. This formulation is useful when describing matched feature points where points from the first
				image are matched to a corresponding point in the second image in order to describe how the feature points moved between the 2 key point sets.

Author:			Frankie (Hoi-Ki) Tong, Jun 22 2015

Changes:		Frankie (Hoi-Ki) Tong, Jun 22 2015
				Initial creation of the file.
*/

#include "KeypointsIO.h"

#ifndef __MatchedPoints_H__
#define __MatchedPoints_H__

template <class FeatureType, unsigned int ImageDimensionality> 
class MatchedPoints {
	public:
		
		//Typedefs to make ITK point sets easier to manipulate
		typedef itk::PointSet< FeatureType, ImageDimensionality,
		itk::DefaultStaticMeshTraits< FeatureType, ImageDimensionality, ImageDimensionality, double > > PointSetType;
				  
		typedef typename PointSetType::PointType PointType;
		typedef typename PointSetType::Pointer PointSetTypePointer;

		typedef typename KeypointsIO<FeatureType, ImageDimensionality> KeypointsIOType;
	  
		//Standard constructors and descructors
		MatchedPoints();                         // constructor; initialize the list to be empty
		~MatchedPoints();						// destructor;

		MatchedPoints(const MatchedPoints &obj);	//copy constuctor
		MatchedPoints& operator=(const MatchedPoints & obj); //= assignment opperator
		
		//Standard get and set commands for data values in the class
		bool setMatchedPoints(PointSetTypePointer * inputmatchedpoints1, PointSetTypePointer * inputmatchedpoints2);
		bool setMatchedPoints(KeypointsIOType * inputmatchedpoints1, KeypointsIOType * inputmatchedpoints2);
		bool getMatchedPoints(PointSetTypePointer * inputmatchedpoints1, PointSetTypePointer * inputmatchedpoints2) const;
		bool getMatchedPoints(KeypointsIOType * inputmatchedpoints1, KeypointsIOType * inputmatchedpoints2) const;

		//Read and write matched key point sets to/from a text file
		bool readMatchedPoints(std::string matchedpoints1filename, std::string matchedpoints2filename);
		bool writeMatchedPoints(std::string matchedpoints1filename, std::string matchedpoints2filename);

		//Export matched key point set data into an Amira and MatLab friendly format
		bool writeAmiraMatchedPoints(std::string amiraLandmarksASCIIFileName);
		bool writeMatlabMatchedPoints(std::string matlabScriptFileName);

		//Print out key points and data onto screen
		void printMatchedPoints();
		void printMatchedPointsWithDistance();

		//Calculate and return data on distances between matched key point sets
		float * getDistance();
		int getDistanceIndex();

		//Print out point distance data to screen and to file
		void printDistance();
		bool writeDistance(std::string distancefilename);
		
	private:
		
		//Private data classes
		KeypointsIOType matchedpoints1;
		KeypointsIOType matchedpoints2;
  
};

#include "MatchedPoints.txx"

#endif