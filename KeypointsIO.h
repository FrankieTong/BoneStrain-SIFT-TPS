/*
File:			KeypointsIO.h

Description:	Header file for KeypointsIO.hxx. Class that defines methods for reading and writing key point sets from ITK to and from text files. 
				Also provides functionality for exporting point data to MatLab and Amira.

Author:			Frankie (Hoi-Ki) Tong, Jun 22 2015

Changes:		Frankie (Hoi-Ki) Tong, Jun 22 2015
				Initial creation of the file.
*/


#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include "itkPointSet.h"
#include "itkVector.h"
#include <limits>

#ifndef __KeypointsIO_H__
#define __KeypointsIO_H__

template <class FeatureType, unsigned int ImageDimensionality> 
class KeypointsIO {
	public:

		//Typedefs to make ITK point sets easier to manipulate
		typedef itk::PointSet< FeatureType, ImageDimensionality,
		itk::DefaultStaticMeshTraits< FeatureType, ImageDimensionality, ImageDimensionality, double > > PointSetType;
					  
		typedef typename PointSetType::PointType PointType;
		typedef typename PointSetType::Pointer PointSetTypePointer;
		typedef typename PointSetType::Pointer * PointSetTypePointerPtr;
		
		//Standard constructors and descructors
		KeypointsIO();                         // constructor; initialize the list to be empty
		KeypointsIO(PointSetTypePointer * inputkeypointSetPtr, std::string inputfileName); // constructor; initialize with preset values
		~KeypointsIO();						// destructor;

		KeypointsIO(const KeypointsIO &obj);	//copy constuctor
		KeypointsIO & operator=(const KeypointsIO & obj ); //= assignment opperator

		//Standard get and set commands for data values in the class
		bool setKeypointSetPtr(PointSetTypePointer * inputkeypointSetPtr);
		PointSetTypePointerPtr getKeypointSetPtr() const;

		bool setFileName(std::string inputfileName);
		std::string getFileName() const;

		//Print out key points onto screen
		void printKeypointSet();

		//Read and write key point sets to/from a text file
		bool readKeypointSet();
		bool writeKeypointSet();

		//Export key point set data into an Amira and MatLab friendly format
		bool writeAmiraKeypointSet();
		bool writeMatlabKeypointSet(std::string variableName);

	private:

		//Private data classes
		PointSetTypePointer * keypointSetPtr;	//Pointer to store ITK style key point sets
		std::string fileName;					//File name for importing/export key point set data
  
};

#include "KeypointsIO.txx"

#endif