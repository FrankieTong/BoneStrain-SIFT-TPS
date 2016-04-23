/*
File:			ImageIOData.h

Description:	Header file for ImageIOData.cxx. Basic class for holding ITK image data for reading and writing to/from raw image data files.

				Consists of the constructor, desctructor, copy constructor and the get and set functions for each private data class member.
				
				DOES NOT DO ANY VERIFICATION OF CONTENTS. THIS IS JUST A BASIC CLASS TO HOLD DATA.

Author:			Frankie (Hoi-Ki) Tong, Sep 6 2014

Changes:		Frankie (Hoi-Ki) Tong, Sep 6 2014
				Initial creation of the file.
*/

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

#include <limits>

#ifndef __ImageIOData_H__
#define __ImageIOData_H__

class ImageIOData {
  public:
  
	//Standard constructors and descructors
	ImageIOData();                      // constructor; initialize the list to be empty
	~ImageIOData();						// destructor;

	ImageIOData(const ImageIOData &obj);	//copy constuctor
	ImageIOData & operator=(const ImageIOData & obj ); //= assignment opperator
	
	
	//Get and set commands for each private data member of this class
	std::string getPixelType() const;
	bool setPixelType(std::string inputPixelType);
	
	int getPixelDimensionality() const;
	bool setPixelDimensionality(int inputPixelDimensionality);
	
	int getImageDimensionality() const;
	bool setImageDimensionality(int inputImageDimensionality);
	
	int * getDimensionSize() const;
	bool setDimensionSize(int * inputDimensionSize);
	
	float * getSpacingSize() const;
	bool setSpacingSize(float * inputSpacingSize);
	
	float * getOrigin() const;
	bool setOrigin(float * inputOrigin);
	
	int getHeaderSize() const;
	bool setHeaderSize(int inputHeaderSize);
	
	std::string getByteOrder() const;
	bool setByteOrder(std::string inputByteOrder);
	
	
	//Print out contents of the member class to screen easily using this function.
	void printImageIOData();

  private:
  
	//Private data classes
	std::string PixelType;		//pixel type of the image
	int PixelDimensionality;	//pixel dimensionality (ie: scalor = 1, vector = 2 or above)
	int ImageDimensionality;	//image dimensionality
	int * DimensionSize;		//image size
	float * SpacingSize;		//image pixel spacing
	float * Origin;				//image origin location
	int HeaderSize;				//size of the data header before first byte of image data is found (usually set to 0)
	std::string ByteOrder;		//byte order of how the image data is stored. Amira likes to use "little" endian.
  
};

#endif