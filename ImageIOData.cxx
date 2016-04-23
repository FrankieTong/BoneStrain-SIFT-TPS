/*
File:			ImageIOData.hxx

Description:	Basic class for holding ITK image data for reading and writing to/from raw image data files.

				Consists of the constructor, desctructor, copy constructor and the get and set functions for each private data class member.
				
				DOES NOT DO ANY VERIFICATION OF CONTENTS. THIS IS JUST A BASIC CLASS TO HOLD DATA.

Author:			Frankie (Hoi-Ki) Tong, Sep 6 2014

Changes:		Frankie (Hoi-Ki) Tong, Sep 6 2014
				Initial creation of the file.
*/

#include "ImageIoData.h"

//Description:	Basic constructor with default values set to either "0" or "NULL"
ImageIOData::ImageIOData(): 
	PixelType(""), PixelDimensionality(0), ImageDimensionality(0), DimensionSize(NULL), SpacingSize(NULL), 
	Origin(NULL), HeaderSize(0), ByteOrder("") {
}

//Description: 	Basic deconstructor which will attempt to clear any allocated space for pointers in this data structure	
ImageIOData::~ImageIOData(){
	if (DimensionSize != NULL) { delete[] DimensionSize; }
	if (SpacingSize != NULL) { delete[] SpacingSize; }
	if (Origin != NULL) { delete[] Origin; }
}

//Description: 	Copy constructor utilizing the get commands defined below
ImageIOData::ImageIOData(const ImageIOData &obj){
	PixelType = obj.getPixelType();
	PixelDimensionality = obj.getPixelDimensionality();
	ImageDimensionality = obj.getImageDimensionality();
	DimensionSize = obj.getDimensionSize();
	SpacingSize = obj.getSpacingSize();
	Origin = obj.getOrigin();
	HeaderSize = obj.getHeaderSize();
	ByteOrder = obj.getByteOrder();
}

//Description: 	Assignment constructor utilizing the get commands defined below
ImageIOData & ImageIOData::operator=(const ImageIOData & obj ) {
	PixelType = obj.getPixelType();
	PixelDimensionality = obj.getPixelDimensionality();
	ImageDimensionality = obj.getImageDimensionality();
	DimensionSize = obj.getDimensionSize();
	SpacingSize = obj.getSpacingSize();
	Origin = obj.getOrigin();
	HeaderSize = obj.getHeaderSize();
	ByteOrder = obj.getByteOrder();

	return *this;
}


//Description: 	Get command for the PixelType private data set	
//Return:		(std::string) - String with pixel type
std::string ImageIOData::getPixelType() const { 
	return PixelType; 
}

//Description: 	Set command for the PixelType private data set
//Input:		inputPixelType - (std::string) String containing pixel type to set data member as
//Return:		(bool) 0
bool ImageIOData::setPixelType(std::string inputPixelType){
	PixelType = inputPixelType;
	return 0;
}


//Description: 	Get command for the PixelDimensionality private data set
//Return:		(std::string) - String with pixel dimensionality	
int ImageIOData::getPixelDimensionality() const { 
	return PixelDimensionality; 
}

//Description: 	Set command for the PixelDimensionality private data set
//Input:		inputPixelDimensionality - (int) Pixel dimensionality to set data member as
//Return:		(bool) 0
bool ImageIOData::setPixelDimensionality(int inputPixelDimensionality){
	PixelDimensionality = inputPixelDimensionality;
	return 0;
}


//Description: 	Get command for the ImageDimensionality private data set		
//Return:		(int) Image dimensionlity value
int ImageIOData::getImageDimensionality() const { 
	return ImageDimensionality; 
}

//Description: 	Set command for the ImageDimensionality private data set
//Input:		inputImageDimensionality - (int) Image dimensionality to set data member as
//Return:		(bool) 0
bool ImageIOData::setImageDimensionality(int inputImageDimensionality) {
	ImageDimensionality = inputImageDimensionality;
	return 0;
}


//Description: 	Get command for the DimensionSize private data set.	
//Return:		(int *) Copy of DimensionSize[ImageDimensionality] if allocated
//				(int *) NULL otherwise
int * ImageIOData::getDimensionSize() const {
	//Check if ImageDimensionality is assigned and DimensionSize[] is allocated and return copy of DimensionSize[]
	if (ImageDimensionality > 0 && DimensionSize != NULL) {
		int * cpyPtr = new int[ImageDimensionality];
		for (int i = 0; i < ImageDimensionality; i++) {
			cpyPtr[i] = DimensionSize[i];
		}
		return cpyPtr;
	//Failed check, return NULL pointer
	} else { return NULL; }
}

//Description: 	Set command for the DimensionSize private data set.
//Input:		inputDimensionSize - (int *) Pointer to DimensionSize[] array to set data member as
//Return:		(bool) 0 if successful
//				(bool) 1 otherwise
bool ImageIOData::setDimensionSize(int * inputDimensionSize){
	//Check if DimensionSize can be set properly. Return 1 if not possible
	if (ImageDimensionality <= 0 || inputDimensionSize == NULL) { 
		return 1; 
	}
	//Delete current DimensionSize if not empty
	if (DimensionSize != NULL) {
		delete[] DimensionSize;
		DimensionSize = NULL;
	}
	//Set new DimensionSize values
	DimensionSize = new int[ImageDimensionality];
	for (int i = 0; i < ImageDimensionality; i++) {
		DimensionSize[i] = inputDimensionSize[i];
	}
	return 0;
}


//Description: 	Get command for the SpacingSize private data set.	
//Return:		(int *) Copy of SpacingSize[ImageDimensionality] if allocated
//				(int *) NULL otherwise
float * ImageIOData::getSpacingSize() const { 	
	//Check if ImageDimensionality is assigned and SpacingSize[] is allocated and return copy of SpacingSize[]
	if (ImageDimensionality > 0 && SpacingSize != NULL) {
		float * cpyPtr = new float[ImageDimensionality];
		for (int i = 0; i < ImageDimensionality; i++) {
			cpyPtr[i] = SpacingSize[i];
		}
		return cpyPtr;
	//Failed check, return NULL pointer
	} else { return NULL; }
}

//Description: 	Set command for the SpacingSize private data set.
//Input:		inputSpacingSize - (float *) Pointer to SpacingSize[] array to set data member as	
//Return:		(bool) 0 if successful
//				(bool) 1 otherwise
bool ImageIOData::setSpacingSize(float * inputSpacingSize){
	//Check if SpacingSize can be set properly. Return 1 if not possible
	if (ImageDimensionality <= 0 || inputSpacingSize == NULL) {
		return 1; 
	}
	//Delete current SpacingSize if not empty
	if (SpacingSize != NULL) {
		delete[] SpacingSize;
		SpacingSize = NULL;
	}
	//Set new SpacingSize values
	SpacingSize = new float[ImageDimensionality];
	for (int i = 0; i < ImageDimensionality; i++) {
		SpacingSize[i] = inputSpacingSize[i];
	}
	return 0;
}


//Description: 	Get command for the Origin private data set.	
//Return:		(int *) Copy of Origin[ImageDimensionality] if allocated
//				(int *) NULL otherwise
float * ImageIOData::getOrigin() const { 	
	//Check if ImageDimensionality is assigned and Origin[] is allocated and return copy of SpacingSize[]
	if (ImageDimensionality > 0 && Origin != NULL) {
		float * cpyPtr = new float[ImageDimensionality];
		for (int i = 0; i < ImageDimensionality; i++) {
			cpyPtr[i] = Origin[i];
		}
		return cpyPtr;
	//Failed check, return NULL pointer
	} else { return NULL; }
}

//Description: 	Set command for the Origin private data set.
//Input:		inputOrigin - (float *) Pointer to Origin[] array to set data member as	
//Return:		(bool) 0 if successful
//				(bool) 1 otherwise
bool ImageIOData::setOrigin(float * inputOrigin){
	//Check if Origin can be set properly. Return 1 if not possible
	if (ImageDimensionality <= 0 || inputOrigin == NULL) { 
		return 1; 
	}
	//Delete current Origin if not empty
	if (Origin != NULL) {
		delete[] Origin;
		Origin = NULL;
	}
	//Set new Origin values
	Origin = new float[ImageDimensionality];
	for (int i = 0; i < ImageDimensionality; i++) {
		Origin[i] = inputOrigin[i];
	}
	return 0;
}


//Description: 	Get command for the HeaderSize private data set		
//Return:		(int) Image dimensionlity value
int ImageIOData::getHeaderSize() const { 
	return HeaderSize; 
}

//Description: 	Set command for the HeaderSize private data set
//Input:		inputHeaderSize - (int) Header size to set data member as
//Return:		(bool) 0
bool ImageIOData::setHeaderSize(int inputHeaderSize){
	HeaderSize = inputHeaderSize;
	return 0;
}


//Description: 	Get command for the ByteOrder private data set		
//Return:		(int) Image dimensionlity value
std::string ImageIOData::getByteOrder()  const { 
	return ByteOrder; 
}

//Description: 	Set command for the ByteOrder private data set
//Input:		inputByteOrder - (int) Header size to set data member as
//Return:		(bool) 0
bool ImageIOData::setByteOrder(std::string inputByteOrder){
	ByteOrder = inputByteOrder;
	return 0;
}

//Description: 	Print out contents of the member class to screen easily using this function.
void ImageIOData::printImageIOData() {
	std::cout<<"Printing ImageIOData:"<<std::endl;
	std::cout<<"PixelDimensionality: "<< PixelDimensionality <<std::endl;
	std::cout<<"ImageDimensionality: "<< ImageDimensionality <<std::endl;

	std::cout<<"DimensionSize: ";
	if (ImageDimensionality <=0) { std::cout<<"N/A"; }
	else { 
		for(int i=0; i<ImageDimensionality; i++) {
			std::cout<<DimensionSize[i];
			if (i<ImageDimensionality - 1) {std::cout<<" ";};
		}
	}
	std::cout<<std::endl;

	std::cout<<"SpacingSize: ";
	if (ImageDimensionality <=0) { std::cout<<"N/A"; }
	else { 
		for(int i=0; i<ImageDimensionality; i++) {
			std::cout<<SpacingSize[i];
			if (i<ImageDimensionality - 1) {std::cout<<" ";};
		}
	}
	std::cout<<std::endl;

	std::cout<<"Origin: ";
	if (ImageDimensionality <=0) { std::cout<<"N/A"; }
	else { 
		for(int i=0; i<ImageDimensionality; i++) {
			std::cout<<Origin[i];
			if (i<ImageDimensionality - 1) {std::cout<<" ";};
		}
	}
	std::cout<<std::endl;

	std::cout<<"HeaderSize: "<< HeaderSize <<std::endl;
	std::cout<<"ByteOrder: "<< ByteOrder <<std::endl;
}