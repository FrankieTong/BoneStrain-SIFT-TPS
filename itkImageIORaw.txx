/*
File:			itkImageIORaw.txx

Description:	Class that defines methods for reading and writing ITK images to/from a raw data format and a easily readable image configuration file.

				Templates in C++ neccessitates the need to include the function defitions inside its header file.

Author:			Frankie (Hoi-Ki) Tong, Oct 4 2014

Changes:		Frankie (Hoi-Ki) Tong, Oct 4 2014
				Initial creation of the file.
*/
#include "itkImageIORaw.h"

//Description:	Basic constructor with default values set to empty 
template <class PixelType, unsigned int ImageDimensionality>
itkImageIORaw<PixelType, ImageDimensionality>::itkImageIORaw(): 
	ImageFileName(""), itkImage(NULL){
}

//Description:	Basic constructor with default values set to parameters defined on construction by the user 
template <class PixelType, unsigned int ImageDimensionality>
itkImageIORaw<PixelType, ImageDimensionality>::itkImageIORaw(std::string inputImageFileName, std::string inputConfigFileName, void * inputitkImage): 
	ImageFileName(inputImageFileName){
	
	this->setConfigFileName(inputConfigFileName);
	this->setitkImage(inputitkImage);
	
}

//Description: Basic deconstructor
template <class PixelType, unsigned int ImageDimensionality>
itkImageIORaw<PixelType, ImageDimensionality>::~itkImageIORaw(){
}

//Description: Copy constructor
template <class PixelType, unsigned int ImageDimensionality>
itkImageIORaw<PixelType, ImageDimensionality>::itkImageIORaw(const itkImageIORaw &obj) {
	ImageFileName = obj.getImageFileName();
	itkImage = obj.getitkImage();
}

//Description: Assignment constructor
template <class PixelType, unsigned int ImageDimensionality>
itkImageIORaw<PixelType, ImageDimensionality>& itkImageIORaw<PixelType, ImageDimensionality>::operator=(const itkImageIORaw & obj){
	ImageFileName = obj.getImageFileName();
	itkImage = obj.getitkImage();

	return *this;
}

//Description: 	Get command for the ImageFileName private data set	
//Return:		(std::string) - String with configuration file and path name
template <class PixelType, unsigned int ImageDimensionality>
std::string itkImageIORaw<PixelType, ImageDimensionality>::getImageFileName() const { return ImageFileName; }

//Description: 	Set command for the inputImageFileName private data set
//Input:		inputImageFileName - (std::string) String with configuration file and path name to image file
//Return:		(bool) 0
template <class PixelType, unsigned int ImageDimensionality>
bool itkImageIORaw<PixelType, ImageDimensionality> ::setImageFileName(std::string inputImageFileName){
	ImageFileName = inputImageFileName;
	return 0;
}

//Description: 	Get command for the inputitkImage private data set. 	
//Input:		inputitkImage - (void *) Pointer to the place where image is to be stored in when returning to main function
//Return:		(bool) 0 if successful
//				(bool) 1 if ImageIOData and/or itkImage is not defined
template <class PixelType, unsigned int ImageDimensionality>
bool itkImageIORaw<PixelType, ImageDimensionality>::getitkImage(void * inputitkImage){
	//check if image meta data is found
	if (this->getImageIOData() == NULL){ inputitkImage = NULL; return 1 };
	
	//see if itkImage is NULL and return it if it is
	if (itkImage == NULL) { inputitkImage = NULL; return 1 };
	
	//give direct copy of pointer to inputitkImage (not a smart copy to keep built in data management from itk)
	typedef itk::Image<PixelType, ImageDimensionality>		ImageType;
	inputitkImage = itkImage;
	
	return 0;
}

//Description: 	Set command for the inputitkImage private data set. 	
//Return:		(bool) 0 if successful
template <class PixelType, unsigned int ImageDimensionality>
bool itkImageIORaw<PixelType, ImageDimensionality>::setitkImage(void * inputitkImage) {
	//set itkImage pointer as inputitkImage pointer (not a smart copy to keep built in data management from itk)
	typedef itk::Image<PixelType, ImageDimensionality>		ImageType;
	itkImage = inputitkImage;
	
	return 0;
}


//Description: 	Attempts to read the raw data image with the path name defined in ImageFileName and the parameters set in ImageIOData
//Return:		(bool) 0 if successful
//				(bool) 1 if there was an error reading in the image
template <class PixelType, unsigned int ImageDimensionality>
bool itkImageIORaw<PixelType, ImageDimensionality>::readRawImageKnown(){
	typedef itk::Image<PixelType, ImageDimensionality>		ImageType;
	typedef itk::RawImageIO<PixelType, ImageDimensionality> ImageIOType;
	
	//Create new raw image data format
	ImageIOType::Pointer imageio = ImageIOType::New();
	
	//Set up raw data format parameters
	imageio->SetFileDimensionality(this->getImageIOData().getImageDimensionality());
	for (int i = 0; i<this->getImageIOData().getImageDimensionality(); i++) {
		imageio->SetDimensions(i,(int)(this->getImageIOData().getDimensionSize())[i]);
		imageio->SetSpacing(i,(float)(this->getImageIOData().getSpacingSize())[i]);
		imageio->SetOrigin(i,(float)(this->getImageIOData().getOrigin())[i]);
	}
	imageio->SetHeaderSize(this->getImageIOData().getHeaderSize());
	
	if (this->getImageIOData().getByteOrder() == "little") {
		imageio->SetByteOrderToLittleEndian();
	} else {
		imageio->SetByteOrderToBigEndian();
	}
	
	//Attempt to read in raw image
	typedef itk::ImageFileReader<ImageType> FileReaderType;
	FileReaderType::Pointer reader = FileReaderType::New();
	reader->SetFileName(this->getImageFileName());
	reader->SetImageIO	(imageio);

	try
	{
		reader->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		//Error reading in raw image data type
		std::cerr << "Exception thrown while reading the image" << std::endl;
		std::cerr << excp << std::endl;
		return 1;
	}
	
	//Store data from image reader in itkImage
	ImageType::Pointer * readImage = (ImageType::Pointer *) itkImage;
	*readImage = reader->GetOutput();
	
	return 0;
}

//Description: 	Attempts to write the raw data image with the path name defined in ImageFileName with the parameters set in ImageIOData
//Return:		(bool) 0 if successful
//				(bool) 1 if there was an error reading in the image
template <class PixelType, unsigned int ImageDimensionality>
bool itkImageIORaw<PixelType, ImageDimensionality>::writeRawImageKnown(){
	typedef itk::Image<PixelType, ImageDimensionality>		ImageType;
	typedef itk::RawImageIO<PixelType, ImageDimensionality> ImageIOType;
	
	//Create new raw image data format
	ImageIOType::Pointer imageio = ImageIOType::New();
	
	//Set up raw data format parameters
	imageio->SetFileDimensionality(this->getImageIOData().getImageDimensionality());
	for (int i = 0; i<this->getImageIOData().getImageDimensionality(); i++) {
		imageio->SetDimensions(i,(int)(this->getImageIOData().getDimensionSize())[i]);
		imageio->SetSpacing(i,(float)(this->getImageIOData().getSpacingSize())[i]);
		imageio->SetOrigin(i,(float)(this->getImageIOData().getOrigin())[i]);
	}
	imageio->SetHeaderSize(this->getImageIOData().getHeaderSize());
	
	if (this->getImageIOData().getByteOrder() == "little") {
		imageio->SetByteOrderToLittleEndian();
	} else {
		imageio->SetByteOrderToBigEndian();
	}
	
	//Attempt to write raw image data
	typedef itk::ImageFileWriter<ImageType> FileWriterType;
	FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetFileName(this->getImageFileName());
	writer->SetImageIO( imageio );
	ImageType::Pointer *writeImage = (ImageType::Pointer *) itkImage;
	writer->SetInput(*writeImage);
	try
	{
		writer->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		//Error in writing raw image file
		std::cerr << "Exception thrown while writing the image" << std::endl;
		std::cerr << excp << std::endl;
		return 1;
	}
	
	return 0;
}

//Description: 	Failed function prototype to make a generalized image reader for the raw image data and configuration file input/output format.
//				More specifically, pixel type and iamge dimensions are defined on compile and can not be modified on execution time which makes this idea not work well.
//Return:		(bool) 1 for error
template <class PixelType, unsigned int ImageDimensionality>
bool itkImageIORaw<PixelType, ImageDimensionality>::readRawImage(){
	std::cout << "readRawImage currently not programmed. Please use readRawImageKnown()" << std::endl;
	return 1;

	/* Will follow the line of something like this (used to be itkImageIORaw::readRawImageKnown
	typedef itk::Image<PixelType, ImageDimensionality>		ImageType;
	typedef itk::RawImageIO<PixelType, ImageDimensionality> ImageIOType;
	
	ImageIOType::Pointer imageio = ImageIOType::New();
	
	imageio->SetFileDimensionality(this->getImageIOData().getImageDimensionality());
	for (int i = 0; i<this->getImageIOData().getImageDimensionality(); i++) {
		imageio->SetDimensions(i,(int)(this->getImageIOData().getDimensionSize())[i]);
		imageio->SetSpacing(i,(float)(this->getImageIOData().getSpacingSize())[i]);
		imageio->SetOrigin(i,(float)(this->getImageIOData().getOrigin())[i]);
	}
	imageio->SetHeaderSize(this->getImageIOData().getHeaderSize());
	
	if (this->getImageIOData().getByteOrder() == "little") {
		imageio->SetByteOrderToLittleEndian();
	} else {
		imageio->SetByteOrderToBigEndian();
	}
	
	typedef itk::ImageFileWriter<ImageType> FileWriterType;
	FileWriterType::Pointer writer = FileWriterType::New();
	writer->SetFileName(this->getImageFileName());
	writer->SetImageIO( imageio );
	ImageType::Pointer *writeImage = (ImageType::Pointer *) itkImage;
	writer->SetInput(*writeImage);
	try
	{
		writer->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << "Exception thrown while writing the image" << std::endl;
		std::cerr << excp << std::endl;
		return 1;
	}
	*/
}

//Description: 	Failed function prototype to make a generalized image writer for the raw image data and configuration file input/output format.
//				More specifically, pixel type and iamge dimensions are defined on compile and can not be modified on execution time which makes this idea not work well.
//Return:		(bool) 1 for error
template <class PixelType, unsigned int ImageDimensionality>
bool itkImageIORaw<PixelType, ImageDimensionality>::writeRawImage(){
	std::cout << "writeRawImage currently not programmed. Please use writeRawImageKnown()" << std::endl;
	return 1;
}
