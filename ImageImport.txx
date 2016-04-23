/*
File:			ImageImport.txx

Description:	Import and export of ITK image data in a raw data format

Author:			Frankie (Hoi-Ki) Tong, Sep 6 2014

Changes:		Frankie (Hoi-Ki) Tong, Sep 6 2014
				Initial creation of the file.
				
				Frankei (Hoi-Ki) Tong, Oct 4 2014
				Depreciated in exchange for ImageIOConfig structure due to awkard method to store configuration parameters
*/

//#include "ImageImport.h"

//Description: 	Hash parameter to be used to determine if the input data type used in PixelType is correct
//Input:		inString - (std::string) input string to be used for comparison with all allowable pixel types
//Return:		eUnknownType::(valid PixelType) if inString matches one of the defined PixelTypes
//				eUnknownType::eUnknownType if inString does not match any defined PixelTypes	
EnumImageType hashitType (std::string inString) {
	if (inString == "char") {return eChar;}
	else if (inString == "unsigned char") {return eUnsignedChar;}
	else if (inString == "int") {return eInt;}
	else if (inString == "unsigned int") {return eUnsignedInt;}
	else if (inString == "short") {return eShort;}
	else if (inString == "unsigned short") {return eUnsignedShort;}
	else if (inString == "float") {return eFloat;}
	else if (inString == "double") {return eDouble;}
	else {return eUnknownType;}
}

//Description: 	Attempts to read the raw data image with the path name defined in InputFile and the parameters set in ConfigFileParam and ConfigFileContent
//
//				NOTE: Does not work properly with float converstions. Abbandonned in favor of using ImageImportReadImage
//
//Input:		ConfigFileParam - (std::list<std::string>) String list where image configuration parameter names are stored
//				ConfigFileContent - (std::list<std::string>) String list where image configuration parameter contents are stored
//				InputFile - (std::string) input string to be used for comparison with all allowable pixel types
//				OutputImage - (void *) Pointer holding ITK image object where image read is to be stored
template <class OutputImageType, class OutputPixelType> 
void ImageImportScalor(std::list<std::string> * ConfigFileParam,  std::list<std::string> * ConfigFileContent,  std::string InputFile, void * OutputImage) {

	//Setup input filter to properly read in the data from the input file
	
	//typedef float							PixelType;
	//typedef itk::Image<PixelType, 3>		ImageType;

	//First, find the input pixel diemsion to determine pixel type.
	std::string ConfigContentString;
	std::istringstream stringConvert;

	//Find pixel dimensionsionality
	ConfigContentString = ReturnConfigIndex (FindParamIndex("PixelDimensionality",ConfigFileParam),ConfigFileContent);
	stringConvert.str(ConfigContentString);
	unsigned int PixelDimensionality;
	stringConvert >> PixelDimensionality;
	stringConvert.str("");
	stringConvert.clear();

	//Find image dimensions
	ConfigContentString = ReturnConfigIndex (FindParamIndex("ImageDimensionality",ConfigFileParam),ConfigFileContent);
	stringConvert.str(ConfigContentString);
	unsigned int ImageDimensionality;
	stringConvert >> ImageDimensionality;
	stringConvert.str("");
	stringConvert.clear();

	//Find pixel type and create RawImageIO object and ImageFileReader object
	ConfigContentString = ReturnConfigIndex (FindParamIndex("PixelType",ConfigFileParam),ConfigFileContent);
	std::string PixelType = ConfigContentString;
	stringConvert.str("");
	stringConvert.clear();

	switch(hashitType(PixelType)) {
		case eFloat:
			if (PixelDimensionality == 1 && ImageDimensionality == 2){
				itk::Image<float,2>::Pointer ImageIn = itk::Image<float,2>::New();
				void * readerImage = &ImageIn;
				ImageImportReadImage<float, 2>(ConfigFileParam, ConfigFileContent, InputFile, readerImage);

				//Cast it to float or whatever we need then pass that back to (void *) imageio
				typedef itk::CastImageFilter< itk::Image<float, 2>, itk::Image<OutputPixelType,2> > CastFilterType;
				CastFilterType::Pointer castFilter = CastFilterType::New();
				castFilter->SetInput(ImageIn);

				itk::Image<OutputPixelType,2>::Pointer * convertedImage = (itk::Image<OutputPixelType,2>::Pointer *) OutputImage;
				*convertedImage = castFilter->GetOutput();
				break;

			} else if (PixelDimensionality == 1 && ImageDimensionality == 3) {
				itk::Image<float,3>::Pointer ImageIn = itk::Image<float,3>::New();
				void * readerImage = &ImageIn;
				ImageImportReadImage<float, 3>(ConfigFileParam, ConfigFileContent, InputFile, readerImage);

				//Cast it to float or whatever we need then pass that back to (void *) imageio
				typedef itk::CastImageFilter< itk::Image<float, 3>, itk::Image<OutputPixelType,3> > CastFilterType;
				CastFilterType::Pointer castFilter = CastFilterType::New();
				castFilter->SetInput(ImageIn);

				itk::Image<OutputPixelType,3>::Pointer * convertedImage = (itk::Image<OutputPixelType,3>::Pointer *) OutputImage;
				*convertedImage = castFilter->GetOutput();
				//*convertedImage = ImageIn;
				break;

			} else {
				//error...
				break;
			}
			break;
		case eChar:
			if (PixelDimensionality == 1 && ImageDimensionality == 2){
				//(itk::RawImageIO<char, 2>::Pointer) imageio = itk::RawImageIO<char, 2>::New();
				//itk::ImageFileReader<itk::Image<char, 2>>::Pointer reader = itk::ImageFileReader<itk::Image<char, 2>>::New();
			} else if (PixelDimensionality == 1 && ImageDimensionality == 3) {
				//(itk::RawImageIO<char, 3>::Pointer) imageio = itk::RawImageIO<char, 3>::New();
				//itk::ImageFileReader<itk::Image<char, 3>>::Pointer reader = itk::ImageFileReader<itk::Image<char, 3>>::New();
			} else {
				//error...
				break;
			}
			break;
		case eUnsignedChar:
			if (PixelDimensionality == 1 && ImageDimensionality == 2){
				//(itk::RawImageIO<char, 2>::Pointer) imageio = itk::RawImageIO<unsigned char, 2>::New();
				//itk::ImageFileReader<itk::Image<char, 2>>::Pointer reader = itk::ImageFileReader<itk::Image<char, 2>>::New();
			} else if (PixelDimensionality == 1 && ImageDimensionality == 3) {
				//(itk::RawImageIO<char, 3>::Pointer) imageio = itk::RawImageIO<char, 3>::New();
				//itk::ImageFileReader<itk::Image<char, 3>>::Pointer reader = itk::ImageFileReader<itk::Image<char, 3>>::New();
			} else {
				//error...
				break;
			}
			break;
		default:
			break;
	}

	//Messing around
	//imageio->SetFileDimensionality(ImageDimensionality);

	//std::list<std::string>::iterator it = std::find(whichList.begin(), whichlist.end(), "");
	

	//OutputImageType

	return;
}

//Description: 	Attempts to read the raw data image with the path name defined in InputFile and the parameters set in ConfigFileParam and ConfigFileContent
//				Image pixel type and dimensions are defined by the template parameters and can not be changed using configuration file settings
//Input:		ConfigFileParam - (std::list<std::string>) String list where image configuration parameter names are stored
//				ConfigFileContent - (std::list<std::string>) String list where image configuration parameter contents are stored
//				InputFile - (std::string) input string to be used for comparison with all allowable pixel types
//				OutputImage - (void *) Pointer holding ITK image object where image read is to be stored
template <class PixelType, unsigned int ImageDimensionality> 
void ImageImportReadImage (std::list<std::string> * ConfigFileParam,  std::list<std::string> * ConfigFileContent, std::string InputFile, void * Image) {

	typedef itk::Image<PixelType, ImageDimensionality>		ImageType;
	typedef itk::RawImageIO<PixelType, ImageDimensionality> ImageIOType;

	//Set up raw data image type for image I/O
	ImageIOType::Pointer imageio = ImageIOType::New();

	std::string  ConfigContentString;
	std::istringstream stringConvert;

	//Find and set image dimensionsionality
	ConfigContentString = ReturnConfigIndex (FindParamIndex("ImageDimensionality",ConfigFileParam),ConfigFileContent);
	stringConvert.str(ConfigContentString);
	unsigned int ImageDimensionality;
	stringConvert >> ImageDimensionality;
	imageio->SetFileDimensionality(ImageDimensionality);
	stringConvert.str("");
	stringConvert.clear();

	//Find and set dimension size
	ConfigContentString = ReturnConfigIndex (FindParamIndex("DimensionSize",ConfigFileParam),ConfigFileContent);
	stringConvert.str(ConfigContentString);
	int DimensionSize = 0;
	for (int i = 0; i < ImageDimensionality; i++ ) {
		stringConvert >> DimensionSize;
		imageio->SetDimensions(i,DimensionSize);
	}
	stringConvert.str("");
	stringConvert.clear();

	//Find and set spacing size
	ConfigContentString = ReturnConfigIndex (FindParamIndex("SpacingSize",ConfigFileParam),ConfigFileContent);
	stringConvert.str(ConfigContentString);
	float SpacingSize = 0;
	for (int i = 0; i < ImageDimensionality; i++ ) {
		stringConvert >> SpacingSize;
		imageio->SetSpacing(i,SpacingSize);
	}
	stringConvert.str("");
	stringConvert.clear();

	//Find and set origin
	ConfigContentString = ReturnConfigIndex (FindParamIndex("Origin",ConfigFileParam),ConfigFileContent);
	stringConvert.str(ConfigContentString);
	float Origin = 0;
	for (int i = 0; i < ImageDimensionality; i++ ) {
		stringConvert >> Origin;
		imageio->SetOrigin(i,Origin);
	}
	stringConvert.str("");
	stringConvert.clear();

	//Find and set header size
	ConfigContentString = ReturnConfigIndex (FindParamIndex("HeaderSize",ConfigFileParam),ConfigFileContent);
	stringConvert.str(ConfigContentString);
	unsigned int HeaderSize;
	stringConvert >> HeaderSize;
	imageio->SetHeaderSize(HeaderSize);
	stringConvert.str("");
	stringConvert.clear();

	//Find and set byte order
	ConfigContentString = ReturnConfigIndex (FindParamIndex("ByteOrder",ConfigFileParam),ConfigFileContent);
	if (ConfigContentString == "little") {
		imageio->SetByteOrderToLittleEndian();
	} else {
		imageio->SetByteOrderToBigEndian();
	}

	//Try to read in raw data file as an image
	typedef itk::ImageFileReader<ImageType> FileReaderType;
	FileReaderType::Pointer reader = FileReaderType::New();
	reader->SetFileName(InputFile);
	reader->SetImageIO( imageio );

	try
	{
		reader->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << "Exception thrown while reading the image" << std::endl;
		std::cerr << excp << std::endl;
	}
	
	//Store image in pointer before return
	ImageType::Pointer *readImage = (ImageType::Pointer *) Image;
	*readImage = reader->GetOutput();

	return;
}

//Description: 	Attempts to write the ITK image as a raw data file to the path name defined in InputFile using the parameters set in ConfigFileParam and ConfigFileContent
//				Image pixel type and dimensions are defined by the template parameters and can not be changed using configuration file settings
//Input:		ConfigFileParam - (std::list<std::string>) String list where image configuration parameter names are stored
//				ConfigFileContent - (std::list<std::string>) String list where image configuration parameter contents are stored
//				OutputFile - (std::string) input string to be used for comparison with all allowable pixel types
//				Image - (void *) Pointer holding ITK image object where image read is to be stored
template <class PixelType, unsigned int ImageDimensionality> 
void ImageImportWriteImage (std::list<std::string> * ConfigFileParam,  std::list<std::string> * ConfigFileContent, std::string OutputFile, void * Image) {

	typedef itk::Image<PixelType, ImageDimensionality>		ImageType;
	typedef itk::RawImageIO<PixelType, ImageDimensionality> ImageIOType;

	//Set up raw data image type for image I/O
	ImageIOType::Pointer imageio = ImageIOType::New();

	std::string  ConfigContentString;
	std::istringstream stringConvert;

	//Find and set image dimensionsionality
	ConfigContentString = ReturnConfigIndex (FindParamIndex("ImageDimensionality",ConfigFileParam),ConfigFileContent);
	stringConvert.str(ConfigContentString);
	unsigned int ImageDimensionality;
	stringConvert >> ImageDimensionality;
	imageio->SetFileDimensionality(ImageDimensionality);
	stringConvert.str("");
	stringConvert.clear();

	//Find and set dimension size
	ConfigContentString = ReturnConfigIndex (FindParamIndex("DimensionSize",ConfigFileParam),ConfigFileContent);
	stringConvert.str(ConfigContentString);
	int DimensionSize = 0;
	for (int i = 0; i < ImageDimensionality; i++ ) {
		stringConvert >> DimensionSize;
		imageio->SetDimensions(i,DimensionSize);
	}
	stringConvert.str("");
	stringConvert.clear();

	//Find and set spacing size
	ConfigContentString = ReturnConfigIndex (FindParamIndex("SpacingSize",ConfigFileParam),ConfigFileContent);
	stringConvert.str(ConfigContentString);
	float SpacingSize = 0;
	for (int i = 0; i < ImageDimensionality; i++ ) {
		stringConvert >> SpacingSize;
		imageio->SetSpacing(i,SpacingSize);
	}
	stringConvert.str("");
	stringConvert.clear();

	//Find and set origin
	ConfigContentString = ReturnConfigIndex (FindParamIndex("Origin",ConfigFileParam),ConfigFileContent);
	stringConvert.str(ConfigContentString);
	float Origin = 0;
	for (int i = 0; i < ImageDimensionality; i++ ) {
		stringConvert >> Origin;
		imageio->SetOrigin(i,Origin);
	}
	stringConvert.str("");
	stringConvert.clear();

	//Find and set header size
	ConfigContentString = ReturnConfigIndex (FindParamIndex("HeaderSize",ConfigFileParam),ConfigFileContent);
	stringConvert.str(ConfigContentString);
	unsigned int HeaderSize;
	stringConvert >> HeaderSize;
	imageio->SetHeaderSize(HeaderSize);
	stringConvert.str("");
	stringConvert.clear();

	//Find and set byte order
	ConfigContentString = ReturnConfigIndex (FindParamIndex("ByteOrder",ConfigFileParam),ConfigFileContent);
	if (ConfigContentString == "little") {
		imageio->SetByteOrderToLittleEndian();
	} else {
		imageio->SetByteOrderToBigEndian();
	}

	//Write ITK image in a raw data format to file
	typedef itk::ImageFileWriter<ImageType> FileWriterType;
	FileWriterType::Pointer writer = FileWriterType::New();

	writer->SetFileName(OutputFile);
	writer->SetImageIO( imageio );

	ImageType::Pointer *writeImage = (ImageType::Pointer *) Image;
	writer->SetInput(*writeImage);

	try
	{
		writer->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << "Exception thrown while writing the image" << std::endl;
		std::cerr << excp << std::endl;
	}

	return;
}

//Description: 	Finds the matching parameter string in the string list and returns the index value of the match
//Input:		Param - (std::string) String storing the parameter name to be searched
//				ConfigFileParam - (std::list<std::string>) String list where image configuration parameters are stored
//Return:		(int) Index number of the matching string found in the string list. Otherwise will return a negative number.
int FindParamIndex (std::string Param, std::list<std::string> * ConfigFileParam ) {
	std::list<std::string>::iterator itParam = std::find(ConfigFileParam->begin(),ConfigFileParam->end(), Param);
	int index = std::distance(ConfigFileParam->begin(),itParam);
	return index;
}

//Description: 	Returns the value in the string list given its index number in the list
//Input:		Index - (int) Index number of the string in the string list to be returned
//				ConfigFileContent - (std::list<std::string>) String list where the return values are stored
//Return:		(std::string) Returns the string located at the index value of the string list.
std::string ReturnConfigIndex (int Index, std::list<std::string> * ConfigFileContent) {
	std::list<std::string>::iterator itContent = ConfigFileContent->begin();
	std::advance(itContent,Index);
	return *itContent;
}
