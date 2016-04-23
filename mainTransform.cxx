/*
File:			mainTransform.cxx

Description:	Subfunction called from imageioraw.cxx applies a transform to an image given a transformation field represented as a displacement vector at every given point in the input image. Input parameters based on input processing from calling mainTransform as a command line function, ignoring the first 3 input parameters. Refer below for input process parameters.

Inputs:			argc - (int) Number of input arguments stored in argv. Default number is 8.
				argv[2]: inputImageFile1 - (std::string) Path name of the raw data file of the input image. Pixels are in 32 bit float format.
				argv[3]: inConfigFile1 - (std::string) Path name of the configuration text file of the input image. Pixels are in 32 bit float format.
				
				argv[4]: inputVectorField1 - (std::string) Path name of the raw data file of the displacement field. Pixels are in 32 bit float format.
				argv[5]: inVectorConfigFile1 - (std::string) Path name of the configuration text file of the displacement field. Pixels are in 32 bit float format.
				
				argv[6]: outputImageFile - (std::string) Path name and file name of the warped image using the displacement field. Pixels are in 32 bit float format.
				argv[7]: outConfigFile - (std::string) Path name of the configuration text file of the output image
				
Author:			Frankie (Hoi-Ki) Tong, Jan 12 2015

Changes:		Frankie (Hoi-Ki) Tong, Jan 12 2015
				Initial creation of the file.
*/


#include "mainTransform.h"

bool mainTransform (int argc, char *argv[]) {

	std::string inputImageFile1 = "";
	std::string inConfigFile1 = "";
	std::string inputVectorField1 = "";
	std::string inVectorConfigFile1 = "";
	std::string outputImageFile = "";
	std::string outConfigFile = "";
	
	if (argc == 8) {

		inputImageFile1 = argv[2];
		inConfigFile1 = argv[3];
		inputVectorField1 = argv[4];
		inVectorConfigFile1 = argv[5];
		outputImageFile = argv[6];
		outConfigFile = argv[7];
	
	}

	else
	{
		std::cerr << "Incorrect number of parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " --transform inputImageFile1 inConfigFile1 inputVectorField1 inVectorConfigFile1 outputImageFile outConfigFile" << std::endl;
		return EXIT_FAILURE;
	}

	
	
	

	//Read in first image
	typedef float									PixelType;
	const unsigned int Dimension = 2;
	typedef itk::Image<PixelType, Dimension>		ImageType;
	
	ImageType::Pointer normalImage = ImageType::New();
	
	itkImageIORaw<PixelType, Dimension>* inputImageImport = new itkImageIORaw<PixelType, Dimension>;

	inputImageImport->setConfigFileName(inConfigFile1);
	if (inputImageImport->readConfigFile()) { 
		std::cerr << "Exception thrown while reading the in configuration file" << std::endl;
		return EXIT_FAILURE;
	}
	inputImageImport->setitkImage((void *) &normalImage);
	inputImageImport->setImageFileName(inputImageFile1);
	if (inputImageImport->readRawImageKnown()) {
		std::cerr << "Image " << inputImageFile1 << " using config file " << inConfigFile1 << " did not read in properly!"<< std::endl;
		return EXIT_FAILURE;
	}

	delete inputImageImport;

	std::cout << "Supposedly, image is read in properly." << std::endl;	

	
	
	
	//Read in vector field
	typedef itk::Vector< PixelType, Dimension > VectorPixelType;  
	typedef itk::Image< VectorPixelType, Dimension > VectorFieldImageType;  
	
	VectorFieldImageType::Pointer vectorField = VectorFieldImageType::New();
	
	itkImageIORaw<VectorPixelType, Dimension>* inputVectorFieldImport = new itkImageIORaw<VectorPixelType, Dimension>;
	
	inputVectorFieldImport->setConfigFileName(inVectorConfigFile1);
	if (inputVectorFieldImport->readConfigFile()) {
		std::cerr << "Exception thrown while reading the in configuration file" << std::endl;
		return EXIT_FAILURE;
	}
	inputVectorFieldImport->setitkImage((void *) &vectorField);
	inputVectorFieldImport->setImageFileName(inputVectorField1);
	if (inputVectorFieldImport->readRawImageKnown()) {
		std::cerr << "Image " << inputVectorField1 << " using config file " << inVectorConfigFile1 << " did not read in properly!"<< std::endl;
		return EXIT_FAILURE;
	}
	
	delete inputVectorFieldImport;
	
	//Transform the normal image using the displacement field given as input to the program


	typedef itk::WarpImageFilter< ImageType, ImageType, VectorFieldImageType >		WarperType;
	typedef	itk::BSplineInterpolateImageFunction< ImageType, double >			InterpolatorType;


	//instantiation and setting the spline order of the interpolator.  
	int splineOrder = 3;
	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	interpolator->SetSplineOrder(splineOrder);


	WarperType::Pointer imageWarper = WarperType::New();

	imageWarper->SetInput( normalImage );
	imageWarper->SetInterpolator( interpolator );
	imageWarper->SetOutputSpacing( normalImage->GetSpacing() );
	imageWarper->SetOutputOrigin( normalImage->GetOrigin() );
	imageWarper->SetDeformationField( vectorField );
  
	imageWarper->Update();

	ImageType::Pointer transformedImage = imageWarper->GetOutput();


	//Output all appropriate outputs

	//Write output image and config file
	itkImageIORaw<PixelType, Dimension>* inputImageImport2 = new itkImageIORaw<PixelType, Dimension>;

	inputImageImport2->setConfigFileName(outConfigFile);

	//Write ConfigFile
	//The way to build ImageIOData structure differs from implementation to implementation which means I can't put this in itkImageIORaw for now...
	ImageIOData outputImageIOData;

	outputImageIOData.setPixelType("float");
	outputImageIOData.setPixelDimensionality(1);
	outputImageIOData.setImageDimensionality(Dimension);

	const ImageType::RegionType& OutputRegion = transformedImage->GetLargestPossibleRegion();
	const ImageType::SizeType& OutputSize = OutputRegion.GetSize();
	unsigned int tmpSize[Dimension];
	for (int i = 0; i < Dimension; i++) {
		tmpSize[i] = OutputSize[i];
	}

	outputImageIOData.setDimensionSize((int *)tmpSize);

	const ImageType::SpacingType OutputSpacing = transformedImage->GetSpacing();
	float tmpSpacing[Dimension];
	for (int i = 0; i < Dimension; i++) {
		tmpSpacing[i] = (float)OutputSpacing[i];
	}
	
	outputImageIOData.setSpacingSize((float *)tmpSpacing);

	const ImageType::PointType& OutputOrigin  = transformedImage->GetOrigin();
	float tmpOrigin[Dimension];
	for (int i=0; i< Dimension; i++) {
		tmpOrigin[i] = (float)OutputOrigin[i];
	}

	outputImageIOData.setOrigin((float *)tmpOrigin);

	outputImageIOData.setHeaderSize(0);
	outputImageIOData.setByteOrder("little");

	inputImageImport2->setImageIOData(outputImageIOData);

	inputImageImport2->writeConfigFile();

	//Write Image
	inputImageImport2->setitkImage((void *) &transformedImage);
	inputImageImport2->setImageFileName(outputImageFile);
	inputImageImport2->writeRawImageKnown();

	delete inputImageImport2;

	std::cout << "Supposedly, output image is written in properly." << std::endl;

	
	return EXIT_SUCCESS;

}