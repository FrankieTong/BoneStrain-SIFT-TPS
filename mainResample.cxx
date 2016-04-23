/*
File:			mainResample.cxx

Description:	Subfunction called from imageioraw.cxx resamples the input image to a different resolution using a B-Spline interpolation scheme. Input parameters based on input processing from calling mainResample as a command line function, ignoring the first 3 input parameters. Refer below for input process parameters.

Inputs:			argc - (int) Number of input arguments stored in argv. Default number is 8.
				argv[2]: inputImageFile1 - (std::string) Path name of the raw data file of the input image. Pixels are in 32 bit float format.
				argv[3]: inConfigFile1 - (std::string) Path name of the configuration text file of the input image. Pixels are in 32 bit float format.
				
				argv[4]: outputImageFile - (std::string) Path name and file name of the warped image using the displacement field. Pixels are in 32 bit float format.
				argv[5]: outConfigFile - (std::string) Path name of the configuration text file of the output image
				
				argv[6]: resampleFactorStr - (std::string) Downsample factor of the output image. Input is a float number in a string format.
				argv[7]: relativePosition - (std::string) Determines where the origin of the output image/field files are placed. They can either be placed at the "origin" or at the "center" of the image/field. Input is a string containing "origin" or "center.
				
Author:			Frankie (Hoi-Ki) Tong, Feb 16 2015

Changes:		Frankie (Hoi-Ki) Tong, Feb 16 2015
				Initial creation of the file.
*/

#include "mainResample.h"

bool mainResample (int argc, char *argv[]) {

	std::string inputImageFile1 = "";
	std::string inConfigFile1 = "";
	std::string outputImageFile = "";
	std::string outConfigFile = "";
	std::string resampleFactorStr = "";
	std::string relativePosition = "";
	
	if (argc == 8) {

		inputImageFile1 = argv[2];
		inConfigFile1 = argv[3];
		outputImageFile = argv[4];
		outConfigFile = argv[5];
		resampleFactorStr = argv[6];
		relativePosition = argv[7];
	
	}

	else
	{
		std::cerr << "Incorrect number of parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " --Resample inputImageFile1 inConfigFile1 outputImageFile outConfigFile resampleFactor relativePosition(origin/center)" << std::endl;
		return EXIT_FAILURE;
	}

	double resampleFactor = 0;
	std::istringstream(resampleFactorStr) >> resampleFactor;	
	
	

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

	
	
	//Resample image here...

	typedef itk::IdentityTransform<double, Dimension> upSampleIdentityTransformType;
	//typedef itk::LinearInterpolateImageFunction<ImageType,double> upSampleInterpolatorType;
	typedef	itk::BSplineInterpolateImageFunction< ImageType, double > upSampleInterpolatorType;
	typedef itk::ResampleImageFilter<ImageType,ImageType,double,double> upSampleResampleImaageFilterType;

	upSampleIdentityTransformType::Pointer upSampleIdentityTransform = upSampleIdentityTransformType::New();
	upSampleIdentityTransform->SetIdentity();

	upSampleInterpolatorType::Pointer upSampleInterpolator = upSampleInterpolatorType::New();

	upSampleResampleImaageFilterType::Pointer upSampleResampleImageFilter = upSampleResampleImaageFilterType::New();
	upSampleResampleImageFilter->SetTransform(upSampleIdentityTransform);
	upSampleResampleImageFilter->SetInterpolator(upSampleInterpolator);

	const ImageType::RegionType& ResampleInputRegion = normalImage->GetLargestPossibleRegion();
	const ImageType::SizeType& ResampleInputSize = ResampleInputRegion.GetSize();
	unsigned int OldSize[Dimension];
	unsigned int NewSize[Dimension];
	for (int i = 0; i < Dimension; i++) {
		OldSize[i] = ResampleInputSize[i];
		NewSize[i] = std::floor((double) OldSize[i]/resampleFactor );
	}

	const ImageType::SpacingType OldSpacing = normalImage->GetSpacing();
	double NewSpacing[Dimension];
	for (int i = 0; i < Dimension; i++) {
		NewSpacing[i] = OldSpacing[i] * (double) OldSize[i]/NewSize[i];
	}
	upSampleResampleImageFilter->SetOutputSpacing(NewSpacing);

	const ImageType::PointType& OldOrigin  = normalImage->GetOrigin();
	double ResampleOrigin[Dimension];
	for (int i=0; i< Dimension; i++) {

		if (relativePosition == "origin") {
			ResampleOrigin[i] = OldOrigin[i]-(OldSpacing[i] - NewSpacing[i])/2.0;
		} else if (relativePosition == "center") {
			ResampleOrigin[i] = OldOrigin[i] - ((OldSpacing[i] * (OldSize[i] - 1)) - (NewSpacing[i] * (NewSize[i] - 1)))/2.0;
		} else {
		
			std::cerr << "Value \"" << relativePosition << "\" for relative position is unknown. Please use either \"origin\" or \"center\"." << std::endl;
			return EXIT_FAILURE;
		}
		
	}
	upSampleResampleImageFilter->SetOutputOrigin(ResampleOrigin);

	itk::Size<Dimension> upSampleSize;
	for (int i = 0; i < Dimension; i++) {
		upSampleSize[i] = NewSize[i];
	}

	upSampleResampleImageFilter->SetSize(upSampleSize);

	upSampleResampleImageFilter->SetInput(normalImage);
	upSampleResampleImageFilter->Update();

	normalImage = upSampleResampleImageFilter->GetOutput();


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

	const ImageType::RegionType& OutputRegion = normalImage->GetLargestPossibleRegion();
	const ImageType::SizeType& OutputSize = OutputRegion.GetSize();
	unsigned int tmpSize[Dimension];
	for (int i = 0; i < Dimension; i++) {
		tmpSize[i] = OutputSize[i];
	}

	outputImageIOData.setDimensionSize((int *)tmpSize);

	const ImageType::SpacingType OutputSpacing = normalImage->GetSpacing();
	float tmpSpacing[Dimension];
	for (int i = 0; i < Dimension; i++) {
		tmpSpacing[i] = (float)OutputSpacing[i];
	}
	
	outputImageIOData.setSpacingSize((float *)tmpSpacing);

	const ImageType::PointType& OutputOrigin  = normalImage->GetOrigin();
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
	inputImageImport2->setitkImage((void *) &normalImage);
	inputImageImport2->setImageFileName(outputImageFile);
	inputImageImport2->writeRawImageKnown();

	delete inputImageImport2;

	std::cout << "Supposedly, output image is written in properly." << std::endl;
	
	return EXIT_SUCCESS;

}