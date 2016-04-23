/*
File:			main.cxx

Description:	Starting main function for utilizing the SIFT feature registration tools. When called from the command line, the input format is dictated as follows:
				
				"main.exe (--mode) (mode parameter 1) (mode parameter 2) ... (mode parameter n)"
				
				(--mode) stands for the type of processing mode you wish to invoke. The possible modes are:
				
				--standard
				--standardResample
				--standardResampleScaleAdjust
				--standardResampleScaleControl
				--standardResampleScaleMatch
				--standardResampleScaleMatchSIFT
				--rotate
				--warp
				--transform
				--resample
				--gradient
				--TPSPerfectMatch
				--TPS
				--BSplinePoints
				
				The input parameters following (--mode) differ depending on the mode called. Refer to the individual function files for more details on the inputs for each mode.

Inputs:			argc - (int) Number of input arguments into the function
				argv - (char *) Character pointer containing the inputs of this function

Author:			Frankie (Hoi-Ki) Tong, September 26 2013

Changes:		Frankie (Hoi-Ki) Tong, September 26 2013
				Initial creation of the file.
				
				Frankie (Hoi-Ki) Tong, April 12 2015
				Start of header documentation.
*/

#include "itkRawImageIO.h"
#include "itkNumericSeriesFileNames.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkSimpleFilterWatcher.h"
#include "itkShrinkImageFilter.h"
#include "itkComposeImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkShrinkImageFilter.h"
#include "itkMeanImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkPermuteAxesImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkBinaryThinningImageFilter.h"

#include "itkBinaryThinningImageFilter3D/itkBinaryThinningImageFilter3D.h"

#include "itkSkeleton/itkSkeletonizeImageFilter.h"
#include "itkSkeleton/itkChamferDistanceTransformImageFilter.h"

#include "itknSift/itkScaleInvariantFeatureImageFilter.h"

#include "ImageStats.h"
#include "itkImageIORaw.h"

#include "makeACube.h"

#include <iostream>
#include <new>
#include <list>
#include <string>
#include <ctime>

#include "mainStandard.h"
#include "mainStandardResample.h"
#include "mainStandardResampleScaleAdjust.h"
#include "mainStandardResampleScaleControl.h"
#include "mainStandardResampleScaleMatch.h"
#include "mainStandardResampleScaleMatchSIFT.h"
#include "mainRotate.h"
#include "mainWarp.h"
#include "mainTransform.h"
#include "mainResample.h"
#include "mainGradient.h"
#include "mainTPSPerfectMatch.h"
#include "mainTPS.h"
#include "mainBSplinePoints.h"

//#define BINARY_TEST

#ifdef BINARY_TEST
int main(int argc, char *argv[])
{
  if ( argc < 5 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImageFile  outputImageFile inConfigFile outConfigFile" << std::endl;
    return EXIT_FAILURE;
    }
	
	std::string inputImageFile = argv[1];
	std::string outputImageFile = argv[2];
	std::string inConfigFile = argv[3];
	std::string outConfigFile = argv[4];

	//Read in inConfigFile and set up import filter
	std::list<std::string> * inConfigFileParam = new std::list<std::string>;
	std::list<std::string> * inConfigFileContent = new std::list<std::string>;

	ImageConfigFileImport(inConfigFile,inConfigFileParam,inConfigFileContent);
	if (!ImageConfigFileCheck(inConfigFileParam,inConfigFileContent)) {
		std::cerr << "Exception thrown while reading the in configuration file" << std::endl;
		//needs to kill program. Find how
	}

	//Subfunction to find a specific parameter
	//Template function that casts input image type into one specific type
	//Want to cast specifically to float

	typedef float							PixelType;
	typedef itk::Image<PixelType, 3>		ImageType;
	ImageType::Pointer image = ImageType::New();

	void * vImage = & image;

	ImageImportReadImage<PixelType, 3>(inConfigFileParam,  inConfigFileContent,  inputImageFile, vImage);

	std::cout << "Supposedly, image is read in properly." << std::endl;

	ImageStats<ImageType>(image);

	//Binary the input image...
	//Limit range of image first...

	typedef itk::RescaleIntensityImageFilter< ImageType, ImageType > RescaleFilterType;
	RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
	rescaleFilter->SetInput(image);
	rescaleFilter->SetOutputMinimum(0);
	rescaleFilter->SetOutputMaximum(255);

	ImageType::Pointer rescaleImage = ImageType::New();
	rescaleImage = rescaleFilter->GetOutput();

	std::cout << "Supposedly, image is rescaled properly." << std::endl;
	ImageStats<ImageType>(rescaleImage);

	//Now cast it to char...
	typedef itk::Image<unsigned char, 3>	CastImageType;
	typedef itk::CastImageFilter< ImageType, CastImageType >	CastFilterType;
	CastFilterType::Pointer castFilter = CastFilterType::New();
	castFilter->SetInput(rescaleImage);

	CastImageType::Pointer castImage = CastImageType::New();
	castImage = castFilter->GetOutput();
 
	std::cout << "Supposedly, image is cast to unsigned char properly." << std::endl;
	ImageStats<CastImageType>(castImage);

	//Binary, lets go!

	int lowerThreshold = 100;
	int upperThreshold = 255;

	typedef itk::BinaryThresholdImageFilter <CastImageType, CastImageType> BinaryThresholdImageFilterType;
	BinaryThresholdImageFilterType::Pointer thresholdFilter = BinaryThresholdImageFilterType::New();

	thresholdFilter->SetInput(castImage);
	thresholdFilter->SetLowerThreshold(lowerThreshold);
	thresholdFilter->SetUpperThreshold(upperThreshold);
	thresholdFilter->SetInsideValue(255);
	thresholdFilter->SetOutsideValue(0);
 
	CastImageType::Pointer binaryImage = CastImageType::New();
	binaryImage = thresholdFilter->GetOutput();

	std::cout << "Supposedly, image is thresholded properly." << std::endl;
	ImageStats<CastImageType>(binaryImage);

	//Thinning algorithm...
	/*typedef itk::BinaryThinningImageFilter <CastImageType, CastImageType> BinaryThinningImageFilterType;
	BinaryThinningImageFilterType::Pointer binaryThinningImageFilter = BinaryThinningImageFilterType::New();
	binaryThinningImageFilter->SetInput(binaryImage);
	binaryThinningImageFilter->Update();
 
	// Rescale the image so that it can be seen (the output is 0 and 1, we want 0 and 255)
	typedef itk::RescaleIntensityImageFilter< CastImageType, CastImageType > RescaleType;
	RescaleType::Pointer rescaler = RescaleType::New();
	rescaler->SetInput( binaryThinningImageFilter->GetOutput() );
	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255);
	rescaler->Update();

	CastImageType::Pointer skeletonImage = CastImageType::New();
	skeletonImage = rescaler->GetOutput();

	std::cout << "Supposedly, image is skeletonized now." << std::endl;
	ImageStats<CastImageType>(skeletonImage);*/

	//3D Thinning Algorithm
	typedef itk::BinaryThinningImageFilter3D <CastImageType, CastImageType> BinaryThinningImageFilterType;
	BinaryThinningImageFilterType::Pointer binaryThinningImageFilter = BinaryThinningImageFilterType::New();
	binaryThinningImageFilter->SetInput(binaryImage);
	binaryThinningImageFilter->Update();
 
	// Rescale the image so that it can be seen (the output is 0 and 1, we want 0 and 255)
	typedef itk::RescaleIntensityImageFilter< CastImageType, CastImageType > RescaleType;
	RescaleType::Pointer rescaler = RescaleType::New();
	rescaler->SetInput( binaryThinningImageFilter->GetOutput() );
	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(255);
	rescaler->Update();

	CastImageType::Pointer skeletonImage = CastImageType::New();
	skeletonImage = rescaler->GetOutput();

	//3D Thinning Algorithm with itkSkeleton
	/*
	unsigned int const Dimension = 3;
	typedef itk::SkeletonizeImageFilter<CastImageType, itk::Connectivity<Dimension, 0> > Skeletonizer;

    typedef itk::ChamferDistanceTransformImageFilter<CastImageType, Skeletonizer::OrderingImageType> DistanceMapFilterType;
    DistanceMapFilterType::Pointer distanceMapFilter = DistanceMapFilterType::New();
    unsigned int weights[] = { 3, 4, 5 };
    distanceMapFilter->SetDistanceFromObject(false);
    distanceMapFilter->SetWeights(weights, weights+3);
    distanceMapFilter->SetInput(binaryImage);
    distanceMapFilter->Update();
    std::clog << "Distance map generated" << std::endl;

    Skeletonizer::Pointer skeletonizer = Skeletonizer::New();
    skeletonizer->SetInput(binaryImage);
    skeletonizer->SetOrderingImage(distanceMapFilter->GetOutput());
    skeletonizer->Update();
    CastImageType::Pointer skeletonImage = skeletonizer->GetOutput();
    std::clog << "Skeleton generated" << std::endl;*/
	
	std::cout << "Supposedly, image is 3D skeletonized now." << std::endl;
	ImageStats<CastImageType>(binaryImage);

	//Read in outConfigFile and set up export filter
	std::list<std::string> * outConfigFileParam = new std::list<std::string>;
	std::list<std::string> * outConfigFileContent = new std::list<std::string>;

	ImageConfigFileImport(outConfigFile,outConfigFileParam,outConfigFileContent);
	if (!ImageConfigFileCheck(outConfigFileParam,outConfigFileContent)) {
		std::cerr << "Exception thrown while reading the out configuration file" << std::endl;
		//needs to kill program. Find how
	}

	std::cout << "Printing Output File" << std::endl;
	vImage = & skeletonImage;

	ImageImportWriteImage<unsigned char, 3>(outConfigFileParam,  outConfigFileContent,  outputImageFile, vImage);

	delete inConfigFileParam;
	delete inConfigFileContent;
	delete outConfigFileParam;
	delete outConfigFileContent;
	return EXIT_SUCCESS;

}

#else

int main(int argc, char *argv[])
{

	typedef itk::Image< double, 2 >         ImageType;
	typedef itk::ImageFileReader<ImageType> ReaderType;
 
	ReaderType::Pointer reader = ReaderType::New();

	std::string errorMsg ("Please use either --standard/--standardResample/--standardResampleScaleAdjust/--standardResampleScaleControl/--standardResampleScaleMatch/--rotate/--warp/--transform/--resample/--gradient/--TPSPerfectMatch/--TPS/--BSplinePoints modes");	


	if (argc < 2) {
		std::cout<< errorMsg <<std::endl;	
		return EXIT_SUCCESS;
	}

	std::string mode = argv[1];

	if (mode == "--standard") {

		//Call the main function line
		if (mainStandard(argc, argv)) {
			return EXIT_SUCCESS;
		}
		return EXIT_SUCCESS;

	}

	else if (mode == "--standardResample") {

		//Call the main function line
		if (mainStandardResample(argc, argv)) {
			return EXIT_SUCCESS;
		}
		return EXIT_SUCCESS;

	}

	else if (mode == "--standardResampleScaleAdjust") {

		//Call the main function line
		if (mainStandardResampleScaleAdjust(argc, argv)) {
			return EXIT_SUCCESS;
		}
		return EXIT_SUCCESS;

	}
	
	else if (mode == "--standardResampleScaleControl") {

		//Call the main function line
		if (mainStandardResampleScaleControl(argc, argv)) {
			return EXIT_SUCCESS;
		}
		return EXIT_SUCCESS;

	}

	else if (mode == "--standardResampleScaleMatch") {

		//Call the main function line
		if (mainStandardResampleScaleMatch(argc, argv)) {
			return EXIT_SUCCESS;
		}
		return EXIT_SUCCESS;

	}

	else if (mode == "--standardResampleScaleMatchSIFT") {

		//Call the main function line
		if (mainStandardResampleScaleMatchSIFT(argc, argv)) {
			return EXIT_SUCCESS;
		}
		return EXIT_SUCCESS;

	}

	else if (mode == "--rotate") {
		if (mainRotate(argc, argv)) {
			return EXIT_SUCCESS;
		}
		return EXIT_SUCCESS;

	}

	else if (mode == "--warp") {
		if (mainWarp(argc, argv)) {
			return EXIT_SUCCESS;
		}
		return EXIT_SUCCESS;
	}

	else if (mode == "--transform") {
		if (mainTransform(argc, argv)) {
			return EXIT_SUCCESS;
		}
		return EXIT_SUCCESS;
	}

	else if (mode == "--resample") {
		if (mainResample(argc, argv)) {
			return EXIT_SUCCESS;
		}
		return EXIT_SUCCESS;
	}

	else if (mode == "--gradient") {
		if (mainGradient(argc, argv)) {
			return EXIT_SUCCESS;
		}
		return EXIT_SUCCESS;
	}

	else if (mode == "--TPSPerfectMatch") {
		if (mainTPSPerfectMatch(argc, argv)) {
			return EXIT_SUCCESS;
		} return EXIT_SUCCESS;
	}

	else if (mode == "--TPS") {
		if (mainTPS(argc, argv)) {
			return EXIT_SUCCESS;
		} return EXIT_SUCCESS;
	}

	else if (mode == "--BSplinePoints") {
		if (mainBSplinePoints(argc, argv)) {
			return EXIT_SUCCESS;
		} return EXIT_SUCCESS;
	}

	/*else if (mode == "--warp") {
		if (mainWarp(argc, argv)) {
			return EXIT_FAILURE;
		}
		return EXIT_SUCCESS;
	}*/

	else {
		std::cout<< errorMsg <<std::endl;	
		return EXIT_SUCCESS;
	}

}

#endif