/*
File:			mainStandardResampleScaleMatchSIFT.cxx

Description:	Subfunction called from imageioraw.cxx Performs SIFT feature registration on 2 input images. Input parameters based on input processing from calling mainStandardResampleScaleMatchSIFT as a command line function, ignoring the first 3 input parameters. Refer below for input process parameters.
				NOTE: This function operates exactly the same as mainStandardResampleScaleMatch() except it does not attempt to calculate and output a displacement field.

Inputs:			argc - (int) Number of input arguments stored in argv. Default number is 28+2*Dimension.
				argv[2]: inputImageFile1 - (std::string) Path name of the raw data file of the fixed image. Pixels are in 32 bit float format. 
				argv[3]: inputImageFile2 - (std::string) Path name of the raw data file of the moving image. Pixels are in 32 bit float format. 
				
				argv[4]: inImageFile1Min - (std::string) Minimum intensity file of the fixed image
				argv[5]: inImageFile1Max - (std::string) Maximum intensity file of the fixed image
				argv[6]: inImageFile2Min - (std::string) Minimum intensity file of the moving image
				argv[7]: inImageFile2Max - (std::string) Maximum intensity file of the moving image
				
				argv[8]: inputImageFile1 - (std::string) Path name of the raw data file of the binarized fixed image. Pixels are in 32 bit float format. 
				argv[9]: inputImageFile2 - (std::string) Path name of the raw data file of the binarized moving image. Pixels are in 32 bit float format. 
				
				argv[10]: inImageFile1Min - (std::string) Minimum intensity file of the binarized fixed image
				argv[11]: inImageFile1Max - (std::string) Maximum intensity file of the binarized fixed image
				argv[12]: inImageFile2Min - (std::string) Minimum intensity file of the binarized moving image
				argv[13]: inImageFile2Max - (std::string) Maximum intensity file of the binarized moving image
				
				argv[14]: inputConfigFile1 - (std::string) Path name of the configuration text file of the fixed image
				argv[15]: inputConfigFile2 - (std::string)  Path name of the configuration text file of the moving image
				
				argv[16]: outName - (std::string) Path name and file name prefix of the output files containing information on SIFT feature point positions
				argv[17]: resampleFactorStr - (std::string) Downsample factor of the input images before using SIFT. Input is a float number in a string format.
				argv[18]: relativePosition - (std::string) Determines where the origin of the output image/field files are placed. They can either be placed at the "origin" or at the "center" of the image/field. Input is a string containing "origin" or "center.
				argv[19]: outputResolution - (std::string) Downsample factor of the output images/fields after using SIFT. Resample based on original input image resolution. Input is a float number in a string format.
				argv[20]: inRemoveEdgePointsPercent - (std::string) Percentage of image removed starting from the image edge and moving inward based on toral image size in each dimension. Input is a float number [0,100] in string format.
				
				argv[21]: inScaleStart - (std::string) Starting octave level to attempt to find features. Input is a int number starting from 0. If you are unsure about this value, leave it at 0.
				argv[22]: inScaleEnd - (std::string) Ending octave level to attempt to find features. Input is a int number starting from , but greater than inScaleStart. If you are unsure about this value, leave it at 2.
				
				argv[23]: inEdgePointRatio - (std::string) Starting edge point ratio value for scale level 0. The higher the number, the more feature points will be accepted for matching. Input is a float number greater than 0. Input can also be a negative number if you want to disregard edge point checking for feature points. If you are unsure about this value, leave it at 10.
				argv[24]: inEdgePointRatioScale - (std::string) Edge point ratio multiplier that multiplies the edge point ratio by this number each time SIFT searches up a higher scale. Input is a float number. Default is 1 for no adjustment in the edge point ratio between scales.

				argv[25]: inMatchScaleOnly - (std::string) Boolean to dictate whether feature point matches can occur only when the features are from the same scale level. Input is either "true" or "false".
				argv[26]: inMatchExtremaOnly - (std::string) Boolean to dictate whether feature point matches can occur only when the features are the same extrema direction (max to max, min to min). Input is either "true" or "false".
				argv[27]: inTurnOffNearestNeighbourInterpolator - (std::string) Boolean to dictate whether the resample step of the binary images before SIFT uses the nearest neighbour filter or the B-Spline filter. Input is either "true" or "false".
				
				argv[28+2n]: inDisplacementLimitMin - (std::string) A series of floats dictating the minimum measured displacement allowable in each image direction. Input is a series of floats in image dimension units (not pixel units)
				argv[29+2n]: inDisplacementLimitMax - (std::string) A series of floats dictating the maximum measured displacement allowable in each image direction. Input is a series of floats in image dimension units (not pixel units)
			
Author:			Frankie (Hoi-Ki) Tong, Aug 20 2015

Changes:		Frankie (Hoi-Ki) Tong, Aug 20 2015
				Initial creation of the file.
*/

#include "mainStandardResampleScaleMatchSIFT.h"

bool mainStandardResampleScaleMatchSIFT (int argc, char *argv[]) {

	const unsigned int Dimension = 2;

	std::string inputImageFile1 = "";
	std::string inputImageFile2 = "";
	std::string inConfigFile1 = "";
	std::string inConfigFile2 = "";

	std::string inImageFile1Min = "";
	std::string inImageFile1Max = "";
	std::string inImageFile2Min = "";
	std::string inImageFile2Max = "";

	std::string inputImageFileBinary1 = "";
	std::string inputImageFileBinary2 = "";
	std::string inConfigFileBinary1 = "";
	std::string inConfigFileBinary2 = "";

	std::string inImageFileBinary1Min = "";
	std::string inImageFileBinary1Max = "";
	std::string inImageFileBinary2Min = "";
	std::string inImageFileBinary2Max = "";

	std::string outName = "";
	std::string resampleFactorStr = "";
	std::string relativePosition = "";
	std::string outputResolution = "";
	std::string inRemoveEdgePointsPixel = "";
	
	std::string inScaleStart = "";
	std::string inScaleEnd = "";

	std::string inEdgePointRatio = "";
	std::string inEdgePointRatioScale = "";

	std::string inMatchScaleOnly = "false";
	std::string inMatchExtremaOnly = "false";
	std::string inTurnOffNearestNeighbourInterpolator = "false";

	std::string * inDisplacementLimitMin = new std::string[Dimension];
	std::string * inDisplacementLimitMax = new std::string[Dimension];

	if (argc == 25) {

		inputImageFile1 = argv[2];
		inputImageFile2 = argv[3];

		inImageFile1Min = argv[4];
		inImageFile1Max = argv[5];
		inImageFile2Min = argv[6];
		inImageFile2Max = argv[7];

		inputImageFileBinary1 = argv[8];
		inputImageFileBinary2 = argv[9];

		inImageFileBinary1Min = argv[10];
		inImageFileBinary1Max = argv[11];
		inImageFileBinary2Min = argv[12];
		inImageFileBinary2Max = argv[13];

		inConfigFile1 = argv[14];
		inConfigFile2 = argv[15];

		outName = argv[16];
		resampleFactorStr = argv[17];
		relativePosition = argv[18];
		outputResolution = argv[19];
		inRemoveEdgePointsPixel = argv[20];
		inScaleStart = argv[21];
		inScaleEnd = argv[22];
		inEdgePointRatio = argv[23];
		inEdgePointRatioScale = argv[24];
		
	} else if (argc == 28) {

		inputImageFile1 = argv[2];
		inputImageFile2 = argv[3];

		inImageFile1Min = argv[4];
		inImageFile1Max = argv[5];
		inImageFile2Min = argv[6];
		inImageFile2Max = argv[7];

		inputImageFileBinary1 = argv[8];
		inputImageFileBinary2 = argv[9];

		inImageFileBinary1Min = argv[10];
		inImageFileBinary1Max = argv[11];
		inImageFileBinary2Min = argv[12];
		inImageFileBinary2Max = argv[13];

		inConfigFile1 = argv[14];
		inConfigFile2 = argv[15];

		outName = argv[16];
		resampleFactorStr = argv[17];
		relativePosition = argv[18];
		outputResolution = argv[19];
		inRemoveEdgePointsPixel = argv[20];
		inScaleStart = argv[21];
		inScaleEnd = argv[22];
		inEdgePointRatio = argv[23];
		inEdgePointRatioScale = argv[24];

		inMatchScaleOnly = argv[25];
		inMatchExtremaOnly = argv[26];
		inTurnOffNearestNeighbourInterpolator = argv[27];
		
	} else if (argc == 28 + 2*Dimension) {

		inputImageFile1 = argv[2];
		inputImageFile2 = argv[3];

		inImageFile1Min = argv[4];
		inImageFile1Max = argv[5];
		inImageFile2Min = argv[6];
		inImageFile2Max = argv[7];

		inputImageFileBinary1 = argv[8];
		inputImageFileBinary2 = argv[9];

		inImageFileBinary1Min = argv[10];
		inImageFileBinary1Max = argv[11];
		inImageFileBinary2Min = argv[12];
		inImageFileBinary2Max = argv[13];

		inConfigFile1 = argv[14];
		inConfigFile2 = argv[15];

		outName = argv[16];
		resampleFactorStr = argv[17];
		relativePosition = argv[18];
		outputResolution = argv[19];
		inRemoveEdgePointsPixel = argv[20];
		inScaleStart = argv[21];
		inScaleEnd = argv[22];
		inEdgePointRatio = argv[23];
		inEdgePointRatioScale = argv[24];

		inMatchScaleOnly = argv[25];
		inMatchExtremaOnly = argv[26];
		inTurnOffNearestNeighbourInterpolator = argv[27];

		for (int i = 0; i < Dimension; i++) {
			inDisplacementLimitMin[i] = argv[28+i];
			inDisplacementLimitMax[i] = argv[28+Dimension+i];
		}
	}

	else
	{
		std::cerr << "Incorrect number of parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " --standardResampleScaleMatch inputImageFile1 inputImageFile2 inImageFile1Min inImageFile1Max inImageFile2Min inImageFile2Max inputImageFileBinary1 inputImageFileBinary2 inImageFileBinary1Min inImageFileBinary1Max inImageFileBinary2Min inImageFileBinary2Max iConfigFile1 iConfigFile2 outName resampleSpacingFactor relativePosition(origin/center) outputDisplacementResolution(original/resampled/both/#) inRemoveEdgePointsPixel inScaleStart inScaleEnd inEdgePointRatio inEdgePointRatioScale [optional input for matching on same scale only] [optional input for matching on same extremas only] [optional input for turning off nearest neighbour intepolator] [optional input for displacement vector limits]"<< std::endl;
		return EXIT_FAILURE;
	}

	float resampleFactor = 0; //defaults
	std::istringstream(resampleFactorStr) >> resampleFactor;
	
	//Get percentage of the image where points found in the area will be removed
	float removeEdgePointsPixel = 0; //defaults
	std::istringstream(inRemoveEdgePointsPixel) >>  removeEdgePointsPixel;
	
	int scaleStart = 0; //defaults
	std::istringstream(inScaleStart) >>  scaleStart;
	
	int scaleEnd = 3; //defaults
	std::istringstream(inScaleEnd) >>  scaleEnd;

	float edgePointRatio = 10; //defaults
	std::istringstream(inEdgePointRatio) >>  edgePointRatio;

	float edgePointRatioScale = 1; //defaults
	std::istringstream(inEdgePointRatioScale) >>  edgePointRatioScale;

	bool matchSameScale = true;
	bool matchSameExtrema = true;

	if (inMatchScaleOnly == "false") {
		matchSameScale = false;
	}
	if (inMatchExtremaOnly == "false") {
		matchSameExtrema = false;
	}

	float * displacementLimitMin = new float[Dimension];
	float * displacementLimitMax = new float[Dimension];
	
	if (argc >= 26 + 2*Dimension) {

		for (int i = 0; i < Dimension; i++) {
			std::istringstream(inDisplacementLimitMin[i]) >>  displacementLimitMin[i];
			std::istringstream(inDisplacementLimitMax[i]) >>  displacementLimitMax[i];
		}

	}

	//Initialize Imagetypes and image pointers
	typedef float									PixelType;
	typedef itk::Image<PixelType, Dimension>		ImageType;
	
	ImageType::Pointer fixedImage = ImageType::New();
	ImageType::Pointer fixedImageBinary = ImageType::New();
	ImageType::Pointer movingImage = ImageType::New();
	ImageType::Pointer movingImageBinary = ImageType::New();
	
	itkImageIORaw<PixelType, Dimension>* inputImageImport = new itkImageIORaw<PixelType, Dimension>;

	//Read in first image
	inputImageImport->setConfigFileName(inConfigFile1);
	if (inputImageImport->readConfigFile()) { 
		std::cerr << "Exception thrown while reading the in configuration file" << std::endl;
		return EXIT_FAILURE;
	}
	inputImageImport->setitkImage((void *) &fixedImage);
	inputImageImport->setImageFileName(inputImageFile1);
	if (inputImageImport->readRawImageKnown()) {
		std::cerr << "Image " << inputImageFile1 << " using config file " << inConfigFile1 << " did not read in properly!"<< std::endl;
		return EXIT_FAILURE;
	}

	//Read in first image binary
	inputImageImport->setConfigFileName(inConfigFile1);
	if (inputImageImport->readConfigFile()) { 
		std::cerr << "Exception thrown while reading the in configuration file" << std::endl;
		return EXIT_FAILURE;
	}
	inputImageImport->setitkImage((void *) &fixedImageBinary);
	inputImageImport->setImageFileName(inputImageFileBinary1);
	if (inputImageImport->readRawImageKnown()) {
		std::cerr << "Image " << inputImageFileBinary1 << " using config file " << inConfigFile1 << " did not read in properly!"<< std::endl;
		return EXIT_FAILURE;
	}

	//Read in second image
	inputImageImport->setConfigFileName(inConfigFile2);
	if (inputImageImport->readConfigFile()) { 
		std::cerr << "Exception thrown while reading the in configuration file" << std::endl;
		return EXIT_FAILURE;
	}
	inputImageImport->setitkImage((void *) &movingImage);
	inputImageImport->setImageFileName(inputImageFile2);
	if (inputImageImport->readRawImageKnown()) {
		std::cerr << "Image " << inputImageFile2 << " using config file " << inConfigFile2 << " did not read in properly!"<< std::endl;
		return EXIT_FAILURE;
	}

	//Read in first image binary
	inputImageImport->setConfigFileName(inConfigFile2);
	if (inputImageImport->readConfigFile()) { 
		std::cerr << "Exception thrown while reading the in configuration file" << std::endl;
		return EXIT_FAILURE;
	}
	inputImageImport->setitkImage((void *) &movingImageBinary);
	inputImageImport->setImageFileName(inputImageFileBinary2);
	if (inputImageImport->readRawImageKnown()) {
		std::cerr << "Image " << inputImageFileBinary2 << " using config file " << inConfigFile2 << " did not read in properly!"<< std::endl;
		return EXIT_FAILURE;
	}

	

	delete inputImageImport;

	std::cout << "Supposedly, images are read in properly." << std::endl;	

#ifdef TIMER
	//Timer functionality
	std::clock_t start;
	double resample_duration;
	double intensity_rescale_duration;
	double find_points_duration;
	double remove_edge_points_duration;
	double match_points_duration;
	double TPS_duration;
	double total_duration;

	start = std::clock();

	total_duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
#endif

	//Upsample the input images if needed
	ImageType::Pointer resampledFixedImage = ImageType::New();
	ImageType::Pointer resampledFixedImageBinary = ImageType::New();
	ImageType::Pointer resampledMovingImage = ImageType::New();
	ImageType::Pointer resampledMovingImageBinary = ImageType::New();
	
	if (resampleFactor == 1) {
		resampledFixedImage = fixedImage;
		resampledFixedImageBinary = fixedImageBinary;
		resampledMovingImage = movingImage;
		resampledMovingImageBinary = movingImageBinary;
	} else {
		
		//Need to resample both images to new resolution
		
		//Resampleing fixedImage and fixedImageBinary
		typedef itk::IdentityTransform<double, Dimension> upSampleIdentityTransformType;
		//typedef itk::LinearInterpolateImageFunction<ImageType,double> upSampleInterpolatorType;
		typedef	itk::BSplineInterpolateImageFunction< ImageType, double > upSampleInterpolatorType;
		typedef itk::ResampleImageFilter<ImageType,ImageType,double,double> upSampleResampleImageFilterType;

		typedef itk::NearestNeighborInterpolateImageFunction< ImageType, double> upSampleBinaryInterpolatorType;
		typedef itk::ResampleImageFilter<ImageType,ImageType,double,double> upSampleBinaryResampleImageFilterType;

		upSampleIdentityTransformType::Pointer upSampleIdentityTransformFixedImage = upSampleIdentityTransformType::New();
		upSampleIdentityTransformFixedImage->SetIdentity();


		upSampleInterpolatorType::Pointer upSampleInterpolatorFixedImage = upSampleInterpolatorType::New();
		upSampleBinaryInterpolatorType::Pointer upSampleBinaryInterpolatorFixedImage = upSampleBinaryInterpolatorType::New();

		upSampleResampleImageFilterType::Pointer upSampleResampleImageFilterFixedImage = upSampleResampleImageFilterType::New();
		upSampleResampleImageFilterFixedImage->SetTransform(upSampleIdentityTransformFixedImage);
		upSampleResampleImageFilterFixedImage->SetInterpolator(upSampleInterpolatorFixedImage);

		upSampleBinaryResampleImageFilterType::Pointer upSampleBinaryResampleImageFilterFixedImage = upSampleBinaryResampleImageFilterType::New();
		upSampleBinaryResampleImageFilterFixedImage->SetTransform(upSampleIdentityTransformFixedImage);
		if (inTurnOffNearestNeighbourInterpolator == "false") {
			upSampleBinaryResampleImageFilterFixedImage->SetInterpolator(upSampleBinaryInterpolatorFixedImage);
		} else {
			upSampleBinaryResampleImageFilterFixedImage->SetInterpolator(upSampleInterpolatorFixedImage);
		}

		const ImageType::RegionType& ResampleInputRegionFixedImage = fixedImage->GetLargestPossibleRegion();
		const ImageType::SizeType& ResampleInputSizeFixedImage = ResampleInputRegionFixedImage.GetSize();
		unsigned int OldSizeFixedImage[Dimension];
		unsigned int NewSizeFixedImage[Dimension];
		for (int i = 0; i < Dimension; i++) {
			OldSizeFixedImage[i] = ResampleInputSizeFixedImage[i];
			NewSizeFixedImage[i] = std::floor((double) OldSizeFixedImage[i]/resampleFactor );
		}

		const ImageType::SpacingType OldSpacingFixedImage = fixedImage->GetSpacing();
		double NewSpacingFixedImage[Dimension];
		for (int i = 0; i < Dimension; i++) {
			NewSpacingFixedImage[i] = OldSpacingFixedImage[i] * (double) OldSizeFixedImage[i]/NewSizeFixedImage[i];
		}
		upSampleResampleImageFilterFixedImage->SetOutputSpacing(NewSpacingFixedImage);
		upSampleBinaryResampleImageFilterFixedImage->SetOutputSpacing(NewSpacingFixedImage);

		const ImageType::PointType& OldOriginFixedImage  = fixedImage->GetOrigin();
		double ResampleOriginFixedImage[Dimension];
		for (int i=0; i< Dimension; i++) {

			if (relativePosition == "origin") {
				ResampleOriginFixedImage[i] = OldOriginFixedImage[i]-(OldSpacingFixedImage[i] - NewSpacingFixedImage[i])/2.0;
			} else if (relativePosition == "center") {
				ResampleOriginFixedImage[i] = OldOriginFixedImage[i] - ((OldSpacingFixedImage[i] * (OldSizeFixedImage[i] - 1)) - (NewSpacingFixedImage[i] * (NewSizeFixedImage[i] - 1)))/2.0;
			} else {
			
				std::cerr << "Value \"" << relativePosition << "\" for relative position is unknown. Please use either \"origin\" or \"center\"." << std::endl;
				return EXIT_FAILURE;
			}
			
		}
		upSampleResampleImageFilterFixedImage->SetOutputOrigin(ResampleOriginFixedImage);
		upSampleBinaryResampleImageFilterFixedImage->SetOutputOrigin(ResampleOriginFixedImage);

		itk::Size<Dimension> upSampleSizeFixedImage;
		for (int i = 0; i < Dimension; i++) {
			upSampleSizeFixedImage[i] = NewSizeFixedImage[i];
		}

		upSampleResampleImageFilterFixedImage->SetSize(upSampleSizeFixedImage);
		upSampleBinaryResampleImageFilterFixedImage->SetSize(upSampleSizeFixedImage);

		upSampleResampleImageFilterFixedImage->SetInput(fixedImage);
		upSampleBinaryResampleImageFilterFixedImage->SetInput(fixedImageBinary);

		//upSampleResampleImageFilterFixedImage->ReleaseDataFlagOn();
		upSampleResampleImageFilterFixedImage->Update();
		upSampleBinaryResampleImageFilterFixedImage->Update();

		resampledFixedImage = upSampleResampleImageFilterFixedImage->GetOutput();
		resampledFixedImageBinary = upSampleBinaryResampleImageFilterFixedImage->GetOutput();
		
		//upSampleResampleImageFilterFixedImage = NULL;
			
		
		
		//Resampleing movingImage

		upSampleIdentityTransformType::Pointer upSampleIdentityTransformMovingImage = upSampleIdentityTransformType::New();
		upSampleIdentityTransformMovingImage->SetIdentity();

		upSampleInterpolatorType::Pointer upSampleInterpolatorMovingImage = upSampleInterpolatorType::New();
		upSampleBinaryInterpolatorType::Pointer upSampleBinaryInterpolatorMovingImage = upSampleBinaryInterpolatorType::New();

		upSampleResampleImageFilterType::Pointer upSampleResampleImageFilterMovingImage = upSampleResampleImageFilterType::New();
		upSampleResampleImageFilterMovingImage->SetTransform(upSampleIdentityTransformMovingImage);
		upSampleResampleImageFilterMovingImage->SetInterpolator(upSampleInterpolatorMovingImage);

		upSampleBinaryResampleImageFilterType::Pointer upSampleBinaryResampleImageFilterMovingImage = upSampleBinaryResampleImageFilterType::New();
		upSampleBinaryResampleImageFilterMovingImage->SetTransform(upSampleIdentityTransformMovingImage);
		upSampleBinaryResampleImageFilterMovingImage->SetInterpolator(upSampleBinaryInterpolatorMovingImage);

		const ImageType::RegionType& ResampleInputRegionMovingImage = movingImage->GetLargestPossibleRegion();
		const ImageType::SizeType& ResampleInputSizeMovingImage = ResampleInputRegionMovingImage.GetSize();
		unsigned int OldSizeMovingImage[Dimension];
		unsigned int NewSizeMovingImage[Dimension];
		for (int i = 0; i < Dimension; i++) {
			OldSizeMovingImage[i] = ResampleInputSizeMovingImage[i];
			NewSizeMovingImage[i] = std::floor((double) OldSizeMovingImage[i]/resampleFactor );
		}

		const ImageType::SpacingType OldSpacingMovingImage = movingImage->GetSpacing();
		double NewSpacingMovingImage[Dimension];
		for (int i = 0; i < Dimension; i++) {
			NewSpacingMovingImage[i] = OldSpacingMovingImage[i] * (double) OldSizeMovingImage[i]/NewSizeMovingImage[i];
		}
		upSampleResampleImageFilterMovingImage->SetOutputSpacing(NewSpacingMovingImage);
		upSampleBinaryResampleImageFilterMovingImage->SetOutputSpacing(NewSpacingMovingImage);

		const ImageType::PointType& OldOriginMovingImage  = movingImage->GetOrigin();
		double ResampleOriginMovingImage[Dimension];
		for (int i=0; i< Dimension; i++) {

			if (relativePosition == "origin") {
				ResampleOriginMovingImage[i] = OldOriginMovingImage[i]-(OldSpacingMovingImage[i] - NewSpacingMovingImage[i])/2.0;
			} else if (relativePosition == "center") {
				ResampleOriginMovingImage[i] = OldOriginMovingImage[i] - ((OldSpacingMovingImage[i] * (OldSizeMovingImage[i] - 1)) - (NewSpacingMovingImage[i] * (NewSizeMovingImage[i] - 1)))/2.0;
			} else {
			
				std::cerr << "Value \"" << relativePosition << "\" for relative position is unknown. Please use either \"origin\" or \"center\"." << std::endl;
				return EXIT_FAILURE;
			}
			
		}
		upSampleResampleImageFilterMovingImage->SetOutputOrigin(ResampleOriginMovingImage);
		upSampleBinaryResampleImageFilterMovingImage->SetOutputOrigin(ResampleOriginMovingImage);

		itk::Size<Dimension> upSampleSizeMovingImage;
		for (int i = 0; i < Dimension; i++) {
			upSampleSizeMovingImage[i] = NewSizeMovingImage[i];
		}

		upSampleResampleImageFilterMovingImage->SetSize(upSampleSizeMovingImage);
		upSampleBinaryResampleImageFilterMovingImage->SetSize(upSampleSizeMovingImage);

		upSampleResampleImageFilterMovingImage->SetInput(movingImage);
		upSampleBinaryResampleImageFilterMovingImage->SetInput(movingImageBinary);

		//upSampleResampleImageFilterMovingImage->ReleaseDataFlagOn();
		upSampleResampleImageFilterMovingImage->Update();
		upSampleBinaryResampleImageFilterMovingImage->Update();

		resampledMovingImage = upSampleResampleImageFilterMovingImage->GetOutput();
		resampledMovingImageBinary = upSampleBinaryResampleImageFilterMovingImage->GetOutput();

		//upSampleResampleImageFilterMovingImage = NULL;
	}
		
#ifdef TIMER
	resample_duration = (( std::clock() - start ) / (double) CLOCKS_PER_SEC) - total_duration;
	total_duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
#endif
			
	//Rescale image intensities

	double imageFile1Min = 0;
	double imageFile1Max = 0;
	double imageFile2Min = 0;
	double imageFile2Max = 0;

	std::istringstream(inImageFile1Min) >> imageFile1Min;
	std::istringstream(inImageFile1Max) >> imageFile1Max;
	std::istringstream(inImageFile2Min) >> imageFile2Min;
	std::istringstream(inImageFile2Max) >> imageFile2Max;

	double imageFileBinary1Min = 0;
	double imageFileBinary1Max = 0;
	double imageFileBinary2Min = 0;
	double imageFileBinary2Max = 0;

	std::istringstream(inImageFileBinary1Min) >> imageFileBinary1Min;
	std::istringstream(inImageFileBinary1Max) >> imageFileBinary1Max;
	std::istringstream(inImageFileBinary2Min) >> imageFileBinary2Min;
	std::istringstream(inImageFileBinary2Max) >> imageFileBinary2Max;


	typedef itk::IntensityWindowingImageFilter <ImageType, ImageType> IntensityWindowingImageFilterType;
	IntensityWindowingImageFilterType::Pointer IntenistyWindowImageFilterFixedImage = IntensityWindowingImageFilterType::New();

	IntenistyWindowImageFilterFixedImage->SetInput(resampledFixedImage);
	IntenistyWindowImageFilterFixedImage->SetWindowMinimum(imageFile1Min);
	IntenistyWindowImageFilterFixedImage->SetWindowMaximum(imageFile1Max);
	IntenistyWindowImageFilterFixedImage->SetOutputMinimum(0);
	IntenistyWindowImageFilterFixedImage->SetOutputMaximum(1);
	IntenistyWindowImageFilterFixedImage->Update();
	resampledFixedImage = IntenistyWindowImageFilterFixedImage->GetOutput();

	IntensityWindowingImageFilterType::Pointer IntenistyWindowImageFilterFixedBinaryImage = IntensityWindowingImageFilterType::New();

	IntenistyWindowImageFilterFixedBinaryImage->SetInput(resampledFixedImageBinary);
	IntenistyWindowImageFilterFixedBinaryImage->SetWindowMinimum(imageFileBinary1Min);
	IntenistyWindowImageFilterFixedBinaryImage->SetWindowMaximum(imageFileBinary1Max);
	IntenistyWindowImageFilterFixedBinaryImage->SetOutputMinimum(0);
	IntenistyWindowImageFilterFixedBinaryImage->SetOutputMaximum(1);
	IntenistyWindowImageFilterFixedBinaryImage->Update();
	resampledFixedImageBinary = IntenistyWindowImageFilterFixedBinaryImage->GetOutput();

	IntensityWindowingImageFilterType::Pointer IntenistyWindowImageFilterMovingImage = IntensityWindowingImageFilterType::New();

	IntenistyWindowImageFilterMovingImage->SetInput(resampledMovingImage);
	IntenistyWindowImageFilterMovingImage->SetWindowMinimum(imageFile2Min);
	IntenistyWindowImageFilterMovingImage->SetWindowMaximum(imageFile2Max);
	IntenistyWindowImageFilterMovingImage->SetOutputMinimum(0);
	IntenistyWindowImageFilterMovingImage->SetOutputMaximum(1);
	IntenistyWindowImageFilterMovingImage->Update();
	resampledMovingImage = IntenistyWindowImageFilterMovingImage->GetOutput();

	IntensityWindowingImageFilterType::Pointer IntenistyWindowImageFilterMovingBinaryImage = IntensityWindowingImageFilterType::New();

	IntenistyWindowImageFilterMovingBinaryImage->SetInput(resampledMovingImageBinary);
	IntenistyWindowImageFilterMovingBinaryImage->SetWindowMinimum(imageFileBinary2Min);
	IntenistyWindowImageFilterMovingBinaryImage->SetWindowMaximum(imageFileBinary2Max);
	IntenistyWindowImageFilterMovingBinaryImage->SetOutputMinimum(0);
	IntenistyWindowImageFilterMovingBinaryImage->SetOutputMaximum(1);
	IntenistyWindowImageFilterMovingBinaryImage->Update();
	resampledMovingImageBinary = IntenistyWindowImageFilterMovingBinaryImage->GetOutput();

	/*std::cout<<std::endl<<"Resampled and intensity adjusted fixed image stats:"<<std::endl;
	ImageStats<ImageType>(resampledFixedImage);
	std::cout<<std::endl<<"Resampled and intensity adjusted moving image stats:"<<std::endl;
	ImageStats<ImageType>(resampledMovingImage);

	std::cout<<std::endl<<"Resampled and intensity adjusted fixed image binary stats:"<<std::endl;
	ImageStats<ImageType>(resampledFixedImageBinary);
	std::cout<<std::endl<<"Resampled and intensity adjusted moving image binary stats:"<<std::endl;
	ImageStats<ImageType>(resampledMovingImageBinary);*/

#ifdef TIMER
	intensity_rescale_duration = (( std::clock() - start ) / (double) CLOCKS_PER_SEC) - total_duration;
	total_duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
#endif


	//Find Keypoints and Matchedpoints betweent the 2 images
	typedef itk::ScaleInvariantFeatureImageFilter<ImageType, Dimension> SiftFilterType;

	SiftFilterType::PointSetTypePointer fixedImageKeyPoints, movingImageKeyPoints;

	//Get Keypoints from fixed image and moving image
	SiftFilterType * siftFilter1 = new SiftFilterType;

	

	if (matchSameScale || matchSameExtrema) {
		
		fixedImageKeyPoints = siftFilter1->getSiftFeaturesRescalingMatch(resampledFixedImage, resampledFixedImageBinary, scaleStart, scaleEnd, edgePointRatio, edgePointRatioScale, true);
		//delete siftFilter1;

		//SiftFilterType * siftFilter2 = new SiftFilterType;
		movingImageKeyPoints = siftFilter1->getSiftFeaturesRescalingMatch(resampledMovingImage, resampledMovingImageBinary, scaleStart, scaleEnd, edgePointRatio, edgePointRatioScale, true);
		//delete siftFilter2;
		

	} else {

		fixedImageKeyPoints = siftFilter1->getSiftFeaturesRescalingMatch(resampledFixedImage, resampledFixedImageBinary, scaleStart, scaleEnd, edgePointRatio, edgePointRatioScale);
		//delete siftFilter1;

		//SiftFilterType * siftFilter2 = new SiftFilterType;
		movingImageKeyPoints = siftFilter1->getSiftFeaturesRescalingMatch(resampledMovingImage, resampledMovingImageBinary, scaleStart, scaleEnd, edgePointRatio, edgePointRatioScale);
		//delete siftFilter2;

	}

#ifdef TIMER
	find_points_duration = (( std::clock() - start ) / (double) CLOCKS_PER_SEC) - total_duration;
	total_duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
#endif

	//Find and remove edge points if requested by user
	if (removeEdgePointsPixel > 0) {

		//Find the bounding box of the fixed image where we keep the points

		const ImageType::RegionType& fixedImageRegion = resampledFixedImage->GetLargestPossibleRegion();
		const ImageType::SizeType& fixedImageSize = fixedImageRegion.GetSize();

		const ImageType::SpacingType fixedImageSpacing = resampledFixedImage->GetSpacing();

		const ImageType::PointType& fixedImageOrigin  = resampledFixedImage->GetOrigin();

		//Find the bounding box of the moving image where we keep the points

		const ImageType::RegionType& movingImageRegion = resampledMovingImage->GetLargestPossibleRegion();
		const ImageType::SizeType& movingImageSize = movingImageRegion.GetSize();

		const ImageType::SpacingType movingImageSpacing = resampledMovingImage->GetSpacing();

		const ImageType::PointType& movingImageOrigin  = resampledMovingImage->GetOrigin();

		//Calculate the minimum and maximum bounding box area to remove points from
		float fixedImageBoundingBoxMin[Dimension];
		float fixedImageBoundingBoxMax[Dimension];

		float movingImageBoundingBoxMin[Dimension];
		float movingImageBoundingBoxMax[Dimension];

		for (int i = 0; i < Dimension; i++) {
			//fixedImageBoundingBoxMin[i] = fixedImageSize[i]*fixedImageSpacing[i]*(removeEdgePointsPercent/100) + fixedImageOrigin[i];
			//fixedImageBoundingBoxMax[i] = fixedImageSize[i]*fixedImageSpacing[i]*(1-(removeEdgePointsPercent/100)) + fixedImageOrigin[i];

			//movingImageBoundingBoxMin[i] = movingImageSize[i]*movingImageSpacing[i]*(removeEdgePointsPercent/100) + movingImageOrigin[i];
			//movingImageBoundingBoxMax[i] = movingImageSize[i]*movingImageSpacing[i]*(1-(removeEdgePointsPercent/100)) + movingImageOrigin[i];

			fixedImageBoundingBoxMin[i] = fixedImageSpacing[i]*removeEdgePointsPixel + fixedImageOrigin[i];
			fixedImageBoundingBoxMax[i] = fixedImageSpacing[i]*(fixedImageSize[i]-removeEdgePointsPixel) + fixedImageOrigin[i];

			movingImageBoundingBoxMin[i] = movingImageSpacing[i]*removeEdgePointsPixel + movingImageOrigin[i];
			movingImageBoundingBoxMax[i] = movingImageSpacing[i]*(movingImageSize[i]-removeEdgePointsPixel) + movingImageOrigin[i];
		}

		//Make new list of points that are inside the minimum and maximum of this bounding box


		//Find if point in fixedImageKeypoints is inside the boundary box. If it is, copy over point and point data
		SiftFilterType::PointSetTypePointer tmpFixedImageKeyPoints;
		tmpFixedImageKeyPoints = SiftFilterType::PointSetType::New();
		int tmpFixedImageKeyPointsCount = 0;

		//Find the point locations inside fixedImageKeyPointsContainer
		for (int i = 0; i < fixedImageKeyPoints->GetNumberOfPoints(); i++) {
			SiftFilterType::PointType tmpp;
			fixedImageKeyPoints->GetPoint(i,&tmpp);

			SiftFilterType::FeatureType tmpf;
			fixedImageKeyPoints->GetPointData(i,&tmpf);

			//Check if point fits inside the bounding box
			bool pointInsideBoundingBox = true;
			for (int j = 0; j < Dimension; j++) {
				if (fixedImageBoundingBoxMin[j] > (float)tmpp[j] || (float)tmpp[j] > fixedImageBoundingBoxMax[j]) {
					pointInsideBoundingBox = false;
				}
			}

			//If point is inside bounding box, then add it to new list
			if (pointInsideBoundingBox) {
				tmpFixedImageKeyPoints->SetPoint(tmpFixedImageKeyPointsCount,tmpp);
				tmpFixedImageKeyPoints->SetPointData(tmpFixedImageKeyPointsCount,tmpf);
				tmpFixedImageKeyPointsCount++;
			}
		}

		fixedImageKeyPoints = tmpFixedImageKeyPoints;
		
		//Find if point in movingImageKeypoints is inside the boundary box. If it is, copy over point and point data
		SiftFilterType::PointSetTypePointer tmpMovingImageKeyPoints;
		tmpMovingImageKeyPoints = SiftFilterType::PointSetType::New();
		int tmpMovingImageKeyPointsCount = 0;

		//Find the point locations inside fixedImageKeyPointsContainer
		for (int i = 0; i < movingImageKeyPoints->GetNumberOfPoints(); i++) {
			SiftFilterType::PointType tmpp;
			movingImageKeyPoints->GetPoint(i,&tmpp);

			SiftFilterType::FeatureType tmpf;
			movingImageKeyPoints->GetPointData(i,&tmpf);

			//Check if point fits inside the bounding box
			bool pointInsideBoundingBox = true;
			for (int j = 0; j < Dimension; j++) {
				if (movingImageBoundingBoxMin[j] > (float)tmpp[j] || (float)tmpp[j] > movingImageBoundingBoxMax[j]) {
					pointInsideBoundingBox = false;
				}
			}

			//If point is inside bounding box, then add it to new list
			if (pointInsideBoundingBox) {
				tmpMovingImageKeyPoints->SetPoint(tmpMovingImageKeyPointsCount,tmpp);
				tmpMovingImageKeyPoints->SetPointData(tmpMovingImageKeyPointsCount,tmpf);
				tmpMovingImageKeyPointsCount++;
			}
		}

		movingImageKeyPoints = tmpMovingImageKeyPoints;
	}

#ifdef TIMER
	remove_edge_points_duration = (( std::clock() - start ) / (double) CLOCKS_PER_SEC) - total_duration;
	total_duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
#endif


	//Setup for getting the matched keypoint lists
	typedef itk::Array<float> FeatureType;
	typedef double CoordinateRepType;
	typedef itk::PointSet< FeatureType, Dimension,
	itk::DefaultStaticMeshTraits< FeatureType, Dimension, Dimension, CoordinateRepType > > PointSetType;
					  
	typedef PointSetType::PointType PointType;
	typedef PointSetType::Pointer PointSetTypePointer;

	PointSetTypePointer fixedImageMatchedPoints;
	PointSetTypePointer movingImageMatchedPoints;

	//Find matched SIFT feature points
	//SiftFilterType * siftFilter3 = new SiftFilterType;
	//siftFilter1->MatchKeypointsFeatures(&fixedImageKeyPoints,&movingImageKeyPoints, &fixedImageMatchedPoints, &movingImageMatchedPoints, NULL);
	//delete siftFilter3;

	if (matchSameScale||matchSameExtrema) {

		siftFilter1->MatchKeypointsFeaturesWithExtraData(&fixedImageKeyPoints,&movingImageKeyPoints, &fixedImageMatchedPoints, &movingImageMatchedPoints, NULL, matchSameScale, matchSameExtrema);
		
	} else {
		
		siftFilter1->MatchKeypointsFeatures(&fixedImageKeyPoints,&movingImageKeyPoints, &fixedImageMatchedPoints, &movingImageMatchedPoints, NULL);
		
	}
#ifdef TIMER
	match_points_duration = (( std::clock() - start ) / (double) CLOCKS_PER_SEC) - total_duration;
	total_duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
#endif

	//Remove feature point matches that produce vectors outside the boundary box
	if (argc == 28 + 2*Dimension) {

		SiftFilterType::PointSetTypePointer tmpFixedImageMatchedPoints = SiftFilterType::PointSetType::New();
		SiftFilterType::PointSetTypePointer tmpMovingImageMatchedPoints = SiftFilterType::PointSetType::New();
		int tmpMatchedPointsCount = 0;

		for (int i = 0; i < fixedImageMatchedPoints->GetNumberOfPoints(); i++) {

			SiftFilterType::PointType tmpp1;
			SiftFilterType::FeatureType tmpf1;
			fixedImageMatchedPoints->GetPoint(i,&tmpp1);
			fixedImageMatchedPoints->GetPointData(i,&tmpf1);

			SiftFilterType::PointType tmpp2;
			SiftFilterType::FeatureType tmpf2;
			movingImageMatchedPoints->GetPoint(i,&tmpp2);
			movingImageMatchedPoints->GetPointData(i,&tmpf2);

			SiftFilterType::PointType p2p1 = tmpp2 - tmpp1;

			bool exceedLimit = false;

			for (int dim = 0; dim < Dimension; dim++) {
				if (p2p1[dim]<=displacementLimitMin[dim] || p2p1[dim]>=displacementLimitMax[dim]) {
					exceedLimit = true;
				}
			}
			
			if (!exceedLimit) {
				
				tmpFixedImageMatchedPoints->SetPoint(tmpMatchedPointsCount,tmpp1);
				tmpFixedImageMatchedPoints->SetPointData(tmpMatchedPointsCount,tmpf1);
				tmpMovingImageMatchedPoints->SetPoint(tmpMatchedPointsCount,tmpp2);
				tmpMovingImageMatchedPoints->SetPointData(tmpMatchedPointsCount,tmpf2);
				
				tmpMatchedPointsCount++;

			}
		
		}

		fixedImageMatchedPoints = tmpFixedImageMatchedPoints;
		movingImageMatchedPoints = tmpMovingImageMatchedPoints;

	}






#ifdef TIMER
	TPS_duration = (( std::clock() - start ) / (double) CLOCKS_PER_SEC) - total_duration;
	total_duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
#endif


#ifdef TIMER

	std::cout<<"\n\nResample Duration"<<std::endl;
	std::cout<<resample_duration<<"\n\n"<<std::endl;

	std::cout<<"\n\nIntensity Rescale Duration"<<std::endl;
	std::cout<<intensity_rescale_duration<<"\n\n"<<std::endl;

	std::cout<<"\n\nFind Points Duration"<<std::endl;
	std::cout<<find_points_duration<<"\n\n"<<std::endl;

	std::cout<<"\n\nRemove Edge Points Duration"<<std::endl;
	std::cout<<remove_edge_points_duration<<"\n\n"<<std::endl;
	
	std::cout<<"\n\nMatch Points Duration"<<std::endl;
	std::cout<<match_points_duration<<"\n\n"<<std::endl;
	
	std::cout<<"\n\nTPS Duration"<<std::endl;
	std::cout<<TPS_duration<<"\n\n"<<std::endl;

	std::cout<<"\n\nTotal Duration"<<std::endl;
	std::cout<<total_duration<<"\n\n"<<std::endl;
#endif



	//Print out the keypoints and matchedpoints files
	KeypointsIO<FeatureType, Dimension> outputKeypoints1;
	outputKeypoints1.setKeypointSetPtr(&fixedImageKeyPoints);
	outputKeypoints1.setFileName(outName + "_fixedImageKeyPoints.txt");
	outputKeypoints1.writeKeypointSet();

	outputKeypoints1.setFileName(outName + "_fixedImageKeyPoints.landmarkASCII");
	outputKeypoints1.writeAmiraKeypointSet();

	outputKeypoints1.setFileName(outName + "_fixedImageKeyPoints.m");
	outputKeypoints1.writeMatlabKeypointSet("SIFT_fixed_points");

	KeypointsIO<FeatureType, Dimension> outputKeypoints2;
	outputKeypoints2.setKeypointSetPtr(&movingImageKeyPoints);
	outputKeypoints2.setFileName(outName + "_movingImageKeyPoints.txt");
	outputKeypoints2.writeKeypointSet();

	outputKeypoints2.setFileName(outName + "_movingImageKeyPoints.landmarkASCII");
	outputKeypoints2.writeAmiraKeypointSet();

	outputKeypoints2.setFileName(outName + "_movingImageKeyPoints.m");
	outputKeypoints2.writeMatlabKeypointSet("SIFT_moving_points");

	KeypointsIO<FeatureType, Dimension> outputKeypoints3;
	outputKeypoints3.setKeypointSetPtr(&fixedImageMatchedPoints);
	outputKeypoints3.setFileName(outName + "_fixedImageMatchedPoints.txt");
	outputKeypoints3.writeKeypointSet();

	KeypointsIO<FeatureType, Dimension> outputKeypoints4;
	outputKeypoints4.setKeypointSetPtr(&movingImageMatchedPoints);
	outputKeypoints4.setFileName(outName + "_movingImageMatchedPoints.txt");
	outputKeypoints4.writeKeypointSet();

	std::cout<<"Supposedly, keypoint lists is written properly."<<std::endl;

	//Want to print out distances as well. Need some fancy footwork here...
	MatchedPoints<FeatureType, Dimension> outputDistance;
	outputDistance.setMatchedPoints(&outputKeypoints3, &outputKeypoints4);
	outputDistance.writeDistance(outName + "_matchedPointsDistance.txt");

	std::cout<<"Supposedly, distance between matched points is written properly."<<std::endl;

	//Print out Amira LandmarksASCII file for matched points
	outputDistance.writeAmiraMatchedPoints(outName + ".landmarkASCII");
	outputDistance.writeMatlabMatchedPoints(outName + ".m");

	delete [] displacementLimitMin;
	delete [] displacementLimitMax;

	delete [] inDisplacementLimitMin;
	delete [] inDisplacementLimitMax;

	return EXIT_SUCCESS;

}