/*
File:			mainStandardResampleScaleControl.cxx

Description:	Subfunction called from imageioraw.cxx Performs SIFT feature registration on 2 input images. Input parameters based on input processing from calling mainStandardResampleScaleControl as a command line function, ignoring the first 3 input parameters. Refer below for input process parameters.

Inputs:			argc - (int) Number of input arguments stored in argv. Default number is 21.
				argv[2]: inputImageFile1 - (std::string) Path name of the raw data file of the fixed image. Pixels are in 32 bit float format.
				argv[3]: inputImageFile2 - (std::string) Path name of the raw data file of the moving image. Pixels are in 32 bit float format.
				argv[4]: inputConfigFile1 - (std::string) Path name of the configuration text file of the fixed image. Pixels are in 32 bit float format.
				argv[5]: inputConfigFile2 - (std::string)  Path name of the configuration text file of the moving image. Pixels are in 32 bit float format. 
				
				argv[6]: inImageFile1Min - (std::string) Minimum intensity file of the fixed image. Pixels are in 32 bit float format. 
				argv[7]: inImageFile1Max - (std::string) Maximum intensity file of the fixed image. Pixels are in 32 bit float format. 
				argv[8]: inImageFile2Min - (std::string) Minimum intensity file of the moving image. Pixels are in 32 bit float format. 
				argv[9]: inImageFile2Max - (std::string) Maximum intensity file of the moving image. Pixels are in 32 bit float format. 
				
				argv[10]: outPointsFile - (std::string) Path name and file name prefix of the output files containing information on SIFT feature point positions
				argv[11]: outputDisplacementField - (std::string) Path name and file name prefix of the output displacement field from SIFT + TPS in raw data format
				argv[12]: outDisplacementConfig - (std::string) Path name of the configuration text file of the displacement field from SIFT + TPS
				argv[13]: resampleFactorStr - (std::string) Downsample factor of the input images before using SIFT. Input is a float number in a string format.
				argv[14]: relativePosition - (std::string) Determines where the origin of the output image/field files are placed. They can either be placed at the "origin" or at the "center" of the image/field. Input is a string containing "origin" or "center.
				argv[15]: outputResolution - (std::string) Downsample factor of the output images/fields after using SIFT. Resample based on original input image resolution. Input is a float number in a string format.
				argv[16]: inRemoveEdgePointsPercent - (std::string) Percentage of image removed starting from the image edge and moving inward based on toral image size in each dimension. Input is a float number [0,100] in string format.
				
				argv[17]: inScaleStart - (std::string)
				argv[18]: inScaleEnd - (std::string)

				argv[19]: inEdgePointRatio - (std::string) Starting edge point ratio value for scale level 0. The higher the number, the more feature points will be accepted for matching. Input is a float number greater than 0. Input can also be a negative number if you want to disregard edge point checking for feature points. If you are unsure about this value, leave it at 10.
				argv[20]: inEdgePointRatioScale - (std::string) Edge point ratio multiplier that multiplies the edge point ratio by this number each time SIFT searches up a higher scale. Input is a float number. Default is 1 for no adjustment in the edge point ratio between scales.

Author:			Frankie (Hoi-Ki) Tong, May 13 2015

Changes:		Frankie (Hoi-Ki) Tong, May 13 2015
				Initial creation of the file.
*/

#include "mainStandardResampleScaleControl.h"

bool mainStandardResampleScaleControl (int argc, char *argv[]) {

	std::string inputImageFile1 = "";
	std::string inputImageFile2 = "";
	std::string inConfigFile1 = "";
	std::string inConfigFile2 = "";

	std::string inImageFile1Min = "";
	std::string inImageFile1Max = "";
	std::string inImageFile2Min = "";
	std::string inImageFile2Max = "";

	std::string outPointsFile = "";
	std::string outputDisplacementField = "";
	std::string outDisplacementConfig = "";
	std::string resampleFactorStr = "";
	std::string relativePosition = "";
	std::string outputResolution = "";
	std::string inRemoveEdgePointsPercent = "";
	
	std::string inScaleStart = "";
	std::string inScaleEnd = "";

	std::string inEdgePointRatio = "";
	std::string inEdgePointRatioScale = "";

	if (argc == 21) {

		inputImageFile1 = argv[2];
		inputImageFile2 = argv[3];
		inConfigFile1 = argv[4];
		inConfigFile2 = argv[5];
		inImageFile1Min = argv[6];
		inImageFile1Max = argv[7];
		inImageFile2Min = argv[8];
		inImageFile2Max = argv[9];
		outPointsFile = argv[10];
		outputDisplacementField = argv[11];
		outDisplacementConfig = argv[12];
		resampleFactorStr = argv[13];
		relativePosition = argv[14];
		outputResolution = argv[15];
		inRemoveEdgePointsPercent = argv[16];
		inScaleStart = argv[17];
		inScaleEnd = argv[18];
		inEdgePointRatio = argv[19];
		inEdgePointRatioScale = argv[20];
		
	} else
	{
		std::cerr << "Incorrect number of parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " --standardResampleScaleControl inputImageFile1 inputImageFile2 inConfigFile1 inConfigFile2 inImageFile1Min inImageFile1Max inImageFile2Min inImageFile2Max outPointsFile outputDisplacementField outDisplacementConfig resampleSpacingFactor relativePosition(origin/center) outputDisplacementResolution(original/resampled/both/#) inRemoveEdgePointsPercent inScaleStart inScaleEnd inEdgePointRatio inEdgePointRatioScale"<< std::endl;
		return EXIT_FAILURE;
	}

	float resampleFactor = 0; //defaults
	std::istringstream(resampleFactorStr) >> resampleFactor;
	
	//Get percentage of the image where points found in the area will be removed
	float removeEdgePointsPercent = 0; //defaults
	std::istringstream(inRemoveEdgePointsPercent) >>  removeEdgePointsPercent;
	
	int scaleStart = 0; //defaults
	std::istringstream(inScaleStart) >>  scaleStart;
	
	int scaleEnd = 3; //defaults
	std::istringstream(inScaleEnd) >>  scaleEnd;

	float edgePointRatio = 10; //defaults
	std::istringstream(inEdgePointRatio) >>  edgePointRatio;

	float edgePointRatioScale = 1; //defaults
	std::istringstream(inEdgePointRatioScale) >>  edgePointRatioScale;

	//Initialize Imagetypes and image pointers
	typedef float									PixelType;
	const unsigned int Dimension = 2;
	typedef itk::Image<PixelType, Dimension>		ImageType;
	
	ImageType::Pointer fixedImage = ImageType::New();
	ImageType::Pointer movingImage = ImageType::New();
	
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
	ImageType::Pointer resampledMovingImage = ImageType::New();
	
	if (resampleFactor == 1) {
		resampledFixedImage = fixedImage;
		resampledMovingImage = movingImage;
	} else {
		
		//Need to resample both images to new resolution
		
		//Resampleing fixedImage
		typedef itk::IdentityTransform<double, Dimension> upSampleIdentityTransformType;
		//typedef itk::LinearInterpolateImageFunction<ImageType,double> upSampleInterpolatorType;
		typedef	itk::BSplineInterpolateImageFunction< ImageType, double > upSampleInterpolatorType;
		typedef itk::ResampleImageFilter<ImageType,ImageType,double,double> upSampleResampleImaageFilterType;

		upSampleIdentityTransformType::Pointer upSampleIdentityTransformFixedImage = upSampleIdentityTransformType::New();
		upSampleIdentityTransformFixedImage->SetIdentity();

		upSampleInterpolatorType::Pointer upSampleInterpolatorFixedImage = upSampleInterpolatorType::New();

		upSampleResampleImaageFilterType::Pointer upSampleResampleImageFilterFixedImage = upSampleResampleImaageFilterType::New();
		upSampleResampleImageFilterFixedImage->SetTransform(upSampleIdentityTransformFixedImage);
		upSampleResampleImageFilterFixedImage->SetInterpolator(upSampleInterpolatorFixedImage);

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

		itk::Size<Dimension> upSampleSizeFixedImage;
		for (int i = 0; i < Dimension; i++) {
			upSampleSizeFixedImage[i] = NewSizeFixedImage[i];
		}

		upSampleResampleImageFilterFixedImage->SetSize(upSampleSizeFixedImage);

		upSampleResampleImageFilterFixedImage->SetInput(fixedImage);

		//upSampleResampleImageFilterFixedImage->ReleaseDataFlagOn();
		upSampleResampleImageFilterFixedImage->Update();

		resampledFixedImage = upSampleResampleImageFilterFixedImage->GetOutput();
		
		//upSampleResampleImageFilterFixedImage = NULL;
			
		
		
		//Resampleing movingImage

		upSampleIdentityTransformType::Pointer upSampleIdentityTransformMovingImage = upSampleIdentityTransformType::New();
		upSampleIdentityTransformMovingImage->SetIdentity();

		upSampleInterpolatorType::Pointer upSampleInterpolatorMovingImage = upSampleInterpolatorType::New();

		upSampleResampleImaageFilterType::Pointer upSampleResampleImageFilterMovingImage = upSampleResampleImaageFilterType::New();
		upSampleResampleImageFilterMovingImage->SetTransform(upSampleIdentityTransformMovingImage);
		upSampleResampleImageFilterMovingImage->SetInterpolator(upSampleInterpolatorMovingImage);

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

		itk::Size<Dimension> upSampleSizeMovingImage;
		for (int i = 0; i < Dimension; i++) {
			upSampleSizeMovingImage[i] = NewSizeMovingImage[i];
		}

		upSampleResampleImageFilterMovingImage->SetSize(upSampleSizeMovingImage);

		upSampleResampleImageFilterMovingImage->SetInput(movingImage);

		//upSampleResampleImageFilterMovingImage->ReleaseDataFlagOn();
		upSampleResampleImageFilterMovingImage->Update();

		resampledMovingImage = upSampleResampleImageFilterMovingImage->GetOutput();

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

	typedef itk::IntensityWindowingImageFilter <ImageType, ImageType> IntensityWindowingImageFilterType;
	IntensityWindowingImageFilterType::Pointer IntenistyWindowImageFilterFixedImage = IntensityWindowingImageFilterType::New();

	IntenistyWindowImageFilterFixedImage->SetInput(resampledFixedImage);
	IntenistyWindowImageFilterFixedImage->SetWindowMinimum(imageFile1Min);
	IntenistyWindowImageFilterFixedImage->SetWindowMaximum(imageFile1Max);
	IntenistyWindowImageFilterFixedImage->SetOutputMinimum(0);
	IntenistyWindowImageFilterFixedImage->SetOutputMaximum(1);
	IntenistyWindowImageFilterFixedImage->Update();
	resampledFixedImage = IntenistyWindowImageFilterFixedImage->GetOutput();

	IntensityWindowingImageFilterType::Pointer IntenistyWindowImageFilterMovingImage = IntensityWindowingImageFilterType::New();
	IntenistyWindowImageFilterMovingImage->SetInput(resampledMovingImage);
	IntenistyWindowImageFilterMovingImage->SetWindowMinimum(imageFile2Min);
	IntenistyWindowImageFilterMovingImage->SetWindowMaximum(imageFile2Max);
	IntenistyWindowImageFilterMovingImage->SetOutputMinimum(0);
	IntenistyWindowImageFilterMovingImage->SetOutputMaximum(1);
	IntenistyWindowImageFilterMovingImage->Update();
	resampledMovingImage = IntenistyWindowImageFilterMovingImage->GetOutput();

	//std::cout<<std::endl<<"Resampled and intensity adjusted fixed image stats:"<<std::endl;
	//ImageStats<ImageType>(resampledFixedImage);
	//std::cout<<std::endl<<"Resampled and intensity adjusted moving image stats:"<<std::endl;
	//ImageStats<ImageType>(resampledMovingImage);

#ifdef TIMER
	intensity_rescale_duration = (( std::clock() - start ) / (double) CLOCKS_PER_SEC) - total_duration;
	total_duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
#endif


	//Find Keypoints and Matchedpoints betweent the 2 images
	typedef itk::ScaleInvariantFeatureImageFilter<ImageType, Dimension> SiftFilterType;

	SiftFilterType::PointSetTypePointer fixedImageKeyPoints, movingImageKeyPoints;

	//Get Keypoints from fixed image and moving image
	SiftFilterType * siftFilter1 = new SiftFilterType;
	fixedImageKeyPoints = siftFilter1->getSiftFeaturesRescaling(resampledFixedImage, scaleStart, scaleEnd, edgePointRatio, edgePointRatioScale);
	//delete siftFilter1;

	//SiftFilterType * siftFilter2 = new SiftFilterType;
	movingImageKeyPoints = siftFilter1->getSiftFeaturesRescaling(resampledMovingImage, scaleStart, scaleEnd, edgePointRatio, edgePointRatioScale);
	//delete siftFilter2;

#ifdef TIMER
	find_points_duration = (( std::clock() - start ) / (double) CLOCKS_PER_SEC) - total_duration;
	total_duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
#endif

	//Find and remove edge points if requested by user
	if (removeEdgePointsPercent > 0) {

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
			fixedImageBoundingBoxMin[i] = fixedImageSize[i]*fixedImageSpacing[i]*(removeEdgePointsPercent/100) + fixedImageOrigin[i];
			fixedImageBoundingBoxMax[i] = fixedImageSize[i]*fixedImageSpacing[i]*(1-(removeEdgePointsPercent/100)) + fixedImageOrigin[i];

			movingImageBoundingBoxMin[i] = movingImageSize[i]*movingImageSpacing[i]*(removeEdgePointsPercent/100) + movingImageOrigin[i];
			movingImageBoundingBoxMax[i] = movingImageSize[i]*movingImageSpacing[i]*(1-(removeEdgePointsPercent/100)) + movingImageOrigin[i];
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

	siftFilter1->MatchKeypointsFeatures(&fixedImageKeyPoints,&movingImageKeyPoints, &fixedImageMatchedPoints, &movingImageMatchedPoints, NULL);

#ifdef TIMER
	match_points_duration = (( std::clock() - start ) / (double) CLOCKS_PER_SEC) - total_duration;
	total_duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
#endif


	std::cout<<"Finding transform from keypoints"<<std::endl;

	//Find transform with keypoints goes here...
	typedef itk::ThinPlateSplineKernelTransform< CoordinateRepType, Dimension> TransformType;

	//Need to convert point set type so transform will accept the points
	//fixed -> source, moving -> target

	//Setup to get points container for SIFT matched points
	typedef PointSetType::PointsContainer::Pointer PointSetTypeContainerPointer;
	typedef PointSetType::PointIdentifier PointIdType;
	PointIdType ID = itk::NumericTraits<PointIdType>::ZeroValue();

	PointSetTypeContainerPointer fixedImageMatchedPointsContainer = fixedImageMatchedPoints->GetPoints();
	PointSetTypeContainerPointer movingImageMatchedPointsContainer = movingImageMatchedPoints->GetPoints();
	
	//Setup to get points container for ITK TPS points container
	typedef itk::Point<CoordinateRepType, Dimension> TransformPointType;
	typedef TransformType::PointSetType TransformPointSetType;
	typedef TransformPointSetType::PointsContainer::Pointer TransformPointSetTypeContainerPointer;
	typedef TransformPointSetType::PointIdentifier TransformPointIdType;
	TransformPointIdType TransformID = itk::NumericTraits<TransformPointIdType>::ZeroValue();

	TransformPointSetType::Pointer sourceLandMarks = TransformPointSetType::New();
	TransformPointSetType::Pointer targetLandMarks = TransformPointSetType::New();

	TransformPointSetTypeContainerPointer sourceLandMarkContainer = sourceLandMarks->GetPoints();
	TransformPointSetTypeContainerPointer targetLandMarkContainer = targetLandMarks->GetPoints();

	//Convert each SIFT found matched point to ITK TPS input points
	for (int i = 0; i < fixedImageMatchedPoints->GetNumberOfPoints(); i++) {

		itk::Point<CoordinateRepType, Dimension> tmpp1;
		itk::Point<CoordinateRepType, Dimension> tmpp2;
		fixedImageMatchedPointsContainer->GetElementIfIndexExists(ID,&tmpp1);
		movingImageMatchedPointsContainer->GetElementIfIndexExists(ID, &tmpp2);

		sourceLandMarkContainer->InsertElement( TransformID, tmpp1 );
		targetLandMarkContainer->InsertElement( TransformID, tmpp2 );

		ID++;	
		TransformID++;

	}

	//Create new TPS transform and input points to the transform 
	TransformType::Pointer transform = TransformType::New();
	transform->SetSourceLandmarks(sourceLandMarks);
	transform->SetTargetLandmarks(targetLandMarks);

	//transform->SetStiffness(0.1);

	transform->ComputeWMatrix();

	std::cout<<"Image transform has been calculated"<<std::endl;




	//Create displacement field here...

	std::cout<<"Making displacement fields"<<std::endl;

	//Make displacement field based off of original image resolution
	//typedef itk::Vector< PixelType, Dimension > VectorPixelType;  
	//typedef itk::Image< VectorPixelType, Dimension > DisplacementFieldImageType;  
	//typedef itk::TransformToDisplacementFieldFilter< DisplacementFieldImageType, CoordinateRepType > DisplacementFieldGeneratorType;

	if (outputResolution == "original" || outputResolution == "both"){

		//Make displacement field based off of original image resolution
		typedef itk::Vector< PixelType, Dimension > VectorPixelType;  
		typedef itk::Image< VectorPixelType, Dimension > DisplacementFieldImageType;  
		typedef itk::TransformToDisplacementFieldFilter< DisplacementFieldImageType, CoordinateRepType > DisplacementFieldGeneratorType;

		DisplacementFieldGeneratorType::Pointer dispfieldGeneratorOriginal = DisplacementFieldGeneratorType::New();  
		DisplacementFieldImageType::Pointer displacementFieldOriginal = DisplacementFieldImageType::New();
	
		dispfieldGeneratorOriginal->UseReferenceImageOn(); 
		dispfieldGeneratorOriginal->SetReferenceImage( fixedImage );  //sets the input parameters for the field generated to be the same as resampledFixedImage
		dispfieldGeneratorOriginal->SetTransform( transform ); 

		displacementFieldOriginal->ReleaseDataFlagOn();

		try {    
			dispfieldGeneratorOriginal->Update();
		}  catch ( itk::ExceptionObject & err )    {
			std::cerr << "Exception detected while generating deformation field";
			std::cerr << " : "  << err << std::endl;    
			return EXIT_FAILURE;    
		}

		displacementFieldOriginal = dispfieldGeneratorOriginal->GetOutput();

		//dispfieldGeneratorOriginal = NULL;

		//Output Displacement Field and config file
		itkImageIORaw<VectorPixelType, Dimension>* inputImageImport3 = new itkImageIORaw<VectorPixelType, Dimension>;

		inputImageImport3->setConfigFileName(outDisplacementConfig + "_original.txt");

		//Write ConfigFile
		ImageIOData outputDisplacementImageIOData;

		outputDisplacementImageIOData.setPixelType("float");
		outputDisplacementImageIOData.setPixelDimensionality(Dimension);
		outputDisplacementImageIOData.setImageDimensionality(Dimension);

		const DisplacementFieldImageType::RegionType& OutputDisplacementRegion = displacementFieldOriginal->GetLargestPossibleRegion();
		const DisplacementFieldImageType::SizeType& OutputDisplacementSize = OutputDisplacementRegion.GetSize();
		unsigned int tmpDisplacementSize[Dimension];
		for (int i = 0; i < Dimension; i++) {
			tmpDisplacementSize[i] = OutputDisplacementSize[i];
		}

		outputDisplacementImageIOData.setDimensionSize((int *)tmpDisplacementSize);

		const DisplacementFieldImageType::SpacingType OutputDisplacementSpacing = displacementFieldOriginal->GetSpacing();
		float tmpDisplacementSpacing[Dimension];
		for (int i = 0; i < Dimension; i++) {
			tmpDisplacementSpacing[i] = (float)OutputDisplacementSpacing[i];
		}
		
		outputDisplacementImageIOData.setSpacingSize((float *)tmpDisplacementSpacing);

		const DisplacementFieldImageType::PointType& OutputDisplacementOrigin  = displacementFieldOriginal->GetOrigin();
		float tmpDisplacementOrigin[Dimension];
		for (int i=0; i< Dimension; i++) {
			tmpDisplacementOrigin[i] = (float)OutputDisplacementOrigin[i];
		}

		outputDisplacementImageIOData.setOrigin((float *)tmpDisplacementOrigin);

		outputDisplacementImageIOData.setHeaderSize(0);
		outputDisplacementImageIOData.setByteOrder("little");

		inputImageImport3->setImageIOData(outputDisplacementImageIOData);

		inputImageImport3->writeConfigFile();

		//Write Displacement Field
		inputImageImport3->setitkImage((void *) &displacementFieldOriginal);
		inputImageImport3->setImageFileName(outputDisplacementField + "_original.raw");
		inputImageImport3->writeRawImageKnown();

		delete inputImageImport3;

		//displacementFieldOriginal = NULL;

		std::cout << "Supposedly, displacement field original is written in properly." << std::endl;
	}

	//Make displacement field based off of resampled resolution
	if (outputResolution == "resampled" || outputResolution == "both") {
		
		//Make displacement field based off of original image resolution
		typedef itk::Vector< PixelType, Dimension > VectorPixelType;  
		typedef itk::Image< VectorPixelType, Dimension > DisplacementFieldImageType;  
		typedef itk::TransformToDisplacementFieldFilter< DisplacementFieldImageType, CoordinateRepType > DisplacementFieldGeneratorType;

		DisplacementFieldGeneratorType::Pointer dispfieldGeneratorResampled = DisplacementFieldGeneratorType::New(); 
		DisplacementFieldImageType::Pointer displacementFieldResampled = DisplacementFieldImageType::New();

		dispfieldGeneratorResampled->UseReferenceImageOn(); 
		dispfieldGeneratorResampled->SetReferenceImage( resampledFixedImage );  //sets the input parameters for the field generated to be the same as resampledFixedImage
		dispfieldGeneratorResampled->SetTransform( transform );  
		dispfieldGeneratorResampled->ReleaseDataFlagOn();

		try {    
			dispfieldGeneratorResampled->Update();
		}  catch ( itk::ExceptionObject & err )    {
			std::cerr << "Exception detected while generating deformation field";
			std::cerr << " : "  << err << std::endl;    
			return EXIT_FAILURE;    
		}

		displacementFieldResampled = dispfieldGeneratorResampled->GetOutput();

		//dispfieldGeneratorResampled = NULL;

		//Output Displacement Field and config file

		itkImageIORaw<VectorPixelType, Dimension>* inputImageImport4 = new itkImageIORaw<VectorPixelType, Dimension>;

		inputImageImport4->setConfigFileName(outDisplacementConfig + "_resampled.txt");

		//Write ConfigFile
		ImageIOData outputDisplacementImageIOData;

		outputDisplacementImageIOData.setPixelType("float");
		outputDisplacementImageIOData.setPixelDimensionality(Dimension);
		outputDisplacementImageIOData.setImageDimensionality(Dimension);

		const DisplacementFieldImageType::RegionType& OutputDisplacementRegion = displacementFieldResampled->GetLargestPossibleRegion();
		const DisplacementFieldImageType::SizeType& OutputDisplacementSize = OutputDisplacementRegion.GetSize();
		unsigned int tmpDisplacementSize[Dimension];
		for (int i = 0; i < Dimension; i++) {
			tmpDisplacementSize[i] = OutputDisplacementSize[i];
		}

		outputDisplacementImageIOData.setDimensionSize((int *)tmpDisplacementSize);

		const DisplacementFieldImageType::SpacingType OutputDisplacementSpacing = displacementFieldResampled->GetSpacing();
		float tmpDisplacementSpacing[Dimension];
		for (int i = 0; i < Dimension; i++) {
			tmpDisplacementSpacing[i] = (float)OutputDisplacementSpacing[i];
		}
		
		outputDisplacementImageIOData.setSpacingSize((float *)tmpDisplacementSpacing);

		const DisplacementFieldImageType::PointType& OutputDisplacementOrigin  = displacementFieldResampled->GetOrigin();
		float tmpDisplacementOrigin[Dimension];
		for (int i=0; i< Dimension; i++) {
			tmpDisplacementOrigin[i] = (float)OutputDisplacementOrigin[i];
		}

		outputDisplacementImageIOData.setOrigin((float *)tmpDisplacementOrigin);

		outputDisplacementImageIOData.setHeaderSize(0);
		outputDisplacementImageIOData.setByteOrder("little");

		inputImageImport4->setImageIOData(outputDisplacementImageIOData);

		inputImageImport4->writeConfigFile();

		//Write Displacement Field
		inputImageImport4->setitkImage((void *) &displacementFieldResampled);
		inputImageImport4->setImageFileName(outputDisplacementField + "_resampled.raw");
		inputImageImport4->writeRawImageKnown();

		delete inputImageImport4;
		
		//displacementFieldResampled = NULL;

		std::cout << "Supposedly, displacement field resampled is written in properly." << std::endl;
		
	}

	if (outputResolution != "original" && outputResolution != "resampled" && outputResolution != "both") {

		//Check if input is a floating number. If it is, then resample the output displacement field to that resample factor
		std::istringstream iss(outputResolution);
		float displacementResampleFactor = 1;
		
		if (!(iss >> displacementResampleFactor)) {

			std::cerr << "Value \"" << outputResolution << "\" for output resolution is unknown. Please use either \"original\" or \"resampled\" or \"both\" or a resampleFactor Number (float)." << std::endl;
			return EXIT_FAILURE;

		} else {

			//Make displacement field based off of original image resolution
			typedef itk::Vector< PixelType, Dimension > VectorPixelType;  
			typedef itk::Image< VectorPixelType, Dimension > DisplacementFieldImageType;  
			typedef itk::TransformToDisplacementFieldFilter< DisplacementFieldImageType, CoordinateRepType > DisplacementFieldGeneratorType;

			//We have a specified resampleFactor that we would like th displacement field to be at. Resize region to allow for this.
			DisplacementFieldGeneratorType::Pointer dispfieldGenerator = DisplacementFieldGeneratorType::New(); 
			DisplacementFieldImageType::Pointer displacementField = DisplacementFieldImageType::New();


			//Make a new image of the same size of parameters of the fixedImage, but downsample it to the size I want

			const ImageType::RegionType& DisplacementInputRegionFixedImage = fixedImage->GetLargestPossibleRegion();
			const ImageType::SizeType& DisplacementInputSizeFixedImage = DisplacementInputRegionFixedImage.GetSize();
			unsigned int OldSizeImage[Dimension];
			unsigned int NewSizeImage[Dimension];
			for (int i = 0; i < Dimension; i++) {
				OldSizeImage[i] = DisplacementInputSizeFixedImage[i];
				NewSizeImage[i] = std::floor((double) OldSizeImage[i]/displacementResampleFactor );
			}

			const ImageType::SpacingType OldSpacingImage = fixedImage->GetSpacing();
			double NewSpacingImage[Dimension];
			for (int i = 0; i < Dimension; i++) {
				NewSpacingImage[i] = OldSpacingImage[i] * (double) OldSizeImage[i]/NewSizeImage[i];
			}
			//dispfieldGenerator->SetOutputSpacing(NewSpacingImage);
			
			itk::Size<Dimension> upSampleSizeImage;
			for (int i = 0; i < Dimension; i++) {
				upSampleSizeImage[i] = NewSizeImage[i];
			}

			//dispfieldGenerator->SetSize(upSampleSizeImage);

			const ImageType::PointType& OldOriginImage  = fixedImage->GetOrigin();
			double ResampleOriginImage[Dimension];
			for (int i=0; i< Dimension; i++) {

				if (relativePosition == "origin") {
					ResampleOriginImage[i] = OldOriginImage[i]-(OldSpacingImage[i] - NewSpacingImage[i])/2.0;
				} else if (relativePosition == "center") {
					ResampleOriginImage[i] = OldOriginImage[i] - ((OldSpacingImage[i] * (OldSizeImage[i] - 1)) - (NewSpacingImage[i] * (NewSizeImage[i] - 1)))/2.0;
				} else {
				
					std::cerr << "Value \"" << relativePosition << "\" for relative position is unknown. Please use either \"origin\" or \"center\"." << std::endl;
					return EXIT_FAILURE;
				}
				
			}
			//dispfieldGenerator->SetOutputOrigin(ResampleOriginImage);

			
			ImageType::Pointer tmpImage = ImageType::New();
			ImageType::RegionType tmpRegion;

			tmpRegion.SetSize(upSampleSizeImage);
			
			ImageType::IndexType upSampleIndexImage;
			for (int i = 0; i < Dimension; i++) {
				upSampleIndexImage[i] = 0;
			}

			tmpRegion.SetIndex(upSampleIndexImage);

			tmpImage->SetRegions(tmpRegion);
			tmpImage->SetSpacing(NewSpacingImage);
			tmpImage->SetOrigin(ResampleOriginImage);
			tmpImage->SetDirection(fixedImage->GetDirection());
			//tmpImage->Allocate();

			dispfieldGenerator->UseReferenceImageOn(); 
			dispfieldGenerator->SetReferenceImage( tmpImage );
			dispfieldGenerator->SetTransform( transform );
			dispfieldGenerator->ReleaseDataFlagOn();

			try {    
				dispfieldGenerator->Update();
			}  catch ( itk::ExceptionObject & err )    {
				std::cerr << "Exception detected while generating deformation field";
				std::cerr << " : "  << err << std::endl;    
				return EXIT_FAILURE;    
			}

			displacementField = dispfieldGenerator->GetOutput();

			//dispfieldGenerator = NULL;
			

			//Output Displacement Field and config file

			itkImageIORaw<VectorPixelType, Dimension>* inputImageImport5 = new itkImageIORaw<VectorPixelType, Dimension>;

			inputImageImport5->setConfigFileName(outDisplacementConfig + ".txt");

			//Write ConfigFile
			ImageIOData outputDisplacementImageIOData;

			outputDisplacementImageIOData.setPixelType("float");
			outputDisplacementImageIOData.setPixelDimensionality(Dimension);
			outputDisplacementImageIOData.setImageDimensionality(Dimension);

			const DisplacementFieldImageType::RegionType& OutputDisplacementRegion = displacementField->GetLargestPossibleRegion();
			const DisplacementFieldImageType::SizeType& OutputDisplacementSize = OutputDisplacementRegion.GetSize();
			unsigned int tmpDisplacementSize[Dimension];
			for (int i = 0; i < Dimension; i++) {
				tmpDisplacementSize[i] = OutputDisplacementSize[i];
			}

			outputDisplacementImageIOData.setDimensionSize((int *)tmpDisplacementSize);

			const DisplacementFieldImageType::SpacingType OutputDisplacementSpacing = displacementField->GetSpacing();
			float tmpDisplacementSpacing[Dimension];
			for (int i = 0; i < Dimension; i++) {
				tmpDisplacementSpacing[i] = (float)OutputDisplacementSpacing[i];
			}
			
			outputDisplacementImageIOData.setSpacingSize((float *)tmpDisplacementSpacing);

			const DisplacementFieldImageType::PointType& OutputDisplacementOrigin  = displacementField->GetOrigin();
			float tmpDisplacementOrigin[Dimension];
			for (int i=0; i< Dimension; i++) {
				tmpDisplacementOrigin[i] = (float)OutputDisplacementOrigin[i];
			}

			outputDisplacementImageIOData.setOrigin((float *)tmpDisplacementOrigin);

			outputDisplacementImageIOData.setHeaderSize(0);
			outputDisplacementImageIOData.setByteOrder("little");

			inputImageImport5->setImageIOData(outputDisplacementImageIOData);

			inputImageImport5->writeConfigFile();

			//Write Displacement Field
			inputImageImport5->setitkImage((void *) &displacementField);
			inputImageImport5->setImageFileName(outputDisplacementField + ".raw");
			inputImageImport5->writeRawImageKnown();

			delete inputImageImport5;
			
			//displacementField = NULL;
			//tmpImage = NULL;

			std::cout << "Supposedly, displacement field resampled is written in properly." << std::endl;
			
		}
	}

	
	




#ifdef TIMER
	TPS_duration = (( std::clock() - start ) / (double) CLOCKS_PER_SEC) - total_duration;
	total_duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
#endif
	
	std::cout<<"Displacement field has been generated."<<std::endl;

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


//Frankie: Commenting this out to save a bit of time at the end...
/*
	//Transform fixed image to match moving image (for show mostly)

	//Set up resample image filter with the proper parameters
	typedef itk::ResampleImageFilter<ImageType, ImageType>  FilterType;
	FilterType::Pointer filter = FilterType::New();

	const ImageType::SpacingType& spacing = fixedImage->GetSpacing();
	const ImageType::PointType& origin  = fixedImage->GetOrigin();
	const ImageType::DirectionType& direction  = fixedImage->GetDirection();
	const ImageType::SizeType size = fixedImage->GetLargestPossibleRegion().GetSize();

	filter->SetSize(size);
	filter->SetOutputSpacing(spacing);
	filter->SetOutputOrigin(origin);
	filter->SetOutputDirection(direction);
	filter->SetDefaultPixelValue( 0 );

	//Default Interpolation Filter is linear interpolate which produces suboptimal interpolation values. Changing it to B-Spline interpolation
	//typedef itk::ConstantBoundaryCondition< ImageType > BoundaryConditionType;
	//const unsigned int WindowRadius = 5;
	//typedef itk::Function::LanczosWindowFunction<WindowRadius> WindowFunctionType;
	//typedef itk::WindowedSincInterpolateImageFunction< ImageType, WindowRadius, WindowFunctionType, BoundaryConditionType, CoordinateRepType > InterpolatorType;
	
	typedef	itk::BSplineInterpolateImageFunction< ImageType, double > InterpolatorType;

	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	filter->SetInterpolator( interpolator );

	filter->SetInput(fixedImage);
	filter->SetTransform(transform);
	filter->Update();

	ImageType::Pointer fixedImageTransformed = ImageType::New();
	fixedImageTransformed = filter->GetOutput();
*/


	//Output all appropriate outputs

	//Frankie: Commenting out to save some time in the end
	/*
	//Write output image and config file
	itkImageIORaw<PixelType, Dimension>* inputImageImport2 = new itkImageIORaw<PixelType, Dimension>;

	inputImageImport2->setConfigFileName(outConfigFile);

	//Write ConfigFile
	//The way to build ImageIOData structure differs from implementation to implementation which means I can't put this in itkImageIORaw for now...
	ImageIOData outputImageIOData;

	outputImageIOData.setPixelType("float");
	outputImageIOData.setPixelDimensionality(1);
	outputImageIOData.setImageDimensionality(Dimension);

	const ImageType::RegionType& OutputRegion = fixedImageTransformed->GetLargestPossibleRegion();
	const ImageType::SizeType& OutputSize = OutputRegion.GetSize();
	unsigned int tmpSize[Dimension];
	for (int i = 0; i < Dimension; i++) {
		tmpSize[i] = OutputSize[i];
	}

	outputImageIOData.setDimensionSize((int *)tmpSize);

	const ImageType::SpacingType OutputSpacing = fixedImageTransformed->GetSpacing();
	float tmpSpacing[Dimension];
	for (int i = 0; i < Dimension; i++) {
		tmpSpacing[i] = (float)OutputSpacing[i];
	}
	
	outputImageIOData.setSpacingSize((float *)tmpSpacing);

	const ImageType::PointType& OutputOrigin  = fixedImageTransformed->GetOrigin();
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
	inputImageImport2->setitkImage((void *) &fixedImageTransformed);
	inputImageImport2->setImageFileName(outputImageFile);
	inputImageImport2->writeRawImageKnown();

	delete inputImageImport2;

	std::cout << "Supposedly, output image is written in properly." << std::endl;
	*/



	//Print out the keypoints and matchedpoints files
	KeypointsIO<FeatureType, Dimension> outputKeypoints1;
	outputKeypoints1.setKeypointSetPtr(&fixedImageKeyPoints);
	outputKeypoints1.setFileName(outPointsFile + "_fixedImageKeyPoints.txt");
	outputKeypoints1.writeKeypointSet();

	KeypointsIO<FeatureType, Dimension> outputKeypoints2;
	outputKeypoints2.setKeypointSetPtr(&movingImageKeyPoints);
	outputKeypoints2.setFileName(outPointsFile + "_movingImageKeyPoints.txt");
	outputKeypoints2.writeKeypointSet();

	KeypointsIO<FeatureType, Dimension> outputKeypoints3;
	outputKeypoints3.setKeypointSetPtr(&fixedImageMatchedPoints);
	outputKeypoints3.setFileName(outPointsFile + "_fixedImageMatchedPoints.txt");
	outputKeypoints3.writeKeypointSet();

	KeypointsIO<FeatureType, Dimension> outputKeypoints4;
	outputKeypoints4.setKeypointSetPtr(&movingImageMatchedPoints);
	outputKeypoints4.setFileName(outPointsFile + "_movingImageMatchedPoints.txt");
	outputKeypoints4.writeKeypointSet();

	std::cout<<"Supposedly, keypoint lists is written properly."<<std::endl;

	//Want to print out distances as well. Need some fancy footwork here...
	MatchedPoints<FeatureType, Dimension> outputDistance;
	outputDistance.setMatchedPoints(&outputKeypoints3, &outputKeypoints4);
	outputDistance.writeDistance(outPointsFile + "_matchedPointsDistance.txt");

	std::cout<<"Supposedly, distance between matched points is written properly."<<std::endl;

	//Print out Amira LandmarksASCII file for matched points
	outputDistance.writeAmiraMatchedPoints(outPointsFile + ".landmarkASCII");


	return EXIT_SUCCESS;

}