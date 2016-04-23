/*
File:			mainStandard.cxx

Description:	Subfunction called from imageioraw.cxx Performs SIFT feature registration on 2 input images. Input parameters based on input processing from calling mainStandard as a command line function, ignoring the first 3 input parameters. Refer below for input process parameters.

Inputs:			argc - (int) Number of input arguments stored in argv. Default number is 12.
				argv[2]: inputImageFile1 - (std::string) Path name of the raw data file of the fixed image. Pixels are in 32 bit float format.
				argv[3]: inputImageFile2 - (std::string) Path name of the raw data file of the moving image. Pixels are in 32 bit float format.
				argv[4]: inputConfigFile1 - (std::string) Path name of the configuration text file of the fixed image. Pixels are in 32 bit float format.
				argv[5]: inputConfigFile2 - (std::string)  Path name of the configuration text file of the moving image. Pixels are in 32 bit float format. 
				
				argv[6]: outputImageFile - (std::string) Path name and file name of the output output displacement field from SIFT + TPS in raw data format
				argv[7]: outConfigFile - (std::string) Path name of the configuration text file of the displacement field from SIFT + TPS
				
				argv[8]: outPointsFile - (std::string) Path name and file name prefix of the output files containing information on SIFT feature point positions
				argv[9]: outputDisplacementField - (std::string) Path name and file name prefix of the output displacement field from SIFT + TPS in raw data format
				argv[10]: outDisplacementConfig - (std::string) Path name of the configuration text file of the displacement field from SIFT + TPS
				argv[11]: inRemoveEdgePointsPercent - (std::string) Percentage of image removed starting from the image edge and moving inward based on toral image size in each dimension. Input is a float number [0,100] in string format.
				
Author:			Frankie (Hoi-Ki) Tong, Feb 10 2015

Changes:		Frankie (Hoi-Ki) Tong, Feb 10 2015
				Initial creation of the file.
*/

#include "mainStandard.h"

bool mainStandard (int argc, char *argv[]) {

	std::string inputImageFile1 = "";
	std::string inputImageFile2 = "";
	std::string inConfigFile1 = "";
	std::string inConfigFile2 = "";
	std::string outputImageFile = "";
	std::string outConfigFile = "";
	std::string outPointsFile = "";
	std::string outputDisplacementField = "";
	std::string outDisplacementConfig = "";
	std::string inRemoveEdgePointsPercent = "";

	if (argc == 11) {

		inputImageFile1 = argv[2];
		inputImageFile2 = argv[3];
		inConfigFile1 = argv[4];
		inConfigFile2 = argv[5];
		outputImageFile = argv[6];
		outConfigFile = argv[7];
		outPointsFile = argv[8];
		outputDisplacementField = argv[9];
		outDisplacementConfig = argv[10];

	}

	else if (argc == 12) {
		inputImageFile1 = argv[2];
		inputImageFile2 = argv[3];
		inConfigFile1 = argv[4];
		inConfigFile2 = argv[5];
		outputImageFile = argv[6];
		outConfigFile = argv[7];
		outPointsFile = argv[8];
		outputDisplacementField = argv[9];
		outDisplacementConfig = argv[10];
		inRemoveEdgePointsPercent = argv[11];
	}

	else
	{
		std::cerr << "Incorrect number of parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " --standard inputImageFile1 inputImageFile2 inConfigFile1 inConfigFIle2 outputImageFile outConfgFile outPointsFile outputDisplacementField outDisplacementConfig (optional)inRemoveEdgePointsPercent" << std::endl;
		return EXIT_FAILURE;
	}




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




	//Find Keypoints and Matchedpoints betweent the 2 images
#ifdef TIMER
	//Timer functionality
	std::clock_t start;
	double duration;

	start = std::clock();
#endif

	typedef itk::ScaleInvariantFeatureImageFilter<ImageType, Dimension> SiftFilterType;

	SiftFilterType::PointSetTypePointer fixedImageKeyPoints, movingImageKeyPoints;
	SiftFilterType siftFilter1, siftFilter2;

	//Get Keypoints from fixed image and moving image
	fixedImageKeyPoints = siftFilter1.getSiftFeatures(fixedImage);
	movingImageKeyPoints = siftFilter1.getSiftFeatures(movingImage);


	//Find and remove edge points if requested by user
	if (argc == 12) {

		//Get percentage of the image where points found in the area will be removed
		std::stringstream ss(inRemoveEdgePointsPercent);
		double removeEdgePointsPercent;
		if (!(ss >> removeEdgePointsPercent)) removeEdgePointsPercent = 0;


		//Find the bounding box of the fixed image where we keep the points

		const ImageType::RegionType& fixedImageRegion = fixedImage->GetLargestPossibleRegion();
		const ImageType::SizeType& fixedImageSize = fixedImageRegion.GetSize();

		const ImageType::SpacingType fixedImageSpacing = fixedImage->GetSpacing();

		const ImageType::PointType& fixedImageOrigin  = fixedImage->GetOrigin();

		//Find the bounding box of the moving image where we keep the points

		const ImageType::RegionType& movingImageRegion = movingImage->GetLargestPossibleRegion();
		const ImageType::SizeType& movingImageSize = movingImageRegion.GetSize();

		const ImageType::SpacingType movingImageSpacing = movingImage->GetSpacing();

		const ImageType::PointType& movingImageOrigin  = movingImage->GetOrigin();

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

		//Swap tmpFixedImageKeyPoints list with fixedImageKeyPoints list
		/*fixedImageKeyPoints = SiftFilterType::PointSetType::New();
		for (int i = 0; i < tmpFixedImageKeyPoints->GetNumberOfPoints(); i++) {
			SiftFilterType::PointType tmpp;
			tmpFixedImageKeyPoints->GetPoint(i,&tmpp);

			SiftFilterType::FeatureType tmpf;
			tmpFixedImageKeyPoints->GetPointData(i,&tmpf);

			fixedImageKeyPoints->SetPoint(i,tmpp);
			fixedImageKeyPoints->SetPointData(i,tmpf);

		}*/

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

		//Swap tmpMovingImageKeyPoints list with movingImageKeyPoints list
		/*movingImageKeyPoints = SiftFilterType::PointSetType::New();
		for (int i = 0; i < tmpMovingImageKeyPoints->GetNumberOfPoints(); i++) {
			SiftFilterType::PointType tmpp;
			tmpMovingImageKeyPoints->GetPoint(i,&tmpp);

			SiftFilterType::FeatureType tmpf;
			tmpMovingImageKeyPoints->GetPointData(i,&tmpf);

			movingImageKeyPoints->SetPoint(i,tmpp);
			movingImageKeyPoints->SetPointData(i,tmpf);

		}*/

		movingImageKeyPoints = tmpMovingImageKeyPoints;
	}



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
	siftFilter1.MatchKeypointsFeatures(&fixedImageKeyPoints,&movingImageKeyPoints, &fixedImageMatchedPoints, &movingImageMatchedPoints, NULL);

	
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

	//Setup to get points conatiner for ITK TPS points contaienr
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
	std::cout<<"Making displacement field"<<std::endl;
	
	typedef itk::Vector< PixelType, Dimension > VectorPixelType;  
	typedef itk::Image< VectorPixelType, Dimension > DisplacementFieldImageType;  
	typedef itk::TransformToDisplacementFieldFilter< DisplacementFieldImageType, CoordinateRepType > DisplacementFieldGeneratorType;

	DisplacementFieldGeneratorType::Pointer dispfieldGenerator = DisplacementFieldGeneratorType::New();  
	dispfieldGenerator->UseReferenceImageOn();  
	dispfieldGenerator->SetReferenceImage( fixedImage );  //sets the input parameters for hte field generated to be the same as fixedImage
	dispfieldGenerator->SetTransform( transform );  

	try    {    
		dispfieldGenerator->Update();
	}  catch ( itk::ExceptionObject & err )    {
		std::cerr << "Exception detected while generating deformation field";
		std::cerr << " : "  << err << std::endl;    
		return EXIT_FAILURE;    
	}

	DisplacementFieldImageType::Pointer displacementField = DisplacementFieldImageType::New();
	displacementField = dispfieldGenerator->GetOutput();




#ifdef TIMER
	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
#endif
	
	std::cout<<"Displacement field has been generated."<<std::endl;

#ifdef TIMER
	std::cout<<"\n\nDuration"<<std::endl;
	std::cout<<duration<<"\n\n"<<std::endl;
#endif




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

	//Default Interpolation Filter is linear interpolate which produces suboptimal interpolation values. Changing it to lanczos interpolation
	typedef itk::ConstantBoundaryCondition< ImageType > BoundaryConditionType;
	const unsigned int WindowRadius = 5;
	typedef itk::Function::LanczosWindowFunction<WindowRadius> WindowFunctionType;
	typedef itk::WindowedSincInterpolateImageFunction< ImageType, WindowRadius, WindowFunctionType, BoundaryConditionType, CoordinateRepType > InterpolatorType;

	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	filter->SetInterpolator( interpolator );

	filter->SetInput(fixedImage);
	filter->SetTransform(transform);
	filter->Update();

	ImageType::Pointer fixedImageTransformed = ImageType::New();
	fixedImageTransformed = filter->GetOutput();



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

	std::cout<<"Supposedly, distance between matched poitns is written properly."<<std::endl;

	//Print out Amira LandmarksASCII file for matched points
	outputDistance.writeAmiraMatchedPoints(outPointsFile + ".landmarkASCII");

	//Write out displacement field and its config file
	//Write displacement field
	itkImageIORaw<VectorPixelType, Dimension>* inputImageImport3 = new itkImageIORaw<VectorPixelType, Dimension>;

	inputImageImport3->setConfigFileName(outDisplacementConfig);

	//Write ConfigFile
	ImageIOData outputDisplacementImageIOData;

	outputDisplacementImageIOData.setPixelType("float");
	outputDisplacementImageIOData.setPixelDimensionality(Dimension);
	outputDisplacementImageIOData.setImageDimensionality(Dimension);

	const DisplacementFieldImageType::RegionType& OutputDisplacementRegion = fixedImage->GetLargestPossibleRegion();
	const DisplacementFieldImageType::SizeType& OutputDisplacementSize = OutputRegion.GetSize();
	unsigned int tmpDisplacementSize[Dimension];
	for (int i = 0; i < Dimension; i++) {
		tmpDisplacementSize[i] = OutputDisplacementSize[i];
	}

	outputDisplacementImageIOData.setDimensionSize((int *)tmpDisplacementSize);

	const DisplacementFieldImageType::SpacingType OutputDisplacementSpacing = fixedImage->GetSpacing();
	float tmpDisplacementSpacing[Dimension];
	for (int i = 0; i < Dimension; i++) {
		tmpDisplacementSpacing[i] = (float)OutputDisplacementSpacing[i];
	}
	
	outputDisplacementImageIOData.setSpacingSize((float *)tmpDisplacementSpacing);

	const DisplacementFieldImageType::PointType& OutputDisplacementOrigin  = fixedImage->GetOrigin();
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
	inputImageImport3->setitkImage((void *) &displacementField);
	inputImageImport3->setImageFileName(outputDisplacementField);
	inputImageImport3->writeRawImageKnown();

	delete inputImageImport3;

	std::cout << "Supposedly, displacement field is written in properly." << std::endl;
	



	return EXIT_SUCCESS;

}