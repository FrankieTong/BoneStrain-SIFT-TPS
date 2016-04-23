/*
File:			mainRotate.cxx

Description:	Subfunction called from imageioraw.cxx applies a rotation onto an input image and applies SIFT to analyze SIFT functionality under rotation. Input parameters based on input processing from calling mainRotate as a command line function, ignoring the first 3 input parameters. Refer below for input process parameters.

Inputs:			argc - (int) Number of input arguments stored in argv. Default number is 8.
				argv[2]: rotateMode - (std::string) Determines the matching mode of feature poitns found by SIFT. Can either be "--euclidean" which will perform feature matching based on closest ecudlian distance after rotational transform or "--sift" which will perform feature matching using SIFT feature descriptors found in the original and rotated image.
				argv[3]: inputImageFile - (std::string) Path name of the raw data file of the input image. Pixels are in 32 bit float format.
				argv[4]: outputImageFile - (std::string) Path name and file name of the rotated imaged image as well as the text files containing statistics for matched feature points. Pixels are in 32 bit float format.
				argv[5]: inConfigFile - (std::string) Path name of the configuration text file of the input image. Pixels are in 32 bit float format.
				argv[6]: outConfigFile - (std::string) Path name of the configuration text file of the output image.
				argv[7]: rotationAngle - (std::string) Rotation angle along the dimension 3 axis in degrees. Rotates image counterclockwise.
				
Author:			Frankie (Hoi-Ki) Tong, Nov 1 2014

Changes:		Frankie (Hoi-Ki) Tong, Nov 1 2014
				Initial creation of the file.
*/

#include "mainRotate.h"

bool mainRotate(int argc, char *argv[]) {
	
	std::string rotateMode = "";
	std::string inputImageFile = "";
	std::string outputImageFile = "";
	std::string inConfigFile = "";
	std::string outConfigFile = "";
	std::string rotationAngle = "";

	if (argc == 8) {

		rotateMode = argv[2];
		inputImageFile = argv[3];
		outputImageFile = argv[4];
		inConfigFile = argv[5];
		outConfigFile = argv[6];
		rotationAngle = argv[7];

	}

	else
	{
		std::cerr << "Incorrect number of parameters " << std::endl;

		std::cerr << "Usage: " << argv[0];
		std::cerr << " --rotate --euclidean/--sift inputImageFile outputImageFile inConfigFile outConfigFile rotationAngle" << std::endl;
		return EXIT_FAILURE;
	}

	if (!(rotateMode == "--euclidean" || rotateMode == "--sift"))
	{
		std::cerr << "Incorrect test mode " << std::endl;

		std::cerr << "Usage: " << argv[0];
		std::cerr << " --rotate --euclidean/--sift inputImageFile outputImageFile inConfigFile outConfigFile rotationAngle" << std::endl;
		return EXIT_FAILURE;
	}
	
	//Initialize Imagetypes and image pointers
	typedef float							PixelType;
	const unsigned int Dimension = 2;
	typedef itk::Image<PixelType, Dimension>		ImageType;
	
	ImageType::Pointer fixedImage = ImageType::New();
	
	itkImageIORaw<PixelType, Dimension>* inputImageImport = new itkImageIORaw<PixelType, Dimension>;
	inputImageImport->setConfigFileName(inConfigFile);
	if (inputImageImport->readConfigFile()) { 
		std::cerr << "Exception thrown while reading the in configuration file" << std::endl;
		return EXIT_FAILURE;
	}

	void * vImage = & fixedImage;
	inputImageImport->setitkImage(vImage);
	inputImageImport->setImageFileName(inputImageFile);
	if (inputImageImport->readRawImageKnown()) {
		std::cerr << "Image " << inputImageFile << " using config file " << inConfigFile << " did not read in properly!"<< std::endl;
		return 1;
	}

	delete inputImageImport;

	std::cout << "Supposedly, image is read in properly." << std::endl;
	
	//Create a second image that is a rotated version of the first.
	
	//Resample here
#ifdef RESAMPLE

	typedef itk::IdentityTransform<double, Dimension> upSampleIdentityTransformType;
	typedef itk::LinearInterpolateImageFunction<ImageType,double> upSampleInterpolatorType;
	typedef itk::ResampleImageFilter<ImageType,ImageType,double,double> upSampleResampleImaageFilterType;

	upSampleIdentityTransformType::Pointer upSampleIdentityTransform = upSampleIdentityTransformType::New();
	upSampleIdentityTransform->SetIdentity();

	upSampleInterpolatorType::Pointer upSampleInterpolator = upSampleInterpolatorType::New();

	upSampleResampleImaageFilterType::Pointer upSampleResampleImageFilter = upSampleResampleImaageFilterType::New();
	upSampleResampleImageFilter->SetTransform(upSampleIdentityTransform);
	upSampleResampleImageFilter->SetInterpolator(upSampleInterpolator);

	const ImageType::RegionType& ResampleInputRegion = fixedImage->GetLargestPossibleRegion();
	const ImageType::SizeType& ResampleInputSize = ResampleInputRegion.GetSize();
	unsigned int OldSize[Dimension];
	unsigned int NewSize[Dimension];
	for (int i = 0; i < Dimension; i++) {
		OldSize[i] = ResampleInputSize[i];
		NewSize[i] = (int)( RESAMPLE_FACTOR * (double) OldSize[i] );
	}

	const ImageType::SpacingType OldSpacing = fixedImage->GetSpacing();
	double NewSpacing[Dimension];
	for (int i = 0; i < Dimension; i++) {
		NewSpacing[i] = OldSpacing[i] * (double) OldSize[i]/NewSize[i];
	}
	upSampleResampleImageFilter->SetOutputSpacing(NewSpacing);

	const ImageType::PointType& OldOrigin  = fixedImage->GetOrigin();
	double ResampleOrigin[Dimension];
	for (int i=0; i< Dimension; i++) {
		ResampleOrigin[i] = OldOrigin[i]-(OldSpacing[i] - NewSpacing[i])/2.0;
	}
	upSampleResampleImageFilter->SetOutputOrigin(ResampleOrigin);

	itk::Size<Dimension> upSampleSize;
	for (int i = 0; i < Dimension; i++) {
		upSampleSize[i] = NewSize[i];
	}

	upSampleResampleImageFilter->SetSize(upSampleSize);

	upSampleResampleImageFilter->SetInput(fixedImage);
	upSampleResampleImageFilter->Update();

	fixedImage = upSampleResampleImageFilter->GetOutput();
#endif

	//First thing to do is to find the rotation angle. (Check for failure, add later...)
	double angleDeg;
	std::istringstream(rotationAngle) >> angleDeg;
	
	//Create duplicate of input image
	typedef itk::ImageDuplicator< ImageType > DuplicatorType;
	DuplicatorType::Pointer duplicator = DuplicatorType::New();
	duplicator->SetInputImage(fixedImage);
	duplicator->Update();
	ImageType::Pointer movingImage = duplicator->GetOutput();
	
	//Rotate movingImage by the rotation angle about the center and find the inverse transform as well
	const ImageType::SpacingType& spacing = movingImage->GetSpacing();
	const ImageType::PointType& origin  = movingImage->GetOrigin();
	const ImageType::DirectionType& direction  = movingImage->GetDirection();
	const ImageType::SizeType size = movingImage->GetLargestPossibleRegion().GetSize();

	//First make the transform
	typedef itk::ScalableAffineTransform< double, Dimension > TransformType;

	TransformType::Pointer transform = TransformType::New();
	TransformType::Pointer inv_transform = TransformType::New();

	TransformType::OutputVectorType translation;

	transform->SetIdentity();
	transform->Scale( 1.0 );

	//Cycle through each dimension and shift by half to center orgin
	for (int k = 0; k < Dimension; ++k) {
		translation[k] = -(origin[k] + spacing[k] * size[k]/2.0);
	}
    transform->Translate(translation);

	//Now rotate the image about the z axis for amount specified by the use
	//Define axis through origin
	itk::Vector<double,Dimension> axis;
    //axis[0] = 0;
    //axis[1] = 0;
    //axis[2] = 1;

	for (int i = 0; i < Dimension; i++) {
		axis[i] = 0;
	}

	axis[Dimension-1] = 1;

	//Change degrees to radians
	double angleRad = 3.1415*angleDeg/180;

	//Rotate transform
	//transform->Rotate3D(axis,angleRad);
	transform->Rotate(0,1,angleRad,true);

	//Change origin back to the original origin
	for (int k = 0; k < Dimension; ++k) {
		translation[k] = -translation[k];
	}
	transform->Translate(translation);

	//Get inverse transform as well;
	transform->GetInverse(inv_transform);

	//Set up resample image filter with the proper parameters
	typedef itk::ResampleImageFilter<ImageType, ImageType >  FilterType;
	FilterType::Pointer filter = FilterType::New();

	filter->SetSize(size);
	filter->SetOutputSpacing(spacing);
	filter->SetOutputOrigin(origin);
	filter->SetOutputDirection(direction);
	filter->SetDefaultPixelValue( 0 );

	//Set extrapolator to take nearest neighbor instead of zero padding
	typedef itk::NearestNeighborExtrapolateImageFunction< ImageType, double > ExtrapolatorType;
	ExtrapolatorType::Pointer extrapolator = ExtrapolatorType::New();
	filter->SetExtrapolator( extrapolator );

	//Default Interpolation Filter is linear interpolate which produces suboptimal interpolation values. Changing it to lanczos interpolation
	typedef itk::ConstantBoundaryCondition< ImageType > BoundaryConditionType;
	const unsigned int WindowRadius = 5;
	typedef itk::Function::LanczosWindowFunction<WindowRadius> WindowFunctionType;
	typedef itk::WindowedSincInterpolateImageFunction< ImageType, WindowRadius, WindowFunctionType, BoundaryConditionType, double > InterpolatorType;

	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	filter->SetInterpolator( interpolator );

	filter->SetInput(movingImage);
	filter->SetTransform(transform);
	filter->Update();

	movingImage = filter->GetOutput();
	
	//Find SIFT features for fixedImage and Moving Image

#ifdef TIMER
	//Timer functionality
	std::clock_t start;
	double duration;

	start = std::clock();
#endif

	typedef itk::ScaleInvariantFeatureImageFilter<ImageType, Dimension> SiftFilterType;

	SiftFilterType::PointSetTypePointer keypoints1, keypoints2;
	SiftFilterType siftFilter1, siftFilter2;

	keypoints1 = siftFilter1.getSiftFeatures(fixedImage);
	keypoints2 = siftFilter1.getSiftFeatures(movingImage);

	//Setup for getting the keypoint list back from testing function.
	typedef itk::Array<float> FeatureType;
	typedef itk::PointSet< FeatureType, Dimension,
	itk::DefaultStaticMeshTraits< FeatureType, Dimension, Dimension, double > > PointSetType;
					  
	typedef PointSetType::PointType PointType;
	typedef PointSetType::Pointer PointSetTypePointer;

	PointSetTypePointer matchedkeypoints1;
	PointSetTypePointer matchedkeypoints2;

	PointSetTypePointer siftmatchedkeypoints1;
	PointSetTypePointer siftmatchedkeypoints2;

	float * distanceBetweenPoints = NULL;
	mainRotateSupportClass<PixelType, Dimension> testPositionDifference;

	//if mode is "--euclidean", match keypoints based on closest distance
	if (rotateMode == "--euclidean") {
		//Transform the feature points back into their original position and match by Ecludian distance
		testPositionDifference.setImageKeyPoints(&keypoints1, &keypoints2);
		testPositionDifference.setInverseTransform(transform);
		testPositionDifference.setMatchedKeyPoints(&matchedkeypoints1, &matchedkeypoints2);
		
		testPositionDifference.updateDistance();
		distanceBetweenPoints=testPositionDifference.getDistanceBetweenPoints();
	}	

	//if mode is "--sift", continue with SIFT matching to find matched keypoints
	else if (rotateMode == "--sift") {

		siftFilter1.MatchKeypointsFeatures(&keypoints1,&keypoints2, &siftmatchedkeypoints1, &siftmatchedkeypoints2, NULL);

		testPositionDifference.setImageKeyPoints(&siftmatchedkeypoints1, &siftmatchedkeypoints2);
		testPositionDifference.setInverseTransform(transform);
		testPositionDifference.setMatchedKeyPoints(&matchedkeypoints1, &matchedkeypoints2);
		
		testPositionDifference.updateMatched();
		distanceBetweenPoints=testPositionDifference.getDistanceBetweenPoints();
	}

	else
	{
		//It should never get here!
		return EXIT_FAILURE;
	}

#ifdef TIMER
	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
#endif

	std::cout<<"Matched points and distance between matched points:"<<std::endl;

	for(int i = 0; i <testPositionDifference.getDistanceBetweenPointsIndex(); i++) {
		PointType p1,p2;
		matchedkeypoints1->GetPoint(i,&p1);
		matchedkeypoints2->GetPoint(i,&p2);

		std::cout<<p1<<" => "<<p2<<" Distance: "<<distanceBetweenPoints[i]<<std::endl;
	}

#ifdef TIMER
	std::cout<<"\n\nDistance (repeated for ease in copying):"<<std::endl;
	for(int i = 0; i <testPositionDifference.getDistanceBetweenPointsIndex(); i++) {
		std::cout<<distanceBetweenPoints[i]<<std::endl;
	}
	std::cout<<"\n\nDuration"<<std::endl;
	std::cout<<duration<<"\n\n"<<std::endl;
#endif

	//Write matched and sift keypoints into a file
	KeypointsIO<FeatureType, Dimension> outputKeypoints1;
	outputKeypoints1.setKeypointSetPtr(&matchedkeypoints1);
	outputKeypoints1.setFileName(outputImageFile + "_matchedkeypoints1.txt");
	outputKeypoints1.writeKeypointSet();

	KeypointsIO<FeatureType, Dimension> outputKeypoints2;
	outputKeypoints2.setKeypointSetPtr(&matchedkeypoints2);
	outputKeypoints2.setFileName(outputImageFile + "_matchedkeypoints2.txt");
	outputKeypoints2.writeKeypointSet();

	KeypointsIO<FeatureType, Dimension> outputKeypoints3;
	outputKeypoints3.setKeypointSetPtr(&siftmatchedkeypoints1);
	outputKeypoints3.setFileName(outputImageFile + "_siftkeypoints1.txt");
	outputKeypoints3.writeKeypointSet();

	KeypointsIO<FeatureType, Dimension> outputKeypoints4;
	outputKeypoints4.setKeypointSetPtr(&siftmatchedkeypoints2);
	outputKeypoints4.setFileName(outputImageFile + "_siftkeypoints2.txt");
	outputKeypoints4.writeKeypointSet();

	//Personal test to see if I can read back properly...
	/*KeypointsIO<FeatureType, Dimension> outputKeypointsTest;
	outputKeypointsTest.setKeypointSetPtr(&matchedkeypoints1);
	outputKeypointsTest.setFileName(outConfigFile.substr(0,outputImageFile.find(".txt")) + "_matchedkeypoints1.txt");
	outputKeypointsTest.readKeypointSet();
	outputKeypointsTest.printKeypointSet();*/

	//Write output image
	itkImageIORaw<PixelType, Dimension>* inputImageImport2 = new itkImageIORaw<PixelType, Dimension>;

	//Write ConfigFile
	//The way to build ImageIOData structure differs from implementation to implementation which means I can't put this in itkImageIORaw for now...
	ImageIOData outputImageIOData;

	outputImageIOData.setPixelType("float");
	outputImageIOData.setPixelDimensionality(1);
	outputImageIOData.setImageDimensionality(Dimension);

	const ImageType::RegionType& OutputRegion = movingImage->GetLargestPossibleRegion();
	const ImageType::SizeType& OutputSize = OutputRegion.GetSize();
	unsigned int tmpSize[Dimension];
	for (int i = 0; i < Dimension; i++) {
		tmpSize[i] = OutputSize[i];
	}

	outputImageIOData.setDimensionSize((int *)tmpSize);

	const ImageType::SpacingType OutputSpacing = movingImage->GetSpacing();
	float tmpSpacing[Dimension];
	for (int i = 0; i < Dimension; i++) {
		tmpSpacing[i] = (float)OutputSpacing[i];
	}
	
	outputImageIOData.setSpacingSize((float *)tmpSpacing);

	const ImageType::PointType& OutputOrigin  = movingImage->GetOrigin();
	float tmpOrigin[Dimension];
	for (int i=0; i< Dimension; i++) {
		tmpOrigin[i] = (float)OutputOrigin[i];
	}

	outputImageIOData.setOrigin((float *)tmpOrigin);

	outputImageIOData.setHeaderSize(0);
	outputImageIOData.setByteOrder("little");

	inputImageImport2->setImageIOData(outputImageIOData);

	inputImageImport2->writeConfigFile();

	void * vImage2 = & movingImage;

	//Write Image
	inputImageImport2->setitkImage(vImage2);
	inputImageImport2->setImageFileName(outputImageFile);
	inputImageImport2->writeRawImageKnown();

	inputImageImport2->setConfigFileName(outConfigFile);

	delete inputImageImport2;
	delete distanceBetweenPoints;

	std::cout << "Supposedly, image is written in properly." << std::endl;

	return EXIT_SUCCESS;















}