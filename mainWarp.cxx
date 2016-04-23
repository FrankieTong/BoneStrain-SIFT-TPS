/*
File:			mainWarp.cxx

Description:	Subfunction called from imageioraw.cxx applies a Thin Plate Spline transform to an image given a set of matched points representing the transformation between the fixed and moving iamges. Input parameters based on input processing from calling mainWarp as a command line function, ignoring the first 3 input parameters. Refer below for input process parameters.

Inputs:			argc - (int) Number of input arguments stored in argv. Default number is 12.
				argv[2]: inputImageFile1 - (std::string) Path name of the raw data file of the fixed image. Pixels are in 32 bit float format.
				argv[3]: inputImageFile2 - (std::string) Path name of the raw data file of the moving image. Pixels are in 32 bit float format.
				
				argv[4]: inKeypoints1 - (std::string) Path name of the matched point set for the fixed image.
				argv[5]: inKeypoints2 - (std::string) Path name of the matched point set for the moving image.
				
					Input keypoint format is a text file with each new line representing a new point and coordinates of the point given in a comma seperated format.
					
					ex)		sample_keypoint_list.txt with points (0.24, 0.487), (0.54, 0.7), and (-4.78, 1.55)
					
							<start text file>
							0.24, 0.487,
							0.54, 0.7,
							-4.78, 1.55,
							<start text file>
				
				argv[6]: inConfigFile1 - (std::string) Path name of the configuration text file of the fixed image. Pixels are in 32 bit float format.
				argv[7]: inConfigFile2 - (std::string) Path name of the configuration text file of the moving image. Pixels are in 32 bit float format.
				
				argv[8]: outputImageFile - (std::string) Path name and file name of the warped image using TPS. Pixels are in 32 bit float format.
				argv[9]: outConfigFile - (std::string) Path name of the configuration text file of the output image
				
				argv[10]: outputDisplacementField - (std::string) Path name and file name prefix of the output displacement field from TPS in raw data format
				argv[11]: outDisplacementConfig - (std::string) Path name of the configuration text file of the displacement field from TPS
				
Author:			Frankie (Hoi-Ki) Tong, Nov 9 2014

Changes:		Frankie (Hoi-Ki) Tong, Nov 9 2014
				Initial creation of the file.
*/


#include "mainWarp.h"

bool mainWarp(int argc, char *argv[]) {
	
	std::string inputImageFile1 = "";
	std::string inputImageFile2 = "";
	std::string inKeypoints1 = "";
	std::string inKeypoints2 = "";
	std::string inConfigFile1 = "";
	std::string inConfigFile2 = "";
	std::string outputImageFile = "";
	std::string outConfigFile = "";
	std::string outputDisplacementField = "";
	std::string outputDisplacementConfig = "";

	if (argc == 12) {

		inputImageFile1 = argv[2];
		inputImageFile2 = argv[3];
		inKeypoints1 = argv[4];
		inKeypoints2 = argv[5];
		inConfigFile1 = argv[6];
		inConfigFile2 = argv[7];
		outputImageFile = argv[8];
		outConfigFile = argv[9];
		outputDisplacementField = argv[10];
		outputDisplacementConfig = argv[11];

	}

	else
	{
		std::cerr << "Incorrect number of parameters " << std::endl;

		std::cerr << "Usage: " << argv[0];
		std::cerr << "--warp inputImageFile1 inputImageFile2 inKeypoints1 inKeypoints2 inConfigFile1 inConfigFile2 outputImageFile outConfigFile outputDisplacementField outputDisplacementConfig" << std::endl;
		return EXIT_FAILURE;
	}
	
	//Initialize Imagetypes and image pointers
	typedef float							PixelType;
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

	//Setup for getting the keypoint lists
	typedef itk::Array<float> FeatureType;
	typedef double CoordinateRepType;
	typedef itk::PointSet< FeatureType, Dimension,
	itk::DefaultStaticMeshTraits< FeatureType, Dimension, Dimension, CoordinateRepType > > PointSetType;
					  
	typedef PointSetType::PointType PointType;
	typedef PointSetType::Pointer PointSetTypePointer;

	//Read in keypoints for corresponding images
	PointSetTypePointer fixedImageKeypoints;
	PointSetTypePointer movingImageKeypoints;

	KeypointsIO<FeatureType, Dimension> inputKeypoints;

	//Read in first set of keypoints
	inputKeypoints.setKeypointSetPtr(&fixedImageKeypoints);
	inputKeypoints.setFileName(inKeypoints1);
	inputKeypoints.readKeypointSet();
	inputKeypoints.printKeypointSet();

	//Read in second set of keypoints
	inputKeypoints.setKeypointSetPtr(&movingImageKeypoints);
	inputKeypoints.setFileName(inKeypoints2);
	inputKeypoints.readKeypointSet();
	inputKeypoints.printKeypointSet();

	//Images and keypoints are now read in as follows: fixedImage -> fixedImageKeypoints, movingImage ->mocingImageKeypoints
	
	//We now want to find the transform that will map fixedImageKeypoints to movingImageKeypoints

#ifdef TIMER
	//Timer functionality
	std::clock_t start;
	double duration;

	start = std::clock();
#endif

	//Find transform with keypoints goes here...
	std::cout<<"Finding transform from keypoints"<<std::endl;
	
	//typedef itk::IdentityTransform< CoordinateRepType, Dimension >  TransformType;

	typedef itk::ThinPlateSplineKernelTransform< CoordinateRepType, Dimension> TransformType;

	//Need to convert point set type so transform will accept the points
	//fixed -> source, moving -> target
	typedef PointSetType::PointsContainer::Pointer PointSetTypeContainerPointer;
	typedef PointSetType::PointIdentifier PointIdType;
	PointIdType ID = itk::NumericTraits<PointIdType>::ZeroValue();

	PointSetTypeContainerPointer fixedImageKeypointsContainer = fixedImageKeypoints->GetPoints();
	PointSetTypeContainerPointer movingImageKeypointsContainer = movingImageKeypoints->GetPoints();


	typedef itk::Point<CoordinateRepType, Dimension> TransformPointType;
	typedef TransformType::PointSetType TransformPointSetType;
	typedef TransformPointSetType::PointsContainer::Pointer TransformPointSetTypeContainerPointer;
	typedef TransformPointSetType::PointIdentifier TransformPointIdType;
	TransformPointIdType TransformID = itk::NumericTraits<TransformPointIdType>::ZeroValue();

	TransformPointSetType::Pointer sourceLandMarks = TransformPointSetType::New();
	TransformPointSetType::Pointer targetLandMarks = TransformPointSetType::New();

	TransformPointSetTypeContainerPointer sourceLandMarkContainer = targetLandMarks->GetPoints();
	TransformPointSetTypeContainerPointer targetLandMarkContainer = sourceLandMarks->GetPoints();

	for (int i = 0; i < fixedImageKeypoints->GetNumberOfPoints(); i++) {

		itk::Point<CoordinateRepType, Dimension> tmpp1;
		itk::Point<CoordinateRepType, Dimension> tmpp2;
		fixedImageKeypointsContainer->GetElementIfIndexExists(ID,&tmpp1);
		movingImageKeypointsContainer->GetElementIfIndexExists(ID, &tmpp2);

		sourceLandMarkContainer->InsertElement( TransformID, tmpp1 );
		targetLandMarkContainer->InsertElement( TransformID, tmpp2 );

		ID++;	
		TransformID++;

	}

	TransformType::Pointer transform = TransformType::New();
	transform->SetSourceLandmarks(sourceLandMarks);
	transform->SetTargetLandmarks(targetLandMarks);

	//transform->SetStiffness(0.1);

	transform->ComputeWMatrix();




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

	//Trasform fixedImage goes here

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

	//Set extrapolator to take nearest neighbor instead of zero padding
	/*typedef itk::NearestNeighborExtrapolateImageFunction< ImageType, double > ExtrapolatorType;
	ExtrapolatorType::Pointer extrapolator = ExtrapolatorType::New();
	filter->SetExtrapolator( extrapolator );*/

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

	fixedImage = filter->GetOutput();
	



	//Write output image
	itkImageIORaw<PixelType, Dimension>* inputImageImport2 = new itkImageIORaw<PixelType, Dimension>;

	inputImageImport2->setConfigFileName(outConfigFile);

	//Write ConfigFile
	//The way to build ImageIOData structure differs from implementation to implementation which means I can't put this in itkImageIORaw for now...
	ImageIOData outputImageIOData;

	outputImageIOData.setPixelType("float");
	outputImageIOData.setPixelDimensionality(1);
	outputImageIOData.setImageDimensionality(Dimension);

	const ImageType::RegionType& OutputRegion = fixedImage->GetLargestPossibleRegion();
	const ImageType::SizeType& OutputSize = OutputRegion.GetSize();
	unsigned int tmpSize[Dimension];
	for (int i = 0; i < Dimension; i++) {
		tmpSize[i] = OutputSize[i];
	}

	outputImageIOData.setDimensionSize((int *)tmpSize);

	const ImageType::SpacingType OutputSpacing = fixedImage->GetSpacing();
	float tmpSpacing[Dimension];
	for (int i = 0; i < Dimension; i++) {
		tmpSpacing[i] = (float)OutputSpacing[i];
	}
	
	outputImageIOData.setSpacingSize((float *)tmpSpacing);

	const ImageType::PointType& OutputOrigin  = fixedImage->GetOrigin();
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
	inputImageImport2->setitkImage((void *) &fixedImage);
	inputImageImport2->setImageFileName(outputImageFile);
	inputImageImport2->writeRawImageKnown();

	delete inputImageImport2;

	std::cout << "Supposedly, image is written in properly." << std::endl;




	//Write displacement field
	itkImageIORaw<VectorPixelType, Dimension>* inputImageImport3 = new itkImageIORaw<VectorPixelType, Dimension>;

	inputImageImport3->setConfigFileName(outputDisplacementConfig);

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