/*
File:			mainTPS.cxx

Description:	Subfunction called from imageioraw.cxx generates a displacement field using thin plate spline interpolation given a set of matched points in 2 text files and an input configuration file for the inital displacement field. Input parameters based on input processing from calling mainTPS as a command line function, ignoring the first 3 input parameters. Refer below for input process parameters.

Inputs:			argc - (int) Number of input arguments stored in argv. Default number is 10.

				argv[2]: inPointsFile1 - (std::string) Path name of the matched point set for the fixed image.
				argv[3]: inPointsFile2 - (std::string) Path name of the matched point set for the moving image.
				
					Input keypoint format is a text file with each new line representing a new point and coordinates of the point given in a comma seperated format.
					
					ex)		sample_keypoint_list.txt with points (0.24, 0.487), (0.54, 0.7), and (-4.78, 1.55)
					
							<start text file>
							0.24, 0.487,
							0.54, 0.7,
							-4.78, 1.55,
							<start text file>
				
				argv[4]: inVectorConfigFile1 - (std::string) Path name of the configuration text file of the displacement field to be generated initally. Pixels are in 32 bit float format.
				
				argv[5]: outPointsFile - (std::string) Path name and prefix of the matched point set for displacement field.
				
				argv[6]: outputDisplacementField - (std::string) Path name and file name prefix of the output displacement field in raw data format
				argv[7]: outDisplacementConfig - (std::string) Path name of the configuration text file of the displacement field from TPS
				
				argv[8]: relativePosition - (std::string) Determines where the origin of the output image/field files are placed. They can either be placed at the "origin" or at the "center" of the image/field. Input is a string containing "origin" or "center.
				argv[9]: outputResolution - (std::string) Downsample factor of the output image. Input is a float number in a string format.
				
Author:			Frankie (Hoi-Ki) Tong, Aug 13 2015

Changes:		Frankie (Hoi-Ki) Tong, Aug 13 2015
				Initial creation of the file.
*/

#include "mainTPS.h"

bool mainTPS (int argc, char *argv[]) {

	std::string inPointsFile1 = "";
	std::string inPointsFile2 = "";

	std::string inVectorConfigFile1 = "";		
	
	std::string outPointsFile = "";
	std::string outputDisplacementField = "";
	std::string outDisplacementConfig = "";
	std::string relativePosition = "";
	std::string outputResolution = "";


	if (argc == 10) {

		inPointsFile1 = argv[2];
		inPointsFile2 = argv[3];
		inVectorConfigFile1 = argv[4];			
		outPointsFile = argv[5];
		outputDisplacementField = argv[6];
		outDisplacementConfig = argv[7];
		relativePosition = argv[8];
		outputResolution = argv[9];

	} else
	{
		std::cerr << "Incorrect number of parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << "--TPS inPointsFile1 inPointsFile2 inVectorConfigFile1 outPointsFile outputDisplacementField outDisplacementConfig relativePosition(origin/center) outputDisplacementResolution(original/#)"<< std::endl;
		return EXIT_FAILURE;
	}

	//Initialize Imagetypes and image pointers
	typedef float									PixelType;
	const unsigned int Dimension = 2;
	typedef itk::Image<PixelType, Dimension>		ImageType;
	
	ImageType::Pointer fixedImage = ImageType::New();
	ImageType::Pointer movingImage = ImageType::New();
	
	//Read in vector field
	typedef itk::Vector< PixelType, Dimension > VectorPixelType;  
	typedef itk::Image< VectorPixelType, Dimension > VectorFieldImageType;  
	
	itkImageIORaw<VectorPixelType, Dimension>* inputVectorFieldImport = new itkImageIORaw<VectorPixelType, Dimension>;
	
	inputVectorFieldImport->setConfigFileName(inVectorConfigFile1);
	if (inputVectorFieldImport->readConfigFile()) {
		std::cerr << "Exception thrown while reading the in configuration file" << std::endl;
		return EXIT_FAILURE;
	}
	
	std::cout << "Supposedly, vector field configuration file is read in properly." << std::endl;	
	
	//Generate template vector field to specify TPS field dimensions and size
	//Make a new image of the same size of parameters of the fixedImage, but downsample it to the size I want

	itk::Size<Dimension> vectorFieldSize;
	float vectorFieldSpacing[Dimension];
	float vectorFieldOrigin[Dimension];

	for (int i = 0; i < Dimension; i++) {
		vectorFieldSize[i] = (int)(inputVectorFieldImport->getImageIOData().getDimensionSize())[i];
		vectorFieldSpacing[i] = (float)(inputVectorFieldImport->getImageIOData().getSpacingSize())[i];
		vectorFieldOrigin[i] = (float)(inputVectorFieldImport->getImageIOData().getOrigin()[i]);
	}
	
	ImageType::Pointer vectorField = ImageType::New();
	ImageType::RegionType vectorRegion;

	vectorRegion.SetSize(vectorFieldSize);
	
	ImageType::IndexType vectorIndex;
	for (int i = 0; i < Dimension; i++) {
		vectorIndex[i] = 0;
	}

	vectorRegion.SetIndex(vectorIndex);

	vectorField->SetRegions(vectorRegion);
	vectorField->SetSpacing(vectorFieldSpacing);
	vectorField->SetOrigin(vectorFieldOrigin);
	

	
	//Setup for getting the matched keypoint lists
	typedef itk::Array<float> FeatureType;
	typedef double CoordinateRepType;
	typedef itk::PointSet< FeatureType, Dimension,
	itk::DefaultStaticMeshTraits< FeatureType, Dimension, Dimension, CoordinateRepType > > PointSetType;
	
	typedef PointSetType::PointType PointType;
	typedef PointSetType::Pointer PointSetTypePointer;
	
	PointSetType::Pointer fixedImageMatchedPoints = PointSetType::New();
	PointSetType::Pointer movingImageMatchedPoints = PointSetType::New();
	
	//Read in the points file 1 and 2
	KeypointsIO<FeatureType, Dimension> inputKeypoints;
	inputKeypoints.setKeypointSetPtr(&fixedImageMatchedPoints);
	inputKeypoints.setFileName(inPointsFile1);
	inputKeypoints.readKeypointSet();

	inputKeypoints.setKeypointSetPtr(&movingImageMatchedPoints);
	inputKeypoints.setFileName(inPointsFile2);
	inputKeypoints.readKeypointSet();



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
		movingImageMatchedPointsContainer->GetElementIfIndexExists(ID,&tmpp2);

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

	if (outputResolution == "original"){

		//Make displacement field based off of original image resolution
		//typedef itk::Vector< PixelType, Dimension > VectorPixelType;  
		typedef itk::Image< VectorPixelType, Dimension > DisplacementFieldImageType;  
		typedef itk::TransformToDisplacementFieldFilter< DisplacementFieldImageType, CoordinateRepType > DisplacementFieldGeneratorType;

		DisplacementFieldGeneratorType::Pointer dispfieldGeneratorOriginal = DisplacementFieldGeneratorType::New();  
		DisplacementFieldImageType::Pointer displacementFieldOriginal = DisplacementFieldImageType::New();
	
		dispfieldGeneratorOriginal->UseReferenceImageOn(); 
		dispfieldGeneratorOriginal->SetReferenceImage( vectorField );  //sets the input parameters for the field generated to be the same as vectorField
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

	
	if (outputResolution != "original") {

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

			const ImageType::RegionType& DisplacementInputRegionFixedImage = vectorField->GetLargestPossibleRegion();
			const ImageType::SizeType& DisplacementInputSizeFixedImage = DisplacementInputRegionFixedImage.GetSize();
			unsigned int OldSizeImage[Dimension];
			unsigned int NewSizeImage[Dimension];
			for (int i = 0; i < Dimension; i++) {
				OldSizeImage[i] = DisplacementInputSizeFixedImage[i];
				NewSizeImage[i] = std::floor((double) OldSizeImage[i]/displacementResampleFactor );
			}

			const ImageType::SpacingType OldSpacingImage = vectorField->GetSpacing();
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

			const ImageType::PointType& OldOriginImage  = vectorField->GetOrigin();
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
			tmpImage->SetDirection(vectorField->GetDirection());
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

	//Applying transform to initial starting points and writing them out for easier input into matlab

	//First convert fixedImageMatchedPoints into a itk::Mesh format
	typedef itk::Mesh<CoordinateRepType,Dimension> MeshType;
	typedef MeshType::PointIdentifier MeshIDType;
	typedef MeshType::PointType MeshPointType;

	ID = itk::NumericTraits<PointIdType>::ZeroValue();
	MeshIDType MeshID = itk::NumericTraits<MeshIDType>::ZeroValue();

	MeshType::Pointer MeshPoints = MeshType::New();
	
	
	//Convert each fixedImageMatchedPoint into mesh points
	for (int i = 0; i < fixedImageMatchedPoints->GetNumberOfPoints(); i++) {

		itk::Point<CoordinateRepType, Dimension> tmpp1;
		fixedImageMatchedPointsContainer->GetElementIfIndexExists(ID,&tmpp1);
		
		MeshPoints->SetPoint( MeshID, tmpp1 );

		ID++;	
		MeshID++;

	}

	//Apply transform onto the mesh points

	typedef itk::TransformMeshFilter< MeshType, MeshType, TransformType > TransformMeshFilterType;
	TransformMeshFilterType::Pointer transformMeshFilter = TransformMeshFilterType::New();
	transformMeshFilter->SetInput( MeshPoints );
	transformMeshFilter->SetTransform( transform );
	transformMeshFilter->Update();
	MeshPoints = transformMeshFilter->GetOutput();

	//Convert the mesh points back into SIFT points

	PointSetType::Pointer fixedImageTransformedPoints = PointSetType::New();
	ID = itk::NumericTraits<PointIdType>::ZeroValue();
	MeshID = itk::NumericTraits<MeshIDType>::ZeroValue();

	for (int i = 0; i < MeshPoints->GetNumberOfPoints(); i++) {

		MeshPointType tmpp1;
		itk::Point<CoordinateRepType, Dimension> tmpp2;
		MeshPoints->GetPoint(MeshID,&tmpp1);

		for (int j = 0; j < Dimension; j++) {
			tmpp2[j] = tmpp1[j];
		}
		
		fixedImageTransformedPoints->SetPoint( ID, tmpp2 );

		ID++;	
		MeshID++;

	}


	//Print out the keypoints and matchedpoints files

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

	//Print out Matlab file for matched points
	outputDistance.writeMatlabMatchedPoints(outPointsFile + ".m");


	//Print out transformed points
	outputKeypoints3.setKeypointSetPtr(&fixedImageMatchedPoints);
	outputKeypoints4.setKeypointSetPtr(&fixedImageTransformedPoints);

	outputDistance.setMatchedPoints(&outputKeypoints3, &outputKeypoints4);

	//Print out Amira LandmarksASCII file for matched points
	outputDistance.writeAmiraMatchedPoints(outPointsFile + "_TransformedPoints.landmarkASCII");

	//Print out Matlab file for matched points
	outputDistance.writeMatlabMatchedPoints(outPointsFile + "_TransformedPoints.m");


	return EXIT_SUCCESS;

}