/*
File:			mainBSplinePoints.cxx

Description:	Subfunction called from imageioraw.cxx generates a displacement field using B-spline interpolation given a set of matched points in 2 text files and an input configuration file for the inital displacement field. Input parameters based on input processing from calling mainBSplinePoints as a command line function, ignoring the first 3 input parameters. Refer below for input process parameters.

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

#include "mainBSplinePoints.h"

bool mainBSplinePoints (int argc, char *argv[]) {

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
		std::cerr << "Usage: " << argv[0] <<" "<< argv[1]<< " ";
		std::cerr << "inPointsFile1 inPointsFile2 inVectorConfigFile1 outPointsFile outputDisplacementField outDisplacementConfig relativePosition(origin/center) outputDisplacementResolution(original/#)"<< std::endl;
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
	
	//Generate template vector field to specify BSpline field dimensions and size
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
	
	typedef itk::PointSet< VectorPixelType, Dimension,
	itk::DefaultStaticMeshTraits< VectorPixelType, Dimension, Dimension, double >> PointSetType;
	
	typedef PointSetType::PointType PointType;
	typedef PointSetType::Pointer PointSetTypePointer;
	
	PointSetType::Pointer fixedImageMatchedPoints = PointSetType::New();
	PointSetType::Pointer movingImageMatchedPoints = PointSetType::New();
	
	//Read in the points file 1 and 2
	KeypointsIO<VectorPixelType, Dimension> inputKeypoints;
	inputKeypoints.setKeypointSetPtr(&fixedImageMatchedPoints);
	inputKeypoints.setFileName(inPointsFile1);
	inputKeypoints.readKeypointSet();

	inputKeypoints.setKeypointSetPtr(&movingImageMatchedPoints);
	inputKeypoints.setFileName(inPointsFile2);
	inputKeypoints.readKeypointSet();

	std::cout<<"Finding transform from keypoints"<<std::endl;

	//PointSetType::PointType::Dimension;

	//Generate displacement vector data for each fixed point in fixedImageMatchedPoints
	PointType fixedPoint;
	PointType movingPoint;

	typedef PointSetType::PointsContainer::Pointer PointSetTypeContainerPointer;
	typedef PointSetType::PointIdentifier PointIdType;
	PointIdType ID = itk::NumericTraits<PointIdType>::ZeroValue();

	PointSetTypeContainerPointer fixedImageMatchedPointsContainer = fixedImageMatchedPoints->GetPoints();
	PointSetTypeContainerPointer movingImageMatchedPointsContainer = movingImageMatchedPoints->GetPoints();

	//Flip point and pointdata such that the point value is in PixelType and the coordinates is in PointType...
	for (int i = 0; i < fixedImageMatchedPoints->GetNumberOfPoints(); i++) {

		fixedImageMatchedPointsContainer->GetElementIfIndexExists(ID,&fixedPoint);
		movingImageMatchedPointsContainer->GetElementIfIndexExists(ID,&movingPoint);

		PointSetType::PixelType Displacement;

		for (int j = 0; j < Dimension; j++) {
			Displacement[j] = movingPoint[j]-fixedPoint[j];
		}

		fixedImageMatchedPoints->SetPointData(i,Displacement);

		ID++;	

	}

	//Instantiate B Spline filter and set the desired parameters
	typedef itk::BSplineScatteredDataPointSetToImageFilter<PointSetType,VectorFieldImageType> FilterType;
	FilterType::Pointer filter = FilterType::New();
	filter->SetSplineOrder(5);
	FilterType::ArrayType ncps;
	ncps.Fill(6);
	filter->SetNumberOfControlPoints(ncps);
	filter->SetNumberOfLevels(5);

	filter->SetInput(fixedImageMatchedPoints);

	std::cout<<"Image transform has been calculated"<<std::endl;


	//Create displacement field here...

	std::cout<<"Making displacement fields"<<std::endl;

	if (outputResolution == "original"){

		//Make displacement field based off of original image resolution
		
		filter->SetOrigin(vectorField->GetOrigin());
		filter->SetSpacing(vectorField->GetSpacing());
		filter->SetSize(vectorField->GetLargestPossibleRegion().GetSize());

		typedef itk::Image< VectorPixelType, Dimension > DisplacementFieldImageType;  
		DisplacementFieldImageType::Pointer displacementFieldOriginal = DisplacementFieldImageType::New();
	
		try {    
			filter->Update();
		}  catch ( itk::ExceptionObject & err )    {
			std::cerr << "Exception detected while generating deformation field";
			std::cerr << " : "  << err << std::endl;    
			return EXIT_FAILURE;    
		}

		displacementFieldOriginal = filter->GetOutput();

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
			
			//We have a specified resampleFactor that we would like th displacement field to be at. Resize region to allow for this.
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

			filter->SetOrigin(tmpImage->GetOrigin());
			filter->SetSpacing(tmpImage->GetSpacing());
			filter->SetSize(tmpImage->GetLargestPossibleRegion().GetSize());
			


			try {    
				filter->Update();
			}  catch ( itk::ExceptionObject & err )    {
				std::cerr << "Exception detected while generating deformation field";
				std::cerr << " : "  << err << std::endl;    
				return EXIT_FAILURE;    
			}

			displacementField = filter->GetOutput();

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

	//Print out the keypoints and matchedpoints files

	KeypointsIO<VectorPixelType, Dimension> outputKeypoints3;
	outputKeypoints3.setKeypointSetPtr(&fixedImageMatchedPoints);
	outputKeypoints3.setFileName(outPointsFile + "_fixedImageMatchedPoints.txt");
	outputKeypoints3.writeKeypointSet();

	KeypointsIO<VectorPixelType, Dimension> outputKeypoints4;
	outputKeypoints4.setKeypointSetPtr(&movingImageMatchedPoints);
	outputKeypoints4.setFileName(outPointsFile + "_movingImageMatchedPoints.txt");
	outputKeypoints4.writeKeypointSet();

	std::cout<<"Supposedly, keypoint lists is written properly."<<std::endl;

	//Want to print out distances as well. Need some fancy footwork here...
	MatchedPoints<VectorPixelType, Dimension> outputDistance;
	outputDistance.setMatchedPoints(&outputKeypoints3, &outputKeypoints4);
	outputDistance.writeDistance(outPointsFile + "_matchedPointsDistance.txt");

	std::cout<<"Supposedly, distance between matched points is written properly."<<std::endl;

	//Print out Amira LandmarksASCII file for matched points
	outputDistance.writeAmiraMatchedPoints(outPointsFile + ".landmarkASCII");

	//Print out Matlab file for matched points
	outputDistance.writeMatlabMatchedPoints(outPointsFile + ".m");

	return EXIT_SUCCESS;

}