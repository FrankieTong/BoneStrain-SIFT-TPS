/*
File:			mainGradient.cxx

Description:	Subfunction called from imageioraw.cxx generates the gradient of an input vector field. Splits the input vector field into its individual components and calculates and outputs the gradient for each component. Input parameters based on input processing from calling mainGradient as a command line function, ignoring the first 3 input parameters. Refer below for input process parameters.

Inputs:			argc - (int) Number of input arguments stored in argv. Default number is 8.
				argv[2]: inputVectorField1 - (std::string) Path name of the raw data file of the input vector field. Pixels are in 32 bit float format.
				argv[3]: inVectorConfigFile1 - (std::string) Path name of the configuration text file of the input vector field. Pixels are in 32 bit float format.
				
				argv[4]: outputImageFile - (std::string) Path name and file name prefix of the output gradient fields. Pixels are in 32 bit float format.
				argv[5]: outConfigFile - (std::string) Path name of the configuration text file of the output gradient fields
				
Author:			Frankie (Hoi-Ki) Tong, Apr 11 2015

Changes:		Frankie (Hoi-Ki) Tong, Apr 11 2015
				Initial creation of the file.
*/

#include "mainGradient.h"

bool mainGradient (int argc, char *argv[]) {

	std::string inputVectorField1 = "";
	std::string inVectorConfigFile1 = "";
	std::string outputImageFile1 = "";
	std::string outConfigFile1 = "";


	if (argc == 6) {

		inputVectorField1 = argv[2];
		inVectorConfigFile1 = argv[3];
		outputImageFile1 = argv[4];
		outConfigFile1 = argv[5];
	
	}

	else
	{
		std::cerr << "Incorrect number of parameters " << std::endl;
		std::cerr << "Usage: " << argv[0];
		std::cerr << " --gradient inputVectorField1 inVectorConfigFile1 outputImageFile1 outConfigFile1" << std::endl;
		return EXIT_FAILURE;
	}	
	
	

	//Read in vector field
	typedef float									PixelType;
	const unsigned int Dimension = 2;
	
	typedef itk::Image< PixelType, Dimension > ScalarImageType;
	
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
	
	
	
	
	//Calculate Gradient Fields Here
	
	//Need to run each component of the image filter into their individual components, then convert the covariant vectors back to normal vectors before stacking all of the
	//vectors back up to get a n times n dimensional image containing the gradient value
	
	//First step is to split the input vector image into n scalar images:
	
	typedef itk::VectorIndexSelectionCastImageFilter<VectorFieldImageType, ScalarImageType> IndexSelectionType;
	IndexSelectionType::Pointer indexSelectionFilter = IndexSelectionType::New();
	//indexSelectionFilter->SetIndex(0);
	indexSelectionFilter->SetInput(vectorField);
	
	typedef itk::CovariantVector<PixelType,Dimension> CovariantPixelType;
	typedef itk::Image<CovariantPixelType,Dimension> CovariantVectorFieldImageType;
	typedef CovariantVectorFieldImageType::Pointer CovariantVectorFieldImageTypePtr;
	
	CovariantVectorFieldImageTypePtr CovariantVectorFieldImageStack[Dimension];
	
	typedef itk::GradientImageFilter<ScalarImageType, PixelType, PixelType> GradientImageFilterType;
	
	GradientImageFilterType::Pointer gradientImageFilter = GradientImageFilterType::New();
	
	
	
	for (int i = 0; i < Dimension; i++) {
		
		indexSelectionFilter->SetIndex(i);
		indexSelectionFilter->Update();
		
		

		gradientImageFilter->SetInput(indexSelectionFilter->GetOutput());
		gradientImageFilter->Update();
		
		CovariantVectorFieldImageStack[i] = CovariantVectorFieldImageType::New();
		CovariantVectorFieldImageStack[i] = gradientImageFilter->GetOutput();

		//Remove the edges of the image since they are going to be bad results
		
		CovariantVectorFieldImageType::SizeType desiredSize;

		desiredSize.Fill(1);

		typedef itk::CropImageFilter < CovariantVectorFieldImageType, CovariantVectorFieldImageType > RemoveEdgeFilterType;
		RemoveEdgeFilterType::Pointer removeEdgeFilter = RemoveEdgeFilterType::New();
		removeEdgeFilter->SetInput(CovariantVectorFieldImageStack[i]);
		removeEdgeFilter->SetBoundaryCropSize(desiredSize);
		removeEdgeFilter->Update();

		CovariantVectorFieldImageStack[i] = removeEdgeFilter->GetOutput();

		//Recalculate the origin
		
		const CovariantVectorFieldImageType::PointType& OldOrigin  = CovariantVectorFieldImageStack[i]->GetOrigin();
		const CovariantVectorFieldImageType::SpacingType& OldSpacing = CovariantVectorFieldImageStack[i]->GetSpacing();
		CovariantVectorFieldImageType::SpacingType NewOrigin;
		for (int j=0; j< Dimension; j++) {
			NewOrigin[j] = (float)OldOrigin[j] + (float)OldSpacing[j];
		}

		CovariantVectorFieldImageStack[i]->SetOrigin(NewOrigin);
	}
	
	
	//Output each CovariantVectorFieldImage in the stack as an individual raw file.
	
	
	/*
	//We now have the gradient vector of each component in the vector image. Lets recompose gradient vector matrix here:
	
	//First allocate image space
	typedef itk::Vector<PixelType, Dimension*Dimension> GradientVectorPixelType;
	typedef itk::Image<GradientVectorPixelType, Dimension> GradientVectorFieldImageType;
	
	GradientVectorFieldImageType::Pointer GradientVectorFieldImage = GradientVectorFieldImageType::New();
	
	const CovariantVectorFieldImageType::RegionType& CovariantVectorFieldRegion = CovariantVectorFieldImageStack[0]->GetLargestPossibleRegion();
	const CovariantVectorFieldImageType::SizeType& CovariantVectorFieldSize = CovariantVectorFieldRegion.GetSize();

	GradientVectorFieldImageType::RegionType GradientVectorFieldImageRegion;
	GradientVectorFieldImageType::IndexType GradientVectorFieldImageIndex;
	GradientVectorFieldImageType::SizeType GradientVectorFieldImageSize;

	for (int i = 0; i < Dimension; i++) {
		GradientVectorFieldImageSize[i] = CovariantVectorFieldSize[i];
	}

	for (int i = 0; i < Dimension; i++) {
		GradientVectorFieldImageIndex[i] = 0;
	}

	GradientVectorFieldImageRegion.SetSize(GradientVectorFieldImageSize);
	GradientVectorFieldImageRegion.SetIndex(GradientVectorFieldImageIndex);

	GradientVectorFieldImage->SetRegions(GradientVectorFieldImageRegion);
	GradientVectorFieldImage->Allocate();

	const CovariantVectorFieldImageType::SpacingType CovariantVectorFieldSpacing = CovariantVectorFieldImageStack[0]->GetSpacing();
	const CovariantVectorFieldImageType::PointType& CovariantVectorFieldOrigin  = CovariantVectorFieldImageStack[0]->GetOrigin();

	float GradientVectorFieldImageSpacing[Dimension];
	float GradientVectorFieldImageOrigin[Dimension];

	for (int i = 0; i < Dimension; i++) {
		GradientVectorFieldImageSpacing[i] = CovariantVectorFieldSpacing[i];
		GradientVectorFieldImageOrigin[i] = CovariantVectorFieldOrigin[i];
	}

	GradientVectorFieldImage->SetSpacing(GradientVectorFieldImageSpacing);
	GradientVectorFieldImage->SetOrigin(GradientVectorFieldImageOrigin);
	*/
	
	//Now calculate the infinitesimal tensor strain matrix from the gradient vectors of the individual displacement components at each vector point. We'll do this mannually.
	/*CovariantPixelType CovariantPixel;
	GradientVectorPixelType GradientVectorPixel;
	
	for (int j = 0; j < Dimension; j++) {
		
		//Set up iterators for both GradientVectorFieldImage and CovariantImages
		itk::ImageRegionIterator<GradientVectorFieldImageType> GradientVectorFieldImageIterator(GradientVectorFieldImage,GradientVectorFieldImageRegion);
		itk::ImageRegionIterator<CovariantVectorFieldImageType> CovariantVectorFieldImageIterator(CovariantVectorFieldImageStack[j],CovariantVectorFieldRegion);

		/*Assign each component of the covariant vector into their proper position on the gradient vector. The idea here is that the gradient matrix in the end will be:
		
		Let u = (u_x1,u_x2,...,u_xn) be the displacement vector at every pixel on the input displacement vector image

		Then the output vector matrix would be this matrix:
		
			| d(u_x1)/dx1	d(u_x1)/dx2		...		d(u_x1)/dxn |
			| d(u_x2)/dx1	d(u_x2)/dx2		...		d(u_x2)/dxn |
			|	...				...			...			...		|
			| d(u_xn)/dx1	d(u_xn)/dx2		...		d(u_xn)/dxn |
			
		The indices of this matrix is then converted to a vector using the formula i_new = j_old * Dimension + i_old
		
		*/
		
	/*	while(!GradientVectorFieldImageIterator.IsAtEnd() || !CovariantVectorFieldImageIterator.IsAtEnd()) {
			//Get the value in CovariantVectorFieldImageIterator
			CovariantPixel = CovariantVectorFieldImageIterator.Get();
			
			//Get the value in the Gradient
			GradientVectorPixel = GradientVectorFieldImageIterator.Get();
			
			//Transfer the CovariantPixel coordinates into an itkVector values
			for (int i = 0; i < Dimension; i++) {
				GradientVectorPixel[j*Dimension + i] = CovariantPixel[i];
			}

			//Set the new value into VectorFieldImageIterator
			GradientVectorFieldImageIterator.Set(GradientVectorPixel);
		}
		
	}*/



	/*typedef itk::VectorGradientMagnitudeImageFilter<VectorFieldImageType> VectorGradientMagnitudeImageFilterType;
	
	VectorGradientMagnitudeImageFilterType::Pointer vectorGradientMagnitudeImageFilter = VectorGradientMagnitudeImageFilterType::New();
	
	vectorGradientMagnitudeImageFilter->SetInput(vectorField);
	vectorGradientMagnitudeImageFilter->Update();
	
	ScalarImageType::Pointer vectorGradientMagnitudeImage = ScalarImageType::New();
	vectorGradientMagnitudeImage = vectorGradientMagnitudeImageFilter->GetOutput();*/
	
	
	//Spit out the physical coordinates and the values for both gradient filters... (?)
	
	
	
	
	
	
	
	
	//Output all appropriate outputs

	for (int j = 0; j < Dimension; j++) {	
		//Write out vector gradient field and its config file
		
		//Setup output vector file name and config file name
		std::ostringstream ss;
		ss << j;
		
		std::string currentDimension = ss.str();
		
		std::string outputConfigFileName = outConfigFile1 + "_x" + currentDimension + ".txt";
		std::string outputImageFileName = outputImageFile1 + "_x" + currentDimension + ".raw";
		
		

		//Write vector gradient field
		itkImageIORaw<CovariantPixelType, Dimension>* outputVectorField1 = new itkImageIORaw<CovariantPixelType, Dimension>;

		outputVectorField1->setConfigFileName(outputConfigFileName);

		//Write ConfigFile
		ImageIOData OutputVectorImageIOData;

		OutputVectorImageIOData.setPixelType("float");
		OutputVectorImageIOData.setPixelDimensionality(Dimension);
		OutputVectorImageIOData.setImageDimensionality(Dimension);

		const VectorFieldImageType::RegionType& OutputVectorRegion = CovariantVectorFieldImageStack[j]->GetLargestPossibleRegion();
		const VectorFieldImageType::SizeType& OutputVectorSize = OutputVectorRegion.GetSize();
		unsigned int tmpVectorSize[Dimension];
		for (int i = 0; i < Dimension; i++) {
			tmpVectorSize[i] = OutputVectorSize[i];
		}

		OutputVectorImageIOData.setDimensionSize((int *)tmpVectorSize);

		const VectorFieldImageType::SpacingType OutputVectorSpacing = CovariantVectorFieldImageStack[j]->GetSpacing();
		float tmpVectorSpacing[Dimension];
		for (int i = 0; i < Dimension; i++) {
			tmpVectorSpacing[i] = (float)OutputVectorSpacing[i];
		}
		
		OutputVectorImageIOData.setSpacingSize((float *)tmpVectorSpacing);

		const VectorFieldImageType::PointType& OutputVectorOrigin  = CovariantVectorFieldImageStack[j]->GetOrigin();
		float tmpVectorOrigin[Dimension];
		for (int i=0; i< Dimension; i++) {
			tmpVectorOrigin[i] = (float)OutputVectorOrigin[i];
		}

		OutputVectorImageIOData.setOrigin((float *)tmpVectorOrigin);

		OutputVectorImageIOData.setHeaderSize(0);
		OutputVectorImageIOData.setByteOrder("little");

		outputVectorField1->setImageIOData(OutputVectorImageIOData);

		outputVectorField1->writeConfigFile();

		//Write Displacement Field
		outputVectorField1->setitkImage((void *) &CovariantVectorFieldImageStack[j]);
		outputVectorField1->setImageFileName(outputImageFileName);
		outputVectorField1->writeRawImageKnown();

		delete outputVectorField1;
		
		std::cout << "Supposedly, vector gradient field is written in properly." << std::endl;
	}
		
	return EXIT_SUCCESS;

}