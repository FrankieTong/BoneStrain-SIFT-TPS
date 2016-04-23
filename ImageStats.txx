/*
File:			ImageStats.hxx

Description:	Function prints out some image statistics onto the command line that is useful as a sanity to see if the image has been modified properly.

				Will display the following image statistics:
				
				Image: 			Starting Index
								Size 
								Spacing 
								Origin 
				Image Inteisty:	Mean 
								Standard Deviation 
								Minimum 
								Maximum 

Inputs:			image - (ImageType *) Pointer that holds the ITK::Image that is to be analyzed

Author:			Frankie (Hoi-Ki) Tong, Sep 6 2014

Changes:		Frankie (Hoi-Ki) Tong, Sep 6 2014
				Initial creation of the file.
				
				Frankie (Hoi-Ki) Tong, Jun 8 2015
				Added functionality to print out image mean, std, min and max pixel values
*/

//#include "ImageStats.h"

template <typename ImageType>
void ImageStats(ImageType * image) {

	ImageType::RegionType region = image->GetLargestPossibleRegion();
	ImageType::IndexType start = region.GetIndex();
	ImageType::SizeType size = region.GetSize();
	ImageType::SpacingType spacing = image->GetSpacing();
	ImageType::PointType origin = image->GetOrigin();

	std::cout << "Start: " << start[0] << ", " << start[1] << ", " << start[2] << std::endl;
	std::cout << "Size: " << size[0] << " by " << size[1] << " by " << size[2] << std::endl;
	std::cout << "Spacing: x=" << spacing[0] << ", y=" << spacing[1] << ", z=" << spacing[2] << std::endl;
	std::cout << "Origin: x=" << origin[0] << ", y=" << origin[1] << ", z=" << origin[2] << std::endl;

	typedef itk::StatisticsImageFilter<ImageType> StatisticsImageFilterType;
	StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New();

	statisticsImageFilter->SetInput(image);

	try
	{
		statisticsImageFilter->Update();
	}
	catch( itk::ExceptionObject & excp )
	{
		std::cerr << "Exception thrown while executing statistics on image" << std::endl;
		std::cerr << excp << std::endl;
		return;
	}

	std::cout << "Mean: " << statisticsImageFilter->GetMean() << std::endl;
	std::cout << "Std.: " << statisticsImageFilter->GetSigma() << std::endl;
	std::cout << "Min: " << statisticsImageFilter->GetMinimum() << std::endl;
	std::cout << "Max: " << statisticsImageFilter->GetMaximum() << std::endl;

	return;
}