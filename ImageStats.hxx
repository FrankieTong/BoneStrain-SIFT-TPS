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