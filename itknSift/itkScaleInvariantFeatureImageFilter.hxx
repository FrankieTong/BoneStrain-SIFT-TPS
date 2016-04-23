/*=========================================================================
Program: ITK nSIFT Implemention (Template Source)
Module: $RCSfile: itkScaleInvariantFeatureImageFilter.txx,v $
Language: C++
Date: $Date: 2007/11/25 15:51:48 $
Version: $Revision: 1.0 $
Copyright (c) 2005,2006,2007 Warren Cheung
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
* Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.
* The name of the Insight Consortium, nor the names of any consortium members,
nor of any contributors, may be used to endorse or promote products derived
from this software without specific prior written permission.
* Modified source versions must be plainly marked as such, and must not be
misrepresented as being the original software.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
=========================================================================*/

#define VERBOSE
//#define DEBUG
//#define DEBUG_VERBOSE

//#define GENERATE_KEYS 
//#define REORIENT

//#define SIFT_FEATURE

//TODO:  Migrate the histogram reorientation code to the SIFT feature
// Make histogram feature special case of the sift feature with
// a single histogram.
//REORIENT NSIFT:  Keep a "global histogram" in addition to the specific ones.
// Or maybe do the global histogram/reorient phase first?

// Granularity of the Histogram causing cycles in reorientation accuracy?
// Or maybe we need to shift all the angles by half a bin width?

//Generate histograms
//compare gradient histograms
//

// Min DoG value for keypoints
// Fix Gaussian Scale
// Only iterate through the region up to 3 sigma from the point

//More advanced features for better matching?
//Simple quadrant system?


/* Based on example code from the ITK Toolkit
* TODO:  Use resampler+identity transform+spacing change instead of scaling
* arbitrary downscale between pyramid levels
* may need to have a threshold
* Gaussian filtration may fail if input image format is integer
* Generate Pointset
* Get Orientation
* Generate Features

* Test vs Noise (no need to for orientation)
* vs Stretch 
* vs Scale
*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#include "itkScaleInvariantFeatureImageFilter.h"
#include <cstdio>

#ifndef SIFTKEY_CLASS
#define SIFTKEY_CLASS


const float PI=3.14159265;

namespace itk
{  

	template <class TFixedImageType, unsigned int VDimension> 

	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::ScaleInvariantFeatureImageFilter() 
	{
#ifdef DO_DOUBLE
		m_ImageScalesTestedNumber = 4;
#else
		m_ImageScalesTestedNumber = 3;
#endif
		m_ScalingFactor = 2.0;
		m_DifferenceOfGaussianTestsNumber = 3;
#ifdef DO_DOUBLE
		m_DoubleOriginalImage = true;
#else
		m_DoubleOriginalImage = false;
#endif
		m_HistogramBinsNumber = 36;      
		m_ErrorThreshold = 0.0;
		m_MaxFeatureDistanceRatio = 0.8;
		m_GaussianSigma = 1.5;  
		m_MinKeypointValue = 0.0075;
		m_SIFTHalfWidth = 8;  // This MUST be a multiple of m_SIFTSubfeatureWidth
		m_SIFTSubfeatureWidth = 4;
		m_SIFTSubfeatureBins = 8;

		// Derived from above
		m_DifferenceOfGaussianImagesNumber = m_DifferenceOfGaussianTestsNumber+2;
		m_GaussianImagesNumber = m_DifferenceOfGaussianImagesNumber+1;
		m_IdentityTransform = IdentityTransformType::New();

#ifdef BAD_POINT_REJECTION
		m_DoGThreshold = 0.03; //Value should be adjusted based on image intensity [0-1]
		m_EdgeThreshold = 10;
		m_RatioEdgeThreshold = 0;
#endif
	}

	
	template <class TFixedImageType, unsigned int VDimension> 
	double
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::GetGaussianScale( int j ) 
	{
		return (pow(2, (double) j / (double) m_DifferenceOfGaussianTestsNumber) * m_GaussianSigma);
	}

	template <class TFixedImageType, unsigned int VDimension> 
	unsigned int
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::HistFeatureSize() 
	{
		unsigned int size = 1;
		// Have m_HistogramBinsNumber for each of the (VDimension-1) Orientation dimensions
		for (unsigned int i = 1; i < VDimension; ++i)
		{
			size *= m_HistogramBinsNumber;
		}
		return size;
	}

	template <class TFixedImageType, unsigned int VDimension> 
	typename ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>::GradientImageType::Pointer
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::GetHypersphericalCoordinates(typename GradientImageType::Pointer inputImg) 
	{
		typedef itk::ImageRegionIteratorWithIndex< GradientImageType > ImageIteratorType;
		typedef itk::ImageRegionConstIteratorWithIndex< GradientImageType > ConstImageIteratorType;
		typedef itk::ImageRegionConstIteratorWithIndex< TFixedImageType > ConstFixedImageIteratorType;
		typename GradientImageType::Pointer outputImg = GradientImageType::New();
		// Copy attributes
		outputImg->SetRegions(inputImg->GetLargestPossibleRegion());
		outputImg->CopyInformation( inputImg );
		outputImg->Allocate();

		ConstImageIteratorType inputIt(inputImg, inputImg->GetLargestPossibleRegion());
		ImageIteratorType outputIt(outputImg, inputImg->GetLargestPossibleRegion());
		
		for ( inputIt.GoToBegin(), outputIt.GoToBegin(); !inputIt.IsAtEnd();
		++inputIt, ++outputIt)
		{
			typename GradientImageType::PixelType x =  inputIt.Get();
			typename GradientImageType::PixelType p;

			// position 0 is the norm
			p[0] = x.GetNorm();

			// position 1 is arctan (x0 / x1)
			p[1] = atan2( x[0],x[1] );

			// Iterate over all the positions
			// position k  is arctan (x_k-1 / (x_k * cos p_k))	  
			for (unsigned int k = 2; k < x.Size(); ++k)
			{
				p[k] = atan2( x[k-1], x[k] * cos(p[k-1]));
			}
			outputIt.Set(p);
		}

		return outputImg;
	}

	// Generates a vector with positions 1 to VDimension-1 filled with
	// histogram bin numbers from a single bin number

	template <class TFixedImageType, unsigned int VDimension> 
	typename ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>::GradientImageType::PixelType
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::BinToVector (unsigned int maxbin) 
	{
		// convert maxpos to orientation bin vector
		typename GradientImageType::PixelType maxp;
		for (unsigned int i = 1; i < VDimension; ++i) {
			maxp[i] = maxbin % m_HistogramBinsNumber;
			maxbin /= m_HistogramBinsNumber;
		}
		return maxp;
	}           

	template <class TFixedImageType, unsigned int VDimension> 
	unsigned int
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::SiftFeatureSize() 
	{
		unsigned int size = 1;
		// Have m_HistogramBinsNumber for each of the (VDimension-1) Orientation dimensions
		for (unsigned int i = 0; i < VDimension; ++i)
		{
			size *= (m_SIFTHalfWidth * 2 / m_SIFTSubfeatureWidth );
			if (i > 0)
			size *= m_SIFTSubfeatureBins;
		}
		
		return size;
	}

	// Convert the delta iterator into index to the
	// start of the SIFT histogram
	template <class TFixedImageType, unsigned int VDimension> 
	unsigned int
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::DeltaToSiftIndex (int delta[]) 
	{
		unsigned int bin = 0;
		unsigned int binpos = 1;
#ifdef DEBUG_VERBOSE
		std::cerr << "Converting delta: ";
#endif
		for (unsigned int i = 0; i < VDimension; ++i) {
#ifdef DEBUG_VERBOSE
			std::cerr << delta[i];
#endif
			unsigned int tmp =  (delta[i] + m_SIFTHalfWidth) / m_SIFTSubfeatureWidth;
			
			bin += tmp * binpos;
			binpos *= (m_SIFTHalfWidth * 2 / m_SIFTSubfeatureWidth );
		}
		for (unsigned int i = 1; i < VDimension; ++i)
		bin *= m_SIFTSubfeatureBins;
		
#ifdef DEBUG_VERBOSE
		std::cerr << "\n";
#endif
		return bin;
	}


#ifdef GENERATE_KEYS
	template <class TFixedImageType, unsigned int VDimension> 
	typename ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>::FeatureType
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::GetSiftKey(typename GradientImageType::Pointer inputImg,
	FixedImagePointer multImg,
	IndexType pixelIndex) 
	{
#ifdef DEBUG_VERBOSE
		std::cerr << "GetSiftKey..." << std::endl; 
#endif
		FeatureType sifthistogram(this->SiftFeatureSize());
		sifthistogram.Fill(0);

		// delta iterates from  -m_SIFTHalfWidth to m_SIFTHalfWidth-1 
		// in each dimensions
		int delta[VDimension];
		for (int k = 0; k < VDimension; ++k) {
			delta[k] = -m_SIFTHalfWidth;
		}

		typename GradientImageType::SizeType regionSize = 
		inputImg->GetLargestPossibleRegion().GetSize();	
		
		while(1) {
			unsigned int siftbin = this->DeltaToSiftIndex(delta);
			
#ifdef DEBUG_VERBOSE
			std::cerr << "Siftbin:" << siftbin << std::endl; 
#endif
			
			// Get pixel index
			// Clamp to image edges
			IndexType tmpIndex;
			for (int k=0; k < VDimension; ++k) {
				if ((pixelIndex[k] + delta[k]) < 0) {
					tmpIndex[k] = 0;
				} else {
					tmpIndex[k] = pixelIndex[k] + delta[k];
					if (tmpIndex[k] >= regionSize[k])
					tmpIndex[k] = regionSize[k]-1;
				}
			}
			
#ifdef DEBUG_VERBOSE
			std::cerr << "Pixel:" << tmpIndex << std::endl; 
#endif
			typename GradientImageType::PixelType x = 
			inputImg->GetPixel(tmpIndex);
			
			// Get histogram bin
			// Iterate over all the positions
			unsigned int bin = 0;
			unsigned int binpos = 1;
			for (unsigned int k = 1; k < x.Size(); ++k)
			{
				// Rescale from -PI to PI ->  0 to m_HistogramBinsNumber-1
				float p;
				p = (x[k] + PI)  * (float) m_SIFTSubfeatureBins / (2.0 * PI);
				
				
				if (p < 0 || p >= m_SIFTSubfeatureBins) 
				p = 0;
				bin += (unsigned int) p * binpos;
				
#ifdef DEBUG_VERBOSE
				std::cout << " " << p;
#endif
				binpos *= m_SIFTSubfeatureBins;
			}
			
			bin += siftbin;
			
			// Fill Sift Index bin
			if (bin > this->SiftFeatureSize()) {
				// VERY BAD
				std::cerr << bin << " > " << this->SiftFeatureSize() << " Warning -- Overload2\n";
			}
			sifthistogram[bin] += x[0] * multImg->GetPixel(tmpIndex);
			
#ifdef DEBUG_VERBOSE
			std::cerr << "Incrementing\n";
#endif	  
			// Increment delta
			bool resetdelta=false;
			for(int k=0; k <= VDimension; ++k) {
#ifdef DEBUG_VERBOSE
				std::cerr << delta[k];
#endif
				if (k == VDimension) {
#ifdef DEBUG_VERBOSE
					std::cerr << "done\n";
#endif
					resetdelta = true;
					break; // done
				}
				// Don't want to go past m_SIFTHalfWidth-1
				if (++delta[k] < (int) m_SIFTHalfWidth) {
					break;
				}
				delta[k] = -m_SIFTHalfWidth; // reset and increment the next pos
			}
			if(resetdelta) break;	  
		}
		
#ifdef DEBUG_VERBOSE
		std::cerr << "SIFT key: " << sifthistogram << "\n";
#endif
		return(sifthistogram);
	}


	// Takes a hyperspherical coordinate gradient and gaussian weights
	// returns a histogram 
	// Each orientation divides the 2PI angles into m_HistogramBinsNumber 
	// The value in each bin is the weighted magnitude
	template <class TFixedImageType, unsigned int VDimension> 
	typename ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>::FeatureType
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::GetHistogram(typename GradientImageType::Pointer inputImg,
	FixedImagePointer multImg) 
	{
		
#ifdef ERROR_CHECK
		std::cerr << "GetHistogram ... ";
#endif
		
		FeatureType histogram(this->HistFeatureSize());

		histogram.Fill(0);
		
		typedef itk::ImageRegionConstIteratorWithIndex< GradientImageType > ConstImageIteratorType;
		typedef itk::ImageRegionConstIteratorWithIndex< TFixedImageType > ConstFixedImageIteratorType;
		
		ConstImageIteratorType inputIt(inputImg, inputImg->GetLargestPossibleRegion());
		ConstFixedImageIteratorType multIt( multImg, inputImg->GetLargestPossibleRegion());
		
		for ( inputIt.GoToBegin(), multIt.GoToBegin(); !inputIt.IsAtEnd();
		++inputIt, ++multIt)
		{
			typename GradientImageType::PixelType x =  inputIt.Get();
			typename TFixedImageType::PixelType m = multIt.Get();
			unsigned int bin = 0;
			unsigned int binpos = 1;
			typename GradientImageType::PixelType p;
			
			// position 0 is the norm
			p[0] = x[0];
			
			//if (std::isnan(p[0]) || (p[0] == 0.0))
			if (p[0] = NAN || (p[0] == 0.0))
			continue;
			
			// multiply by m
			p[0] *= m;
			
#ifdef DEBUG_VERBOSE
			std::cout << "Bin: ";
#endif	  
			// Iterate over all the positions
			for (unsigned int k = 1; k < x.Size(); ++k)
			{
				// Rescale from -PI to PI ->  0 to m_HistogramBinsNumber-1
				p[k] = (x[k] + PI)  * (float) m_HistogramBinsNumber / (2.0 * PI);
				
				
				if (p[k] < 0 || p[k] >= m_HistogramBinsNumber) 
				p[k] = 0;
				bin += (unsigned int) p[k] * binpos;
				
#ifdef DEBUG_VERBOSE
				std::cout << " " << p[k];
#endif
				binpos *= m_HistogramBinsNumber;
			}
#ifdef DEBUG_VERBOSE
			std::cout << " Value: " << p[0] << std::endl;
#endif
			if (bin > this->HistFeatureSize()) {
				// VERY BAD
				std::cerr << x << " -> " << p << "\n";
				std::cerr << bin << " > " << this->HistFeatureSize() << " Warning -- Overload2\n";
			}
			histogram[bin] += p[0];
		}
#ifdef DEBUG
		// Print the Histogram
		std::cout << histogram << std::endl;
#endif
		
		// Since we are going to use this as a feature
		// Normalise
		float hmag = 0.0;
#ifdef REORIENT
		float maxmag = -1;
		unsigned int maxbin;
#endif
		for (unsigned int i = 0; i < this->HistFeatureSize(); ++i) {
			float mag = histogram[i]*histogram[i];
			hmag += mag;
#ifdef REORIENT
			if (maxmag < 0 || mag > maxmag) {
				maxmag = mag;
				maxbin = i;
			}
#endif
		}
		hmag = sqrt(hmag);
		
#ifdef REORIENT
		typename GradientImageType::PixelType maxp = this->BinToVector(maxbin);
		
		FeatureType histogram2(this->HistFeatureSize());
		histogram2.Fill(0);
#endif
		
		for (unsigned int i = 0; i < this->HistFeatureSize(); ++i) {
			histogram[i] /= hmag;
			
#ifdef REORIENT
			typename GradientImageType::PixelType bini = this->BinToVector(i);
			
			unsigned int bin = 0;
			unsigned int binpos = 1;
			for (unsigned int k = 1; k < VDimension; ++k) {
				bini[k] = ((int) (bini[k] - maxp[k] + m_HistogramBinsNumber)) % m_HistogramBinsNumber;
				if (bini[k] < 0 || bini[k] >= m_HistogramBinsNumber)
				bini[k] = 0;
				bin += (int) bini[k] * binpos;
				binpos *= m_HistogramBinsNumber;
			}
			histogram2[bin] = histogram[i];
#endif
		}
		
		
#ifdef DEBUG
		std::cout << histogram << std::endl;
#ifdef REORIENT
		std::cout << "Max Bin: " << maxbin << " Max Mag: " << maxmag << std::endl;
		std::cout << "Reoriented: " << histogram2 << std::endl;
#endif
#endif
#ifdef ERROR_CHECK
		std::cerr << "OK\n ";
#endif
		
#ifdef REORIENT
		return histogram2;
#else
		return histogram;
#endif
	}
#endif

	template <class TFixedImageType, unsigned int VDimension> 
	void
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::CheckLocalExtrema(FixedImagePointer image, IndexType pixelIndex,
	PixelType pixelValue,
	bool &isMax, bool &isMin,
	bool checkCentre)
	{
#ifdef ERROR_CHECK
		std::cerr << "CheckLocalExtrema ... ";
#endif
		int delta[VDimension];
		for (int k = 0; k < VDimension; ++k) {
			delta[k] = -1;
		}
		while(1) {
			bool isZero=true;
			if (!checkCentre) {
				for (int k=0; k < VDimension; ++k) {
					if(delta[k] != 0) {
						isZero = false;
						break;
					}
				}
			}
			
			if (checkCentre || !isZero) {
				// Check if not the centre
				IndexType tmpIndex;
				for (int k=0; k < VDimension; ++k) {
					tmpIndex[k] = pixelIndex[k] + delta[k];
				}
				
				typename TFixedImageType::PixelType tmpValue = 
				image->GetPixel(tmpIndex);
				
#ifdef DEBUG_VERBOSE
				std::cout << "...Comparing to ( ";
				for (int k = 0; k < VDimension; ++k)
				std::cout << tmpIndex[k] << " ";
				std::cout << ") = " << tmpValue << "\n";
#endif  
				// Treat as equality if within the error bound
				if (((tmpValue - pixelValue) <= m_ErrorThreshold) && 
						((tmpValue - pixelValue) >= -m_ErrorThreshold)) {
					isMax = false;
					isMin = false;
				} else	
				if (tmpValue > pixelValue) {
					isMax = false;
				} else
				if (tmpValue < pixelValue) {
					isMin = false;
				}
				if (!isMax && !isMin) break;
			}
			// Increment delta
			bool resetdelta=false;
			for(int k=0; k <= VDimension; ++k) {
				if (k == VDimension) {
					resetdelta = true;
					break; // done
				}
				if (delta[k] < 1) {
					++delta[k];
					break;
				}
				delta[k] = -1; // reset and increment the next pos
			}
			if(resetdelta) break;	  
		}  
#ifdef ERROR_CHECK
		std::cerr << "OK\n";
#endif
	}

#ifdef BAD_POINT_REJECTION

	
	/*template <class TFixedImageType, unsigned int VDimension> 
	void
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::CorrectBadFeaturePointScaleSpace(
			FixedImagePointer previousScaleImage,
			FixedImagePointer currentScaleImage,
			FixedImagePointer nextScaleImage,
			IndexType pixelIndex,
			PointType* pointPtr,
			bool & isGoodPointPtr)
	{

		//Might want to move all this up to main function anyway. The hessian and vector calculations can stay in seperate functions.

		PointType point = *pointPtr;
		bool isGoodPoint = *isGoodPoint;

		//Create Hessian Function and gradient vector
		m_ScaleSpaceHessianMatrixType hessianMatrix;
		m_ScaleSpaceVectorType gradientVector;


		//Calculate local gradient vector values. Values in this vector is arrages as such (Dx0, Dx1, ..., Dxn, Dsg) where Dsg = Dsigma
		this->CalculateFirstOrderGradientScaleSpace(previousScaleImage, currentScaleImage, nextScaleImage, pixelIndex, &gradientVector);

		//Calculate local Hessian Matrix. Matrix is defined in a similar fashion to the gradient vector:

		//| Dx0x0 Dx0x1 ... Dx0xn Dx0sg |
		//| Dx1x0 Dx1x1 ... Dx1xn Dx1sg |
		//|             ...             |
		//| Dxnx0 Dxnx1 ... Dxnxn Dxnsg |
		//| Dsgx0 Dsgx1 ... Dsgxn Dsgsg |

		this->CalculateHessianMatrixScaleSpace(previousScaleImage, currentScaleImage, nextScaleImage, pixelIndex, &hessianMatrix);

		//Get inverse of Hessian Matrix
		m_ScaleSpaceHessianMatrixType inverseHessianMatrix = hessianMatrix.GetInverse();
		
		//Calculate displacement error
		m_ScaleSpaceVectorType displacementError = -(inverseHessianMatrix*gradientVector);

		//Determine if the displacementError is greater than 0.5 in each direction. If it is, shift the point from current index to new index.
		//IndexType newPixelIndex = pixelIndex;
		bool changeInIndex = false;

		for (int i = 0; i < VDimension) {

			if  (displacementError[i] > 0.5) {
				pixelIndex[i] = pixelIndex[i] + 1;
				changeInIndex = true;
			}

			if (displacementError[i] < 0.5) {
				pixelIndex[i] = pixelIndex[i] - 1;
				changeInIndex = true;
			}
		}

		//Calculate change in scale space values if scale space index is greater than 0.5 in each direction

		//Recalculate gradient and hessian and adjust point if there is a change in index

		if (changeInIndex) {
			this->CalculateFirstOrderGradientScaleSpace(previousScaleImage, currentScaleImage, nextScaleImage, pixelIndex, &gradientVector);
			this->CalculateHessianMatrixScaleSpace(previousScaleImage, currentScaleImage, nextScaleImage, pixelIndex, &hessianMatrix);
			currentScaleImage->TransformIndexToPhysicalPoint (pixelIndex, point);
		}

		//Calculate scale space value at extremum to determine if we should keep the point
		
		double scaleSpaceValue = (double)(currentScaleImage->GetPixel(pixelIndex) + 0.5*(gradientVector*point.GetVnlVector()));

		isGoodPoint = fabs(scaleSpaceValue) >= fabs(m_DoGThreshold);


		//Remove bad edge responses
		

		PixelType lHessianTrace2 = (dxx + dyy) * (dxx + dyy);
		PixelType lHessianDet = dxx * dyy - dxy * dxy;

	}*/

	/*template <class TFixedImageType, unsigned int VDimension> 
	void
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::CalculateFirstOrderGradientScaleSpace(
	FixedImagePointer previousScaleImage,
	FixedImagePointer currentScaleImage,
	FixedImagePointer nextScaleImage,
	IndexType pixelIndex,
	m_ScaleSpaceVectorType * gradientVectorPtr){

		//Calculte local gradient vector values for current image scale
		for (int i = 0; i < VDimension; i++) {
			
			//Get index value of current point +1/-1 in current scaleSpaceDimension
			IndexType positiveIndex = pixelIndex;
			positiveIndex[i] = pixelIndex[i] + 1;

			IndexType negativeIndex = pixelIndex;
			negativeIndex[i] = pixelIndex[i] - 1;

			//Get the values of the points inside the indexes
			typename TFixedImageType::PixelType positiveValue = currentScaleImage->GetPixel(positiveIndex);
			typename TFixedImageType::PixelType negativeValue = currentScaleImage->GetPixel(negativeIndex);

			//Get the gradient value for current scaleSpaceDimension
			(*gradientVectorPtr)[i] = (double)(0.5 * (positiveValue - negativeValue));
		}

		//Calculate local gradient vector between image scales
		typename TFixedImageType::PixelType positiveValue = nextScaleImage->GetPixel(pixelIndex);
		typename TFixedImageType::PixelType negativeValue = previousScaleImage->GetPixel(pixelIndex);
		(*gradientVectorPtr)[VDimension] = (double)(0.5 * (positiveValue - negativeValue));

		return

	}*/

	/*template <class TFixedImageType, unsigned int VDimension> 
	void
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::CalculateHessianMatrixScaleSpace(
	FixedImagePointer previousScaleImage,
	FixedImagePointer currentScaleImage,
	FixedImagePointer nextScaleImage,
	IndexType pixelIndex,
	m_ScaleSpaceHessianMatrixType * hessianMatrixPtr){

		//First calculate the Hessian values for current scale space only:
		for (int i = 0; i < VDimension+1; i++) {

			for (int j = 0; j < VDimension+1; j++) {

				if (i == j) {
					//Calculate local second derivative when derivative is taken to the same variable

					if (i == VDimension) {

						//Calculaitng second derivative when the derivative is the scale space variable (twice)
						typename TFixedImageType::PixelType positiveValue = nextScaleImage->GetPixel(pixelIndex);
						typename TFixedImageType::PixelType negativeValue = previousScaleImage->GetPixel(pixelIndex);
						typename TFixedImageType::PixelType currentValue = currentScaleImage->GetPixel(pixelIndex);

						//std::cout<<positiveValue<<", "<<currentValue<<", "<<negativeValue<<std::endl;

						(*hessianMatrixPtr)(i,i) = (double)(0.25 *  positiveValue - 2 * currentValue + negativeValue);
					
					} else {

						//Get index value of current point +1/-1 in current scaleSpaceDimension
						IndexType positiveIndex = pixelIndex;
						positiveIndex[i] = pixelIndex[i] + 1;

						IndexType negativeIndex = pixelIndex;
						negativeIndex[i] = pixelIndex[i] - 1;

						//Get the values of the points inside the indexes
						typename TFixedImageType::PixelType positiveValue = currentScaleImage->GetPixel(positiveIndex);
						typename TFixedImageType::PixelType negativeValue = currentScaleImage->GetPixel(negativeIndex);
						typename TFixedImageType::PixelType currentValue = currentScaleImage->GetPixel(pixelIndex);

						//std::cout<<positiveValue<<", "<<currentValue<<", "<<negativeValue<<std::endl;

						(*hessianMatrixPtr)(i,j) = positiveValue - 2 * currentValue + negativeValue;
					}

				} else {
					//Calculate local second derivative when derivative is taken to different variables
					
					typename TFixedImageType::PixelType ipositive_jpositive_Value;
					typename TFixedImageType::PixelType inegative_jnegative_Value;
					typename TFixedImageType::PixelType ipositive_jnegative_Value;
					typename TFixedImageType::PixelType inegative_jpositive_Value;

					if (i == VDimension) {

						IndexType ipositive_jpositive_Index = pixelIndex;
						//ipositive_jpositive_Index[i] = pixelIndex[i];
						ipositive_jpositive_Index[j] = pixelIndex[j] + 1;

						IndexType inegative_jnegative_Index = pixelIndex;
						//inegative_jnegative_Index[i] = pixelIndex[i];
						inegative_jnegative_Index[j] = pixelIndex[j] - 1;

						IndexType ipositive_jnegative_Index = pixelIndex;
						//ipositive_jnegative_Index[i] = pixelIndex[i];
						ipositive_jnegative_Index[j] = pixelIndex[j] - 1;

						IndexType inegative_jpositive_Index = pixelIndex;
						//inegative_jpositive_Index[i] = pixelIndex[i];
						inegative_jpositive_Index[j] = pixelIndex[j] + 1;

						//Get the values of the points inside the indexes
						ipositive_jpositive_Value = nextScaleImage->GetPixel(ipositive_jpositive_Index);
						inegative_jnegative_Value = previousScaleImage->GetPixel(inegative_jnegative_Index);
						ipositive_jnegative_Value = nextScaleImage->GetPixel(ipositive_jnegative_Index);
						inegative_jpositive_Value = previousScaleImage->GetPixel(inegative_jpositive_Index);

					}

					else if (j == VDimension) {

						IndexType ipositive_jpositive_Index = pixelIndex;
						ipositive_jpositive_Index[i] = pixelIndex[i] + 1;
						//ipositive_jpositive_Index[j] = pixelIndex[j];

						IndexType inegative_jnegative_Index = pixelIndex;
						inegative_jnegative_Index[i] = pixelIndex[i] - 1;
						//inegative_jnegative_Index[j] = pixelIndex[j];

						IndexType ipositive_jnegative_Index = pixelIndex;
						ipositive_jnegative_Index[i] = pixelIndex[i] + 1;
						//ipositive_jnegative_Index[j] = pixelIndex[j];

						IndexType inegative_jpositive_Index = pixelIndex;
						inegative_jpositive_Index[i] = pixelIndex[i] - 1;
						//inegative_jpositive_Index[j] = pixelIndex[j];

						//Get the values of the points inside the indexes
						ipositive_jpositive_Value = nextScaleImage->GetPixel(ipositive_jpositive_Index);
						inegative_jnegative_Value = previousScaleImage->GetPixel(inegative_jnegative_Index);
						ipositive_jnegative_Value = previousScaleImage->GetPixel(ipositive_jnegative_Index);
						inegative_jpositive_Value = nextScaleImage->GetPixel(inegative_jpositive_Index);

					} else {

						//Get index value of current point (+1,+1)/(+1,-1)/(-1,+1)/(-1,-1) in current scaleSpaceDimension
						IndexType ipositive_jpositive_Index = pixelIndex;
						ipositive_jpositive_Index[i] = pixelIndex[i] + 1;
						ipositive_jpositive_Index[j] = pixelIndex[j] + 1;

						IndexType inegative_jnegative_Index = pixelIndex;
						inegative_jnegative_Index[i] = pixelIndex[i] - 1;
						inegative_jnegative_Index[j] = pixelIndex[j] - 1;

						IndexType ipositive_jnegative_Index = pixelIndex;
						ipositive_jnegative_Index[i] = pixelIndex[i] + 1;
						ipositive_jnegative_Index[j] = pixelIndex[j] - 1;

						IndexType inegative_jpositive_Index = pixelIndex;
						inegative_jpositive_Index[i] = pixelIndex[i] - 1;
						inegative_jpositive_Index[j] = pixelIndex[j] + 1;

						//Get the values of the points inside the indexes
						ipositive_jpositive_Value = currentScaleImage->GetPixel(ipositive_jpositive_Index);
						inegative_jnegative_Value = currentScaleImage->GetPixel(inegative_jnegative_Index);
						ipositive_jnegative_Value = currentScaleImage->GetPixel(ipositive_jnegative_Index);
						inegative_jpositive_Value = currentScaleImage->GetPixel(inegative_jpositive_Index);
					}
					
					//std::cout<<ipositive_jpositive_Value<<", "<<inegative_jnegative_Value<<", "<<ipositive_jnegative_Value<<", "<<inegative_jpositive_Value<<std::endl;
					(*hessianMatrixPtr)(i,j) = (double)(0.25 * (ipositive_jpositive_Value + inegative_jnegative_Value - ipositive_jnegative_Value - inegative_jpositive_Value));
				}
			}
		}
		//std::cout << "M1: " << (*hessianMatrixPtr) << std::endl;
		return;
	}*/
	








/*
	template <class TFixedImageType, unsigned int VDimension> 
	void
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::CorrectBadFeaturePoint(
			FixedImagePointer currentScaleImage,
			IndexType pixelIndex,
			PointType* pointPtr,
			bool & isGoodPointPtr)
	{


		PointType point = *pointPtr;
		bool isGoodPoint = *isGoodPoint;

		//Create Hessian Function and gradient vector
		m_HessianMatrixType hessianMatrix;
		m_VectorType gradientVector;


		//Calculate local gradient vector values. Values in this vector is arrages as such (Dx0, Dx1, ..., Dxn, Dsg) where Dsg = Dsigma
		this->CalculateFirstOrderGradient(currentScaleImage, pixelIndex, &gradientVector);

		//Calculate local Hessian Matrix. Matrix is defined in a similar fashion to the gradient vector:

		//| Dx0x0 Dx0x1 ... Dx0xn Dx0sg |
		//| Dx1x0 Dx1x1 ... Dx1xn Dx1sg |
		//|             ...             |
		//| Dxnx0 Dxnx1 ... Dxnxn Dxnsg |
		//| Dsgx0 Dsgx1 ... Dsgxn Dsgsg |

		this->CalculateHessianMatrix(currentScaleImage, pixelIndex, &hessianMatrix);

		//Get inverse of Hessian Matrix
		m_HessianMatrixType inverseHessianMatrix = hessianMatrix.GetInverse();
		
		//Calculate displacement error
		m_VectorType displacementError = -(inverseHessianMatrix*gradientVector);

		//Determine if the displacementError is greater than 0.5 in each direction. If it is, shift the point from current index to new index.
		//IndexType newPixelIndex = pixelIndex;
		bool changeInIndex = false;

		for (int i = 0; i < VDimension) {

			if  (displacementError[i] > 0.5) {
				pixelIndex[i] = pixelIndex[i] + 1;
				changeInIndex = true;
			}

			if (displacementError[i] < 0.5) {
				pixelIndex[i] = pixelIndex[i] - 1;
				changeInIndex = true;
			}
		}

		//Recalculate gradient and hessian and adjust point if there is a change in index

		if (changeInIndex) {
			this->CalculateFirstOrderGradient(previousScaleImage, currentScaleImage, nextScaleImage, pixelIndex, &gradientVector);
			this->CalculateHessianMatrix(previousScaleImage, currentScaleImage, nextScaleImage, pixelIndex, &hessianMatrix);
			currentScaleImage->TransformIndexToPhysicalPoint (pixelIndex, point);
		}

		//Calculate scale space value at extremum to determine if we should keep the point
		
		double scaleSpaceValue = (double)(currentScaleImage->GetPixel(pixelIndex) + 0.5*(gradientVector*point.GetVnlVector()));

		isGoodPoint = fabs(scaleSpaceValue) >= fabs(m_DoGThreshold);


		//Remove bad edge responses
		

		
	}*/

	template <class TFixedImageType, unsigned int VDimension> 
	void
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::CalculateFirstOrderGradient(
	FixedImagePointer currentScaleImage,
	IndexType pixelIndex,
	m_VectorType * gradientVectorPtr){

		//Calculate local gradient vector values for current image scale
		for (int i = 0; i < VDimension; i++) {
			
			//Get index value of current point +1/-1 in current scaleSpaceDimension
			IndexType positiveIndex = pixelIndex;
			positiveIndex[i] = pixelIndex[i] + 1;

			IndexType negativeIndex = pixelIndex;
			negativeIndex[i] = pixelIndex[i] - 1;

			//Get the values of the points inside the indexes
			typename TFixedImageType::PixelType positiveValue = currentScaleImage->GetPixel(positiveIndex);
			typename TFixedImageType::PixelType negativeValue = currentScaleImage->GetPixel(negativeIndex);

			//Get the gradient value for current scaleSpaceDimension
			(*gradientVectorPtr)[i] = (double)(0.5 * (positiveValue - negativeValue));
		}

		return;

	}

	template <class TFixedImageType, unsigned int VDimension> 
	void
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::CalculateHessianMatrix(
	FixedImagePointer currentScaleImage,
	IndexType pixelIndex,
	m_HessianMatrixType * hessianMatrixPtr){

		//First calculate the Hessian values for current scale space only:
		for (int i = 0; i < VDimension; i++) {

			for (int j = 0; j < VDimension; j++) {

				if (i == j) {
					//Calculate local second derivative when derivative is taken to the same variable

					//Get index value of current point +1/-1 in current scaleSpaceDimension
					IndexType positiveIndex = pixelIndex;
					positiveIndex[i] = pixelIndex[i] + 1;

					IndexType negativeIndex = pixelIndex;
					negativeIndex[i] = pixelIndex[i] - 1;

					//Get the values of the points inside the indexes
					typename TFixedImageType::PixelType positiveValue = currentScaleImage->GetPixel(positiveIndex);
					typename TFixedImageType::PixelType negativeValue = currentScaleImage->GetPixel(negativeIndex);
					typename TFixedImageType::PixelType currentValue = currentScaleImage->GetPixel(pixelIndex);

					(*hessianMatrixPtr)(i,j) = positiveValue - 2 * currentValue + negativeValue;

				} else {
					//Calculate local second derivative when derivative is taken to different variables
					
					typename TFixedImageType::PixelType ipositive_jpositive_Value;
					typename TFixedImageType::PixelType inegative_jnegative_Value;
					typename TFixedImageType::PixelType ipositive_jnegative_Value;
					typename TFixedImageType::PixelType inegative_jpositive_Value;

				
					//Get index value of current point (+1,+1)/(+1,-1)/(-1,+1)/(-1,-1) in current scaleSpaceDimension
					IndexType ipositive_jpositive_Index = pixelIndex;
					ipositive_jpositive_Index[i] = pixelIndex[i] + 1;
					ipositive_jpositive_Index[j] = pixelIndex[j] + 1;

					IndexType inegative_jnegative_Index = pixelIndex;
					inegative_jnegative_Index[i] = pixelIndex[i] - 1;
					inegative_jnegative_Index[j] = pixelIndex[j] - 1;

					IndexType ipositive_jnegative_Index = pixelIndex;
					ipositive_jnegative_Index[i] = pixelIndex[i] + 1;
					ipositive_jnegative_Index[j] = pixelIndex[j] - 1;

					IndexType inegative_jpositive_Index = pixelIndex;
					inegative_jpositive_Index[i] = pixelIndex[i] - 1;
					inegative_jpositive_Index[j] = pixelIndex[j] + 1;

					//Get the values of the points inside the indexes
					ipositive_jpositive_Value = currentScaleImage->GetPixel(ipositive_jpositive_Index);
					inegative_jnegative_Value = currentScaleImage->GetPixel(inegative_jnegative_Index);
					ipositive_jnegative_Value = currentScaleImage->GetPixel(ipositive_jnegative_Index);
					inegative_jpositive_Value = currentScaleImage->GetPixel(inegative_jpositive_Index);
					
					(*hessianMatrixPtr)(i,j) = (double)(0.25 * (ipositive_jpositive_Value + inegative_jnegative_Value - ipositive_jnegative_Value - inegative_jpositive_Value));
				}
			}
		}
		return;
	}
	

#endif

#ifdef GENERATE_KEYS
	// Input:  Image, Gradient Image, Point
	// Output:  Vector of direction
	template <class TFixedImageType, unsigned int VDimension> 
	typename ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>::FeatureType
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::GetFeatures( FixedImagePointer fixedImage, 
	typename GradientImageType::Pointer hgradImage,
	PointType &point, float currScale) 
	{
#ifdef ERROR_CHECK
		std::cerr << "GetFeatures ... ";
#endif
		
		// Generate the Gaussian
		typedef GaussianImageSource<TFixedImageType> GaussianImageSourceType;
		typename GaussianImageSourceType::Pointer gaussImgSource;
		typename GaussianImageSourceType::ArrayType sigma;

		gaussImgSource = GaussianImageSourceType::New();
		gaussImgSource->SetNormalized(true);
		//gaussImgSource->SetNormalized(false);
		//gaussImgSource->SetScale(255);  
		gaussImgSource->SetSpacing(fixedImage->GetSpacing());
		gaussImgSource->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
		gaussImgSource->SetMean(point);
		
#if 0
		// If we wanted to find the orientation,  we would use this
		for (int i = 0; i < VDimension; ++i)
		sigma[i] = 3.0 * currScale;
		gaussImgSource->SetSigma(sigma);
#endif
		
		// Simulate the 16x16 Gaussian window descriptor
		// use sigma equal to half the descriptor window width
		for (int i = 0; i < VDimension; ++i)
		sigma[i] = 8.0;
		gaussImgSource->SetSigma(sigma);
		
		gaussImgSource->Update();
		
		IndexType pixelIndex;
		fixedImage->TransformPhysicalPointToIndex(point, pixelIndex);
		
#if 0
		// Only iterate through the region that is within 3 sigma of the mean
		
		IndexType regionStart;
		for (int k=0; k < VDimension; ++k)
		{
			if ((pixelIndex[k] - 3*sigma[k]) > 0)
			regionStart[k] = (int) floor(pixelIndex[k] - 3*sigma[k]);
			else 
			regionStart[k] = 0;
		}
		
		typename TFixedImageType::SizeType regionSize = 
		fixedImage->GetLargestPossibleRegion().GetSize();
		
		// Avoid far edge
		for (int k=0; k < VDimension; ++k) {
			if ( ceil(regionStart[k] + 6*sigma[k]) < regionSize[k])
			regionSize[k] =  (int) ceil(6*sigma[k]);
			else
			regionSize[k] -=  regionStart[k];
		}
		
		typename TFixedImageType::RegionType itregion;
		itregion.SetIndex(regionStart);
		itregion.SetSize(regionSize);
#endif
		
		// return the Gaussian weighted Histogram
#ifdef SIFT_FEATURE
		return this->GetSiftKey(hgradImage, gaussImgSource->GetOutput(),
		pixelIndex);
#else
		return this->GetHistogram(hgradImage, gaussImgSource->GetOutput());
#endif
	}



	//Frankie: Edit to allow storing of current scale into data array
	// Input:  Image, Gradient Image, Point
	// Output:  Vector of direction
	template <class TFixedImageType, unsigned int VDimension> 
	typename ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>::FeatureType
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::GetFeaturesWithExtraData( FixedImagePointer fixedImage, 
	typename GradientImageType::Pointer hgradImage,
	PointType &point, float currScale, int featureSizeExtra = 0) 
	{
#ifdef ERROR_CHECK
		std::cerr << "GetFeatures ... ";
#endif
		
		// Generate the Gaussian
		typedef GaussianImageSource<TFixedImageType> GaussianImageSourceType;
		typename GaussianImageSourceType::Pointer gaussImgSource;
		typename GaussianImageSourceType::ArrayType sigma;

		gaussImgSource = GaussianImageSourceType::New();
		gaussImgSource->SetNormalized(true);
		//gaussImgSource->SetNormalized(false);
		//gaussImgSource->SetScale(255);  
		gaussImgSource->SetSpacing(fixedImage->GetSpacing());
		gaussImgSource->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
		gaussImgSource->SetMean(point);
		
#if 0
		// If we wanted to find the orientation,  we would use this
		for (int i = 0; i < VDimension; ++i)
		sigma[i] = 3.0 * currScale;
		gaussImgSource->SetSigma(sigma);
#endif
		
		// Simulate the 16x16 Gaussian window descriptor
		// use sigma equal to half the descriptor window width
		for (int i = 0; i < VDimension; ++i)
		sigma[i] = 8.0;
		gaussImgSource->SetSigma(sigma);
		
		gaussImgSource->Update();
		
		IndexType pixelIndex;
		fixedImage->TransformPhysicalPointToIndex(point, pixelIndex);
		
#if 0
		// Only iterate through the region that is within 3 sigma of the mean
		
		IndexType regionStart;
		for (int k=0; k < VDimension; ++k)
		{
			if ((pixelIndex[k] - 3*sigma[k]) > 0)
			regionStart[k] = (int) floor(pixelIndex[k] - 3*sigma[k]);
			else 
			regionStart[k] = 0;
		}
		
		typename TFixedImageType::SizeType regionSize = 
		fixedImage->GetLargestPossibleRegion().GetSize();
		
		// Avoid far edge
		for (int k=0; k < VDimension; ++k) {
			if ( ceil(regionStart[k] + 6*sigma[k]) < regionSize[k])
			regionSize[k] =  (int) ceil(6*sigma[k]);
			else
			regionSize[k] -=  regionStart[k];
		}
		
		typename TFixedImageType::RegionType itregion;
		itregion.SetIndex(regionStart);
		itregion.SetSize(regionSize);
#endif
		
		// return the Gaussian weighted Histogram

		

#ifdef SIFT_FEATURE
		FeatureType result(this->SiftFeatureSize());
		FeatureType resultWithScale(this->SiftFeatureSize() + featureSizeExtra);
		result = this->GetSiftKey(hgradImage, gaussImgSource->GetOutput(), pixelIndex);

		//Run through array and put copy content from result to resultWithScale then append scale level at the end
		for (int i = 0; i < this->SiftFeatureSize(); i++) {
			resultWithScale[i] = result[i];
		}

		return resultWithScale;
#else
		FeatureType result(this->HistFeatureSize());
		FeatureType resultWithScale(this->HistFeatureSize() + 1);
		result = this->GetHistogram(hgradImage, gaussImgSource->GetOutput());

		//Run through array and put copy content from result to resultWithScale then append scale level at the end
		for (int i = 0; i < this->HistFeatureSize(); i++) {
			resultWithScale[i] = result[i];
		}

		resultWithScale[this->HistFeatureSize()] = currScale;

		return resultWithScale;
#endif
	}






#endif

	template <class TFixedImageType, unsigned int VDimension> 
	void
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::writeImage(FixedImagePointer fixedImage, const char *filename)
	{
		
		//Frankie: Hacking in a print functionaility for 3D images
#ifdef IMAGEIMPORT
		
			
		//Read in outConfigFile and set up export filter
		std::list<std::string> * outConfigFileParam = new std::list<std::string>;
		std::list<std::string> * outConfigFileContent = new std::list<std::string>;

		std::string outputImageFile(filename);
		std::string outConfigFile = "oconfig_" + outputImageFile+ ".txt";
		
		outputImageFile = outputImageFile + ".raw";

		//Write output image and config file
		itkImageIORaw<PixelType, VDimension>* inputImageImport = new itkImageIORaw<PixelType, VDimension>;

		inputImageImport->setConfigFileName(outConfigFile);

		ImageIOData outputImageIOData;

		outputImageIOData.setPixelType("float");
		outputImageIOData.setPixelDimensionality(1);
		outputImageIOData.setImageDimensionality(VDimension);

		const TFixedImageType::RegionType& OutputRegion = fixedImage->GetLargestPossibleRegion();
		const TFixedImageType::SizeType& OutputSize = OutputRegion.GetSize();
		unsigned int tmpSize[VDimension];
		for (int i = 0; i < VDimension; i++) {
			tmpSize[i] = OutputSize[i];
		}

		outputImageIOData.setDimensionSize((int *)tmpSize);

		const TFixedImageType::SpacingType OutputSpacing = fixedImage->GetSpacing();
		float tmpSpacing[VDimension];
		for (int i = 0; i < VDimension; i++) {
			tmpSpacing[i] = (float)OutputSpacing[i];
		}
		
		outputImageIOData.setSpacingSize((float *)tmpSpacing);

		const TFixedImageType::PointType& OutputOrigin  = fixedImage->GetOrigin();
		float tmpOrigin[VDimension];
		for (int i=0; i< VDimension; i++) {
			tmpOrigin[i] = (float)OutputOrigin[i];
		}

		outputImageIOData.setOrigin((float *)tmpOrigin);

		outputImageIOData.setHeaderSize(0);
		outputImageIOData.setByteOrder("little");

		inputImageImport->setImageIOData(outputImageIOData);

		inputImageImport->writeConfigFile();

		//Write Image
		inputImageImport->setitkImage((void *) &fixedImage);
		inputImageImport->setImageFileName(outputImageFile);
		inputImageImport->writeRawImageKnown();

		delete inputImageImport;

		
		return;
		
#endif

		typedef itk::Image< unsigned char, VDimension >  OutImageType;
		typedef typename itk::ImageFileWriter< OutImageType  >  FixedWriterType;
		typedef  itk::ResampleImageFilter< TFixedImageType,  OutImageType    >    
		OutResampleFilterType;
		
		
		typename OutResampleFilterType::Pointer resampler = 
		OutResampleFilterType::New();
		
		resampler->SetSize( fixedImage->GetLargestPossibleRegion().GetSize() );
		resampler->SetOutputSpacing( fixedImage->GetSpacing() );
		resampler->SetTransform(m_IdentityTransform);
		resampler->SetInput(fixedImage);
		
		typename FixedWriterType::Pointer fixedWriter = FixedWriterType::New();
		
		fixedWriter->SetFileName(filename);
		fixedWriter->SetInput( resampler->GetOutput() );
		
		std::cout << "[Writing file << " << filename << "]";
		
		try 
		{
			fixedWriter->Update();
		}
		catch( itk::ExceptionObject & excep )
		{
			std::cerr << "Exception caught !" << std::endl;
			std::cerr << excep << std::endl;
		}

		return;

	}

	// create a filter that resamples the image (scale up or down) 
	template <class TFixedImageType, unsigned int VDimension> 
	typename ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>::ResampleFilterType::Pointer
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::getScaleResampleFilter ( typename TFixedImageType::Pointer fixedImage, float scale )
	{
		typename ResampleFilterType::Pointer scaler = ResampleFilterType::New();
		typedef	itk::BSplineInterpolateImageFunction< TFixedImageType, double > upSampleInterpolatorType;
		
		upSampleInterpolatorType::Pointer upSampleInterpolatorFixedImage = upSampleInterpolatorType::New();

		scaler->SetInput( fixedImage );
		scaler->SetInterpolator(upSampleInterpolatorFixedImage);
		
		
		// Change the size of the image
		typename TFixedImageType::SizeType size = 
		fixedImage->GetLargestPossibleRegion().GetSize();
		for (int k = 0; k < VDimension; ++k)
		size[k] = (unsigned int) floor(size[k] * scale);
		scaler->SetSize( size );

		// Change the spacing
		typename TFixedImageType::SpacingType spacing = 
		fixedImage->GetSpacing();
		for (int k = 0; k < VDimension; ++k)
		spacing[k] = (spacing[k] / scale);
		scaler->SetOutputSpacing( spacing );

		scaler->SetTransform( m_IdentityTransform );
		scaler->SetDefaultPixelValue( 100 );
		scaler->Update();
		
		return scaler;
	}

	template <class TFixedImageType, unsigned int VDimension> 
	void
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::SetDoubling (bool tmp) 
	{
		m_DoubleOriginalImage = tmp;
	}

	template <class TFixedImageType, unsigned int VDimension>
	void
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::SetNumBins( unsigned int tmp) 
	{
		m_HistogramBinsNumber = tmp;
	}

	template <class TFixedImageType, unsigned int VDimension> 
	void
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::SetSigma( double tmp) 
	{
		m_GaussianSigma = tmp;
	}      

	template <class TFixedImageType, unsigned int VDimension> 
	void
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::SetNumScales ( unsigned int tmp) 
	{
		m_ImageScalesTestedNumber = tmp;
	}
	
	template <class TFixedImageType, unsigned int VDimension> 
	void
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::SetMatchRatio ( float tmp) 
	{
		m_MaxFeatureDistanceRatio = tmp;
	}


	template <class TFixedImageType, unsigned int VDimension> 
	typename ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>::PointSetTypePointer
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::getSiftFeatures(FixedImagePointer fixedImage, float edge_point_ratio, float edge_point_ratio_scale) 
	{
		
		unsigned int numMin = 0, numMax = 0, numReject = 0;
		m_KeypointSet = PointSetType::New();

		m_PointsCount = 0;      

		// Declare Gaussian 
		typedef itk::DiscreteGaussianImageFilter<TFixedImageType, TFixedImageType > 
		GaussianFilterType;

		//typename GaussianFilterType::Pointer
		//  gaussianFilter[m_GaussianImagesNumber]; 
		GaussianFilterType::Pointer * gaussianFilter = new GaussianFilterType::Pointer[m_GaussianImagesNumber];


		//typename TFixedImageType::Pointer
		//  gaussianImage[m_GaussianImagesNumber];
		TFixedImageType::Pointer * gaussianImage = new TFixedImageType::Pointer[m_GaussianImagesNumber];

		// Declare DoG 
		typedef itk::SubtractImageFilter<TFixedImageType, TFixedImageType, 
		TFixedImageType> DifferenceFilterType;
		//typename DifferenceFilterType::Pointer dogFilter[m_DifferenceOfGaussianImagesNumber];
		//typename TFixedImageType::Pointer dogImage[m_DifferenceOfGaussianImagesNumber];
		DifferenceFilterType::Pointer * dogFilter = new DifferenceFilterType::Pointer[m_DifferenceOfGaussianImagesNumber];
		TFixedImageType::Pointer * dogImage = new TFixedImageType::Pointer[m_DifferenceOfGaussianImagesNumber];

		// Resampled image filters
		//typename ResampleFilterType::Pointer scaler[m_ImageScalesTestedNumber];
		//typename TFixedImageType::Pointer scaleImage[m_ImageScalesTestedNumber];
		ResampleFilterType::Pointer * scaler = new ResampleFilterType::Pointer[m_ImageScalesTestedNumber];
		TFixedImageType::Pointer * scaleImage = new TFixedImageType::Pointer[m_ImageScalesTestedNumber];

#ifdef GENERATE_KEYS
		// Declare Gradient
		//typename GradientFilterType::Pointer gradFilter[m_ImageScalesTestedNumber];
		//typename GradientImageType::Pointer gradImage[m_ImageScalesTestedNumber];
		//typename GradientImageType::Pointer hgradImage[m_ImageScalesTestedNumber];
		GradientFilterType::Pointer * gradFilter = new GradientFilterType::Pointer[m_ImageScalesTestedNumber];
		GradientImageType::Pointer * gradImage = new GradientImageType::Pointer[m_ImageScalesTestedNumber];
		GradientImageType::Pointer * hgradImage = new GradientImageType::Pointer[m_ImageScalesTestedNumber];

		//typename GradientMagFilterType::Pointer gradMagFilter[m_ImageScalesTestedNumber];
		GradientMagFilterType::Pointer * gradMagFilter = new GradientMagFilterType::Pointer[m_ImageScalesTestedNumber];
		
		//FixedImagePointer gradMagImage[m_ImageScalesTestedNumber];
		FixedImagePointer * gradMagImage = new FixedImagePointer[m_ImageScalesTestedNumber];
#endif


		float currScale = 0.5;

		// For each scale
		for (int i = 0; i < m_ImageScalesTestedNumber; ++i) {
			std::cout << "Computing Scale Level " << i << "... (";

			if (i == 0 && !m_DoubleOriginalImage) {
				scaleImage[0] = fixedImage;
			} else {
				if (i == 0) {
					// Input is the fixed Image.  
					scaler[i] = getScaleResampleFilter ( fixedImage, m_ScalingFactor );
				} else {
					// Input is the 2*sigma smoothed image from the previous octave
					scaler[i] = getScaleResampleFilter ( gaussianImage[m_DifferenceOfGaussianTestsNumber] , 1.0 / m_ScalingFactor );
				}
				scaleImage[i] = scaler[i]->GetOutput();
			}

			{
				typename TFixedImageType::SizeType gsize = 
				scaleImage[i]->GetLargestPossibleRegion().GetSize();
				for (int j = 0; j < VDimension; ++j)
				std::cout << gsize[j] << " ";
			}
			
			std::cout << ") Done\n";

#ifdef DEBUG
			char filename[256];
			sprintf(filename, "gauss-%d-0.png", i);
			this->writeImage(scaleImage[i], filename);
#endif

#ifdef GENERATE_KEYS
			
			std::cout << "...Computing Gradient for Binary Image...";
			gradFilter[i] = GradientFilterType::New();
			gradFilter[i]->SetInput(scaleImage[i]);
			// Do this in pixel space
			gradFilter[i]->SetUseImageSpacing(false);
			gradFilter[i]->Update();
			gradImage[i] = gradFilter[i]->GetOutput();
			hgradImage[i] = this->GetHypersphericalCoordinates(gradImage[i]);

			gradMagFilter[i] = GradientMagFilterType::New();
			gradMagFilter[i]->SetInput(scaleImage[i]);
			// Do this in pixel space
			gradMagFilter[i]->SetUseImageSpacing(false);
			gradMagFilter[i]->Update();
			gradMagImage[i] = gradMagFilter[i]->GetOutput();

			std::cout << "...Done\n";
#endif

			// ...Compute Gaussians
			for (int j = 0; j < m_GaussianImagesNumber; ++j) {
#ifdef VERBOSE
				std::cout << "Setting Up Gaussian Filter " << i << "-" << j << "...";
				std::cout.flush();
#endif
				/* Variance is square of the sigma
	* sigma = (2^(j/s)*sigma)
	*/

				double variance = this->GetGaussianScale(j);
				variance *= variance;

				gaussianFilter[j] = GaussianFilterType::New();

				gaussianFilter[j]->SetVariance(variance);
				gaussianFilter[j]->SetInput( scaleImage[i] );
				// pixel-wise smoothing
				gaussianFilter[j]->SetUseImageSpacing(false); 
				try {
					gaussianFilter[j]->Update();
				}
				catch( itk::ExceptionObject & excep ) {
					std::cerr << "Exception caught !" << std::endl;
					std::cerr << excep << std::endl;
				}

				gaussianImage[j] = gaussianFilter[j]->GetOutput();

#ifdef DEBUG
				char filename[256];
				sprintf(filename, "gauss-%d-%d.png", i, j);
				this->writeImage(gaussianImage[j], filename);
#endif

#ifdef VERBOSE
				std::cout << "Done\n";
				std::cout.flush();
#endif
			}
			
			// ...Compute Difference of Gaussians
			for (int j = 0; j < (m_DifferenceOfGaussianImagesNumber); ++j) {
#ifdef VERBOSE
				std::cout << "Setting Up DoG Binary Filter " << i << "-" << j << "...";
				std::cout.flush();
#endif
				dogFilter[j] = DifferenceFilterType::New();
				dogFilter[j]->SetInput1( gaussianImage[j] );
				dogFilter[j]->SetInput2( gaussianImage[j+1] );
				dogFilter[j]->Update();
				dogImage[j] = dogFilter[j]->GetOutput();

#ifdef DEBUG
				char filename[256];
				sprintf(filename, "dog-%d-%d.png", i, j);
				this->writeImage(dogImage[j], filename);
#endif

#ifdef VERBOSE
				std::cout << "Done\n";
				std::cout.flush();
#endif 
			}

			for (int j=1; j < (m_DifferenceOfGaussianImagesNumber - 1); ++j) {
				// Search the dogImages for local maxima,  w.r.t. corresponding
				// point in the scale above and below
				// level 0 is the "doubled" image
				// Iterate over the various doG filters
				// Only use the middle dogs (ones with both neighbours above and below)
				// Iterate over each position in the dog filter
				typedef itk::ImageRegionIteratorWithIndex< TFixedImageType > 
				ImageIteratorType;

				IndexType regionStart;
				// Avoid the edges
				for (int k=0; k < VDimension; ++k)
				regionStart[k] = 1;

				typename TFixedImageType::SizeType regionSize = 
				dogImage[j]->GetLargestPossibleRegion().GetSize();

#ifdef VERBOSE
				std::cout << "Searching for Extrema in DoG Binary Image " << i << "-" << j;
				std::cout << " ( ";
				for (int k=0; k < VDimension; ++k)
				std::cout << regionSize[k] << " ";
				std::cout << ") Scale " << currScale << "\n";
				std::cout.flush();
#endif

				// Avoid far edge
				for (int k=0; k < VDimension; ++k)
				regionSize[k] -=  2;
				
				typename TFixedImageType::RegionType itregion;
				itregion.SetIndex(regionStart);
				itregion.SetSize(regionSize);
				
				ImageIteratorType pixelIt(dogImage[j],
				itregion);
				
				for ( pixelIt.GoToBegin(); !pixelIt.IsAtEnd(); ++pixelIt) {
					// Make sure to start sufficiently into the image so that all
					// neighbours are present
					IndexType pixelIndex = pixelIt.GetIndex();
					typename TFixedImageType::PixelType pixelValue = pixelIt.Get();

					PointType point;
					dogImage[j]->TransformIndexToPhysicalPoint (pixelIndex, point);

#ifdef ERROR_CHECK
					std::cerr << "Checking ( ";
					for (int k = 0; k < VDimension; ++k)
					std::cerr << pixelIndex[k] << " ";
					std::cerr << ") = " << pixelValue <<"\n";
#endif

					// Compare to the 8 immediate neighbours
					bool isMax=true;
					bool isMin=true;

					this->CheckLocalExtrema(dogImage[j], 
					pixelIndex, pixelValue, isMax, isMin, false);

					if (!isMax && !isMin) continue;
					
					// Compare to scale above
					if (j < (m_GaussianImagesNumber-1)) {
#ifdef DEBUG_VERBOSE
						std::cout << "...Checking Scale Above\n";
#endif
						//dogImage[i+1][j]->TransformPhysicalPointToIndex (point, tmpIndex);

						this->CheckLocalExtrema(dogImage[j+1], 
						pixelIndex, pixelValue, isMax, isMin, true);
					}
					if (!isMax && !isMin) continue;

					// Compare to scale below
					if (j > 0) {
#ifdef DEBUG_VERBOSE
						std::cout << "...Checking Scale Below\n";
#endif
						//dogImage[i-1][j]->TransformPhysicalPointToIndex (point, tmpIndex);

						this->CheckLocalExtrema(dogImage[j-1], 
						pixelIndex, pixelValue, isMax, isMin, true);
					}
					if (!isMax && !isMin) continue;
					

#ifndef BAD_POINT_REJECTION

					// Check if it is sufficiently large (absolute value)
					if (fabs(pixelValue) < m_MinKeypointValue) {
						++numReject;
						continue;
					}

#else	
					//Frankie: Will want to add in the extra bit with the Hessian checks here...
					//Frankie: Maybe the std and abs pixel value from original image check here as well...
					if (DoGAndEdgePointRatioThreshold(dogImage[j], pixelIndex, pixelValue, edge_point_ratio * pow(edge_point_ratio_scale,i))) {
						++numReject;
						continue;
					}


#endif						  
					
					// Passed all checks:

#ifdef DEBUG
					std::cout << point << std::endl;
#endif
					m_KeypointSet->SetPoint( m_PointsCount, point);
#ifdef GENERATE_KEYS
					// Generate features
					// Space used is the (smoothed) original image)
					m_KeypointSet->SetPointData( m_PointsCount, 
					this->GetFeatures( gaussianImage[0], 
					hgradImage[i], point,
					this->GetGaussianScale(j)));
#else
					m_KeypointSet->SetPointData( m_PointsCount, currScale);
#endif
					++m_PointsCount;

					if (isMax) {
						// Maxima detected.  
						++numMax;
#ifdef DEBUG
						std::cout << "Found Maxima! ";
#endif
					}
					if (isMin) {
						// Minima detected.  
						++numMin;
#ifdef DEBUG
						std::cout << "Found Minima! ";
#endif
					}
				}
#ifdef VERBOSE
				std::cout << "Acc. Num Max: " << numMax 
				<< "\nAcc. Num Min: " << numMin 
				<< "\nAcc. Num Reject: " << numReject 
				<< std::endl;
				std::cout.flush();
#endif
			}
			currScale *= m_ScalingFactor;
		}
		
#ifdef VERBOSE
		std::cout << "Total Num Max: " << numMax 
		<< "\nTotal Num Min: " << numMin 
		<< "\nTotal Num Reject: " << numReject
		<< std::endl;	
		std::cout.flush();
#endif
		return m_KeypointSet;

	}
	
	//Frankie: Adaptation of getSiftFeatures that adjusts the scale levels to match the rescaling processes
	template <class TFixedImageType, unsigned int VDimension> 
	typename ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>::PointSetTypePointer
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::getSiftFeaturesRescaling(FixedImagePointer fixedImage, int scaleStart, int scaleEnd, float edge_point_ratio, float edge_point_ratio_scale)
	{
		unsigned int numMin = 0, numMax = 0, numReject = 0;
		m_KeypointSet = PointSetType::New();

		m_PointsCount = 0;   

		//We want to change what scale levels are analysed based off of how much the resample factor adjust the scale of the original image already
		
		float ScalingFactorAdjust = pow(m_ScalingFactor, scaleStart);
		int ImageScalesTestedNumber = scaleEnd - scaleStart;
		
		// Declare Gaussian 
		typedef itk::DiscreteGaussianImageFilter<TFixedImageType, TFixedImageType > 
		GaussianFilterType;

		//typename GaussianFilterType::Pointer
		//  gaussianFilter[m_GaussianImagesNumber]; 
		GaussianFilterType::Pointer * gaussianFilter = new GaussianFilterType::Pointer[m_GaussianImagesNumber];

		//typename TFixedImageType::Pointer
		//  gaussianImage[m_GaussianImagesNumber];
		TFixedImageType::Pointer * gaussianImage = new TFixedImageType::Pointer[m_GaussianImagesNumber];

		// Declare DoG 
		typedef itk::SubtractImageFilter<TFixedImageType, TFixedImageType, 
		TFixedImageType> DifferenceFilterType;
		//typename DifferenceFilterType::Pointer dogFilter[m_DifferenceOfGaussianImagesNumber];
		//typename TFixedImageType::Pointer dogImage[m_DifferenceOfGaussianImagesNumber];
		DifferenceFilterType::Pointer * dogFilter = new DifferenceFilterType::Pointer[m_DifferenceOfGaussianImagesNumber];
		TFixedImageType::Pointer * dogImage = new TFixedImageType::Pointer[m_DifferenceOfGaussianImagesNumber];

		// Resampled image filters
		//typename ResampleFilterType::Pointer scaler[ImageScalesTestedNumber];
		//typename TFixedImageType::Pointer scaleImage[ImageScalesTestedNumber];
		ResampleFilterType::Pointer * scaler = new ResampleFilterType::Pointer[ImageScalesTestedNumber];
		TFixedImageType::Pointer * scaleImage = new TFixedImageType::Pointer[ImageScalesTestedNumber];

#ifdef GENERATE_KEYS
		// Declare Gradient
		//typename GradientFilterType::Pointer gradFilter[ImageScalesTestedNumber];
		//typename GradientImageType::Pointer gradImage[ImageScalesTestedNumber];
		//typename GradientImageType::Pointer hgradImage[ImageScalesTestedNumber];
		GradientFilterType::Pointer * gradFilter = new GradientFilterType::Pointer[ImageScalesTestedNumber];
		GradientImageType::Pointer * gradImage = new GradientImageType::Pointer[ImageScalesTestedNumber];
		GradientImageType::Pointer * hgradImage = new GradientImageType::Pointer[ImageScalesTestedNumber];

		//typename GradientMagFilterType::Pointer gradMagFilter[ImageScalesTestedNumber];
		GradientMagFilterType::Pointer * gradMagFilter = new GradientMagFilterType::Pointer[ImageScalesTestedNumber];
		
		//FixedImagePointer gradMagImage[ImageScalesTestedNumber];
		FixedImagePointer * gradMagImage = new FixedImagePointer[ImageScalesTestedNumber];
#endif

		float currScale = 0.5 * ScalingFactorAdjust;
		
		// For each scale
		for (int i = 0; i < ImageScalesTestedNumber; ++i) {
			std::cout << "Computing Scale Level " << i << "... (";

			if (i == 0 && !m_DoubleOriginalImage) {

				if (ScalingFactorAdjust == 1) {
					scaleImage[0] = fixedImage;
				} else {
					scaler[0] = getScaleResampleFilter( fixedImage, 1.0 / ScalingFactorAdjust);
					scaleImage[0] = scaler[0]->GetOutput();
				}
			} else {
				if (i == 0) {
					// Input is the fixed Image.  
					scaler[i] = getScaleResampleFilter ( fixedImage, m_ScalingFactor / ScalingFactorAdjust);
				} else {
					// Input is the 2*sigma smoothed image from the previous octave
					scaler[i] = getScaleResampleFilter ( gaussianImage[m_DifferenceOfGaussianTestsNumber] , 1.0 / m_ScalingFactor );
				}
				scaleImage[i] = scaler[i]->GetOutput();
			}

			{
				typename TFixedImageType::SizeType gsize = 
				scaleImage[i]->GetLargestPossibleRegion().GetSize();
				for (int j = 0; j < VDimension; ++j)
				std::cout << gsize[j] << " ";
			}
			
			std::cout << ") Done\n";

#ifdef DEBUG
			char filename[256];
			sprintf(filename, "gauss-%d-0.png", i);
			this->writeImage(scaleImage[i], filename);
#endif

#ifdef GENERATE_KEYS
			// ...Compute Gradient
			std::cout << "...Computing Gradient...";
			gradFilter[i] = GradientFilterType::New();
			gradFilter[i]->SetInput(scaleImage[i]);
			// Do this in pixel space
			gradFilter[i]->SetUseImageSpacing(false);
			gradFilter[i]->Update();
			gradImage[i] = gradFilter[i]->GetOutput();
			hgradImage[i] = this->GetHypersphericalCoordinates(gradImage[i]);

			gradMagFilter[i] = GradientMagFilterType::New();
			gradMagFilter[i]->SetInput(scaleImage[i]);
			// Do this in pixel space
			gradMagFilter[i]->SetUseImageSpacing(false);
			gradMagFilter[i]->Update();
			gradMagImage[i] = gradMagFilter[i]->GetOutput();
			std::cout << "...Done\n";
#endif

			// ...Compute Gaussians
			for (int j = 0; j < m_GaussianImagesNumber; ++j) {
#ifdef VERBOSE
				std::cout << "Setting Up Gaussian Filter " << i << "-" << j << "...";
				std::cout.flush();
#endif
				gaussianFilter[j] = GaussianFilterType::New();

				/* Variance is square of the sigma
	* sigma = (2^(j/s)*sigma)
	*/

				double variance = this->GetGaussianScale(j);
				variance *= variance;
				gaussianFilter[j]->SetVariance(variance);
				gaussianFilter[j]->SetInput( scaleImage[i] );
				// pixel-wise smoothing
				gaussianFilter[j]->SetUseImageSpacing(false); 
				try {
					gaussianFilter[j]->Update();
				}
				catch( itk::ExceptionObject & excep ) {
					std::cerr << "Exception caught !" << std::endl;
					std::cerr << excep << std::endl;
				}

				gaussianImage[j] = gaussianFilter[j]->GetOutput();

#ifdef DEBUG
				char filename[256];
				sprintf(filename, "gauss-%d-%d.png", i, j);
				this->writeImage(gaussianImage[j], filename);
#endif

#ifdef VERBOSE
				std::cout << "Done\n";
				std::cout.flush();
#endif
			}
			
			// ...Compute Difference of Gaussians
			for (int j = 0; j < (m_DifferenceOfGaussianImagesNumber); ++j) {
#ifdef VERBOSE
				std::cout << "Setting Up DoG Filter " << i << "-" << j << "...";
				std::cout.flush();
#endif
				dogFilter[j] = DifferenceFilterType::New();
				dogFilter[j]->SetInput1( gaussianImage[j] );
				dogFilter[j]->SetInput2( gaussianImage[j+1] );
				dogFilter[j]->Update();
				dogImage[j] = dogFilter[j]->GetOutput();

#ifdef DEBUG
				char filename[256];
				sprintf(filename, "dog-%d-%d.png", i, j);
				this->writeImage(dogImage[j], filename);
#endif

#ifdef VERBOSE
				std::cout << "Done\n";
				std::cout.flush();
#endif 
			}

			for (int j=1; j < (m_DifferenceOfGaussianImagesNumber - 1); ++j) {
				// Search the dogImages for local maxima,  w.r.t. corresponding
				// point in the scale above and below
				// level 0 is the "doubled" image
				// Iterate over the various doG filters
				// Only use the middle dogs (ones with both neighbours above and below)
				// Iterate over each position in the dog filter
				typedef itk::ImageRegionIteratorWithIndex< TFixedImageType > 
				ImageIteratorType;

				IndexType regionStart;
				// Avoid the edges
				for (int k=0; k < VDimension; ++k)
				regionStart[k] = 1;

				typename TFixedImageType::SizeType regionSize = 
				dogImage[j]->GetLargestPossibleRegion().GetSize();

#ifdef VERBOSE
				std::cout << "Searching for Extrema in DoG Image " << i << "-" << j;
				std::cout << " ( ";
				for (int k=0; k < VDimension; ++k)
				std::cout << regionSize[k] << " ";
				std::cout << ") Scale " << currScale << "\n";
				std::cout.flush();
#endif

				// Avoid far edge
				for (int k=0; k < VDimension; ++k)
				regionSize[k] -=  2;
				
				typename TFixedImageType::RegionType itregion;
				itregion.SetIndex(regionStart);
				itregion.SetSize(regionSize);
				
				ImageIteratorType pixelIt(dogImage[j],
				itregion);
				
				for ( pixelIt.GoToBegin(); !pixelIt.IsAtEnd(); ++pixelIt) {
					// Make sure to start sufficiently into the image so that all
					// neighbours are present
					IndexType pixelIndex = pixelIt.GetIndex();
					typename TFixedImageType::PixelType pixelValue = pixelIt.Get();

					PointType point;
					dogImage[j]->TransformIndexToPhysicalPoint (pixelIndex, point);

#ifdef ERROR_CHECK
					std::cerr << "Checking ( ";
					for (int k = 0; k < VDimension; ++k)
					std::cerr << pixelIndex[k] << " ";
					std::cerr << ") = " << pixelValue <<"\n";
#endif

					// Compare to the 8 immediate neighbours
					bool isMax=true;
					bool isMin=true;

					this->CheckLocalExtrema(dogImage[j], 
					pixelIndex, pixelValue, isMax, isMin, false);

					if (!isMax && !isMin) continue;
					
					// Compare to scale above
					if (j < (m_GaussianImagesNumber-1)) {
#ifdef DEBUG_VERBOSE
						std::cout << "...Checking Scale Above\n";
#endif
						//dogImage[i+1][j]->TransformPhysicalPointToIndex (point, tmpIndex);

						this->CheckLocalExtrema(dogImage[j+1], 
						pixelIndex, pixelValue, isMax, isMin, true);
					}
					if (!isMax && !isMin) continue;

					// Compare to scale below
					if (j > 0) {
#ifdef DEBUG_VERBOSE
						std::cout << "...Checking Scale Below\n";
#endif
						//dogImage[i-1][j]->TransformPhysicalPointToIndex (point, tmpIndex);

						this->CheckLocalExtrema(dogImage[j-1], 
						pixelIndex, pixelValue, isMax, isMin, true);
					}
					if (!isMax && !isMin) continue;
					

#ifndef BAD_POINT_REJECTION

					// Check if it is sufficiently large (absolute value)
					if (fabs(pixelValue) < m_MinKeypointValue) {
						++numReject;
						continue;
					}

#else	
					if (DoGAndEdgePointRatioThreshold(dogImage[j], pixelIndex, pixelValue, edge_point_ratio * pow(edge_point_ratio_scale,i))) {
						++numReject;
						continue;
					}


#endif						  
					
					// Passed all checks:

#ifdef DEBUG
					std::cout << point << std::endl;
#endif
					m_KeypointSet->SetPoint( m_PointsCount, point);
#ifdef GENERATE_KEYS
					// Generate features
					// Space used is the (smoothed) original image)
					m_KeypointSet->SetPointData( m_PointsCount, 
					this->GetFeatures( gaussianImage[0], 
					hgradImage[i], point,
					this->GetGaussianScale(j)));
#else
					m_KeypointSet->SetPointData( m_PointsCount, currScale);
#endif
					++m_PointsCount;

					if (isMax) {
						// Maxima detected.  
						++numMax;
#ifdef DEBUG
						std::cout << "Found Maxima! ";
#endif
					}
					if (isMin) {
						// Minima detected.  
						++numMin;
#ifdef DEBUG
						std::cout << "Found Minima! ";
#endif
					}
				}
#ifdef VERBOSE
				std::cout << "Acc. Num Max: " << numMax 
				<< "\nAcc. Num Min: " << numMin 
				<< "\nAcc. Num Reject: " << numReject 
				<< std::endl;
				std::cout.flush();
#endif
			}
			currScale *= m_ScalingFactor;
		}
		
#ifdef VERBOSE
		std::cout << "Total Num Max: " << numMax 
		<< "\nTotal Num Min: " << numMin 
		<< "\nTotal Num Reject: " << numReject
		<< std::endl;	
		std::cout.flush();
#endif
		return m_KeypointSet;
	}
	
	
	

	template <class TFixedImageType, unsigned int VDimension> 
	void
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::MatchKeypointsPos(PointSetTypePointer keypoints1, PointSetTypePointer keypoints2,
	typename TransformType::Pointer inverse_transform)
	{
		// Compare Keypoints.  Check Coverage and Accuracy
		// This does the comparison based on position of the keypoints
		// Find:  
		// # of points that match which will tell us
		// # of points that did not scale
		// # of points created by the scale
		unsigned int numMatches;
		unsigned int numMatches2;
		unsigned int numMatches5;
		const double MATCH_THRESHOLD = 1.5;
		typename PointSetType::PointsContainer::Pointer keyps1, keyps2;

		unsigned long numpoints1, numpoints2;
		numpoints1 = keypoints1->GetNumberOfPoints();
		std::cout << "Keypoints1 Found: " << numpoints1 << std::endl;
		numpoints2 = keypoints2->GetNumberOfPoints();
		std::cout << "Keypoints2 Found: " << numpoints2 << std::endl;

		if (!inverse_transform)
		return;

		numMatches = 0;
		numMatches2 = 0;
		numMatches5 = 0;
		for (unsigned int i = 0; i < numpoints2; ++i) {
			PointType pp2;
			keypoints2->GetPoint(i, &pp2);	
			
			bool match = false;
			bool match2 = false;
			bool match5 = false;
			for (unsigned int j = 0; j < numpoints1; ++j) {
				PointType tmpp, pp;
				keypoints1->GetPoint(j, &tmpp);
				pp = inverse_transform->TransformPoint(tmpp);
				//	    for (int k = 0; k < VDimension; ++k)
				//	      pp[k] *= scale_test;	      
				if(!match  && pp.EuclideanDistanceTo(pp2) <= MATCH_THRESHOLD)
				{
					++numMatches;
					match = true;
				}
				if(!match2 && pp.EuclideanDistanceTo(pp2) <= 2*MATCH_THRESHOLD)
				{
					++numMatches2;
					match2 = true;
				}
				if(!match5 && pp.EuclideanDistanceTo(pp2) <= 5*MATCH_THRESHOLD)
				{
					++numMatches5;
					match5 = true;
				}
				if (match && match2 && match5)
				break;
			}      
		}

		std::cout << "Keypoints 2 Matching to Keypoints 1 (<" << MATCH_THRESHOLD << "): " << numMatches << std::endl;
		std::cout << "% of Keypoints 2 Matching (<" << MATCH_THRESHOLD << "):  " << (float) numMatches / numpoints2 << std::endl;
		std::cout << "Keypoints 2 Matching to Keypoints 1 (<" << 2*MATCH_THRESHOLD << "): " << numMatches2 << std::endl;
		std::cout << "% of Keypoints 2 Matching (<" << 2*MATCH_THRESHOLD << "):  " << (float) numMatches2 / numpoints2 << std::endl;
		std::cout << "Keypoints 2 Matching to Keypoints 1 (<" << 5*MATCH_THRESHOLD << "): " << numMatches5 << std::endl;
		std::cout << "% of Keypoints 2 Matching (<" << 5*MATCH_THRESHOLD << "):  " << (float) numMatches5 / numpoints2 << std::endl;


		numMatches = 0;
		numMatches2 = 0;
		numMatches5 = 0;
		for (unsigned int j = 0; j < numpoints1; ++j) {
			PointType tmpp, pp;
			keypoints1->GetPoint(j, &tmpp);
			pp = inverse_transform->TransformPoint(tmpp);
			//	  for (int k = 0; k < VDimension; ++k)
			//	    pp[k] *= scale_test;

			bool match = false;
			bool match2 = false;
			bool match5 = false;	
			for (unsigned int i = 0; i < numpoints2; ++i) {
				PointType pp2;
				keypoints2->GetPoint(i, &pp2);	

				if(!match  && pp.EuclideanDistanceTo(pp2) <= MATCH_THRESHOLD)
				{
					++numMatches;
					match = true;
				}
				if(!match2 && pp.EuclideanDistanceTo(pp2) <= 2*MATCH_THRESHOLD)
				{
					++numMatches2;
					match2 = true;
				}
				if(!match5 && pp.EuclideanDistanceTo(pp2) <= 5*MATCH_THRESHOLD)
				{
					++numMatches5;
					match5 = true;
				}
				if (match && match2 && match5)
				break;
			}      
		}

		std::cout << "Keypoints 1 Matching to Keypoints 2 (<" << MATCH_THRESHOLD << "): " << numMatches << std::endl;
		std::cout << "% of Keypoints 1 Matching (<" << MATCH_THRESHOLD << "):  " << (float) numMatches / numpoints1 << std::endl;
		std::cout << "Keypoints 1 Matching to Keypoints 2 (<" << 2*MATCH_THRESHOLD << "): " << numMatches2 << std::endl;
		std::cout << "% of Keypoints 1 Matching (<" << 2*MATCH_THRESHOLD << "):  " << (float) numMatches2 / numpoints1 << std::endl;
		std::cout << "Keypoints 1 Matching to Keypoints 2 (<" << 5*MATCH_THRESHOLD << "): " << numMatches5 << std::endl;
		std::cout << "% of Keypoints 1 Matching (<" << 5*MATCH_THRESHOLD << "):  " << (float) numMatches5 / numpoints1 << std::endl;
	}

#ifdef GENERATE_KEYS
	template <class TFixedImageType, unsigned int VDimension> 
	void
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::MatchKeypointsFeatures(PointSetTypePointer keypoints1, PointSetTypePointer keypoints2,
	typename TransformType::Pointer inverse_transform)
	{
		// Compare Keypoints.  Check Coverage and Accuracy
		// This does the comparison based on position of the keypoints
		// Find:  
		// # of points that match which will tell us
		// # of points that did not scale
		// # of points created by the scale
		unsigned int numMatches;
		unsigned int numMatchesTried;
		unsigned int numMatches2;
		unsigned int numMatches5;
		const double MATCH_THRESHOLD = 1.5;
		typename PointSetType::PointsContainer::Pointer keyps1, keyps2;

		unsigned long numpoints1, numpoints2;
		numpoints1 = keypoints1->GetNumberOfPoints();
		std::cout << "Keypoints1 Found: " << numpoints1 << std::endl;
		numpoints2 = keypoints2->GetNumberOfPoints();
		std::cout << "Keypoints2 Found: " << numpoints2 << std::endl;

		std::cout << "***Keypoint Matches***\n";

		numMatches = 0;
		numMatches2 = 0;
		numMatches5 = 0;
		numMatchesTried = 0;
		for (unsigned int i = 0; i < numpoints2; ++i) {
			PointType pp2;
			keypoints2->GetPoint(i, &pp2);	
			FeatureType ft2;
			keypoints2->GetPointData(i, &ft2);	
			
			FeatureType bestft;
			float bestdist = -1.0;
			float nextbestdist = -1.0;
			unsigned int bestj;
			for (unsigned int j = 0; j < numpoints1; ++j) {
				PointType pp;
				keypoints1->GetPoint(j, &pp);
				FeatureType ft;
				keypoints1->GetPointData(j, &ft);	

				float dist = 0.0;
				for (int k = 0; k < ft.Size(); ++k) {
					dist += (ft[k] - ft2[k])*(ft[k] - ft2[k]);
				} 
				
				if (nextbestdist < 0.0 || dist < bestdist)
				{
					nextbestdist = bestdist;
					bestdist=dist;
					bestft = ft;
					bestj = j;
				}	  
			}

			/* Reject "too close" matches */
			if ((bestdist / nextbestdist) >  m_MaxFeatureDistanceRatio) {
				continue;
			}

			/* NEW IDEA -- look to make sure it is a reciprocal best match */
			/* Take the best feature found,  see if pp2 makes the cut */
			bestdist = -1.0;
			nextbestdist = -1.0;
			FeatureType bestft2;
			unsigned int bestj2;

			for (unsigned int j = 0; j < numpoints2; ++j) {
				PointType pp;
				keypoints2->GetPoint(j, &pp);
				FeatureType ft;
				keypoints2->GetPointData(j, &ft);	

				float dist = 0.0;
				for (int k = 0; k < ft.Size(); ++k)
				dist += (ft[k] - bestft[k])*(ft[k] - bestft[k]);
				
				if (nextbestdist < 0.0 || dist < bestdist)
				{
					nextbestdist = bestdist;
					bestdist=dist;
					bestft2 = ft;
					bestj2 = j;
				}	  
			}

			/* Reject if not reciprocal best hit or "too close" matches */
			if ( bestft2 != ft2 || ((bestdist / nextbestdist) >  m_MaxFeatureDistanceRatio))
			continue;	
			/* END NEW IDEA */

			++numMatchesTried;

			// Check goodness of best feature
			PointType tmpp, pp;
			keypoints1->GetPoint(bestj, &tmpp);

			// Print the match
			std::cout << tmpp << " => " << pp2 << std::endl;

			if (!inverse_transform) {
				continue;
			}
			pp = inverse_transform->TransformPoint(tmpp);

			if(pp.EuclideanDistanceTo(pp2) <= MATCH_THRESHOLD)
			++numMatches;
			if(pp.EuclideanDistanceTo(pp2) <= (2*MATCH_THRESHOLD))
			++numMatches2;
			if(pp.EuclideanDistanceTo(pp2) <= (5*MATCH_THRESHOLD))
			++numMatches5;
		}

		std::cout << "\n***Features 2 Matches Attempted: " << numMatchesTried << std::endl;
		std::cout << "Features 2 Matching to Features 1 (<" << MATCH_THRESHOLD << "): " << numMatches << std::endl;
		std::cout << "% of Features 2 Matching (<" << MATCH_THRESHOLD << "):  " << (float) numMatches / numMatchesTried << std::endl;
		std::cout << "Features 2 Matching to Features 1 (<" << 2*MATCH_THRESHOLD << "): " << numMatches2 << std::endl;
		std::cout << "% of Features 2 Matching (<" << 2*MATCH_THRESHOLD << "):  " << (float) numMatches2 / numMatchesTried << std::endl;
		std::cout << "Features 2 Matching to Features 1 (<" << 5*MATCH_THRESHOLD << "): " << numMatches5 << std::endl;
		std::cout << "% of Features 2 Matching (<" << 5*MATCH_THRESHOLD << "):  " << (float) numMatches5 / numMatchesTried << std::endl;

	}


	//Frankie: Second level hack to match features only when they match on current scale
	template <class TFixedImageType, unsigned int VDimension> 
	void
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::MatchKeypointsFeaturesWithExtraData(PointSetTypePointer * inputkeypoints1, PointSetTypePointer * inputkeypoints2, PointSetTypePointer * matchedKeypoints1,
	PointSetTypePointer * matchedKeypoints2, typename TransformType::Pointer inverse_transform, bool matchSameScale, bool matchSameExtrema, int featureSizeExtra)
	{
		// Compare Keypoints.  Check Coverage and Accuracy
		// This does the comparison based on position of the keypoints
		// Find:  
		// # of points that match which will tell us
		// # of points that did not scale
		// # of points created by the scale
		unsigned int numMatches;
		unsigned int numMatchesTried;
		unsigned int numMatches2;
		unsigned int numMatches5;
		const double MATCH_THRESHOLD = 1.5;
		typename PointSetType::PointsContainer::Pointer keyps1, keyps2;

		unsigned int matchedpointindex = 0;
		PointSetType::Pointer * matchedkeypointsptr1 = (PointSetTypePointer *) matchedKeypoints1;
		PointSetType::Pointer * matchedkeypointsptr2 = (PointSetTypePointer *) matchedKeypoints2;

		PointSetTypePointer keypoints1 = *inputkeypoints1;
		PointSetTypePointer keypoints2 = *inputkeypoints2;
		
		*matchedkeypointsptr1 = PointSetType::New();
		*matchedkeypointsptr2 = PointSetType::New();
		PointSetType::PointsContainerPointer matchedpoints1 = (*matchedkeypointsptr1)->GetPoints();
		PointSetType::PointsContainerPointer matchedpoints2 = (*matchedkeypointsptr2)->GetPoints();

		unsigned long numpoints1, numpoints2;
		numpoints1 = keypoints1->GetNumberOfPoints();
		std::cout << "Keypoints1 Found: " << numpoints1 << std::endl;
		numpoints2 = keypoints2->GetNumberOfPoints();
		std::cout << "Keypoints2 Found: " << numpoints2 << std::endl;

		//std::cout << "***Keypoint Matches***\n";

		numMatches = 0;
		numMatches2 = 0;
		numMatches5 = 0;
		numMatchesTried = 0;
		for (unsigned int i = 0; i < numpoints2; ++i) {
			PointType pp2;
			keypoints2->GetPoint(i, &pp2);	
			FeatureType ft2;
			keypoints2->GetPointData(i, &ft2);	



			//Readjust the feature size such that we seperate the histogram data from the scale and extrema data
			int ftSizeHist = ft2.Size() - featureSizeExtra;

			//Define the location of the scale number and the extrema definition in the feature array
			int ftScale = ft2.Size() - featureSizeExtra;
			int ftExtrema = ft2.Size() - featureSizeExtra + 1;



			
			FeatureType bestft;
			float bestdist = -1.0;
			float nextbestdist = -1.0;
			unsigned int bestj;
			for (unsigned int j = 0; j < numpoints1; ++j) {
				PointType pp;
				keypoints1->GetPoint(j, &pp);
				FeatureType ft;
				keypoints1->GetPointData(j, &ft);	

				float dist = 0.0;
				int size = ft.Size();

				//Check if the scale level (situated at the last index of the feature vector) match
				if (matchSameScale && (ft[ftScale] != ft2[ftScale])) {
					continue;
				}

				//Check if the extrema direction (maximum or minimum) is the same for both points
				if (matchSameExtrema && (ft[ftExtrema] != ft2[ftExtrema])) {
					continue;
				}


				for (int k = 0; k < ftSizeHist; ++k) {
					dist += (ft[k] - ft2[k])*(ft[k] - ft2[k]);
				}
				
				if (nextbestdist < 0.0 || dist < bestdist)
				{
					nextbestdist = bestdist;
					bestdist=dist;
					bestft = ft;
					bestj = j;
				}	  
			}

			/* Reject "too close" matches */
			if ((bestdist / nextbestdist) >  m_MaxFeatureDistanceRatio) {
				continue;
			}

			/* NEW IDEA -- look to make sure it is a reciprocal best match */
			/* Take the best feature found,  see if pp2 makes the cut */
			bestdist = -1.0;
			nextbestdist = -1.0;
			FeatureType bestft2;
			unsigned int bestj2;

			for (unsigned int j = 0; j < numpoints2; ++j) {
				PointType pp;
				keypoints2->GetPoint(j, &pp);
				FeatureType ft;
				keypoints2->GetPointData(j, &ft);	

				float dist = 0.0;

				//Check if the scale level (situated at the last index of the feature vector) match
				if (matchSameScale && (ft[ftScale] != bestft[ftScale])) {
					continue;
				}

				//Check if the extrema direction (maximum or minimum) is the same for both points
				if (matchSameExtrema && (ft[ftExtrema] != bestft[ftExtrema])) {
					continue;
				}

				for (int k = 0; k < ftSizeHist; ++k)
				dist += (ft[k] - bestft[k])*(ft[k] - bestft[k]);
				
				if (nextbestdist < 0.0 || dist < bestdist)
				{
					nextbestdist = bestdist;
					bestdist=dist;
					bestft2 = ft;
					bestj2 = j;
				}	  
			}

			/* Reject if not reciprocal best hit or "too close" matches */
			if ( bestft2 != ft2 || ((bestdist / nextbestdist) >  m_MaxFeatureDistanceRatio)) {
				continue;	
			}



			/* END NEW IDEA */

			++numMatchesTried;

			// Check goodness of best feature
			PointType tmpp, pp;
			keypoints1->GetPoint(bestj, &tmpp);

			// Print the match
			//std::cout << tmpp << " => " << pp2 << std::endl;
			
			//Store matched points into matched points array
			matchedpoints1->InsertElement(matchedpointindex,tmpp);
			matchedpoints2->InsertElement(matchedpointindex,pp2);
			matchedpointindex++;

			if (!inverse_transform) {
				continue;
			}
			
			pp = inverse_transform->TransformPoint(tmpp);

			if(pp.EuclideanDistanceTo(pp2) <= MATCH_THRESHOLD)
			++numMatches;
			if(pp.EuclideanDistanceTo(pp2) <= (2*MATCH_THRESHOLD))
			++numMatches2;
			if(pp.EuclideanDistanceTo(pp2) <= (5*MATCH_THRESHOLD))
			++numMatches5;
		}

		std::cout << "\n***Features 2 Matches Attempted: " << numMatchesTried << std::endl;
		std::cout << "Features 2 Matching to Features 1 (<" << MATCH_THRESHOLD << "): " << numMatches << std::endl;
		std::cout << "% of Features 2 Matching (<" << MATCH_THRESHOLD << "):  " << (float) numMatches / numMatchesTried << std::endl;
		std::cout << "Features 2 Matching to Features 1 (<" << 2*MATCH_THRESHOLD << "): " << numMatches2 << std::endl;
		std::cout << "% of Features 2 Matching (<" << 2*MATCH_THRESHOLD << "):  " << (float) numMatches2 / numMatchesTried << std::endl;
		std::cout << "Features 2 Matching to Features 1 (<" << 5*MATCH_THRESHOLD << "): " << numMatches5 << std::endl;
		std::cout << "% of Features 2 Matching (<" << 5*MATCH_THRESHOLD << "):  " << (float) numMatches5 / numMatchesTried << std::endl;

	}
	

	//Frankie: Hack to allow sending back matched keypoints up a level.
	template <class TFixedImageType, unsigned int VDimension> 
	void
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::MatchKeypointsFeatures(PointSetTypePointer * inputkeypoints1, PointSetTypePointer * inputkeypoints2, PointSetTypePointer * matchedKeypoints1,
	PointSetTypePointer * matchedKeypoints2, typename TransformType::Pointer inverse_transform)
	{
		// Compare Keypoints.  Check Coverage and Accuracy
		// This does the comparison based on position of the keypoints
		// Find:  
		// # of points that match which will tell us
		// # of points that did not scale
		// # of points created by the scale
		unsigned int numMatches;
		unsigned int numMatchesTried;
		unsigned int numMatches2;
		unsigned int numMatches5;
		const double MATCH_THRESHOLD = 1.5;
		typename PointSetType::PointsContainer::Pointer keyps1, keyps2;

		//Frankie: Check if we have vector limit conditions
		bool displacementLimits = false;
		
		unsigned int matchedpointindex = 0;
		PointSetType::Pointer * matchedkeypointsptr1 = (PointSetTypePointer *) matchedKeypoints1;
		PointSetType::Pointer * matchedkeypointsptr2 = (PointSetTypePointer *) matchedKeypoints2;

		PointSetTypePointer keypoints1 = *inputkeypoints1;
		PointSetTypePointer keypoints2 = *inputkeypoints2;
		
		*matchedkeypointsptr1 = PointSetType::New();
		*matchedkeypointsptr2 = PointSetType::New();
		PointSetType::PointsContainerPointer matchedpoints1 = (*matchedkeypointsptr1)->GetPoints();
		PointSetType::PointsContainerPointer matchedpoints2 = (*matchedkeypointsptr2)->GetPoints();

		unsigned long numpoints1, numpoints2;
		numpoints1 = keypoints1->GetNumberOfPoints();
		std::cout << "Keypoints1 Found: " << numpoints1 << std::endl;
		numpoints2 = keypoints2->GetNumberOfPoints();
		std::cout << "Keypoints2 Found: " << numpoints2 << std::endl;

		//std::cout << "***Keypoint Matches***\n";

		numMatches = 0;
		numMatches2 = 0;
		numMatches5 = 0;
		numMatchesTried = 0;
		for (unsigned int i = 0; i < numpoints2; ++i) {
			PointType pp2;
			keypoints2->GetPoint(i, &pp2);	
			FeatureType ft2;
			keypoints2->GetPointData(i, &ft2);	
			
			FeatureType bestft;
			float bestdist = -1.0;
			float nextbestdist = -1.0;
			unsigned int bestj;
			for (unsigned int j = 0; j < numpoints1; ++j) {
				PointType pp;
				keypoints1->GetPoint(j, &pp);
				FeatureType ft;
				keypoints1->GetPointData(j, &ft);	

				float dist = 0.0;
				for (int k = 0; k < ft.Size(); ++k) {
					dist += (ft[k] - ft2[k])*(ft[k] - ft2[k]);
				}
				
				if (nextbestdist < 0.0 || dist < bestdist)
				{

					nextbestdist = bestdist;
					bestdist=dist;
					bestft = ft;
					bestj = j;
				}	  
			}

			/* Reject "too close" matches */
			if ((bestdist / nextbestdist) >  m_MaxFeatureDistanceRatio) {
				continue;
			}

			/* NEW IDEA -- look to make sure it is a reciprocal best match */
			/* Take the best feature found,  see if pp2 makes the cut */
			bestdist = -1.0;
			nextbestdist = -1.0;
			FeatureType bestft2;
			unsigned int bestj2;

			for (unsigned int j = 0; j < numpoints2; ++j) {
				PointType pp;
				keypoints2->GetPoint(j, &pp);
				FeatureType ft;
				keypoints2->GetPointData(j, &ft);	

				float dist = 0.0;
				for (int k = 0; k < ft.Size(); ++k)
				dist += (ft[k] - bestft[k])*(ft[k] - bestft[k]);
				
				if (nextbestdist < 0.0 || dist < bestdist)
				{
					nextbestdist = bestdist;
					bestdist=dist;
					bestft2 = ft;
					bestj2 = j;
				}	  
			}

			/* Reject if not reciprocal best hit or "too close" matches */
			if ( bestft2 != ft2 || ((bestdist / nextbestdist) >  m_MaxFeatureDistanceRatio)) {
				continue;	
			}
			/* END NEW IDEA */

			++numMatchesTried;

			// Check goodness of best feature
			PointType tmpp, pp;
			keypoints1->GetPoint(bestj, &tmpp);

			// Print the match
			//std::cout << tmpp << " => " << pp2 << std::endl;
			
			//Store matched points into matched points array
			matchedpoints1->InsertElement(matchedpointindex,tmpp);
			matchedpoints2->InsertElement(matchedpointindex,pp2);
			matchedpointindex++;

			if (!inverse_transform) {
				continue;
			}
			
			pp = inverse_transform->TransformPoint(tmpp);

			if(pp.EuclideanDistanceTo(pp2) <= MATCH_THRESHOLD)
			++numMatches;
			if(pp.EuclideanDistanceTo(pp2) <= (2*MATCH_THRESHOLD))
			++numMatches2;
			if(pp.EuclideanDistanceTo(pp2) <= (5*MATCH_THRESHOLD))
			++numMatches5;
		}

		std::cout << "\n***Features 2 Matches Attempted: " << numMatchesTried << std::endl;
		std::cout << "Features 2 Matching to Features 1 (<" << MATCH_THRESHOLD << "): " << numMatches << std::endl;
		std::cout << "% of Features 2 Matching (<" << MATCH_THRESHOLD << "):  " << (float) numMatches / numMatchesTried << std::endl;
		std::cout << "Features 2 Matching to Features 1 (<" << 2*MATCH_THRESHOLD << "): " << numMatches2 << std::endl;
		std::cout << "% of Features 2 Matching (<" << 2*MATCH_THRESHOLD << "):  " << (float) numMatches2 / numMatchesTried << std::endl;
		std::cout << "Features 2 Matching to Features 1 (<" << 5*MATCH_THRESHOLD << "): " << numMatches5 << std::endl;
		std::cout << "% of Features 2 Matching (<" << 5*MATCH_THRESHOLD << "):  " << (float) numMatches5 / numMatchesTried << std::endl;

	}

#endif




	//Frankie: Adaptation of getSiftFeaturesRescaling that does matching with feature descriptors found in the original image
	template <class TFixedImageType, unsigned int VDimension> 
	typename ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>::PointSetTypePointer
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::getSiftFeaturesRescalingMatch(FixedImagePointer fixedImage, FixedImagePointer fixedImageBinary, int scaleStart, int scaleEnd, float edge_point_ratio, float edge_point_ratio_scale, bool matchToScaleOrExtrema)
	{
		unsigned int numMin = 0, numMax = 0, numReject = 0;
		m_KeypointSet = PointSetType::New();

		m_PointsCount = 0;   

		//We want to change what scale levels are analysed based off of how much the resample factor adjust the scale of the original image already
		
		float ScalingFactorAdjust = pow(m_ScalingFactor, scaleStart);
		int ImageScalesTestedNumber = scaleEnd - scaleStart;
		
		// Declare Gaussian 
		typedef itk::DiscreteGaussianImageFilter<TFixedImageType, TFixedImageType > 
		GaussianFilterType;

		//typename GaussianFilterType::Pointer
		//  gaussianFilter[m_GaussianImagesNumber]; 
		GaussianFilterType::Pointer gaussianFilterZero = GaussianFilterType::New();
		GaussianFilterType::Pointer gaussianFilterNextOctave = GaussianFilterType::New();
		GaussianFilterType::Pointer * gaussianFilterBinary = new GaussianFilterType::Pointer[m_GaussianImagesNumber];


		//typename TFixedImageType::Pointer
		//  gaussianImage[m_GaussianImagesNumber];
		TFixedImageType::Pointer gaussianImageZero = TFixedImageType::New();
		TFixedImageType::Pointer gaussianImageNextOctave = TFixedImageType::New();
		TFixedImageType::Pointer * gaussianImageBinary = new TFixedImageType::Pointer[m_GaussianImagesNumber];

		// Declare DoG 
		typedef itk::SubtractImageFilter<TFixedImageType, TFixedImageType, 
		TFixedImageType> DifferenceFilterType;
		//typename DifferenceFilterType::Pointer dogFilter[m_DifferenceOfGaussianImagesNumber];
		//typename TFixedImageType::Pointer dogImage[m_DifferenceOfGaussianImagesNumber];
		DifferenceFilterType::Pointer * dogFilterBinary = new DifferenceFilterType::Pointer[m_DifferenceOfGaussianImagesNumber];
		TFixedImageType::Pointer * dogImageBinary = new TFixedImageType::Pointer[m_DifferenceOfGaussianImagesNumber];

		// Resampled image filters
		//typename ResampleFilterType::Pointer scaler[m_ImageScalesTestedNumber];
		//typename TFixedImageType::Pointer scaleImage[m_ImageScalesTestedNumber];
		ResampleFilterType::Pointer * scaler = new ResampleFilterType::Pointer[ImageScalesTestedNumber];
		TFixedImageType::Pointer * scaleImage = new TFixedImageType::Pointer[ImageScalesTestedNumber];

		ResampleFilterType::Pointer * scalerBinary = new ResampleFilterType::Pointer[ImageScalesTestedNumber];
		TFixedImageType::Pointer * scaleImageBinary = new TFixedImageType::Pointer[ImageScalesTestedNumber];

#ifdef GENERATE_KEYS
		// Declare Gradient
		//typename GradientFilterType::Pointer gradFilter[m_ImageScalesTestedNumber];
		//typename GradientImageType::Pointer gradImage[m_ImageScalesTestedNumber];
		//typename GradientImageType::Pointer hgradImage[m_ImageScalesTestedNumber];
		GradientFilterType::Pointer * gradFilter = new GradientFilterType::Pointer[ImageScalesTestedNumber];
		GradientImageType::Pointer * gradImage = new GradientImageType::Pointer[ImageScalesTestedNumber];
		GradientImageType::Pointer * hgradImage = new GradientImageType::Pointer[ImageScalesTestedNumber];

		//typename GradientMagFilterType::Pointer gradMagFilter[m_ImageScalesTestedNumber];
		GradientMagFilterType::Pointer * gradMagFilter = new GradientMagFilterType::Pointer[ImageScalesTestedNumber];
		
		//FixedImagePointer gradMagImage[m_ImageScalesTestedNumber];
		FixedImagePointer * gradMagImage = new FixedImagePointer[ImageScalesTestedNumber];
#endif

		float currScale = 0.5 * ScalingFactorAdjust;
		
		// For each scale
		for (int i = 0; i < ImageScalesTestedNumber; ++i) {
			std::cout << "Computing Scale Level " << i << "... (";

			if (i == 0 && !m_DoubleOriginalImage) {

				if (ScalingFactorAdjust == 1) {
					scaleImage[0] = fixedImage;
					scaleImageBinary[0] = fixedImageBinary;
				} else {
					scaler[0] = getScaleResampleFilter( fixedImage, 1.0 / ScalingFactorAdjust);
					scaleImage[0] = scaler[0]->GetOutput();
					scalerBinary[0] = getScaleResampleFilter( fixedImageBinary, 1.0 / ScalingFactorAdjust);
					scaleImageBinary[0] = scalerBinary[0]->GetOutput();
				}
			} else {
				if (i == 0) {
					// Input is the fixed Image.  
					scaler[i] = getScaleResampleFilter ( fixedImage, m_ScalingFactor / ScalingFactorAdjust);
					scalerBinary[i] = getScaleResampleFilter ( fixedImageBinary, m_ScalingFactor / ScalingFactorAdjust );
				} else {
					// Input is the 2*sigma smoothed image from the previous octave
					scaler[i] = getScaleResampleFilter ( gaussianImageNextOctave , 1.0 / m_ScalingFactor );
					scalerBinary[i] = getScaleResampleFilter ( gaussianImageBinary[m_DifferenceOfGaussianTestsNumber] , 1.0 / m_ScalingFactor );
				}
				scaleImage[i] = scaler[i]->GetOutput();
				scaleImageBinary[i] = scalerBinary[i]->GetOutput();
			}

			{
				typename TFixedImageType::SizeType gsize = 
				scaleImage[i]->GetLargestPossibleRegion().GetSize();
				for (int j = 0; j < VDimension; ++j)
				std::cout << gsize[j] << " ";
			}
			
			std::cout << ") Done\n";


#ifdef DEBUG
			char filename[256];
			sprintf(filename, "gauss-%d-0.png", i);
			this->writeImage(scaleImageBinary[i], filename);
#endif

#ifdef GENERATE_KEYS
			
			std::cout << "...Computing Gradient for Image...";
			gradFilter[i] = GradientFilterType::New();
			gradFilter[i]->SetInput(scaleImage[i]);
			// Do this in pixel space
			gradFilter[i]->SetUseImageSpacing(false);
			gradFilter[i]->Update();
			gradImage[i] = gradFilter[i]->GetOutput();
			hgradImage[i] = this->GetHypersphericalCoordinates(gradImage[i]);

			gradMagFilter[i] = GradientMagFilterType::New();
			gradMagFilter[i]->SetInput(scaleImage[i]);
			// Do this in pixel space
			gradMagFilter[i]->SetUseImageSpacing(false);
			gradMagFilter[i]->Update();
			gradMagImage[i] = gradMagFilter[i]->GetOutput();

			std::cout << "...Done\n";
#endif

			// ...Compute Gaussians
			for (int j = 0; j < m_GaussianImagesNumber; ++j) {
#ifdef VERBOSE
				std::cout << "Setting Up Gaussian Filter " << i << "-" << j << "...";
				std::cout.flush();
#endif
				/* Variance is square of the sigma
	* sigma = (2^(j/s)*sigma)
	*/

				double variance = this->GetGaussianScale(j);
				variance *= variance;

				if (j == 0) {

					gaussianFilterZero = GaussianFilterType::New();

					gaussianFilterZero->SetVariance(variance);
					gaussianFilterZero->SetInput( scaleImage[i] );
					// pixel-wise smoothing
					gaussianFilterZero->SetUseImageSpacing(false); 
					try {
						gaussianFilterZero->Update();
					}
					catch( itk::ExceptionObject & excep ) {
						std::cerr << "Exception caught !" << std::endl;
						std::cerr << excep << std::endl;
					}

					gaussianImageZero = gaussianFilterZero->GetOutput();
				}

				if (j == m_GaussianImagesNumber - 1) {

					gaussianFilterNextOctave = GaussianFilterType::New();

					gaussianFilterNextOctave->SetVariance(variance);
					gaussianFilterNextOctave->SetInput( scaleImage[i] );
					// pixel-wise smoothing
					gaussianFilterNextOctave->SetUseImageSpacing(false); 
					try {
						gaussianFilterNextOctave->Update();
					}
					catch( itk::ExceptionObject & excep ) {
						std::cerr << "Exception caught !" << std::endl;
						std::cerr << excep << std::endl;
					}

					gaussianImageNextOctave = gaussianFilterNextOctave->GetOutput();
				}


				gaussianFilterBinary[j] = GaussianFilterType::New();

				gaussianFilterBinary[j]->SetVariance(variance);
				gaussianFilterBinary[j]->SetInput( scaleImageBinary[i] );
				// pixel-wise smoothing
				gaussianFilterBinary[j]->SetUseImageSpacing(false); 
				try {
					gaussianFilterBinary[j]->Update();
				}
				catch( itk::ExceptionObject & excep ) {
					std::cerr << "Exception caught !" << std::endl;
					std::cerr << excep << std::endl;
				}

				gaussianImageBinary[j] = gaussianFilterBinary[j]->GetOutput();

#ifdef DEBUG
				char filename[256];
				sprintf(filename, "gauss-%d-%d.png", i, j);
				this->writeImage(gaussianImageBinary[j], filename);
#endif

#ifdef VERBOSE
				std::cout << "Done\n";
				std::cout.flush();
#endif
			}
			
			// ...Compute Difference of Gaussians
			for (int j = 0; j < (m_DifferenceOfGaussianImagesNumber); ++j) {
#ifdef VERBOSE
				std::cout << "Setting Up DoG Binary Filter " << i << "-" << j << "...";
				std::cout.flush();
#endif
				dogFilterBinary[j] = DifferenceFilterType::New();
				dogFilterBinary[j]->SetInput1( gaussianImageBinary[j] );
				dogFilterBinary[j]->SetInput2( gaussianImageBinary[j+1] );
				dogFilterBinary[j]->Update();
				dogImageBinary[j] = dogFilterBinary[j]->GetOutput();

#ifdef DEBUG
				char filename[256];
				sprintf(filename, "dog-%d-%d.png", i, j);
				this->writeImage(dogImageBinary[j], filename);
#endif

#ifdef VERBOSE
				std::cout << "Done\n";
				std::cout.flush();
#endif 
			}

			for (int j=1; j < (m_DifferenceOfGaussianImagesNumber - 1); ++j) {
				// Search the dogImages for local maxima,  w.r.t. corresponding
				// point in the scale above and below
				// level 0 is the "doubled" image
				// Iterate over the various doG filters
				// Only use the middle dogs (ones with both neighbours above and below)
				// Iterate over each position in the dog filter
				typedef itk::ImageRegionIteratorWithIndex< TFixedImageType > 
				ImageIteratorType;

				IndexType regionStart;
				// Avoid the edges
				for (int k=0; k < VDimension; ++k)
				regionStart[k] = 1;

				typename TFixedImageType::SizeType regionSize = 
				dogImageBinary[j]->GetLargestPossibleRegion().GetSize();

#ifdef VERBOSE
				std::cout << "Searching for Extrema in DoG Binary Image " << i << "-" << j;
				std::cout << " ( ";
				for (int k=0; k < VDimension; ++k)
				std::cout << regionSize[k] << " ";
				std::cout << ") Scale " << currScale << "\n";
				std::cout.flush();
#endif

				// Avoid far edge
				for (int k=0; k < VDimension; ++k)
				regionSize[k] -=  2;
				
				typename TFixedImageType::RegionType itregion;
				itregion.SetIndex(regionStart);
				itregion.SetSize(regionSize);
				
				ImageIteratorType pixelIt(dogImageBinary[j],
				itregion);
				
				for ( pixelIt.GoToBegin(); !pixelIt.IsAtEnd(); ++pixelIt) {
					// Make sure to start sufficiently into the image so that all
					// neighbours are present
					IndexType pixelIndex = pixelIt.GetIndex();
					typename TFixedImageType::PixelType pixelValue = pixelIt.Get();

					PointType point;
					dogImageBinary[j]->TransformIndexToPhysicalPoint (pixelIndex, point);

#ifdef ERROR_CHECK
					std::cerr << "Checking ( ";
					for (int k = 0; k < VDimension; ++k)
					std::cerr << pixelIndex[k] << " ";
					std::cerr << ") = " << pixelValue <<"\n";
#endif

					// Compare to the 8 immediate neighbours
					bool isMax=true;
					bool isMin=true;

					this->CheckLocalExtrema(dogImageBinary[j], 
					pixelIndex, pixelValue, isMax, isMin, false);

					if (!isMax && !isMin) continue;
					
					// Compare to scale above
					if (j < (m_GaussianImagesNumber-1)) {
#ifdef DEBUG_VERBOSE
						std::cout << "...Checking Scale Above\n";
#endif
						//dogImage[i+1][j]->TransformPhysicalPointToIndex (point, tmpIndex);

						this->CheckLocalExtrema(dogImageBinary[j+1], 
						pixelIndex, pixelValue, isMax, isMin, true);
					}
					if (!isMax && !isMin) continue;

					// Compare to scale below
					if (j > 0) {
#ifdef DEBUG_VERBOSE
						std::cout << "...Checking Scale Below\n";
#endif
						//dogImage[i-1][j]->TransformPhysicalPointToIndex (point, tmpIndex);

						this->CheckLocalExtrema(dogImageBinary[j-1], 
						pixelIndex, pixelValue, isMax, isMin, true);
					}
					if (!isMax && !isMin) continue;
					
					//Setup variables to store reference and HGrad images when getting feature descriptors
					FixedImagePointer referenceImage;
					GradientImageType::Pointer hypersphericalGradientImage;

					referenceImage = gaussianImageZero;
					hypersphericalGradientImage = hgradImage[i];

#ifndef BAD_POINT_REJECTION

					// Check if it is sufficiently large (absolute value)
					if (fabs(pixelValue) < m_MinKeypointValue) {
						++numReject;
						continue;
					}

#else	
					//Frankie: Will want to add in the extra bit with the Hessian checks here...
					//Frankie: Maybe the std and abs pixel value from original image check here as well...

					m_VectorType displacementError;

					if (DoGAndEdgePointRatioThreshold(dogImageBinary[j], pixelIndex, pixelValue, edge_point_ratio * pow(edge_point_ratio_scale,i), &displacementError) == false) {
						++numReject;
						continue;
					}
					
#ifdef INTERPOLATE_FEATURES

					//std::cout << point << std::endl;

					//Adjust point index to new location. (real point = original point - displacement error)
					for (int tmp = 0; tmp < VDimension; tmp++) {
						point[tmp] = point[tmp] - displacementError[tmp];
					}

					//std::cout << point << std::endl <<std::endl;

					//Now we want to resample the scaleImage[i], gradImage[i] and hgradImage[i] images to the actual corners
					
					//Setting up the translation filter
					TranslationFilterType::Pointer translationTransform = TranslationFilterType::New();
					TranslationFilterType::OutputVectorType translation;
					
					for (int tmp = 0; tmp < VDimension; tmp++) {
						translation[tmp] = - displacementError[tmp];
					}

					translationTransform->Translate(translation);

					//Set up resample filter to apply the transform
					typedef	itk::BSplineInterpolateImageFunction< TFixedImageType, double > BSplineInterpolatorType;
					BSplineInterpolatorType::Pointer BSplineInterpolator = BSplineInterpolatorType::New();
					ResampleFilterType::Pointer translationFilter = ResampleFilterType::New();
					IdentityTransformType::Pointer identityTransform = IdentityTransformType::New();

					translationFilter->SetTransform(identityTransform);
					//translationFilter->SetTransform(translationTransform);	//Don't need to shift if I just change original origin before resample
					translationFilter->SetInterpolator(BSplineInterpolator);

					//Setup size and spacing to be the same as the input image
					const TFixedImageType::RegionType& imageRegion = scaleImage[i]->GetLargestPossibleRegion();
					const TFixedImageType::SizeType& imageSize = imageRegion.GetSize(); 
					const TFixedImageType::SpacingType& imageSpacing = scaleImage[i]->GetSpacing();
					const TFixedImageType::DirectionType& imageDirection = scaleImage[i]->GetDirection();
					
					//Adjust output origin based on our translation
					const TFixedImageType::PointType& imageOriginOld  = scaleImage[i]->GetOrigin();
					double imageOriginNew[VDimension];

					for (int tmp = 0; tmp < VDimension; tmp++) {
						imageOriginNew[tmp] = imageOriginOld[tmp] - displacementError[tmp];
					}
					
					




					//Setup image container with the exact same parameters as the output image from the resmple filter in order to
					//calculate the exact index values that we want to calculate to get the feature descriptors

					FixedImagePointer TranslatedImage = TFixedImageType::New();
					
					IndexType TranslatedImageStartIndex;
					TranslatedImageStartIndex.Fill(0);

					SizeType TranslatedImageSize;
					TranslatedImageSize = imageSize;

					RegionType TranslatedImageRegion(TranslatedImageStartIndex, TranslatedImageSize);
					TranslatedImage->SetRegions(TranslatedImageRegion);

					TranslatedImage->SetOrigin(imageOriginNew);
					TranslatedImage->SetSpacing(imageSpacing);
					TranslatedImage->SetDirection(imageDirection);


					//Define the region that we only care about when calculating for interpolated descriptor after transform
					
					//We want to get the region about 1.5 times in size of the actual region used for calculating the descriptor
					//just to be safe.

					SizeType RequestedSize;
					IndexType RequestedStartIndex;

					IndexType CenterIndex;

					TranslatedImage->TransformPhysicalPointToIndex (point, CenterIndex);
					
					//Calculate the starting index and size based off of this initial center index value

					//Requested region will be 1.5*(2m_SIFTHalfWidth - 1) in size centered about CenterIndex
					for (int tmp = 0; tmp < VDimension; tmp++) {
						RequestedStartIndex[tmp] = CenterIndex[tmp] - 1.5*m_SIFTHalfWidth;
						RequestedSize[tmp] = 1.5*(2*m_SIFTHalfWidth - 1);

						//Clamping to image edge if neccessary

						if (RequestedStartIndex[tmp] < 0) {
							RequestedStartIndex[tmp] = 0;
						}

						if (RequestedStartIndex[tmp] + RequestedSize[tmp] >= TranslatedImageSize[tmp]) {
							RequestedSize[tmp] = TranslatedImageSize[tmp] - RequestedStartIndex[tmp] - 1;
						}
					}

					RegionType RequestedRegion(RequestedStartIndex, RequestedSize);

					//Recalculate origin of requested image
					OriginType RequestedOrigin;
					TranslatedImage->TransformIndexToPhysicalPoint(RequestedStartIndex, RequestedOrigin);

					//Now compute the requested region for the translated image
					//std::cout<<"...Recomputing Scaled Image for Region Around Point..."<<std::endl;

					translationFilter->SetSize(RequestedSize);
					translationFilter->SetOutputSpacing(imageSpacing);
					translationFilter->SetOutputOrigin(RequestedOrigin);
					translationFilter->SetOutputDirection(imageDirection);
					translationFilter->SetInput(scaleImage[i]);

					//translationFilter->GetOutput()->SetRequestedRegion(RequestedRegion);
					//translationFilter->Update();
					TranslatedImage = translationFilter->GetOutput();



					/*translationFilter->SetSize(imageSize);
					translationFilter->SetOutputSpacing(imageSpacing);
					translationFilter->SetOutputOrigin(imageOriginNew);
					translationFilter->SetOutputDirection(imageDirection);
					translationFilter->SetInput(scaleImage[i]);

					translationFilter->GetOutput()->SetRequestedRegion(RequestedRegion);
					translationFilter->Update();

					//std::cout<<std::endl<<"Translated Image Stats:"<<std::endl;
					//ImageStats<TFixedImageType>(TranslatedImage);


					//Extracting region of interest from the scale image
					RegionOfInterestFilterType::Pointer ROIFilter = RegionOfInterestFilterType::New();
					ROIFilter->SetRegionOfInterest(RequestedRegion);
					ROIFilter->SetInput(translationFilter->GetOutput());
					TranslatedImage = ROIFilter->GetOutput();*/

					
					
					//Recomputing Gradient about the point 
					//std::cout << "...Recomputing Gradient for Image..."<<std::endl;

					GradientFilterType::Pointer ExtractedTranslatedGradFilter = GradientFilterType::New();
					GradientImageType::Pointer ExtractedTranslatedGradImage = GradientImageType::New();
					GradientImageType::Pointer ExtractedTranslatedHGradImage = GradientImageType::New();


					ExtractedTranslatedGradFilter->SetInput(TranslatedImage);
					// Do this in pixel space
					ExtractedTranslatedGradFilter->SetUseImageSpacing(false);
					ExtractedTranslatedGradFilter->Update();
					ExtractedTranslatedGradImage = ExtractedTranslatedGradFilter->GetOutput();

					//Recompute the gradint in hyperspherical coordinates
					ExtractedTranslatedHGradImage = this->GetHypersphericalCoordinates(ExtractedTranslatedGradImage);
				

					referenceImage = TranslatedImage;
					hypersphericalGradientImage = ExtractedTranslatedHGradImage;

#endif		

#endif
					
					// Passed all checks:

#ifdef DEBUG
					std::cout << point << std::endl;
#endif
					m_KeypointSet->SetPoint( m_PointsCount, point);
#ifdef GENERATE_KEYS

					if (matchToScaleOrExtrema) {
						
						FeatureType featureData = this->GetFeaturesWithExtraData( referenceImage, hypersphericalGradientImage, point, this->GetGaussianScale(j), 2);
						int featureDataSize = featureData.size();


						//Assign values of the scale and extrema of the point into the feature data array
						featureData[featureDataSize - 2] = this->GetGaussianScale(j);

						//0 = minimum, 1 = maximum
						if (isMax) {
							featureData[featureDataSize - 1] = 1;
						} else {
							featureData[featureDataSize - 1] = 0;
						}


						m_KeypointSet->SetPointData( m_PointsCount, featureData);

					} else {

						// Generate features
						// Space used is the (smoothed) original image)
						m_KeypointSet->SetPointData( m_PointsCount, 
						this->GetFeatures( referenceImage, 
						hypersphericalGradientImage, point,
						this->GetGaussianScale(j)));

					}
#else
					m_KeypointSet->SetPointData( m_PointsCount, currScale);
#endif
					++m_PointsCount;

					if (isMax) {
						// Maxima detected.  
						++numMax;
#ifdef DEBUG
						std::cout << "Found Maxima! ";
#endif
					}
					if (isMin) {
						// Minima detected.  
						++numMin;
#ifdef DEBUG
						std::cout << "Found Minima! ";
#endif
					}
				}
#ifdef VERBOSE
				std::cout << "Acc. Num Max: " << numMax 
				<< "\nAcc. Num Min: " << numMin 
				<< "\nAcc. Num Reject: " << numReject 
				<< std::endl;
				std::cout.flush();
#endif
			}
			currScale *= m_ScalingFactor;
		}
		
#ifdef VERBOSE
		std::cout << "Total Num Max: " << numMax 
		<< "\nTotal Num Min: " << numMin 
		<< "\nTotal Num Reject: " << numReject
		<< std::endl;	
		std::cout.flush();
#endif
		return m_KeypointSet;
	}
	
#ifdef BAD_POINT_REJECTION

template <class TFixedImageType, unsigned int VDimension> 
	bool
	ScaleInvariantFeatureImageFilter<TFixedImageType,VDimension>
	::DoGAndEdgePointRatioThreshold(FixedImagePointer image, IndexType pixelIndex, PixelType pixelValue, float edge_point_ratio, m_VectorType * displacementErrorPtr) {

		//Frankie: Will want to add in the extra bit with the Hessian checks here...
		//Frankie: Maybe the std and abs pixel value from original image check here as well...
					
		if (fabs(pixelValue) < m_DoGThreshold) {
			return false;
		}

		//Check if it is a poorly defined edge point.

		//First, need to calculate Hessian.
		m_HessianMatrixType hessianMatrix;
		
		//Do Hessian calculation for only dimensional space.
		this->CalculateHessianMatrix(image, pixelIndex, &hessianMatrix);

		//See if we need to change the edgeThreshold value according to an input value from user. Otherwise, set as default.
		float edgeThreshold = 0;

		if (edge_point_ratio > 0) {
			edgeThreshold = edge_point_ratio;
			//std::cout<<"Garbage: "<<edgeThreshold<<std::endl;
		} else {
			edgeThreshold = m_EdgeThreshold;
		}



		//1) Direct implementation of the ratio (Probably won't work)
		/*//Get VNL Matrix from Hessian Matrix
		HessianDeterminant = vnl_determinant(hessianMatrix.GetVnlMatrix());

		//Calculate Trace of VNL Matrix (manually since vnl_trace() is not included in ITK...)
		for (int traceIndex = 0; traceIndex < VDimension; traceIndex ++) {
			HessianTrace = HessianTrace + hessianMatrix(traceIndex,traceIndex);
		}
		m_RatioEdgeThreshold = pow((m_EdgeThreshold+(VDimension-1)),VDimension)/m_EdgeThreshold;
		if (fabs(pow(HessianTrace,(int)VDimension)) >= fabs(m_RatioEdgeThreshold*HessianDeterminant)) {
			return false;
		}*/




		//2) More loose implementaiton of the ratio.
		/*//Get VNL Matrix from Hessian Matrix
		HessianDeterminant = vnl_determinant(hessianMatrix.GetVnlMatrix());

		//Calculate Trace of VNL Matrix (manually since vnl_trace() is not included in ITK...)
		for (int traceIndex = 0; traceIndex < VDimension; traceIndex ++) {
			HessianTrace = HessianTrace + hessianMatrix(traceIndex,traceIndex);
		}
		
		// Assume ratio between largest and second largest eigenvalues is less than R, second largest to third largest is less than R, third largest and forth largest is less than R, ...
		float m_RatioEdgeThreshold_numerator = 0;
		float m_RatioEdgeThreshold_denominator = 0;
		int m_RatioEdgeThreshold_denominator_pow = 0;
		

		for(int tmp1 = 0; tmp1 < VDimension; tmp1++) {
			m_RatioEdgeThreshold_numerator = m_RatioEdgeThreshold_numerator + pow(edgeThreshold,tmp1);
		}

		m_RatioEdgeThreshold_numerator = pow(m_RatioEdgeThreshold_numerator, (int)VDimension);

		for (int tmp1 = 0; tmp1 < VDimension; tmp1++) {
			m_RatioEdgeThreshold_denominator_pow = m_RatioEdgeThreshold_denominator_pow + tmp1;
		}

		m_RatioEdgeThreshold_denominator = pow(edgeThreshold, m_RatioEdgeThreshold_denominator_pow);
		
		m_RatioEdgeThreshold = m_RatioEdgeThreshold_numerator/m_RatioEdgeThreshold_denominator;

		if (fabs(pow(HessianTrace,(int)VDimension)) >= fabs(m_RatioEdgeThreshold*HessianDeterminant)) {
			return false;
		}*/



		//3) Compare each pairs of directionality independently and get their principal curvature ratios. Feature point is found if at least
		//   one of the direcitonality pair passes the corner test.

		//Run the 2 by 2 ratio test for each direction combination and see if at least 1 of them is valid
		bool validCorner = false;
		double HessianTrace = 0.0;
		double HessianDeterminant = 0.0;
		double m_RatioEdgeThreshold_numerator = 0;
		double m_RatioEdgeThreshold_denominator = 0;

		for (int row = 0; row < VDimension-1; row++) {
			for (int col = row+1; col < VDimension; col++) {

				//Need to break down the big hessian matrix into the individual 2 by 2 hessian matrix comparing 2 directions only
				itk::Matrix<double, 2, 2> hessianMatrix2by2;

				hessianMatrix2by2(0,0) = hessianMatrix(row,row);
				hessianMatrix2by2(0,1) = hessianMatrix(row,col);
				hessianMatrix2by2(1,0) = hessianMatrix(col,row);
				hessianMatrix2by2(1,1) = hessianMatrix(col,col);


				//Get VNL Matrix from Hessian Matrix
				HessianDeterminant = vnl_determinant(hessianMatrix2by2.GetVnlMatrix());

				//Calculate Trace of VNL Matrix (manually since vnl_trace() is not included in ITK...)
				HessianTrace = hessianMatrix2by2(0,0) + hessianMatrix2by2(1,1);

				double m_RatioEdgeThreshold_numerator = pow(edgeThreshold + 1, 2);
				double m_RatioEdgeThreshold_denominator = edgeThreshold;
				
				m_RatioEdgeThreshold = m_RatioEdgeThreshold_numerator/m_RatioEdgeThreshold_denominator;
				
				if (fabs(pow(HessianTrace,2)) < fabs(m_RatioEdgeThreshold*HessianDeterminant)) {
					validCorner = true;
					break;
				}
			}
			if (validCorner) {
				break;
			}
		}

		if (!validCorner) {
				return false;
		}

		//If function reaches this part, then it passed all the tests

		//Return displacement error if we need to
		if (displacementErrorPtr != NULL) {

			//Calculate gradient vector
			m_VectorType gradientVector;
			this->CalculateFirstOrderGradient(image,pixelIndex,&gradientVector);

			//Get spacings of the input image
			const SpacingType imageSpacing = image->GetSpacing();

			//Get inverse of Hessian Matrix
			m_HessianMatrixType inverseHessianMatrix = hessianMatrix.GetInverse();
			
			//Calculate displacement error in index spacing
			m_VectorType displacementErrorIndex = -(inverseHessianMatrix*gradientVector);

			//Calculate displacment in image distance
			for (int tmp = 0; tmp < VDimension; tmp++) {
				(*displacementErrorPtr)[tmp] = displacementErrorIndex[tmp] * imageSpacing[tmp];
			}
		}

		return true;

	}
#endif


} // end namespace itk

#endif
