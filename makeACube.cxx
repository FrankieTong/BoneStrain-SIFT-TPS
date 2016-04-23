/*
File:			makeACube.cxx

Description:	Function effectively sets the square from index value (0,0,0) to (50,50,50) of the input image to 0 and sets the square from index value (10,10,10) to (40,40,40) to 255.

Inputs:			OutputImage - (void *) Void pointer that points to the ITK::image type pointer that holds the image that will be modified by this function. Due to how this function was programmed, the pointer has to point to a "itk::Image< float, 3 >"  image.

Author:			Frankie (Hoi-Ki) Tong, Sep 6 2014

Changes:		Frankie (Hoi-Ki) Tong, Sep 6 2014
				Initial creation of the file.
*/

#include "makeACube.h"
 
int makeACube(void * OutputImage)
{
  typedef itk::Image< float, 3 >  ImageType;
 
  ImageType::Pointer *image = (ImageType::Pointer *) OutputImage;
 
  // Set all pixels to black
  for(unsigned int r = 0; r < 50; r++)
  {
      for(unsigned int c = 0; c < 50; c++)
      {
		for(unsigned int d = 0; d < 50; d++) {
		
			  ImageType::IndexType pixelIndex;
			  pixelIndex[0] = r;
			  pixelIndex[1] = c;
			  pixelIndex[2] = d;

			  (*image)->SetPixel(pixelIndex, 0);
			}
      }
  }
 
  // Set pixels in a square to one value
  for(unsigned int r = 10; r < 40; r++)
  {
      for(unsigned int c = 10; c < 40; c++)
      {
		for(unsigned int d = 10; d < 40; d++) {
		
			  ImageType::IndexType pixelIndex;
			  pixelIndex[0] = r;
			  pixelIndex[1] = c;
			  pixelIndex[2] = d;
	 
			  (*image)->SetPixel(pixelIndex, 255);
			}
      }
  }

  // Set pixels in a square to one value
/*  for(unsigned int r = 15; r < 35; r++)
  {
      for(unsigned int c = 15; c < 35; c++)
      {
		for(unsigned int d = 15; d < 35; d++) {
		
			  ImageType::IndexType pixelIndex;
			  pixelIndex[0] = r;
			  pixelIndex[1] = c;
			  pixelIndex[2] = d;
	 
			  (*image)->SetPixel(pixelIndex, 0);
			}
      }
  }*/

	return 0;

}


  