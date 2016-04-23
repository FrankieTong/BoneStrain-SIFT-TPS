/*
File:			itkImageIORaw.h

Description:	Header file for itlImageIORaw.txx. Class that defines methods for reading and writing ITK images to/from a raw data format and a easily readable image configuration file.

Author:			Frankie (Hoi-Ki) Tong, Oct 4 2014

Changes:		Frankie (Hoi-Ki) Tong, Oct 4 2014
				Initial creation of the file.
*/

#include "ImageIOConfig.h"
#include <string>

#include "itkRawImageIO.h"
#include "itkNumericSeriesFileNames.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"

#ifndef __itkImageIORaw_H__
#define __itkImageIORaw_H__

template <class PixelType, unsigned int ImageDimensionality> 
class itkImageIORaw:public ImageIOConfig {
		public:
		
			//Standard constructors and descructors
			itkImageIORaw();                         // constructor; initialize the list to be empty
			itkImageIORaw(std::string inputImageFileName, std::string inputConfigFileName, void * inputitkImage); // constructor; initialize with preset values
			~itkImageIORaw();						// destructor;

			itkImageIORaw(const itkImageIORaw &obj);	//copy constuctor
			itkImageIORaw& operator=(const itkImageIORaw & obj); //Assignment constructor
			
			//Standard get and set commands for data values in the class
			std::string getImageFileName() const;
			bool setImageFileName(std::string inputImageFileName);
			
			bool getitkImage(void * inputitkImage);
			bool setitkImage(void * inputitkImage);
						
			//Function calls for reading and writing raw images to and from ITK along with reading and writing to/from a easy to read image ocnfiguration file
			bool readRawImageKnown();
			bool writeRawImageKnown();

		private:

		//Private data classes
		std::string ImageFileName;	//file name of the image
		void * itkImage;			//image pointer to store ITK image in
  
		//The functions are meant to work as universal read and write functions. Limitations with template structure of ITK prevented this idea form working and thus are non-functional.
		bool readRawImage();	//Does not function. Use readRawImageKnown()
		bool writeRawImage();	//Does not function. Use writeRawImageKnown()
};

#include "itkImageIORaw.txx"

#endif