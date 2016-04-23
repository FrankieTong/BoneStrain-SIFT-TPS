# README #

This README is a preliminary draft for the instructions to setup and run the algorithm to use SIFT and TPS to measure strain. The current implementation is a bit ugly in terms of structure, but is being documented anyways as both a backup as well as a starting point to get this documentation process started.

-----------------------

### What is this repository for? ###

This repository is for the C++ code used to perform Scale Invarient Feature Transform (SIFT) for feature finding and matching as well as the Thin Plate Splines (TPS) method for strain calculation.

Version: (Last Edited): April 20, 2016

-----------------------

### How do I get set up? ###

This code was originally built using the Micosoft Visual Studio 2008 x64 compiler from the full installation of the program (not trial). Note that the x64 compiler has to be toggled on during the installation process in the custom installation settings in order for it to be installed (do not blindly click standard installation during the instalation process). I have yet to try this with a VS2009 trial version so it will be up to you to trouble shoot.

After VS2009 is installed, the next step is to build ITK on the computer. ITK v4.6.0 was used for building this code although newer version should work as well. In order to install ITK, you will first need to install CMake in order to build ITK on your computer. Follow the instructions to install CMake and ITK on your computer.

Once CMake and ITK have been installed and properly linked to the PATH variable on your computer, we need build this project using CMake while making sure the ITK dependencies are porperly assigned. This is easily done by chaningthe `ITK_DIR` line to point to where the ITK _binaries_ are installed. First open the CMake gui and point the source to where you stored this project's repository. Then assign the place where you want the binaries for this project to build. Finally assign the location of the `ITK_DIR`. Choose VS2009 x64 compiler when prompted (first time only) and hit configure and generate. Once the VS2009 files are generated, you can then open project using VS2009 and compile the code. I tend to try and compile in debug mode when testing things and then in relwithdebug mode when actually using the program. There might be some troubleshooting that needs to be done in order for this to compile properly.  

### Running the Program ###

After the program is compiled, open the directory containing the executable (should be called `main.exe`). From there open a command window to that program (either by typing `cmd` in the windows explorer address of the folder where the executable is stored or by typing `cmd` in the run program and navigating to the applicable folder).

After you have navigated to the executable with the command window, start the program by typing `main.exe` into the command window. If the program is built properly, it should display additional parameters that could be used to have `main.exe` execute what you want.

#### Executable Commands ####

Although there are a lot of parameters and options originally listed for the `main.exe` executable, there are only 2 relevant ones that a user would worry about:

Execution Mode                 | Description   
------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
`--transform`                  | Applies a given transform described by a vector field with the same dimensions and spacing as the input image and outputs the transformed image. Useful for generating the deformed image.
`--standardResampleScaleMatch` | Applies the SIFT algorithm to the input fixed and moving images and then generates a uniform displacement field by applying TPS onto the matched feature points. Outputs the matched SIFT points as well as the uniform displacement field.

An example of an execution of a mode is `main.exe --transform <Nessesary Parameters>` into the command window. For reference on what type of parameters are needed for each mode, just leave the additional parameters section blank after typing in the mode and the program will tell you what are the additional parameters needed. More information on the input parameters for each mode is described below.

#### Trasnform Mode ####

`main.exe --transform <inputImageFile1> <inConfigFile1> <inputVectorField1> <inVectorConfigFile1> <outputImageFile> <outConfigFile>`

##### Input Parameters #####

Parameter                      | Description   
------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
`inputImageFile1`              | File name of the input image file to be transformed using the displacement field. The input image has to be in a `.raw` format with the individual pixel values in float (i.e: no header).
`inConfigFile1`                | File name of the configuration file that describes the parameters of `inputImageFile1`. The configuration file description is described in a later section of this README.
`inputVectorField1`            | File name of the input image file that contains the displacement field used for the transform. The input image has to be in a `.raw` format with the individual pixel values in float (i.e: no header).
`inVectorConfigFile1`          | File name of the configuration file that describes the parameters of `inputVectorField1`. The configuration file description is described in a later section of this README.
`outputImageFile`              | File name of the output image file. The output file is a `.raw` file with float values for pixels that should have the same configuration file as `inConfigFile1`.
`outConfigFile`                | File name of the configuration file that describes the parameters of `outputImageFile`. The configuration file description is described in a later section of this README.

##### Outputs #####

* Output image file in `.raw` format
* Output config file

### Standard Resample Scale Match  Mode ###

`main.exe --standardResampleScaleMatch inputImageFile1 inputImageFile2 inImageFile1Min inImageFile1Max inImageFile2Min inImageFile2Max inputImageFileBinary1 inputImageFileBinary2 inImageFileBinary1Min inImageFileBinary1Max inImageFileBinary2Min inImageFileBinary2Max iConfigFile1 iConfigFile2 outName resampleSpacingFactor relativePosition(origin/center) outputDisplacementResolution(original/resampled/both/#) inRemoveEdgePointsPixel inScaleStart inScaleEnd inEdgePointRatio inEdgePointRatioScale [optional input for matching on same scale only] [optional input for matching on same extremas only] [optional input for turning off nearest neighbour intepolator] [optional input for displacement vector limits]`

Parameter                      | Description   
------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
`inputImageFile1`              | File name of the input image file to be used as the fixed image for SIFT. The input image has to be in a `.raw` format with the individual pixel values in float (i.e: no header).
`inputImageFile2`              | File name of the input image file to be used as the fixed image for SIFT. The input image has to be in a `.raw` format with the individual pixel values in float (i.e: no header).
`inImageFile1Min`              | The minimum intensity threshold of the fixed image.
`inImageFile1Max`              | The maximum intensity threshold of the fixed image.
`inImageFile2Min`              | The minimum intensity threshold of the moving image.
`inImageFile2Max`              | The maximum intensity threshold of the moving image.
`inputImageFileBinary1`        | File name of the input image file to be for the binarized version of the fixed image. The input image has to be in a `.raw` format with the individual pixel values in float (i.e: no header).
`inputImageFileBinary2`        | File name of the input image file to be for the binarized version of the moving image. The input image has to be in a `.raw` format with the individual pixel values in float (i.e: no header).
`inImageFileBinary1Min`        | The minimum intensity threshold of the binarized fixed image.
`inImageFileBinary1Max`        | The maximum intensity threshold of the binarized fixed image.
`inImageFileBinary2Min`        | The minimum intensity threshold of the binarized moving image.
`inImageFileBinary2Min`        | The maximum intensity threshold of the binarized moving image.
`inConfigFile1`                | File name of the configuration file that describes the parameters of the input fixed image. The configuration file description is described in a later section of this README.
`inConfigFile2`                | File name of the configuration file that describes the parameters of the input fixed image. The configuration file description is described in a later section of this README.
`outName`                      | Starting prefix that the output files from this program produces. Do not need to add file type name at the end.
`resampleSpacingFactor`        | Resample the input images with this scaling factor before running the SIFT algorithm. The input image will be downsampled by this factor (i.e: `resampleSpacingFactor = 0.5` means downsample by a factor of 0.5, or upsample by a factor of 2). Should be a value between 0 and 1.
`relativePosition`             | Allows the user to choose how the coordinates of the upsampled image is set relative to the original image. `relativePosition = origin` sets the origin of the resampled image to the origin of the original image while `relativePosition = center` sets the center of the resampled image to be the same as the original image.
`outputDisplacementResolution` | Identifies resolution of the output TPS field. `outputDisplacementResolution = original` sets the output TPS field at the same resolution as the input image resolution. `outputDisplacementResolution = resampled` sets the output TPS field to the resampled image resolution. `outputDisplacementResolution = both` outputs the TPS field in both the original and the resampled resolution. Finally, a number can be used as the input which will make the TPS field resolution that of which would be created by upsampling relative to the original input image resolution. In this case, value should be less than or equal to 1.
`inRemoveEdgePointsPixel`      | Identifies the percentrage of the outer image where found feature points are discarded. (i.e: `inRemoveEdgePointsPixel = 5` means remove the feature points found within 5% of the outer edges of the image from the matching list).
`inScaleStart`                 | Identifies the scale level in which the SIFT algorithm initailly starts looking for feature points. `inScaleStart = 0` would tell SIFT to start lookung for feature points at the base resampled image resolution and can be considered the default.
`inScaleEnd`                   | Identifies the scale level in which the SIFT algorithm stops looking for feature points. `inScaleEnd = 2` means to stop looking for feature points after going through 2 octaves, which would be the default for the normal SIFT algorithm.
`inEdgePointRatio`             | Identifies the edge point ratio value threshold where feature points found to have a corner ratio of more than this value will be rejected. Value should be greater than 1.
`inEdgePointRatioScale`        | Idenfifies the scaling factor of the `inEdgePointRatio` value that changes for each increasing octave. Value should be greater or equal to 1.
`[optional input for matching same scale only]` | Boolean value that determines if the matching of feature points will be conducted only on feature points that have the same scale. `false` is typically the default to ignore this parameter.
`[optional input for matching on same extremas only]` | Boolean value that determines if matching of feature points will be conducted only on feature points that have the same extrema direction. `false` is typically the default to ignore this parameter.
`[optional input for turning off nearest neighbour intepolator]` | Boolean value that determines if the interpolator for the binary images uses a nearest pixel interpolator rather than a BSplineInterpolator. `false` is typically the default to ignore this parameter.
`[optional input for displacement vector limits]` | Identifies the limit range of the physical constraint limit filter. Values are given as a white space seperated list of minimum and maximum boundary range centered about the fixed point for each dimensional direction given in absolute distance (i.e: For a 2D image using `-0.05 0.05 -1 1` for this variable means a displacement vector will be accepted in if it falls within the boundaries of having an x displacement value between -0.5 and 0.5 units and a y displacement value between -1 and 1 units).

##### Outputs #####

* Keypoint list of keypoints found by SIFT in the fixed image in a `.txt` format
* Keypoint list of keypoints found by SIFT in the fixed image in a `.landmarkASCII` format (Amira) for ease of import into Amira
* Keypoint list of keypoints found by SIFT in the moving image in a `.txt` format
* Keypoint list of keypoints found by SIFT in the fixed image in a `.landmarkASCII` format (Amira) for ease of import into Amira
* Matched keypoint list of keypoints found by SIFT in a `.txt` format
* Matched keypoint list of keypoints found by SIFT in a `.landmarkASCII` format (Amira) for ease of import into Amira
* Matched keypoint list of keypoints found by SIFT in a `.m` format (MatLab) for ease of import into MatLab
* Output displacement field in the orignal image resolution/resampled/both/# image resolution in a `.raw` format
* Output config file for each output displacement field in a `.txt` format


### Image Reading and Writing ###

The input image type is a `.raw` type image with a `float` datatype for the individual pixels. The input image dimensionality is currently set to 2 but can be easily changed to support any dimension by changing the `const unsigned int Dimension = 2;` in the source to the number of dimensions that you want to process with this program and recompiling the program.

The image and vector files for this program was orginally generated from Amira. Amira can be used to generate and save images in `.raw` format using the `Raw Data 3D (.raw)` option in the save as option when saving an image in Amira. Create an accompaning configuration file for the saved image that outlines the parameters of the image and use the raw image file along with the config file to import the image into the program.

#### Configuration File Format ####

Refer to `sample_files/image_config_file` for an example of the image configuration file format.

Example file from `iconfig.txt`

    //Contains parameters for configuration of the main.cxx program
    
    PixelType float
    PixelDimensionality 1
    ImageDimensionality 3
    DimensionSize 50 50 50
    SpacingSize 0.002 0.002 0.002
    Origin 0 0 0
    HeaderSize 0
    ByteOrder little

Parameter                      | Description   
------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
`PixelType`                    | The data type of the individual image pixels. Orignally intended to vary depending on what type is given, it is now irrelevant as the code is hard coded to accept float data type only.
`PixelDimensionality`          | The number of values each pixel stores. Typically 1 for a standard image or `ImageDimensionality` for a vector image.
`ImageDimensionality`          | The number of dimensions in the image file. (i.e.: `ImageDimensionality = 2` for a 2D image)
`DimensionSize`                | The space seperated list for the total size of the image in each dimension in pixels.
`SpacingSize`                  | The space seperated list for the total spacing between each pixel each dimension of the image.
`Origin`                       | The space seperated list for the coordinates of the origin of the image.
`HeaderSize`                   | The size of the inital header. Since I was using `.raw` files, this was typically 0
`ByteOrder`                    | The byte order/endianess of the image. Can either be set to `little` or `big`. Amira tends to save `.raw` images in `little` endian.

Other examples for the configuration of keypoint files and vector field configuration files can be found in  the`sample_files` folder as well.

### Who do I talk to? ###

Repo Owner: Frankie (Hoi-Ki) Tong <hoiki.tong@mail.utoronto.ca\>