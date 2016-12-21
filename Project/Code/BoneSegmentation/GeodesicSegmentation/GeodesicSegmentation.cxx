// Libraries to include
#include "itkGeodesicActiveContourLevelSetImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkExtractImageFilter.h"
#include "itkTimeProbe.h"
#include "itkResampleImageFilter.h"
#include "itkFastMarchingImageFilter.h"
#include "itkAddImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkVarianceImageFunction.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkAndImageFilter.h"
#include "itkJoinSeriesImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"

//definitions of used types

// Type for input 3D Images
typedef   float   InputPixelType;
const     unsigned int    InputDimension = 3;
typedef itk::Image< InputPixelType, InputDimension >  InputImageType;

// Type for 2D slices internally extracted from 3D volume 
typedef float InternalPixelType;
const     unsigned int    InternalDimension = 2;
typedef itk::Image< InternalPixelType, InternalDimension >  InternalImageType;

// Type for segmented 2D Images 
typedef unsigned char BinaryPixelType;
const     unsigned int    BinaryDimension = 2;
typedef itk::Image< BinaryPixelType, BinaryDimension >  BinaryImageType;

// Type for output 3D volumes
typedef float OutputPixelType;
const     unsigned int    OutputDimension = 3;
typedef itk::Image< OutputPixelType, OutputDimension > OutputImageType;

typedef itk::VarianceImageFunction<InternalImageType,float> VarianceFilterType;
typedef itk::ExtractImageFilter< InputImageType, InternalImageType > ExtractFilterType;
typedef itk::ImageFileReader< InputImageType > ReaderType;
typedef itk::ImageFileWriter< OutputImageType  > WriterType;
typedef itk::BinaryThresholdImageFilter< InternalImageType, BinaryImageType > ThresholdingFilterType;
typedef itk::CurvatureAnisotropicDiffusionImageFilter< InternalImageType, InternalImageType >  SmoothingFilterType;
typedef itk::GradientMagnitudeRecursiveGaussianImageFilter< InternalImageType, InternalImageType >  GradientFilterType;
typedef itk::GeodesicActiveContourLevelSetImageFilter< InternalImageType, InternalImageType > GeodesicActiveContourFilterType;
typedef itk::SigmoidImageFilter <InternalImageType,InternalImageType> SigmoidFilterType;
typedef itk::BinaryBallStructuringElement< InternalImageType::PixelType,InternalDimension>  StructuringElementType;
typedef itk::BinaryDilateImageFilter <BinaryImageType, BinaryImageType, StructuringElementType> BinaryDilateImageFilterType;
typedef itk::FastMarchingImageFilter< InternalImageType, InternalImageType >  FastMarchingFilterType;
typedef itk::AndImageFilter <BinaryImageType> AndImageFilterType;
typedef itk::AddImageFilter <BinaryImageType, BinaryImageType, InternalImageType> AddImageFilterType;
typedef itk::MultiplyImageFilter<BinaryImageType,BinaryImageType,BinaryImageType> MultiplyFilterType;
typedef itk::JoinSeriesImageFilter<InternalImageType, OutputImageType> JoinSeriesImageFilterType;
typedef itk::CastImageFilter< BinaryImageType, InternalImageType > CastFilterType;
typedef itk::MinimumMaximumImageCalculator <InternalImageType> ImageCalculatorFilterType;
typedef itk::RescaleIntensityImageFilter< InternalImageType, InternalImageType > RescaleType;
 
//compute pixel variance at each pixel in the 2D image
void computeVarianceImage (VarianceFilterType::Pointer varfilter,
    InternalImageType::Pointer inputImage, InternalImageType::Pointer varianceImage)
{
	// Set the varianceImage to have the same size and origin as the input 2D image
    varianceImage->CopyInformation(inputImage);
    varianceImage->SetRegions(inputImage->GetLargestPossibleRegion());
    varianceImage->Allocate();
    varianceImage->FillBuffer(0);
    
    InternalImageType::IndexType pi;
 
	// Loop over all pixels in the image 
    for (unsigned x=0; x<=inputImage->GetLargestPossibleRegion().GetSize(0)-1; x++)
    {
        pi.SetElement(0,x);
        for (unsigned y=0; y<=inputImage->GetLargestPossibleRegion().GetSize(1)-1; y++)
        {
            pi.SetElement(1,y);       
			varianceImage->SetPixel(pi, varfilter->EvaluateAtIndex(pi));  // Calculate variance
		}
            
	}
 }

// main function
int main( int argc, char *argv[] )
{
  if( argc < 23 )   // Check to see if all parameters are provided
    {
    std::cerr << "Parameters not provided" << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImage outputImage";
	std::cerr << " FStartSlice1 FSeedx FSeedy TStartSlice1 TSeedx TSeedy";
	std::cerr << " FStartSlice2 FEndSlice2 FSeedx FSeedy TStartSlice2 TEndSlice2 TSeedx TSeedy";
	std::cerr << " FEndSlice3 FSeedx FSeedy TEndSlice3 TSeedx TSeedy";
    std::cerr << "" << std::endl;
    return EXIT_FAILURE;
    }


  itk::TimeProbe clock;
  // Measure execution time
  clock.Start();

  // Read input arguments
  int FStartSlice1 = atoi(argv[3]);   // Start slice of Femur section 1
  int FSeedx1 = atoi(argv[4]);        // Seed x location 
  int FSeedy1 = atoi(argv[5]);        // Seed y location 
  int TStartSlice1 = atoi(argv[6]);   // Start slice of Tibia section 1
  int TSeedx1 = atoi(argv[7]);        // Seed x location 
  int TSeedy1 = atoi(argv[8]);        // Seed y location

  int FStartSlice2 = atoi(argv[9]);   // Start slice of Femur section 2
  int FEndSlice2 = atoi(argv[10]);    // End slice of Femur section 2
  int FSeedx2 = atoi(argv[11]);       // Seed x location
  int FSeedy2 = atoi(argv[12]);       // Seed y location
  int TStartSlice2 = atoi(argv[13]);  // Start slice of Tibia section 2
  int TEndSlice2 = atoi(argv[14]);    // End slice of Tibia section 2
  int TSeedx2 = atoi(argv[15]);       // Seed x location
  int TSeedy2 = atoi(argv[16]);       // Seed y location

  int FEndSlice3 = atoi(argv[17]);    // End slice of Femur section 3
  int FSeedx3 = atoi(argv[18]);       // Seed x location
  int FSeedy3 = atoi(argv[19]);       // Seed y location 
  int TEndSlice3 = atoi(argv[20]);    // End slice of Tibia section 3
  int TSeedx3 = atoi(argv[21]);       // Seed x location
  int TSeedy3 = atoi(argv[22]);       // Seed y location

  // Initialize reader and writer objects
  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();
  reader->SetFileName( argv[1] );
  writer->SetFileName( argv[2] );
  
  // Read input volume
  try{
	reader->UpdateOutputInformation();
  }
  catch( itk::ExceptionObject & excep )
  {
	std::cerr << "Exception caught while reading Input File!" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }

  InputImageType::Pointer inputImage = reader->GetOutput();    // Input 3D volume

  // Check to see if writer can write to the given location
  try{
	writer->SetInput(inputImage);
	writer->Update();
  }
  catch( itk::ExceptionObject & excep )
  {
	std::cerr << "Exception caught while accessing write location!" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }

  // Object to join 2D slices to produce final 3D volume
  JoinSeriesImageFilterType::Pointer joinFilter = JoinSeriesImageFilterType::New();
  joinFilter->SetOrigin(inputImage->GetOrigin()[2]);       // Set origin in z direction from input Image
  joinFilter->SetSpacing(inputImage->GetSpacing()[2]);     // Set spacing in z direction from input Image

  // Loop over all slices in the Input image
  for (unsigned z=0; z<=reader->GetOutput()->GetLargestPossibleRegion().GetSize(2)-1; z++)
  {
    std::cout << "Processing Slice " << z+1 << std::endl;

	// Extract one slice from the volume
	InputImageType::RegionType inputRegion =
             reader->GetOutput()->GetLargestPossibleRegion();
	InputImageType::SizeType size = inputRegion.GetSize();
	size[2] = 0;
	InputImageType::IndexType start = inputRegion.GetIndex();
	const unsigned int sliceNumber = z;     // gets slice number from value of z
	start[2] = sliceNumber;
	InputImageType::RegionType desiredRegion;
	desiredRegion.SetSize(  size  );
	desiredRegion.SetIndex( start );

	ExtractFilterType::Pointer extractfilter = ExtractFilterType::New();  // object to extract 2D slice
	extractfilter->InPlaceOn();
	extractfilter->SetDirectionCollapseToSubmatrix();
	extractfilter->SetExtractionRegion( desiredRegion );
	extractfilter->SetInput(reader->GetOutput());
	extractfilter->Update();

	InternalImageType::Pointer image=extractfilter->GetOutput();	   // extracted 2D image

	// Get Maximum pixel value from 2D slice
	ImageCalculatorFilterType::Pointer imageCalculatorFilter
          = ImageCalculatorFilterType::New ();
	imageCalculatorFilter->SetImage(image);
	imageCalculatorFilter->Compute();
	float maxVal = imageCalculatorFilter->GetMaximum();

	// Object to rescale intensity values between 0 and 127 (optimal parameters used for that range)
	RescaleType::Pointer rescaler = RescaleType::New();
	rescaler->SetInput( image );
	rescaler->SetOutputMinimum(0);
	rescaler->SetOutputMaximum(127);

	// alpha and beta for sigmoid filter for Variance
	double alpha = -3.0;
	double beta = 15.0;
	
	if (maxVal > 511)         // if intensity range higher
	{
		rescaler->Update();
		image= rescaler->GetOutput();             // rescale image
		alpha = -0.5;
		beta = 3.0;
	}

	// Return zero intensity image for sections where bone does not exist
	InternalImageType::Pointer segmentedImage = InternalImageType::New();
	segmentedImage->CopyInformation(image);               // copy information from 2D image
	segmentedImage->SetRegions(image->GetLargestPossibleRegion());
    segmentedImage->Allocate();
    segmentedImage->FillBuffer(0);                        // All pixels have intensity zero

	// Slices at the beginning and at the end of the volume with no bone
	if (((z<FStartSlice1-1) && (z<TStartSlice1-1)) || ((z>FEndSlice3-1) && (z>TEndSlice3-1)))
	{
		joinFilter->SetInput(z,segmentedImage);            // add slice to joinFilter 
		continue;                                          // go to next slice
	}



	// Object to threshold the output of the level set from Geodesic Active Contours Segmentation 
	// For Femur intensity image output
	ThresholdingFilterType::Pointer thresholder1 = ThresholdingFilterType::New();
	thresholder1->SetLowerThreshold( -1000.0 );
	thresholder1->SetUpperThreshold(  0.0 );
	thresholder1->SetOutsideValue(  0  );
	thresholder1->SetInsideValue(  1 ); 
	// For Femur variance image output
	ThresholdingFilterType::Pointer thresholder2 = ThresholdingFilterType::New();
	thresholder2->SetLowerThreshold( -1000.0 );
	thresholder2->SetUpperThreshold(  0.0 );
	thresholder2->SetOutsideValue(  0  );
	thresholder2->SetInsideValue(  1 );
	// For Tibia intensity image output
	ThresholdingFilterType::Pointer thresholder3 = ThresholdingFilterType::New();
	thresholder3->SetLowerThreshold( -1000.0 );
	thresholder3->SetUpperThreshold(  0.0 );
	thresholder3->SetOutsideValue(  0  );
	thresholder3->SetInsideValue(  1 );
	// For Tibia variance image output
	ThresholdingFilterType::Pointer thresholder4 = ThresholdingFilterType::New();
	thresholder4->SetLowerThreshold( -1000.0 );
	thresholder4->SetUpperThreshold(  0.0 );
	thresholder4->SetOutsideValue(  0  );
	thresholder4->SetInsideValue(  1 );
  
	// Filter to perform Curvature Anisotropic diffusion smoothing
	SmoothingFilterType::Pointer smoothing1 = SmoothingFilterType::New();    // For Intensity image
	SmoothingFilterType::Pointer smoothing2 = SmoothingFilterType::New();    // For Variance image

	// Filter to compute Gradient Magnitudes
	GradientFilterType::Pointer gradientMagnitude1 = GradientFilterType::New();  // For Intensity image
	GradientFilterType::Pointer gradientMagnitude2 = GradientFilterType::New();  // For Variance image

	// Sigmoid Filters to generate edge images
	SigmoidFilterType::Pointer sigmoid1 = SigmoidFilterType::New();      // For Intensity image
	SigmoidFilterType::Pointer sigmoid2 = SigmoidFilterType::New();      // For Variance image

	// Fast Marching Filter to generate initial level set 
	FastMarchingFilterType::Pointer fastMarching1 = FastMarchingFilterType::New();  // For femur 
	FastMarchingFilterType::Pointer fastMarching2 = FastMarchingFilterType::New();  // For tibia

	typedef FastMarchingFilterType::NodeContainer  NodeContainer;   // Node to contain seeds 
	typedef FastMarchingFilterType::NodeType       NodeType;

	NodeContainer::Pointer seeds1 = NodeContainer::New();        // For Femur
	InternalImageType::IndexType  seedPosition1;
	double initialDistance1 = 0;          // Initial contour defined from the seed position

	if (z<FStartSlice2-1)
	{
		seedPosition1[0] = FSeedx1;
		seedPosition1[1] = FSeedy1;
		initialDistance1 = 5;
	}
	else if (z<FEndSlice2-1)
	{
		seedPosition1[0] = FSeedx2;
		seedPosition1[1] = FSeedy2;
		initialDistance1 = 10;      // Larger distance for central slices
	}
	else
	{
		seedPosition1[0] = FSeedx3;
		seedPosition1[1] = FSeedy3;
		initialDistance1 = 5;
	}


	NodeType node1;
	const double seedValue1 = - initialDistance1;
	node1.SetValue( seedValue1 );
	node1.SetIndex( seedPosition1 );

	seeds1->Initialize();
	seeds1->InsertElement( 0, node1 );

	fastMarching1->SetTrialPoints( seeds1 );         // Set seeds
	fastMarching1->SetSpeedConstant( 1.0 );          // Constant speed
	fastMarching1->SetOutputSize(
           extractfilter->GetOutput()->GetBufferedRegion().GetSize() );
	fastMarching1->SetOutputOrigin(extractfilter->GetOutput()->GetOrigin());
	fastMarching1->SetOutputDirection(extractfilter->GetOutput()->GetDirection());
	fastMarching1->SetOutputSpacing(extractfilter->GetOutput()->GetSpacing());


	NodeContainer::Pointer seeds2 = NodeContainer::New();        // For Tibia
	InternalImageType::IndexType  seedPosition2;

	double initialDistance2 = 0;
	if (z<TStartSlice2-1)
	{
		seedPosition2[0] = TSeedx1;
		seedPosition2[1] = TSeedy1;
		initialDistance2 = 5;
	}
	else if (z<TEndSlice2-1)
	{
		seedPosition2[0] = TSeedx2;
		seedPosition2[1] = TSeedy2;
		initialDistance2 = 10;
	}
	else
	{
		seedPosition2[0] = TSeedx3;
		seedPosition2[1] = TSeedy3;
		initialDistance2 = 5;
	}

	 

	NodeType node2;
	const double seedValue2 = - initialDistance2;
	node2.SetValue( seedValue2 );
	node2.SetIndex( seedPosition2 );

	seeds2->Initialize();
	seeds2->InsertElement( 0, node2 );

	fastMarching2->SetTrialPoints( seeds2 );
	fastMarching2->SetSpeedConstant( 1.0 );
	fastMarching2->SetOutputSize(
           extractfilter->GetOutput()->GetBufferedRegion().GetSize() );
	fastMarching2->SetOutputOrigin(extractfilter->GetOutput()->GetOrigin());
	fastMarching2->SetOutputDirection(extractfilter->GetOutput()->GetDirection());
	fastMarching2->SetOutputSpacing(extractfilter->GetOutput()->GetSpacing());

	// Geodesic Active Contour Filter objects
	GeodesicActiveContourFilterType::Pointer geodesicActiveContour1 =
                                       GeodesicActiveContourFilterType::New();  // For femur (intensity)
	GeodesicActiveContourFilterType::Pointer geodesicActiveContour2 =
                                       GeodesicActiveContourFilterType::New();  // For femur (variance)
	GeodesicActiveContourFilterType::Pointer geodesicActiveContour3 =
                                       GeodesicActiveContourFilterType::New();  // For tibia (intensity)
	GeodesicActiveContourFilterType::Pointer geodesicActiveContour4 =
                                       GeodesicActiveContourFilterType::New();  // For tibia (variance)

	
	// Structuring element to dilate the segmentation output from variance
	StructuringElementType structuringElement;
	structuringElement.SetRadius(11);                // Radius of the dilation structuring element
	structuringElement.CreateStructuringElement();
 
	// Objects to dilate the segmentation output from variance 
	BinaryDilateImageFilterType::Pointer dilateFilter1     // For Femur 
          = BinaryDilateImageFilterType::New();
	dilateFilter1->SetKernel(structuringElement);
	dilateFilter1->SetDilateValue(1);
	BinaryDilateImageFilterType::Pointer dilateFilter2     // For Tibia
          = BinaryDilateImageFilterType::New();
	dilateFilter2->SetKernel(structuringElement);
	dilateFilter2->SetDilateValue(1);
  
	
	// Parameters for different modules - Set optimally but might need to be adjusted for images in different intensity ranges
	// For simplicity same parameters have been used for intensity and variance images, femur and tibia

	sigmoid1->SetAlpha(-3.0);      // Alpha and beta of the sigmoid filter transformation      
	sigmoid1->SetBeta(15.0);  
	sigmoid1->SetOutputMinimum(0.0);
	sigmoid1->SetOutputMaximum(1.0);

	smoothing1->SetTimeStep( 0.02 );    
	smoothing1->SetNumberOfIterations(  5 );
	smoothing1->SetConductanceParameter( 9.0 );

	gradientMagnitude1->SetSigma( 0.01 );      // Sigma of Gaussian used in Gradient Magnitude Recursive Gaussian   

	sigmoid2->SetAlpha(alpha);                 // use parameters for variance defined earlier
	sigmoid2->SetBeta(beta);  
	sigmoid2->SetOutputMinimum(0.0);
	sigmoid2->SetOutputMaximum(1.0);

	smoothing2->SetTimeStep( 0.02 );
	smoothing2->SetNumberOfIterations(  5 );
	smoothing2->SetConductanceParameter( 9.0 );

	gradientMagnitude2->SetSigma( 0.01 );


	geodesicActiveContour1->SetPropagationScaling( 10.0 );    // Propagation parameter 
	geodesicActiveContour1->SetCurvatureScaling( 1.0 );       // Curvation parameter
	geodesicActiveContour1->SetAdvectionScaling( 1.0 );       // Advection parameter
	geodesicActiveContour1->SetMaximumRMSError( 0.0002 );     // Maximum RMS Error allowed
	geodesicActiveContour1->SetNumberOfIterations( 400 );     // Maximum number of iterations   

	geodesicActiveContour2->SetPropagationScaling( 10.0 );
	geodesicActiveContour2->SetCurvatureScaling( 1.0 );
	geodesicActiveContour2->SetAdvectionScaling( 1.0 );
	geodesicActiveContour2->SetMaximumRMSError( 0.0002 );          
	geodesicActiveContour2->SetNumberOfIterations( 400 );        

	geodesicActiveContour3->SetPropagationScaling( 10.0 );
	geodesicActiveContour3->SetCurvatureScaling( 1.0 );
	geodesicActiveContour3->SetAdvectionScaling( 1.0 );
	geodesicActiveContour3->SetMaximumRMSError( 0.0002 );          
	geodesicActiveContour3->SetNumberOfIterations( 400 );        


	geodesicActiveContour4->SetPropagationScaling( 10.0 );
	geodesicActiveContour4->SetCurvatureScaling( 1.0 );
	geodesicActiveContour4->SetAdvectionScaling( 1.0 );
	geodesicActiveContour4->SetMaximumRMSError( 0.0002 );          
	geodesicActiveContour4->SetNumberOfIterations( 400 );        
	
	         

	// Filter to compute variance
	VarianceFilterType::Pointer varfilter = VarianceFilterType::New();
	InternalImageType::Pointer varianceImage = InternalImageType::New();
	varfilter->SetInputImage(image);
	varfilter->SetNeighborhoodRadius( 5 );    // Neighborhood of variance calculation
	computeVarianceImage (varfilter,image,varianceImage);  // Function to compute variance


	// Set up the segmentation flow 
	// For Intensity Image 
	smoothing1->SetInput( image );
	gradientMagnitude1->SetInput( smoothing1->GetOutput() );
	sigmoid1->SetInput(gradientMagnitude1->GetOutput() );
	thresholder1->SetInput( geodesicActiveContour1->GetOutput() );
	thresholder3->SetInput( geodesicActiveContour3->GetOutput() );

	// For Variance Image
	smoothing2->SetInput( varianceImage );
	gradientMagnitude2->SetInput( smoothing2->GetOutput() );
	sigmoid2->SetInput(gradientMagnitude2->GetOutput() );
	thresholder2->SetInput( geodesicActiveContour2->GetOutput() );
	dilateFilter1->SetInput(thresholder2->GetOutput());
	thresholder4->SetInput( geodesicActiveContour4->GetOutput() );
	dilateFilter2->SetInput(thresholder4->GetOutput());
  
	// For Intensity Image (Femur)
	geodesicActiveContour1->SetInput(  fastMarching1->GetOutput() );                    
	geodesicActiveContour1->SetFeatureImage( sigmoid1->GetOutput() );    

    // For Variance Image (Femur)
	geodesicActiveContour2->SetInput(  fastMarching1->GetOutput() );                    
	geodesicActiveContour2->SetFeatureImage( sigmoid2->GetOutput() );  

	// For Intensity Image (Tibia)
	geodesicActiveContour3->SetInput(  fastMarching2->GetOutput() );                    
	geodesicActiveContour3->SetFeatureImage( sigmoid1->GetOutput() );  

	// For Intensity Image (Femur)
	geodesicActiveContour4->SetInput( fastMarching2->GetOutput() );                    
	geodesicActiveContour4->SetFeatureImage( sigmoid2->GetOutput() );


	BinaryImageType::Pointer image1 = thresholder1->GetOutput();   // Segmentation output of Intensity (femur)
	BinaryImageType::Pointer image2 = dilateFilter1->GetOutput();  // Segmentation output of Variance (femur)
	BinaryImageType::Pointer image3 = thresholder3->GetOutput();   // Segmentation output of Intensity (tibia)
	BinaryImageType::Pointer image4 = dilateFilter2->GetOutput();  // Segmentation output of Variance (tibia)

    // Object for Logical AND
	AndImageFilterType::Pointer andFilter1
          = AndImageFilterType::New();       // For Femur
	andFilter1->SetInput(0, image1);
	andFilter1->SetInput(1, image2);

	AndImageFilterType::Pointer andFilter2
          = AndImageFilterType::New();       // For Tibia
	andFilter2->SetInput(0, image3);
	andFilter2->SetInput(1, image4);

	MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
	multiplyFilter->SetInput(andFilter2->GetOutput());
	multiplyFilter->SetConstant(3);         // Set intensity for Tibia to 3

	// Object to combine femur and tibia segmented output
	AddImageFilterType::Pointer addFilter
    = AddImageFilterType::New ();

	// Cast from binary to float
	CastFilterType::Pointer castFilter = CastFilterType::New();
  
	try
    {   // Decide based on slice number
		if ((z<FStartSlice1-1) && (z>=TStartSlice1-1))
		{
			castFilter->SetInput(multiplyFilter->GetOutput());     // only tibia
			castFilter->Update();
			joinFilter->SetInput(z,castFilter->GetOutput());
		}
		else if ((z>=FStartSlice1-1) && (z<TStartSlice1-1)) 
		{
			castFilter->SetInput(andFilter1->GetOutput());        // only femur
			castFilter->Update();
			joinFilter->SetInput(z,castFilter->GetOutput());
		}
		else if ((z>FEndSlice3-1) && (z<=TEndSlice3-1))
		{
			castFilter->SetInput(multiplyFilter->GetOutput());    // only tibia
			castFilter->Update();
			joinFilter->SetInput(z,castFilter->GetOutput());
		}
		else if ((z<=FEndSlice3-1) && (z>TEndSlice3-1))
		{
			castFilter->SetInput(andFilter1->GetOutput());        // only femur
			castFilter->Update();
			joinFilter->SetInput(z,castFilter->GetOutput());
		}
		else
		{
			addFilter->SetInput1(andFilter1->GetOutput());
			addFilter->SetInput2(multiplyFilter->GetOutput());
			addFilter->Update();                                    // combine both femur and tibia
			joinFilter->SetInput(z,addFilter->GetOutput());
		}
	
	}
	catch( itk::ExceptionObject & excep )
    {
		std::cerr << "Exception caught !" << std::endl;
		std::cerr << excep << std::endl;
		return EXIT_FAILURE;
    }


  }

  try
  {
	joinFilter->Update();
	writer->SetInput(joinFilter->GetOutput());           // write to output 3D volume
	writer->Update();
  }
  catch( itk::ExceptionObject & excep )
  {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
  }
  
  clock.Stop();
  std::cout << "Time Elapsed: " << clock.GetTotal() << std::endl;         // Display elapsed time
  return EXIT_SUCCESS;
}