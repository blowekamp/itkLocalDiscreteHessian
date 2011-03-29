#ifndef __itkHessianImageFilter_txx
#define __itkHessianImageFilter_txx


#include "itkHessianImageFilter.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkNeighborhoodAlgorithm.h"

#include "itkProgressReporter.h"
#include "itkProgressAccumulator.h"

#include "vnl/vnl_trace.h"

namespace itk
{
namespace Local
{

/**
 *  Constructor
 */

template <typename TInputImage, typename TOutputImage >
HessianImageFilter<TInputImage,TOutputImage>
::HessianImageFilter( void )
{
}

/**
 * Enlarge Input Requested Region
 */
template< class TInputImage, class TOutputImage >
void
HessianImageFilter< TInputImage, TOutputImage >
::GenerateInputRequestedRegion()
throw ( InvalidRequestedRegionError )
{
  // call the superclass' implementation of this method. this should
  // copy the output requested region to the input requested region
  Superclass::GenerateInputRequestedRegion();

  // get pointers to the input and output
  typename InputImageType::Pointer inputPtr =
    const_cast< TInputImage * >( this->GetInput() );

  if ( !inputPtr )
    {
    return;
    }


  // the hessaion just needs a 1 radius neighborhood
  const unsigned int radius = 1;

  // get a copy of the input requested region (should equal the output
  // requested region)
  typename TInputImage::RegionType inputRequestedRegion;
  inputRequestedRegion = inputPtr->GetRequestedRegion();

  // pad the input requested region by the operator radius
  inputRequestedRegion.PadByRadius( radius );

  // crop the input requested region at the input's largest possible region
  if ( inputRequestedRegion.Crop( inputPtr->GetLargestPossibleRegion() ) )
    {
    inputPtr->SetRequestedRegion(inputRequestedRegion);
    return;
    }
  else
    {
    // Couldn't crop the region (requested region is outside the largest
    // possible region).  Throw an exception.

    // store what we tried to request (prior to trying to crop)
    inputPtr->SetRequestedRegion(inputRequestedRegion);

    // build an exception
    InvalidRequestedRegionError e(__FILE__, __LINE__);
    e.SetLocation(ITK_LOCATION);
    e.SetDescription("Requested region is (at least partially) outside the largest possible region.");
    e.SetDataObject(inputPtr);
    throw e;
    }
}

/**
 * Threaded Data Generation
 */
template <typename TInputImage, typename TOutputImage >
void
HessianImageFilter<TInputImage,TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                       int threadId)
{
  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels(), 100, 0.5, 0.5);

  const TInputImage *input = this->GetInput();

  const unsigned int ImageDimension = TInputImage::ImageDimension;

  TOutputImage *output = this->GetOutput();


  typedef TOutputImage                        OutputImageType;
  typedef typename OutputImageType::PixelType HessianType;
  ImageRegionIterator<OutputImageType> oit;

  itk::Size<ImageDimension> radius;
  radius.Fill( 1 );
  unsigned long center;
  unsigned long stride[ImageDimension];

  typename TInputImage::SpacingType spacing = input->GetSpacing();


  // compute the boundary faces of our region
  typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< TInputImage >::FaceListType faceList;
  NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< TInputImage > bC;
  faceList = bC( input, outputRegionForThread, radius );

  typename NeighborhoodAlgorithm::ImageBoundaryFacesCalculator< TInputImage >::FaceListType::iterator fit;

  typedef ConstNeighborhoodIterator< TInputImage > NeighborhoodType;

  // get center and dimension strides for iterator neighborhoods
  NeighborhoodType it( radius, input, *faceList.begin() );
  center = it.Size()/2;
  for ( unsigned int i = 0; i < ImageDimension; ++i )
    {
    stride[i] = it.GetStride(i);
    }

  // process each of the region "faces"
  for ( fit = faceList.begin(); fit != faceList.end(); ++fit )
    {
    // set up the iterator for the "face" and let the automatic
    // boundary condition detection work as needed
    it = NeighborhoodType( radius, input, *fit);

    oit = ImageRegionIterator<OutputImageType>( output, *fit );

    while ( !it.IsAtEnd() )
      {
      // symetric hessian
      HessianType H;


      //Calculate 2nd order derivative on the diaganal
      for ( unsigned int i = 0; i < ImageDimension; i++ )
        {
        H(i,i) = it.GetPixel(center + stride[i]) + it.GetPixel(center - stride[i])
          - 2.0 * it.GetPixel(center);
        H(i,i) /= spacing[i] * spacing[i];
        }

      //Calculate the 2nd derivatives
      for ( unsigned int i = 0; i < ImageDimension - 1; i++ )
        {
        for ( unsigned int j = i + 1; j < ImageDimension; j++ )
          {
          H(i,j) = it.GetPixel(center - stride[i] - stride[j])
            - it.GetPixel(center - stride[i] + stride[j])
            - it.GetPixel(center + stride[i] - stride[j])
            + it.GetPixel(center + stride[i] + stride[j]);
          H(i,j) /= 4.0 * spacing[i] * spacing[j];
          }
        }

      oit.Set( H );

      ++oit;
      ++it;
       progress.CompletedPixel();
      }
    }

}


} // end namespace Local
} // end namespace itk

#endif // __itkHessianImageFilter_txx
