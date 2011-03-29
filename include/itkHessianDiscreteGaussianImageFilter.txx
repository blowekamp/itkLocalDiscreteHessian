#ifndef __itkHessianDiscreteGaussianImageFilter_txx
#define __itkHessianDiscreteGaussianImageFilter_txx

#include "itkHessianDiscreteGaussianImageFilter.h"
#include "itkHessianImageFilter.h"
namespace itk
{
namespace Local
{

/**
 * Constructor
 */
template< typename TInputImage, typename TOutputImage >
HessianDiscreteGaussianImageFilter< TInputImage, TOutputImage >
::HessianDiscreteGaussianImageFilter()
{
  m_NormalizeAcrossScale = false;
  m_Sigma = 1.0;
}

template< typename TInputImage, typename TOutputImage >
void
HessianDiscreteGaussianImageFilter< TInputImage, TOutputImage >
::GenerateInputRequestedRegion()
  throw( InvalidRequestedRegionError )
{
  // call the superclass' implementation of this method. this should
  // copy the output requested region to the input requested region
  Superclass::GenerateInputRequestedRegion();

  // This filter needs all of the input
  typename InputImageType::Pointer image = const_cast< InputImageType * >( this->GetInput() );

  if ( image )
    {
    image->SetRequestedRegion( this->GetInput()->GetLargestPossibleRegion() );
    }
}

template< typename TInputImage, typename TOutputImage >
void
HessianDiscreteGaussianImageFilter< TInputImage, TOutputImage >
::GenerateData(void)
{
  itkDebugMacro(<< "HessianDiscreteGaussianImageFilter generating data ");

  // Create a process accumulator for tracking the progress of this
  // minipipeline
  ProgressAccumulator::Pointer progress = ProgressAccumulator::New();
  progress = ProgressAccumulator::New();
  progress->SetMiniPipelineFilter(this);


  // Get the input and output pointers
  typename InputImageType::ConstPointer  input = this->GetInput();
  typename OutputImageType::Pointer      output = this->GetOutput();


  typedef itk::Image< double, ImageDimension> RealImageType;
  typedef itk::RecursiveGaussianImageFilter< InputImageType, RealImageType > FirstGaussianFilterType;
  typedef itk::RecursiveGaussianImageFilter< RealImageType, RealImageType > RealGaussianFilterType;

  // The first Gaussian filter in the mini pipeline
  //
  // Convert the input to the real pixel type,
  // Do not perform operation in-place as not to steal input data
  typename FirstGaussianFilterType::Pointer firstGaussian = FirstGaussianFilterType::New();
  firstGaussian->SetOrder( FirstGaussianFilterType::ZeroOrder );
  firstGaussian->SetNormalizeAcrossScale( false );
  firstGaussian->SetSigma( this->m_Sigma );
  firstGaussian->SetDirection( ImageDimension - 1 );
  firstGaussian->InPlaceOff();
  firstGaussian->ReleaseDataFlagOn();
  firstGaussian->SetInput( input );

  progress->RegisterInternalFilter( firstGaussian,  1.0/(ImageDimension+1) );


  // Assemble remaining gaussian filters
  // All are inplace and real image to real image
  typename RealGaussianFilterType::Pointer gaussianFilters[ImageDimension - 1];
  for( unsigned int i = 0; i < ImageDimension-1; ++i )
    {
    gaussianFilters[i] = RealGaussianFilterType::New();
    if ( i == 0 )
      {
      gaussianFilters[i]->SetInput( firstGaussian->GetOutput() );
      }
    else
      {
      gaussianFilters[i]->SetInput( gaussianFilters[i-1]->GetOutput() );
      }
    gaussianFilters[i]->SetOrder( RealGaussianFilterType::ZeroOrder );
    gaussianFilters[i]->SetNormalizeAcrossScale( false );
    gaussianFilters[i]->SetSigma( this->m_Sigma );
    gaussianFilters[i]->SetDirection( ImageDimension-i-2 );
    gaussianFilters[i]->InPlaceOn();
    gaussianFilters[i]->ReleaseDataFlagOn();

    progress->RegisterInternalFilter( gaussianFilters[i],  1.0/(ImageDimension+1) );
    }


  typedef itk::Local::HessianImageFilter<RealImageType, OutputImageType> HessianImageFilterType;
  typename HessianImageFilterType::Pointer hessian = HessianImageFilterType::New();
  hessian->SetInput( gaussianFilters[ ImageDimension - 2 ]->GetOutput() );

  progress->RegisterInternalFilter( hessian, 1.0/(ImageDimension+1) );

  // Perform standard graft-update-graft
  //
  // This saves memory, and copies the requested region to the filter
  // so that hessian filter updates the correct region
  hessian->GraftOutput( output );
  hessian->Update();
  this->GraftOutput( hessian->GetOutput() );
}

template< typename TInputImage, typename TOutputImage >
void
HessianDiscreteGaussianImageFilter< TInputImage, TOutputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << "NormalizeAcrossScale: " << m_NormalizeAcrossScale << std::endl;
  os << "Sigma: " << m_Sigma << std::endl;
}

} // end namespace Local
} // end namespace itk

#endif //  __itkHessianDiscreteGaussianImageFilter_txx
