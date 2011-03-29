#include "itkHessianImageFilter.h"
#include "itkGaussianImageSource.h"

int itkHessianImageFilterTest( int argc, char *argv[] )
{
  const unsigned int Dimension = 3;
  typedef itk::Image< float, Dimension > ImageType;

  const unsigned int imageSize = 64;

  ImageType::SizeType size;
  size.Fill( imageSize );

  ImageType::SpacingType spacing;
  spacing.Fill( 1.0 );
  spacing[0] = 1.0;

  typedef itk::GaussianImageSource<ImageType> GaussianSourceType;
  GaussianSourceType::Pointer gaussianSource = GaussianSourceType::New();
  gaussianSource->SetSize( size );
  gaussianSource->SetSpacing( spacing );
  gaussianSource->SetMean( itk::FixedArray< double, Dimension>( imageSize/2 ) );
  gaussianSource->SetSigma( itk::FixedArray< double, Dimension>( 10.0 ) );
  gaussianSource->SetNormalized( false );
  gaussianSource->SetScale( 1.0 ); // dark blob


  typedef itk::Local::HessianImageFilter< ImageType > HessianFilterType;
  HessianFilterType::Pointer hessian = HessianFilterType::New();
  hessian->SetInput( gaussianSource->GetOutput() );
  hessian->Update();

  ImageType::IndexType idx;
  idx.Fill( imageSize/2 );

  std::cout << hessian->GetOutput()->GetPixel( idx );

  return 0;
}
