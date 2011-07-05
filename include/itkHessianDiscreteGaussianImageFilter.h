#ifndef __itkHessianDiscreteGaussianImageFilter_h
#define __itkHessianDiscreteGaussianImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkPixelTraits.h"


namespace itk
{
namespace Local
{


/**
 * \class HessianDiscreteGaussianImageFilter
 * \brief Convolves the image with a discrete Gaussian kernel then
 * performs central differences.
 *
 * This filter is a composite filter of the
 * DiscreteGaussianImageFilter and the HessianImageFilter.
 *
 * The additional feature added is the normalization across scale.
 *
 * THIS FILTER IS NOT IMPLEMENTED WITH DISCRETE GAUSSIAN FILTER YET
 *
 * \ingroup GradientFilters
 * \ingroup Streamed
 */
template< typename TInputImage,
          typename TOutputImage = Image< SymmetricSecondRankTensor<
                                           ITK_TYPENAME NumericTraits< ITK_TYPENAME TInputImage::PixelType >::RealType,
                                           ::itk::GetImageDimension< TInputImage >::ImageDimension >,
                                         ::itk::GetImageDimension< TInputImage >::ImageDimension > >
class ITK_EXPORT HessianDiscreteGaussianImageFilter:
  public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef HessianDiscreteGaussianImageFilter     Self;
  typedef ImageToImageFilter< TInputImage, TOutputImage > Superclass;
  typedef SmartPointer< Self >                            Pointer;
  typedef SmartPointer< const Self >                      ConstPointer;

  /** Pixel Type of the input image */
  typedef TInputImage                                   InputImageType;
  typedef typename TInputImage::PixelType               PixelType;
  typedef typename NumericTraits< PixelType >::RealType RealType;

  /** Image dimension. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      ::itk::GetImageDimension< TInputImage >::ImageDimension);

  /** Type of the output Image */
  typedef TOutputImage                                       OutputImageType;
  typedef typename          OutputImageType::PixelType       OutputPixelType;
  typedef typename PixelTraits< OutputPixelType >::ValueType OutputComponentType;

  /** Run-time type information (and related methods).   */
  itkTypeMacro(HessianDiscreteGaussianImageFilter, ImageToImageFilter);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Set Sigma value. Sigma is measured in the units of image spacing.  */
  itkSetMacro( Sigma, double );
  itkGetConstMacro( Sigma, double );

  /** Define if scale-space normalization factor will be used
   *
   *  \sa  RecursiveGaussianImageFilter::SetNormalizeAcrossScale
   */
  itkSetMacro( NormalizeAcrossScale, bool );
  itkGetConstMacro( NormalizeAcrossScale, bool );
  itkBooleanMacro( NormalizeAcrossScale );

  /** HessianDiscreteGaussianImageFilter needs all of the input to produce an
   * output. Therefore, HessianDiscreteGaussianImageFilter needs to provide
   * an implementation for GenerateInputRequestedRegion in order to inform
   * the pipeline execution model.
   * \sa ImageToImageFilter::GenerateInputRequestedRegion() */
  virtual void GenerateInputRequestedRegion()
  throw( InvalidRequestedRegionError );

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro( InputHasNumericTraitsCheck,
                   ( Concept::HasNumericTraits< PixelType > ) );
  itkConceptMacro( OutputHasPixelTraitsCheck,
                   ( Concept::HasPixelTraits< OutputPixelType > ) );
  /** End concept checking */
#endif
protected:

  HessianDiscreteGaussianImageFilter();
  virtual ~HessianDiscreteGaussianImageFilter() {}
  void PrintSelf(std::ostream & os, Indent indent) const;

  /** Generate Data */
  void GenerateData(void);

private:

  HessianDiscreteGaussianImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);                      //purposely not implemented

  /** Normalize the image across scale space */
  bool m_NormalizeAcrossScale;
  double m_Sigma;
};


} // end namespace Local
} // end namespace itk

#include "itkHessianDiscreteGaussianImageFilter.txx"

#endif // __itkHessianDiscreteGaussianImageFilter_h
