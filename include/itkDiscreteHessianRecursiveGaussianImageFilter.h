#ifndef __itkDiscreteHessianRecursiveGaussianImageFilter_h
#define __itkDiscreteHessianRecursiveGaussianImageFilter_h

#include "itkRecursiveGaussianImageFilter.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkPixelTraits.h"


namespace itk
{
namespace Local
{


/**
 * \class DiscreteHessianRecursiveGaussianImageFilter
 * \brief Preforms an approximation of convolution of the image with a
 * Gaussian kernel by a recursive algorithm then performs central
 * differences.
 *
 * This filter is a composite filter of the
 * RecursiveGaussianImageFilter and the HessianImageFilter.
 *
 * The additional feature added is the normalization across scale.
 *
 * \ingroup GradientFilters
 * \ingroup ITKDiscreteHessian
 */
template< typename TInputImage,
          typename TOutputImage = Image< SymmetricSecondRankTensor<
                                           typename NumericTraits< typename TInputImage::PixelType >::RealType,
                                           ::itk::GetImageDimension< TInputImage >::ImageDimension >,
                                         ::itk::GetImageDimension< TInputImage >::ImageDimension > >
class ITK_EXPORT DiscreteHessianRecursiveGaussianImageFilter:
  public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Standard class typedefs. */
  typedef DiscreteHessianRecursiveGaussianImageFilter     Self;
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
  itkTypeMacro(DiscreteHessianRecursiveGaussianImageFilter, ImageToImageFilter);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Set Sigma value. Sigma is measured in the units of image spacing.  */
  itkSetMacro( Sigma, double );
  itkGetConstMacro( Sigma, double );

  /** Define which normalization factor will be used for the Gaussian */
  itkSetMacro( NormalizeAcrossScale, bool );
  itkGetConstMacro( NormalizeAcrossScale, bool );
  itkBooleanMacro( NormalizeAcrossScale );

  /** DiscreteHessianRecursiveGaussianImageFilter needs all of the input to produce an
   * output. Therefore, DiscreteHessianRecursiveGaussianImageFilter needs to provide
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

  DiscreteHessianRecursiveGaussianImageFilter();
  virtual ~DiscreteHessianRecursiveGaussianImageFilter() {}
  void PrintSelf(std::ostream & os, Indent indent) const;

  /** Generate Data */
  void GenerateData(void);

private:

  DiscreteHessianRecursiveGaussianImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);                      //purposely not implemented

  // special binary functor to perform A+B*ConstValue
  //
  class MultHessianComponentwiseConstFunctor
  {
  public:
    typedef MultHessianComponentwiseConstFunctor Self;

    MultHessianComponentwiseConstFunctor( void ) : m_Value( NumericTraits<OutputPixelType>::One ) {}

    bool operator!=( const Self &other ) const { return !(*this==other); }
    bool operator==( const Self &other ) const { return m_Value == other.m_Value; }

    inline OutputPixelType operator()( const OutputPixelType &a ) const
    {
      OutputPixelType o;
      for ( unsigned int i = 0; i < m_Value.GetNumberOfComponents(); ++i )
        {
        o.SetNthComponent( i, a.GetNthComponent(i) * m_Value.GetNthComponent(i) );
        }
      return o;
    }

    OutputPixelType m_Value;
  };


  /** Normalize the image across scale space */
  bool m_NormalizeAcrossScale;
  double m_Sigma;
};


} // end namespace Local
} // end namespace itk

#include "itkDiscreteHessianRecursiveGaussianImageFilter.txx"

#endif // __itkDiscreteHessianRecursiveGaussianImageFilter_h
