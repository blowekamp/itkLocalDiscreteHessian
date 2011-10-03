#ifndef __itkHessianImageFilter_h
#define __itkHessianImageFilter_h

#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include "itkSymmetricSecondRankTensor.h"

namespace itk
{
namespace Local
{

/**
 * \class HessianImageFilter
 * \brief Computes the Hessian matrix by central differences
 *
 * \ingroup GradientFilters
 * \ingroup Streamed
 * \ingroup ITK-DiscreteHessian
 */
template <typename TInputImage,
          typename TOutputImage = Image< SymmetricSecondRankTensor<
            typename NumericTraits< typename TInputImage::PixelType >::RealType,
            ::itk::GetImageDimension< TInputImage >::ImageDimension >,
                                         ::itk::GetImageDimension< TInputImage >::ImageDimension > >
class ITK_EXPORT HessianImageFilter :
    public ImageToImageFilter< TInputImage, TOutputImage>
{
public:

  /** Standard typedefs */
  typedef HessianImageFilter                           Self;
  typedef ImageToImageFilter<TInputImage,TOutputImage> Superclass;
  typedef SmartPointer<Self>                           Pointer;
  typedef SmartPointer<const Self>                     ConstPointer;


  /** Pixel Type of the input image */
  typedef TInputImage                                    InputImageType;
  typedef typename InputImageType::PixelType             PixelType;

  /** Type of the output Image */
  typedef TOutputImage                                      OutputImageType;
  typedef typename          OutputImageType::PixelType      OutputPixelType;
  typedef typename OutputImageType::RegionType              OutputImageRegionType;


 /** Run-time type information (and related methods).   */
  itkTypeMacro( HessianImageFilter, ImageToImageFilter );

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  virtual void GenerateInputRequestedRegion()
    throw( InvalidRequestedRegionError );


#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputHasNumericTraitsCheck,
                  (Concept::HasNumericTraits<PixelType>));
  itkConceptMacro(OutputHasPixelTraitsCheck,
                  (Concept::HasPixelTraits<OutputPixelType>));
  /** End concept checking */
#endif

protected:

  HessianImageFilter( void );

  void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, ThreadIdType threadId);


private:

  HessianImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};

} // end namepace Local
} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkHessianImageFilter.txx"
#endif


#endif //__itkHessianImageFilter_h
