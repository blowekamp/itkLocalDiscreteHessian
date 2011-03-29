#ifndef __itkHessianImageFilter_h
#define __itkHessianImageFilter_h

#include "itkImage.h"
#include "itkImageToImageFilter.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkGradientRecursiveGaussianImageFilter.h"

namespace itk
{
namespace Local
{

/**
 * \class t
 *
 * \todo fix me
 */
template <typename TInputImage,
          typename TOutputImage = Image< SymmetricSecondRankTensor<
            ITK_TYPENAME NumericTraits< ITK_TYPENAME TInputImage::PixelType >::RealType,
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

  void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread, int threadId);


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
