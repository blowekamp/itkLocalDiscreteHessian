itk_module_test()
set(ITKLocalDiscreteHessianTests
  itkDiscreteHessianHeaderTest.cxx
  itkHessianImageFilterTest.cxx
  itkDiscreteHessianRecursiveGaussianImageFilterTest.cxx
)

CreateTestDriver(ITKLocalDiscreteHessian  "${ITKLocalDiscreteHessian-Test_LIBRARIES}" "${ITKLocalDiscreteHessianTests}")

add_test(NAME itkLocalDiscreteHessianHeaderTest
      COMMAND ITKLocalDiscreteHessianTestDriver itkDiscreteHessianHeaderTest)

add_test(NAME itkLocalHessianImageFilterTest
      COMMAND ITKLocalDiscreteHessianTestDriver itkHessianImageFilterTest )

add_test(NAME itkLocalDiscreteHessianRecursiveGaussianImageFilterTest
      COMMAND ITKLocalDiscreteHessianTestDriver itkDiscreteHessianRecursiveGaussianImageFilterTest )
