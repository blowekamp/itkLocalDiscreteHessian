itk_module_test()
set(ITK-LocalDiscreteHessianTests
  itkDiscreteHessianHeaderTest.cxx
  itkHessianImageFilterTest.cxx
)

CreateTestDriver(ITK-LocalDiscreteHessian  "${ITK-LocalDiscreteHessian-Test_LIBRARIES}" "${ITK-LocalDiscreteHessianTests}")

add_test(NAME itkLocalDiscreteHessianHeaderTest
      COMMAND ITK-LocalDiscreteHessianTestDriver itkDiscreteHessianHeaderTest)


add_test(NAME itkLocalHessianImageFilterTest
      COMMAND ITK-LocalDiscreteHessianTestDriver itkHessianImageFilterTest )
