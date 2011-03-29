itk_module_test()
set(ITK-LocalDiscreteHessianTests
  itkDiscreteHessianHeaderTest.cxx
)

CreateTestDriver(ITK-LocalDiscreteHessian  "${ITK-LocalDiscreteHessian-Test_LIBRARIES}" "${ITK-LocalDiscreteHessianTests}")

add_test(NAME itkLocalDiscreteHessianHeaderTest
      COMMAND ITK-LocalDiscreteHessianTestDriver itkDiscreteHessianHeaderTest)
