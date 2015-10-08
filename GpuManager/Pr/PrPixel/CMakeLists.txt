gaudi_subdir(PrPixel v2r1)

gaudi_depends_on_subdirs(
  GpuManager/GpuService)

gaudi_add_module(
  PrPixel
  src/*.cpp
  INCLUDE_DIRS AIDA GSL Pr/PrKernel GpuService
  LINK_LIBRARIES VPDetLib DAQEventLib TrackEvent GaudiAlgLib GaudiKernel)
