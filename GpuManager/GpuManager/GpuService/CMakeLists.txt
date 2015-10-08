gaudi_subdir(GpuService v0r1)

gaudi_depends_on_subdirs(
  GpuManager/GpuIpc)

gaudi_add_module(
  GpuService
  src/*.cpp
  INCLUDE_DIRS GpuIpc
  LINK_LIBRARIES GpuIpc GaudiKernel)

gaudi_install_headers(GpuService)
