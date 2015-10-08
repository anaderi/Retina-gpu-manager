gaudi_subdir(GpuDriver v0r1)

gaudi_depends_on_subdirs(
  GpuManager/GpuIpc)

find_package(Boost COMPONENTS filesystem program_options system thread)

gaudi_add_executable(
  gpudriver
  src/*.cpp
  INCLUDE_DIRS Boost GpuIpc GpuHandler
  LINK_LIBRARIES Boost GpuIpc GpuHandler GaudiPluginService)
