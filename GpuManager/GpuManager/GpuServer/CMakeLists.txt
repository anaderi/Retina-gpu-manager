gaudi_subdir(GpuServer v0r1)

gaudi_depends_on_subdirs(
  GpuManager/GpuIpc
  GpuManager/GpuHandler)

find_package(Boost COMPONENTS filesystem program_options system thread)

gaudi_add_executable(
  gpuserver
  src/*.cpp
  INCLUDE_DIRS Boost GpuIpc GpuHandler
  LINK_LIBRARIES Boost GpuIpc GpuHandler GaudiPluginService)
