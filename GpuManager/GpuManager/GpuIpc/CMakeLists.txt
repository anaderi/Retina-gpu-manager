gaudi_subdir(GpuIpc v0r1)

find_package(Boost COMPONENTS system thread)

gaudi_add_library(
  GpuIpc
  src/*.cpp
  INCLUDE_DIRS Boost GpuIpc
  PUBLIC_HEADERS GpuIpc
  LINK_LIBRARIES Boost)
