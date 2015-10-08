gaudi_subdir(PrPixelCuda v0r1)

gaudi_depends_on_subdirs(
  GpuManager/GpuHandler)


set(CUDA_PROPAGATE_HOST_FLAGS OFF)

gaudi_add_cuda_module(
  PrPixelCuda
  src/*.cpp #src/*.cu
  INCLUDE_DIRS GpuHandler CUDA
  LINK_LIBRARIES GpuHandler CUDA)
