find_package(CUDA REQUIRED)

macro(gaudi_add_cuda_module library)
  # Process both CUDA and Gaudi arguments

  cmake_parse_arguments(
    ARG
    ""
    "GENCONF_PRELOAD;GENCONF_USER_MODULE"
    "LINK_LIBRARIES;INCLUDE_DIRS"
    ${ARGN})

  cuda_add_cuda_include_once()

  cuda_get_sources_and_options(srcs cmake_options options ${ARG_UNPARSED_ARGUMENTS})

  gaudi_common_add_build(
    ${srcs}
    INCLUDE_DIRS   ${ARG_INCLUDE_DIRS}
    LINK_LIBRARIES ${ARG_LINK_LIBRARIES})

  # Generate CUDA object files

  cuda_wrap_srcs(
    ${library} OBJ
    generated_files
    ${srcs}
    ${cmake_options}
    SHARED
    OPTIONS ${options} -std=c++11)

  # Add a library including CUDA's object files

  add_library(
    ${library} MODULE ${cmake_options}
    ${generated_files}
    ${srcs})

  # Configure as a Gaudi module

  target_link_libraries(${library} GaudiPluginService ${ARG_LINK_LIBRARIES})
  _gaudi_detach_debinfo(${library})

  gaudi_generate_componentslist(${library})
  set(ARG_GENCONF)
  foreach(genconf_sub_opt PRELOAD USER_MODULE)
    if(ARG_GENCONF_${genconf_sub_opt})
      set(ARG_GENCONF ${ARG_GENCONF} ${genconf_sub_opt} ${ARG_GENCONF_${genconf_sub_opt}})
    endif()
  endforeach()
  gaudi_generate_configurables(${library} ${ARG_GENCONF})

  set_property(GLOBAL APPEND PROPERTY COMPONENT_LIBRARIES ${library})

  gaudi_add_genheader_dependencies(${library})

  #----Installation details-------------------------------------------------------
  install(TARGETS ${library} LIBRARY DESTINATION lib OPTIONAL)
  gaudi_export(MODULE ${library})
endmacro(gaudi_add_cuda_module library)
