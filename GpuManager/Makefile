build_tool=cmake
ifeq ($(wildcard build.conf),build.conf)
include build.conf
endif

ifeq ($(build_tool),cmake)
######## CMake case ########

# custom flags
CMAKEFLAGS="-DGAUDI_USE_EXE_SUFFIX=OFF"

include $(LBUTILSROOT)/data/Makefile-cmake.mk

else
######### CMT case #########

CMTPROJECTPATH := $(shell bash -c ". build_env.sh ; printenv CMTPROJECTPATH")

all:
%:
	@env CMTPROJECTPATH="$(CMTPROJECTPATH)" $(MAKE) -f $(LBCONFIGURATIONROOT)/data/Makefile-cmt.mk $*

endif
