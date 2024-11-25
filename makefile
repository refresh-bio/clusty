all: clusty

# *** REFRESH makefile utils
include refresh.mk

$(call INIT_SUBMODULES)
$(call INIT_GLOBALS)
$(call CHECK_OS_ARCH, $(PLATFORM))

# *** Project directories
$(call SET_SRC_OBJ_BIN,src,obj,bin)
3RD_PARTY_DIR := ./libs

# *** Project configuration
$(call ADD_MIMALLOC, $(3RD_PARTY_DIR)/mimalloc)
#$(call ADD_REFRESH_LIB, $(3RD_PARTY_DIR))
$(call SET_STATIC, $(STATIC_LINK))
$(call SET_C_CPP_STANDARDS, c11, c++17)
$(call SET_GIT_COMMIT)

ifeq ($(LEIDEN),true) 
$(call ADD_IGRAPH, $(3RD_PARTY_DIR)/igraph)
else
DEFINE_FLAGS += -DNO_LEIDEN
endif

$(call SET_FLAGS, $(TYPE))

$(call SET_COMPILER_VERSION_ALLOWED, GCC, Linux_x86_64, 10, 20)
$(call SET_COMPILER_VERSION_ALLOWED, GCC, Linux_aarch64, 11, 20)
$(call SET_COMPILER_VERSION_ALLOWED, GCC, Darwin_x86_64, 11, 13)
$(call SET_COMPILER_VERSION_ALLOWED, GCC, Darwin_arm64, 11, 13)

ifneq ($(MAKECMDGOALS),clean)
$(call CHECK_COMPILER_VERSION)
endif

# *** Source files and rules
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,MAIN,.))

# *** Targets
clusty: $(OUT_BIN_DIR)/clusty
$(OUT_BIN_DIR)/clusty: mimalloc_obj \
	$(OBJ_MAIN)
	-mkdir -p $(OUT_BIN_DIR)	
	$(CXX) -o $@  \
	$(MIMALLOC_OBJ) \
	$(OBJ_MAIN) \
	$(LIBRARY_FILES) $(LINKER_FLAGS) $(LINKER_DIRS)

# *** Cleaning
.PHONY: clean init
clean: clean-zlib-ng clean-isa-l clean-mimalloc_obj clean-igraph
	-rm -r $(OBJ_DIR)
	-rm -r $(OUT_BIN_DIR)

init:
	$(call INIT_SUBMODULES)

