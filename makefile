all: clusty


####################

ifdef MSVC     # Avoid the MingW/Cygwin sections
    uname_S := Windows
else                          # If uname not available => 'not'
    uname_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
	uname_M := $(shell sh -c 'uname -m 2>/dev/null || echo not')
endif

ifeq ($(uname_S),Linux)   
	# check if CPU supports AVX2
#	HAVE_AVX2 = $(filter-out 0,$(shell grep avx2 /proc/cpuinfo | wc -l))
#	OMP_FLAGS = -fopenmp
	ABI_FLAGS = -fabi-version=6
	MIMALLOC_OBJ=libs/mimalloc/mimalloc.o
endif
ifeq ($(uname_S),Darwin)
	 # check if CPU supports SSE4.2
#	HAVE_AVX2 = $(filter-out 0,$(shell  sysctl -n machdep.cpu.features machdep.cpu.leaf7_features| grep AVX2 - | wc -l))
#	OMP_FLAGS = -Xpreprocessor -fopenmp 
	ABI_FLAGS = 
	MIMALLOC_OBJ=
endif


ifeq ($(PLATFORM), arm8)
$(info *** ARMv8 with NEON extensions ***)
    ARCH_FLAGS := -march=armv8-a  -DARCH_ARM
else ifeq ($(PLATFORM), m1)
$(info *** Apple M1(or never) with NEON extensions ***)
    ARCH_FLAGS := -march=armv8.4-a  -DARCH_ARM
else ifeq ($(PLATFORM), sse2)
$(info *** x86-64 with SSE2 extensions ***)
    ARCH_FLAGS := -msse2 -m64 -DARCH_X64 
else ifeq ($(PLATFORM), avx)
$(info *** x86-64 with AVX extensions ***)
    ARCH_FLAGS := -mavx -m64  -DARCH_X64
else ifeq ($(PLATFORM), avx2)
$(info *** x86-64 with AVX2 extensions ***)
    ARCH_FLAGS := -mavx2 -m64  -DARCH_X64
else
$(info *** Unspecified platform - use native compilation)
    ifeq ($(uname_M),x86_64)
        ARCH_FLAGS := -march=native -DARCH_X64
    else
        ARCH_FLAGS := -march=native -DARCH_ARM
    endif	
endif


GIT_COMMIT = $(shell git describe --always --dirty) 

#####################
ROOT_DIR = .
MAIN_DIR = src
INCLUDES = -I libs/mimalloc/include
DEFINE_FLAGS := -DGIT_COMMIT=$(GIT_COMMIT)

ifeq ($(LEIDEN), true) 
    INCLUDES += -I ./libs/igraph/include -I ./libs/igraph/build/include
    LIB_IGRAPH = ./libs/igraph/build/src/libigraph.a
else
    DEFINE_FLAGS += -DNO_LEIDEN
endif


ifeq ($(DYNAMIC_LINK), true)
    CFLAGS	= -Wall -O3 $(ARCH_FLAGS) -std=c++17 $(DEFINE_FLAGS) $(INCLUDES) -pthread 
    CLINK	= -lm -O3 -std=c++17 -pthread $(ABI_FLAGS) 
else
ifeq ($(uname_S),Darwin)
    CFLAGS  = -Wall -O3 $(ARCH_FLAGS) -std=c++17 $(DEFINE_FLAGS) $(INCLUDES)
    CLINK	= -lm -O3 -std=c++17 $(ABI_FLAGS) -static-libgcc
else
    CFLAGS  = -Wall -O3 $(ARCH_FLAGS) -std=c++17 $(DEFINE_FLAGS) $(INCLUDES) -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
    CLINK	= -lm -static -O3 -std=c++17 $(ABI_FLAGS) -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
endif  
endif



CMAKE_OSX_SYSROOT_FLAG =
ifeq ($(uname_S),Darwin)
	SDK_PATH := $(shell $(CXX) -v 2>&1 | grep -- '--with-sysroot' | sed -E 's/.*--with-sysroot=([^ ]+).*/\1/')
	CMAKE_OSX_SYSROOT_FLAG := -DCMAKE_OSX_SYSROOT=$(SDK_PATH)
endif

ifeq ($(LEIDEN), true)  
igraph:
	mkdir libs/igraph/build 
	cmake $(CMAKE_OSX_SYSROOT_FLAG) -DIEEE754_DOUBLE_ENDIANNESS_MATCHES=TRUE -DCMAKE_C_COMPILER=$(CC) -DCMAKE_CXX_COMPILER=$(CXX) -S libs/igraph -B libs/igraph/build 
	cmake $(CMAKE_OSX_SYSROOT_FLAG) -DIEEE754_DOUBLE_ENDIANNESS_MATCHES=TRUE -DCMAKE_C_COMPILER=$(CC) -DCMAKE_CXX_COMPILER=$(CXX) -S libs/igraph -B libs/igraph/build  
	cmake --build libs/igraph/build
else
igraph:

endif


$(MIMALLOC_OBJ):
	$(CC) -DMI_MALLOC_OVERRIDE -O3 -DNDEBUG -fPIC -Wall -Wextra -Wno-unknown-pragmas -fvisibility=hidden -Wstrict-prototypes -ftls-model=initial-exec -fno-builtin-malloc -std=gnu11 -c -I libs/mimalloc/include libs/mimalloc/src/static.c -o $(MIMALLOC_OBJ)


  
OBJS := \
	$(MIMALLOC_OBJ) \
	$(MAIN_DIR)/console.o \
	$(MAIN_DIR)/conversion.o \
	$(MAIN_DIR)/graph.o \
	$(MAIN_DIR)/log.o \
	$(MAIN_DIR)/main.o \
	$(MAIN_DIR)/params.o \

%.o: %.cpp igraph
	$(CXX) $(CFLAGS) -c $< -o $@

clusty: $(OBJS)
	$(CXX) $(CLINK) $(LDFLAGS) -o $(ROOT_DIR)/$@ $(OBJS) $(LIB_IGRAPH)

clean:
	-rm $(MAIN_DIR)/*.o
	-rm $(MIMALLOC_OBJ)
	-rm clusty
	-rm -r libs/igraph/build
	
