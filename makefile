# Inital MakeFile for building StableFluids3D-GL On Linux. 
# Tested on Centos 7.8 

src_dir := src
include_dir := include 
# Check if explcitlly Need Devtoolset and OpenMP Include location + GLEW and GLFW Vendor Headers (within /usr)
obj_dir := build
bin_dir := build/bin 

TGT := ${bin_dir}/SFGL3D
SRCS := $(wildcard ${src_dir}/*.cpp)
INCLUDE := -I$(wildcard ${include_dir}/*.h) # I include dir pre-pending each header.h. 

#OBJ := $(SRCS:${src_dir}/%.cpp=${obj_dir}/%.o)
# Replace all .cpp file with .o name for obj targets using patsubst.
OBJ := $(patsubst SRCS:${src_dir}/%.cpp, ${obj_dir}/%.o, ${SRCS})

# Compiler and Linker Flags. 
CXXFLAGS = -Wall -std=c++14 
#-O2 -fopenmp -fpermissive
# Assuming -L Dir is /usr/lib64 can just specify lib name -l**
LDFLAGS += -Lusr/lib64/
LDFLAGS += -lGL
LDFLAGS += -lglfw
LDFLAGS += -lGLEW

# Phony Targets
.PHONY: all clean

# Build All
all: ${TGT}

# Build Target 
${TGT}: ${OBJ} | ${bin_dir}
	${CXX} ${CXXFLAGS} ${SRCS} ${INCLUDE} -o sfgl3d ${LDFLAGS}
	#${CXX} ${CXXFLAGS} ${OBJ} ${INCLUDE} -o sfgl3d ${LDFLAGS}

# Build Clean
clean: 
	rm -r -f ${obj_dir} *.o 
