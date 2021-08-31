 
src_dir := src
include_dir := include 
obj_dir := build
bin_dir := build/bin 

TGT := ${bin_dir}/SFGL3D
SRCS := $(wildcard ${src_dir}/*.cpp)
INCLUDE := -I$(wildcard ${include_dir}/*.h)
OBJ := $(patsubst SRCS:${src_dir}/%.cpp, ${obj_dir}/%.o, ${SRCS})

CXXFLAGS = -Wall -std=c++14 
#-O2 -fopenmp -fpermissive

LDFLAGS += -Lusr/lib64/
LDFLAGS += -lGL
LDFLAGS += -lglfw
LDFLAGS += -lGLEW

# Targets 
.PHONY: all clean

all: ${TGT}

${TGT}: ${OBJ} | ${bin_dir}
	${CXX} ${CXXFLAGS} ${SRCS} ${INCLUDE} -o sfgl3d ${LDFLAGS}
	#${CXX} ${CXXFLAGS} ${OBJ} ${INCLUDE} -o sfgl3d ${LDFLAGS}
clean: 
	rm -r -f ${obj_dir} *.o 
