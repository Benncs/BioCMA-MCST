.PHONY: all clean external_module

BUILD_DIR=./build
INCLUDE_DIR=./../../apps/libs/biomodels/ext/ 
LIB_DIR=/home/benjamin/Documents/code/cpp/DynamicLibrary/public/

INC_MC=./../../apps/libs/mc/public

all: dir main

dir:
	mkdir -p ${BUILD_DIR}

main:
	g++ -shared -fopenmp -fvisibility=hidden -fPIC ./impl.cpp -o ${BUILD_DIR}/impl.so -I${INCLUDE_DIR} -I${LIB_DIR} -I${INC_MC}

clean:
	rm -rf ${BUILD_DIR}
