
include_directories(${CUDA_INCLUDE_DIRS})
include_directories(${THRUST_INCLUDE_DIR})
include_directories(${SRC}) #include for higher source files include

MESSAGE("SRC in srcedge: " ${SRC})
MESSAGE("SRC in CUDA_INCLUDE_DIRS: " ${CUDA_INCLUDE_DIRS})

add_library(edgeswapLib STATIC
	Edgeswap_test.cu
	Edgeswap_test.h)


target_compile_features(edgeswapLib PUBLIC cxx_std_11)