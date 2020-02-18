
include_directories(${CUDA_INCLUDE_DIRS})
include_directories(${THRUST_INCLUDE_DIR})
include_directories(${SRC}) #include for higher source files include

MESSAGE("SRC in srcsprings: " ${SRC})
MESSAGE("SRC in CUDA_INCLUDE_DIRS: " ${CUDA_INCLUDE_DIRS})

add_library(springsLib STATIC
	AreaTriangles.cu
	AreaTriangles.h
	BendingTriangles.cu
	BendingTriangles.h
	MemRepulsionSprings.cu
	MemRepulsionSprings.h
	LinearSprings.cu
	LinearSprings.h
	LJSprings.cu
	LJSprings.h
	LJSprings_LJ.cu
	LJSprings_LJ.h
	VolumeComp.cu
	VolumeComp.h
	VolumeSprings.cu
	VolumeSprings.h
	LineTensionSprings.cu
	LineTensionSprings.h
	NodeAdvance.cu
	NodeAdvance.h
	AreaTrianglesEnergy.cu
	AreaTrianglesEnergy.h
	BendingTrianglesEnergy.cu
	BendingTrianglesEnergy.h
	LinearSpringsEnergy.cu
	LinearSpringsEnergy.h
	MemRepulsionEnergy.cu
	MemRepulsionEnergy.h)


target_compile_features(springsLib PUBLIC cxx_std_11)