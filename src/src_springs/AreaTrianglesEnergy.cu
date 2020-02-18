#include "System.h"
#include "SystemStructures.h"
#include "AreaTrianglesEnergy.h"

double ComputeAreaTrianglesEnergy(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    AreaTriangleInfoVecs& areaTriangleInfoVecs) {

        thrust::counting_iterator<int> elemId(0);
    
        
        areaTriangleInfoVecs.area_triangle_energy = thrust::transform_reduce( 
			thrust::make_zip_iterator(
				thrust::make_tuple(
                    elemId,
					coordInfoVecs.triangles2Nodes_1.begin(),
					coordInfoVecs.triangles2Nodes_2.begin(),
					coordInfoVecs.triangles2Nodes_3.begin())),
			thrust::make_zip_iterator(
				thrust::make_tuple(
                    elemId,
					coordInfoVecs.triangles2Nodes_1.begin(),
					coordInfoVecs.triangles2Nodes_2.begin(),
					coordInfoVecs.triangles2Nodes_3.begin())) + coordInfoVecs.num_triangles,
            AreaEnergyFunctor( 
                areaTriangleInfoVecs.initial_area,
                areaTriangleInfoVecs.spring_constant,
                thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
                thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
                thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()),
				thrust::raw_pointer_cast(areaTriangleInfoVecs.tempNodeIdUnreduced.data())),
            0.0, thrust::plus<double>() );
            
        return areaTriangleInfoVecs.area_triangle_energy;
};