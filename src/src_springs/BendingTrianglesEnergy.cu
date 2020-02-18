
#include "System.h"
#include "SystemStructures.h"
#include "BendingTrianglesEnergy.h"

double ComputeCosTrianglesEnergy(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    BendingTriangleInfoVecs& bendingTriangleInfoVecs) {
	
    
    thrust::counting_iterator<int> elemId(0); 

	//bendingTriangleInfoVecs.initial_angle = 1.5707963267/2.0;
	/*thrust::fill(bendingTriangleInfoVecs.tempNodeForceXReduced.begin(),bendingTriangleInfoVecs.tempNodeForceXReduced.end(),0.0);
	thrust::fill(bendingTriangleInfoVecs.tempNodeForceYReduced.begin(),bendingTriangleInfoVecs.tempNodeForceYReduced.end(),0.0);
	thrust::fill(bendingTriangleInfoVecs.tempNodeForceZReduced.begin(),bendingTriangleInfoVecs.tempNodeForceZReduced.end(),0.0);
	thrust::fill(bendingTriangleInfoVecs.tempNodeForceXUnreduced.begin(),bendingTriangleInfoVecs.tempNodeForceXUnreduced.end(),0.0);
	thrust::fill(bendingTriangleInfoVecs.tempNodeForceYUnreduced.begin(),bendingTriangleInfoVecs.tempNodeForceYUnreduced.end(),0.0);
	thrust::fill(bendingTriangleInfoVecs.tempNodeForceZUnreduced.begin(),bendingTriangleInfoVecs.tempNodeForceZUnreduced.end(),0.0);*/

    //apply force to temporary vectors.
    bendingTriangleInfoVecs.bending_triangle_energy= 
    thrust::transform_reduce(
        thrust::make_zip_iterator(
            thrust::make_tuple(
				elemId,
                coordInfoVecs.edges2Triangles_1.begin(),
                coordInfoVecs.edges2Triangles_2.begin(),
                coordInfoVecs.edges2Nodes_1.begin(),
                coordInfoVecs.edges2Nodes_2.begin())),
        thrust::make_zip_iterator(
            thrust::make_tuple(
				elemId, 
                coordInfoVecs.edges2Triangles_1.begin(),
                coordInfoVecs.edges2Triangles_2.begin(),
                coordInfoVecs.edges2Nodes_1.begin(),
                coordInfoVecs.edges2Nodes_2.begin())) + coordInfoVecs.num_edges,
        CosBendingEnergyFunctor(
            bendingTriangleInfoVecs.spring_constant,
            bendingTriangleInfoVecs.initial_angle,        
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()),
            
            thrust::raw_pointer_cast(bendingTriangleInfoVecs.tempNodeIdUnreduced.data()),
            //thrust::raw_pointer_cast(bendingTriangleInfoVecs.tempNodeForceXUnreduced.data()),
            //thrust::raw_pointer_cast(bendingTriangleInfoVecs.tempNodeForceYUnreduced.data()),
            //thrust::raw_pointer_cast(bendingTriangleInfoVecs.tempNodeForceZUnreduced.data()),
            thrust::raw_pointer_cast(coordInfoVecs.triangles2Nodes_1.data()), 
            thrust::raw_pointer_cast(coordInfoVecs.triangles2Nodes_2.data()), 
            thrust::raw_pointer_cast(coordInfoVecs.triangles2Nodes_3.data())),
		0.0, thrust::plus<double>() );
	
    return bendingTriangleInfoVecs.bending_triangle_energy;
};
