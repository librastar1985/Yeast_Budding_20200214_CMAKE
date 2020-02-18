#ifndef AREATRIANGLESENERGY_H_
#define AREATRIANGLESENERGY_H_ 

struct GeneralParams;
struct CoordInfoVecs;
struct AreaTriangleInfoVecs;
struct BendingTriangleInfoVecs;
struct LinearSpringInfoVecs;

struct HostSetInfoVecs;
struct LinearSpringInfoVecs;


double ComputeAreaTrianglesEnergy(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    AreaTriangleInfoVecs& areaTriangleInfoVecs);


struct AreaEnergyFunctor {
	double area_0;
    double spring_constant;
    double* locXAddr;
    double* locYAddr;
    double* locZAddr;

    int* idKey;
    
	__host__ __device__ AreaEnergyFunctor(
        double& _area_0,
        double& _spring_constant,
        double* _locXAddr,
        double* _locYAddr,
        double* _locZAddr,
        int* _idKey):
        area_0(_area_0),
        spring_constant(_spring_constant),
        locXAddr(_locXAddr),
        locYAddr(_locYAddr),
        locZAddr(_locZAddr),
        idKey(_idKey) {}

	//hand in counting iterator and id of triangle
	__device__ double operator()(const Tuuuu &u4) {
        //test placing the ids of the nodes and then get positions. 
		int counter = thrust::get<0>(u4);
		int place = 3 * counter;//represents location in write to vector.

        int id_i = thrust::get<1>(u4);
        int id_j = thrust::get<2>(u4);
        int id_k = thrust::get<3>(u4);

		
		CVec3 ri = thrust::make_tuple<double>(locXAddr[id_i], locYAddr[id_i], locZAddr[id_i]);
		CVec3 rj = thrust::make_tuple<double>(locXAddr[id_j], locYAddr[id_j], locZAddr[id_j]);
		CVec3 rk = thrust::make_tuple<double>(locXAddr[id_k], locYAddr[id_k], locZAddr[id_k]);

		CVec3 rkj = CVec3_subtract(rk, rj);
		CVec3 rij = CVec3_subtract(ri, rj);

        double area_current = sqrt( CVec3_dot( CVec3_cross(rkj, rij), CVec3_cross(rkj, rij) ) )/2;

        double energy =  spring_constant/2 * (area_current - area_0) * (area_current - area_0) / area_0;
        return energy;

    };
};
#endif //AREATRIANGLESENERGY_H_