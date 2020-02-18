#ifndef CAPSIDESPRINGS_H_
#define CAPSIDESPRINGS_H_ 

struct GeneralParams;
struct CoordInfoVecs;
struct AreaTriangleInfoVecs;
struct BendingTriangleInfoVecs;
struct LinearSpringInfoVecs;

struct HostSetInfoVecs;
struct LinearSpringInfoVecs;


void ComputeCapsideSprings(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    CapsidInfoVecs& capsidInfoVecs,
    AuxVecs& auxVecs);

void AdvanceCapsideParticles(
    GeneralParams& generalParams,
    CapsidInfoVecs& capsidInfoVecs);
    
struct CapsidSpringFunctor {
    double length_cutoff;
    double spring_constant;
    double length_zero;
    int capsidMaxNode;
    int membraneMaxNode;
    
    double* capsidNodeXAddr;
    double* capsidNodeYAddr;
    double* capsidNodeZAddr;
    double* membraneNodeXAddr;
    double* membraneNodeYAddr;
    double* membraneNodeZAddr;

    int* bucketNbrsExp;
    int* keyBegin;
    int* keyEnd;
    
	__host__ __device__ 
    CapsidSpringFunctor(  
    double& _length_cutoff,  
    double& _spring_constant,
    double& _length_zero,
    int& _capsidMaxNode,
    int& _membraneMaxNode,
    
    double* _capsidNodeXAddr,
    double* _capsidNodeYAddr,
    double* _capsidNodeZAddr,
    double* _membraneNodeXAddr,
    double* _membraneNodeYAddr,
    double* _membraneNodeZAddr,

    int* _bucketNbrsExp,
    int* _keyBegin,
    int* _keyEnd) :
    length_cutoff(_length_cutoff),
    spring_constant(_spring_constant),
    length_zero(_length_zero),
    capsidMaxNode(_capsidMaxNode),
    membraneMaxNode(_membraneMaxNode),

    capsidNodeXAddr(_capsidNodeXAddr),
    capsidNodeYAddr(_capsidNodeYAddr),
    capsidNodeZAddr(_capsidNodeZAddr),
    membraneNodeXAddr(_membraneNodeXAddr),
    membraneNodeYAddr(_membraneNodeYAddr),
    membraneNodeZAddr(_membraneNodeZAddr),

    bucketNbrsExp(_bucketNbrsExp),
    keyBegin(_keyBegin),
    keyEnd(_keyEnd) {}

    //each capside point chooses a single membrane point.
	__device__
    U2CVec4 operator() (const Tuu& u2) {
		int bucketId = thrust::get<0>(u2);//bucket containing nodeId
		int capsidId = thrust::get<1>(u2);//node to attempt link from.
	
		//beginning and end of attempted attachment id's in bucketNbrsExp
        //these indices are membrane
		int beginIndex = keyBegin[bucketId];
		int endIndex = keyEnd[bucketId];

        double xLoc_LR;
        double yLoc_LR;
        double zLoc_LR;

        double length_current = 100.0;//default higher than possible lengths
        int memIdFinalChoice = membraneMaxNode+1;//begin id as higher than possible id's

        //iterate through membrane id's and choose closest one under cutoff length
       // for (int i = beginIndex; i < endIndex; i++ ) {
        for (int memId = 0; memId < membraneMaxNode; memId++) {
            //int memId = bucketNbrsExp[i];
            if (memId < membraneMaxNode) {
                xLoc_LR = membraneNodeXAddr[memId] - capsidNodeXAddr[capsidId];
                yLoc_LR = membraneNodeYAddr[memId] - capsidNodeYAddr[capsidId];
                zLoc_LR = membraneNodeZAddr[memId] - capsidNodeZAddr[capsidId];

        	    double length_temp = sqrt( 
                                        (xLoc_LR) * (xLoc_LR) + 
                                        (yLoc_LR) * (yLoc_LR)  + 
                                        (zLoc_LR) * (zLoc_LR) );
                
                //only choose id, if length is lower than cutoff and minimal.
                if ((length_temp < length_current) && (length_temp < length_cutoff) ) {
                    length_current = length_temp;
                    memIdFinalChoice = memId;
                } 
                
            }
        }
        //now capsidId has found a membrane node to bind to.
        //Compute forces
		
        double forceX = 0.0;
        double forceY = 0.0;
        double forceZ = 0.0;

        //only apply force if id that is found is not default site.
        if (memIdFinalChoice < membraneMaxNode) {

            //reset direction vectors and recalculate current length. 
            xLoc_LR = membraneNodeXAddr[memIdFinalChoice] - capsidNodeXAddr[capsidId];
            yLoc_LR = membraneNodeYAddr[memIdFinalChoice] - capsidNodeYAddr[capsidId];
            zLoc_LR = membraneNodeZAddr[memIdFinalChoice] - capsidNodeZAddr[capsidId];

            length_current = sqrt( 
                                (xLoc_LR) * (xLoc_LR) + 
                                (yLoc_LR) * (yLoc_LR) + 
                                (zLoc_LR) * (zLoc_LR) );

            double magnitude = -spring_constant * (length_current - length_zero);
        
            forceX = magnitude * ( xLoc_LR / length_current );
            forceY = magnitude * ( yLoc_LR / length_current );
            forceZ = magnitude * ( zLoc_LR / length_current );
        }

        //uu,dddd
        return thrust::make_tuple( memIdFinalChoice, capsidId, length_current, forceX, forceY, forceZ );
    }
};

struct SaxpyFunctorCapsid {
	double dt;
    double viscosity;

    double forceX;
    double forceY;
    double forceZ;

	__host__ __device__
		//
		SaxpyFunctorCapsid(
			double& _dt, 
			double& _viscosity,
			double& _forceX,
			double& _forceY,
			double& _forceZ) :
		dt(_dt),
        viscosity(_viscosity),
		forceX(_forceX),
		forceY(_forceY),
		forceZ(_forceZ) {}

	__device__
		CVec3 operator()(const CVec3 &p3) {

			double xLocNew = thrust::get<0>(p3) + (dt / viscosity)* forceX;
			double yLocNew = thrust::get<1>(p3) + (dt / viscosity)* forceY;
			double zLocNew = thrust::get<2>(p3) + (dt / viscosity)* forceZ;

			return thrust::make_tuple(xLocNew, yLocNew, zLocNew);
	    }                                 

};


#endif