#ifndef CAPSIDESPRINGSAll_H_
#define CAPSIDESPRINGSAll_H_ 

struct GeneralParams;
struct CoordInfoVecs;
struct AreaTriangleInfoVecs;
struct BendingTriangleInfoVecs;
struct LinearSpringInfoVecs;

struct HostSetInfoVecs;
struct LinearSpringInfoVecs;

void ComputeCapsideSpringsAll(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    CapsidInfoVecs& capsidInfoVecs,
    AuxVecs& auxVecs);
    
struct CapsidSpringFunctorAll {
    int factor;
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

    int* capsidIdUnreduced;
    int* membraneIdUnreduced;
    double* forceXAddr;
    double* forceYAddr;
    double* forceZAddr;
    
	__host__ __device__ 
    CapsidSpringFunctorAll(  
    int& _factor,
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

    int* _capsidIdUnreduced,
    int* _membraneIdUnreduced,
    double* _forceXAddr,
    double* _forceYAddr,
    double* _forceZAddr ) :
    factor(_factor),
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
    
    capsidIdUnreduced(_capsidIdUnreduced),
    membraneIdUnreduced(_membraneIdUnreduced),
    forceXAddr(_forceXAddr),
    forceYAddr(_forceYAddr),
    forceZAddr(_forceZAddr) {}

    //each capside point chooses a single membrane point.
	__device__
    void operator() (const Tuuu& u3) {

        int counter = thrust::get<0>(u3);
        int place = factor * counter;//force and id writing location to unreduced vectors

		int bucketId = thrust::get<1>(u3);//bucket containing nodeId
		int capsidId = thrust::get<2>(u3);//node to attempt link from.
	

        double xLoc_LR;
        double yLoc_LR;
        double zLoc_LR;

        double forceX = 0.0;
        double forceY = 0.0;
        double forceZ = 0.0;


        //iterate through membrane id's and choose closest one under cutoff length
       // for (int i = beginIndex; i < endIndex; i++ ) {
        for (int memId = 0; memId < membraneMaxNode; memId++) {
            //int memId = bucketNbrsExp[i];
            if (memId < membraneMaxNode) {
                xLoc_LR = membraneNodeXAddr[memId] - capsidNodeXAddr[capsidId];
                yLoc_LR = membraneNodeYAddr[memId] - capsidNodeYAddr[capsidId];
                zLoc_LR = membraneNodeZAddr[memId] - capsidNodeZAddr[capsidId];

        	    double length_current = sqrt( (xLoc_LR) * (xLoc_LR) + 
                            (yLoc_LR) * (yLoc_LR)  + 
                            (zLoc_LR) * (zLoc_LR) );
                
                if ( (length_current < length_cutoff) ) {
                    //then we can apply force
                    double magnitude = -spring_constant * (length_current - length_zero);
                    forceX = magnitude * (xLoc_LR/length_current);
                    forceY = magnitude * (yLoc_LR/length_current);
                    forceZ = magnitude * (zLoc_LR/length_current);

                    //store
                        membraneIdUnreduced[place] = memId;
                        capsidIdUnreduced[place] = capsidId;
                        forceXAddr[place] = forceX;
                        forceYAddr[place] = forceY;
                        forceZAddr[place] = forceZ;
                        
                        //after saving data, write to next location.  
                        place+=1;
                    
                } 
                
            }
        }

    }
};



#endif