
#include "System.h"
#include "SystemStructures.h"
#include "CapsidSprings.h"


void ComputeCapsideSprings(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    CapsidInfoVecs& capsidInfoVecs,
    AuxVecs& auxVecs) {

    thrust::fill(capsidInfoVecs.tempNodeForceX.begin(),capsidInfoVecs.tempNodeForceX.end(),0.0);
    thrust::fill(capsidInfoVecs.tempNodeForceY.begin(),capsidInfoVecs.tempNodeForceY.end(),0.0);
    thrust::fill(capsidInfoVecs.tempNodeForceZ.begin(),capsidInfoVecs.tempNodeForceZ.end(),0.0);

    thrust::fill(capsidInfoVecs.tempLengthsPairs.begin(),capsidInfoVecs.tempLengthsPairs.end(),100.0);

    thrust::fill(capsidInfoVecs.tempMembraneId.begin(),capsidInfoVecs.tempMembraneId.end(), generalParams.maxNodeCount + 1);
    thrust::fill(capsidInfoVecs.tempCapsideId.begin(),capsidInfoVecs.tempCapsideId.end(), capsidInfoVecs.maxNodeCount + 1);
    
    //for each capsid node, choose the closest membrane node. 
    //what if I chose the closest capsid node for each membrane node?
    thrust::transform(  
        thrust::make_zip_iterator(
            thrust::make_tuple(
                capsidInfoVecs.bucketKeys.begin(),
                capsidInfoVecs.bucketValues.begin())),
        
        thrust::make_zip_iterator(
            thrust::make_tuple(
                capsidInfoVecs.bucketKeys.begin(),
                capsidInfoVecs.bucketValues.begin())) + capsidInfoVecs.maxNodeCount,

        thrust::make_zip_iterator(//output.
            thrust::make_tuple(//WARNING: Membrande id must come first, followed by capside id.
                capsidInfoVecs.tempMembraneId.begin(),
                capsidInfoVecs.tempCapsideId.begin(),
                capsidInfoVecs.tempLengthsPairs.begin(),
                capsidInfoVecs.tempNodeForceX.begin(),
                capsidInfoVecs.tempNodeForceY.begin(),
                capsidInfoVecs.tempNodeForceZ.begin() )),
                  
            CapsidSpringFunctor(
                capsidInfoVecs.length_cutoff,
                capsidInfoVecs.spring_constant, 
                capsidInfoVecs.length_zero, 
                capsidInfoVecs.maxNodeCount,
                generalParams.maxNodeCount,

                thrust::raw_pointer_cast(capsidInfoVecs.nodeLocX.data()),
                thrust::raw_pointer_cast(capsidInfoVecs.nodeLocY.data()),
                thrust::raw_pointer_cast(capsidInfoVecs.nodeLocZ.data()),
                thrust::raw_pointer_cast(coordInfoVecs.nodeLocX.data()),
                thrust::raw_pointer_cast(coordInfoVecs.nodeLocY.data()),
                thrust::raw_pointer_cast(coordInfoVecs.nodeLocZ.data()),

                thrust::raw_pointer_cast(auxVecs.bucketValuesIncludingNeighbor.data()),
                thrust::raw_pointer_cast(auxVecs.keyBegin.data()),
                thrust::raw_pointer_cast(auxVecs.keyEnd.data())) );
                
  /*  std::cout<<"here"<<std::flush;
    int bucket = capsidInfoVecs.bucketKeys[35];
    int keyBegin = auxVecs.keyBegin[bucket];
    int keyEnd = auxVecs.keyEnd[bucket];
    
    std::cout<<"bucket35: "<< bucket<<std::endl;
    for (int i = keyBegin; i < keyEnd; i++) { 
        
        int memId = auxVecs.bucketValuesIncludingNeighbor[i];
        double xLoc_LR = capsidInfoVecs.nodeLocX[35] - coordInfoVecs.nodeLocX[memId];
        double yLoc_LR = capsidInfoVecs.nodeLocY[35] - coordInfoVecs.nodeLocY[memId];
        double zLoc_LR = capsidInfoVecs.nodeLocZ[35] - coordInfoVecs.nodeLocZ[memId];

         double dist = sqrt( (xLoc_LR) * (xLoc_LR) + 
         (yLoc_LR) * (yLoc_LR)  + 
         (zLoc_LR) * (zLoc_LR) );

         std::cout<<"capsid 35 neighbor: "<< memId <<std::endl;
         std::cout<<"capsid 35 dist: "<< dist <<std::endl;
         
           
    }*/

    //now we have un reduced forces. Sort by id and reduce. 
    //key, then value. Each vector returns sorted		
    //since each capsid point chose as membrane point, to acheieve 1-1 ratio we first sort by membrane node,
    //then reduce by unique membrane id.
    
    //shuffle indices for permutation. This must be redone each time step
    std::vector<int> hostIndices;

    hostIndices.reserve(capsidInfoVecs.tempMembraneId.size());
    for (int i = 0; i < capsidInfoVecs.tempMembraneId.size(); i++)
        hostIndices.push_back(i);
    
    std::random_device seed;
    std::mt19937 rng(seed()); 

    std::shuffle(hostIndices.begin(), hostIndices.end(), rng);


    thrust::device_vector<int> indexes = hostIndices;

    //randomly permute items first so that multiple id's have a chance of binding.
   /* thrust::sort_by_key(
        thrust::make_permutation_iterator(capsidInfoVecs.tempMembraneId.begin(), indexes.begin() ),
        thrust::make_permutation_iterator(capsidInfoVecs.tempMembraneId.begin(), indexes.end()) ,
        
        thrust::make_zip_iterator(
                thrust::make_tuple(
                    thrust::make_permutation_iterator(
                        capsidInfoVecs.tempCapsideId.begin(),
                        indexes.begin()),
                    thrust::make_permutation_iterator(
                        capsidInfoVecs.tempLengthsPairs.begin(),
                        indexes.begin()),
                    thrust::make_permutation_iterator(
                        capsidInfoVecs.tempNodeForceX.begin(),
                        indexes.begin()),
                    thrust::make_permutation_iterator(
                        capsidInfoVecs.tempNodeForceY.begin(),
                        indexes.begin()),
                    thrust::make_permutation_iterator(
                        capsidInfoVecs.tempNodeForceZ.begin(),
                        indexes.begin()) )), 
                thrust::less<int>());
*/
    //capsid id's are already unique, so we need to order mem_id
    //first sort by lengths to take minimal distance, 
    //comment this out if you don't want to sort via distances
    thrust::sort_by_key(capsidInfoVecs.tempLengthsPairs.begin(), capsidInfoVecs.tempLengthsPairs.end(),
        thrust::make_zip_iterator(
            thrust::make_tuple(
                capsidInfoVecs.tempMembraneId.begin(),
                capsidInfoVecs.tempCapsideId.begin(),
                capsidInfoVecs.tempNodeForceX.begin(),
                capsidInfoVecs.tempNodeForceY.begin(),
                capsidInfoVecs.tempNodeForceZ.begin() )), thrust::less<double>());

    thrust::stable_sort_by_key(capsidInfoVecs.tempMembraneId.begin(), capsidInfoVecs.tempMembraneId.end(),
        thrust::make_zip_iterator(
            thrust::make_tuple(
                capsidInfoVecs.tempLengthsPairs.begin(),
                capsidInfoVecs.tempCapsideId.begin(),
                capsidInfoVecs.tempNodeForceX.begin(),
                capsidInfoVecs.tempNodeForceY.begin(),
                capsidInfoVecs.tempNodeForceZ.begin() )), thrust::less<int>());//now sorting int
    
    //now given that the membrane id's are sorted, we sort them by ascending lengths while 
    //keeping them sorted.
    
   //if you want to take minimal distance between capsid and membrane, use a stable sort to keep relative order
   //on id's so that distance is preserved. 
   /*
if (generalParams.iteration % 100 == 0) {
    for (int i = 0; i < capsidInfoVecs.tempLengthsPairs.size(); i++){
        std::cout<< "post stable sort " <<" cap Id: "<< capsidInfoVecs.tempCapsideId[i]<< " << mem Id: "<< capsidInfoVecs.tempMembraneId[i]<< " "<< capsidInfoVecs.tempNodeForceX[i] << " "<<capsidInfoVecs.tempNodeForceY[i] << " "<<capsidInfoVecs.tempNodeForceZ[i] << std::endl;
        std::cout<<"dist: "<< capsidInfoVecs.tempLengthsPairs[i]<< std::endl;
    }
}*/
    int endKey = thrust::get<0>(
        thrust::unique_by_key(
            //input key
            capsidInfoVecs.tempMembraneId.begin(), capsidInfoVecs.tempMembraneId.end(),
        thrust::make_zip_iterator(//input val
            thrust::make_tuple( 
                capsidInfoVecs.tempCapsideId.begin(),
                capsidInfoVecs.tempLengthsPairs.begin(),
                capsidInfoVecs.tempNodeForceX.begin(),
                capsidInfoVecs.tempNodeForceY.begin(),
                capsidInfoVecs.tempNodeForceZ.begin())),
        thrust::equal_to<int>() )) - capsidInfoVecs.tempMembraneId.begin();//binary_pred 
    
    capsidInfoVecs.num_connections=endKey-1;//endkey includes the last bad id which is larger than the maxnodecount

    //apply reduced force to all nodes. 
    
    UCVec3 init(0, 0.0,0.0,0.0);
   
/*
if (generalParams.iteration % 100 == 0) {
    for (int i = 0; i < endKey; i++){
        std::cout<< "post unique " <<" cap Id: "<< capsidInfoVecs.tempCapsideId[i]<< " << mem Id: "<< capsidInfoVecs.tempMembraneId[i]<< " "<< capsidInfoVecs.tempNodeForceX[i] << " "<<capsidInfoVecs.tempNodeForceY[i] << " "<<capsidInfoVecs.tempNodeForceZ[i] << std::endl;
        std::cout<<"dist: "<< capsidInfoVecs.tempLengthsPairs[i]<< std::endl;
    }
}*/
    thrust::for_each(
        thrust::make_zip_iterator(//1st begin
            thrust::make_tuple(
                capsidInfoVecs.tempMembraneId.begin(),
                capsidInfoVecs.tempNodeForceX.begin(),
                capsidInfoVecs.tempNodeForceY.begin(),
                capsidInfoVecs.tempNodeForceZ.begin())),
        thrust::make_zip_iterator(//1st end
            thrust::make_tuple(
                capsidInfoVecs.tempMembraneId.begin(),
                capsidInfoVecs.tempNodeForceX.begin(),
                capsidInfoVecs.tempNodeForceY.begin(),
                capsidInfoVecs.tempNodeForceZ.begin())) + endKey,
        AddForceFunctorAlt(
            capsidInfoVecs.maxNodeCount,
            thrust::raw_pointer_cast(coordInfoVecs.nodeForceX.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeForceY.data()),
            thrust::raw_pointer_cast(coordInfoVecs.nodeForceZ.data())));


   /* if (generalParams.iteration %10000 == 0) {
        for (int i = 0; i < endKey; i++){
            std::cout<< "post add " <<" cap Id: "<< capsidInfoVecs.tempCapsideId[i]<< " << mem Id: "<< capsidInfoVecs.tempMembraneId[i]<< " "<< capsidInfoVecs.tempNodeForceX[i] << " "<<capsidInfoVecs.tempNodeForceY[i] << " "<<capsidInfoVecs.tempNodeForceZ[i] << std::endl;
            std::cout<<"dist: "<< capsidInfoVecs.tempLengthsPairs[i]<< std::endl;
        }
    }*/
        UCVec3 temp = thrust::reduce(        
            thrust::make_zip_iterator(//1st begin
                thrust::make_tuple(
                    capsidInfoVecs.tempMembraneId.begin(),
                    capsidInfoVecs.tempNodeForceX.begin(),
                    capsidInfoVecs.tempNodeForceY.begin(),
                    capsidInfoVecs.tempNodeForceZ.begin())),
            thrust::make_zip_iterator(//1st end
                thrust::make_tuple(
                    capsidInfoVecs.tempMembraneId.begin(),
                    capsidInfoVecs.tempNodeForceX.begin(),
                    capsidInfoVecs.tempNodeForceY.begin(),
                    capsidInfoVecs.tempNodeForceZ.begin())) + endKey,
            init, UCVec3Add());

            
    //if (generalParams.iteration %10000 == 0) {
      /*  if ( (generalParams.iteration % 100) == 0 ){
    
            std::cout<<"linear applied to capsid: "<< -thrust::get<1>(temp) << " "<< -thrust::get<2>(temp) << " "<< -thrust::get<3>(temp) << std::endl;
        }*/
      //  std::cout<<"pos: "<< capsidInfoVecs.nodeLocX[0] << " "<< capsidInfoVecs.nodeLocY[0] << " "<< capsidInfoVecs.nodeLocZ[0] << std::endl;
    //}
        capsidInfoVecs.forceX += -thrust::get<1>(temp);
        capsidInfoVecs.forceY += -thrust::get<2>(temp);
        capsidInfoVecs.forceZ += -thrust::get<3>(temp);
    
};

void AdvanceCapsideParticles(
    GeneralParams& generalParams,
    CapsidInfoVecs& capsidInfoVecs) {
      
      /*  if ( (generalParams.iteration % 100) == 0 ){
    
            std::cout<<"force applied to capsid: "<< capsidInfoVecs.forceX << " "<< capsidInfoVecs.forceY << " "<< capsidInfoVecs.forceZ << std::endl;
        }*/
        //apply force uniformly to each node of capside set
        thrust::transform( 
			thrust::make_zip_iterator(
				thrust::make_tuple(
					capsidInfoVecs.nodeLocX.begin(),
					capsidInfoVecs.nodeLocY.begin(),
					capsidInfoVecs.nodeLocZ.begin() )),
			thrust::make_zip_iterator(
				thrust::make_tuple(
			        capsidInfoVecs.nodeLocX.begin(),
			        capsidInfoVecs.nodeLocY.begin(),
			        capsidInfoVecs.nodeLocZ.begin() )) + capsidInfoVecs.maxNodeCount,
			//save result in third vector to test values
			thrust::make_zip_iterator(
				thrust::make_tuple(
					capsidInfoVecs.nodeLocX.begin(),
					capsidInfoVecs.nodeLocY.begin(),
					capsidInfoVecs.nodeLocZ.begin() )),
			SaxpyFunctorCapsid(
                generalParams.dt,
                capsidInfoVecs.viscosity,
                capsidInfoVecs.forceX,
                capsidInfoVecs.forceY,
                capsidInfoVecs.forceZ));
};
