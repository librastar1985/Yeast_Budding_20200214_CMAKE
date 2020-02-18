
#include "System.h"
#include "SystemStructures.h"
#include "CapsidSpringsAll.h"


void ComputeCapsideSpringsAll(
    GeneralParams& generalParams,
    CoordInfoVecs& coordInfoVecs,
    CapsidInfoVecs& capsidInfoVecs,
    AuxVecs& auxVecs) {

    int factor = 15;//number of interactions a single capside point can have.
    //use these as reduced vectors.    
    
    capsidInfoVecs.tempCapsideId.resize(factor * capsidInfoVecs.maxNodeCount);
    capsidInfoVecs.tempMembraneId.resize(factor * capsidInfoVecs.maxNodeCount);

    thrust::device_vector<int> tempMembraneIdReduced;
    thrust::device_vector<int> tempCapsidIdReduced;
    thrust::device_vector<double> tempNodeForceXReduced;
    thrust::device_vector<double> tempNodeForceYReduced;
    thrust::device_vector<double> tempNodeForceZReduced;

    tempNodeForceXReduced.resize(factor * capsidInfoVecs.maxNodeCount);
    tempNodeForceYReduced.resize(factor * capsidInfoVecs.maxNodeCount);
    tempNodeForceZReduced.resize(factor * capsidInfoVecs.maxNodeCount);
    tempMembraneIdReduced.resize(factor * capsidInfoVecs.maxNodeCount);
    tempCapsidIdReduced.resize(factor * capsidInfoVecs.maxNodeCount);



    thrust::device_vector<int> tempMembraneIdUnreduced;
    thrust::device_vector<int> tempCapsidIdUnreduced;
    thrust::device_vector<double> tempNodeForceXUnreduced;
    thrust::device_vector<double> tempNodeForceYUnreduced;
    thrust::device_vector<double> tempNodeForceZUnreduced;
    
    tempCapsidIdUnreduced.resize(factor * capsidInfoVecs.maxNodeCount);
    tempMembraneIdUnreduced.resize(factor * capsidInfoVecs.maxNodeCount);
    tempNodeForceXUnreduced.resize(factor * capsidInfoVecs.maxNodeCount);
    tempNodeForceYUnreduced.resize(factor * capsidInfoVecs.maxNodeCount);
    tempNodeForceZUnreduced.resize(factor * capsidInfoVecs.maxNodeCount);
    thrust::fill(tempMembraneIdUnreduced.begin(),tempMembraneIdUnreduced.end(), generalParams.maxNodeCount);

    
    thrust::counting_iterator<int> begin(0);
    //for each capsid node, choose the closest membrane node. 
    //what if I chose the closest capsid node for each membrane node?
    thrust::for_each(  
        thrust::make_zip_iterator(
            thrust::make_tuple(
                begin,
                capsidInfoVecs.bucketKeys.begin(),
                capsidInfoVecs.bucketValues.begin())),
        
        thrust::make_zip_iterator(
            thrust::make_tuple(
                begin,
                capsidInfoVecs.bucketKeys.begin(),
                capsidInfoVecs.bucketValues.begin())) + capsidInfoVecs.maxNodeCount,
            CapsidSpringFunctorAll(
                factor,
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
                
                thrust::raw_pointer_cast(tempCapsidIdUnreduced.data()),
                thrust::raw_pointer_cast(tempMembraneIdUnreduced.data()),
                thrust::raw_pointer_cast(tempNodeForceXUnreduced.data()),
                thrust::raw_pointer_cast(tempNodeForceYUnreduced.data()),
                thrust::raw_pointer_cast(tempNodeForceZUnreduced.data()) ));

    //now we have un reduced forces. Sort by id and reduce. 
    //key, then value. Each vector returns sorted	  
    //sort by capsid  
    thrust::sort_by_key(tempCapsidIdUnreduced.begin(), tempCapsidIdUnreduced.end(),
    thrust::make_zip_iterator(
        thrust::make_tuple(
            tempMembraneIdUnreduced.begin(),
            tempNodeForceXUnreduced.begin(),
            tempNodeForceYUnreduced.begin(),
            tempNodeForceZUnreduced.begin())), thrust::less<int>());	

    //stable sort by membrane for reduciton prep
    thrust::stable_sort_by_key(tempMembraneIdUnreduced.begin(), tempMembraneIdUnreduced.end(),
    thrust::make_zip_iterator(
        thrust::make_tuple(
            tempCapsidIdUnreduced.begin(),
            tempNodeForceXUnreduced.begin(),
            tempNodeForceYUnreduced.begin(),
            tempNodeForceZUnreduced.begin())), thrust::less<int>());
    
    //to image, we need the corresponding id's of edges. 
    //copy these to tempCapsidId and tempMembraneId
    //the number of edges made is the amount of id's smaller than generalParams.maxNodeCount
    thrust::copy(tempMembraneIdUnreduced.begin(),tempMembraneIdUnreduced.end(), capsidInfoVecs.tempMembraneId.begin());
    thrust::copy(tempCapsidIdUnreduced.begin(),tempCapsidIdUnreduced.end(), capsidInfoVecs.tempCapsideId.begin());

    capsidInfoVecs.num_connections= thrust::count_if(
        capsidInfoVecs.tempMembraneId.begin(),
        capsidInfoVecs.tempMembraneId.end(), is_less_than(generalParams.maxNodeCount) );
    

    //Reduce and apply force
    int endKey = thrust::get<0>(
       thrust::reduce_by_key(
            tempMembraneIdUnreduced.begin(), 
            tempMembraneIdUnreduced.end(),
       thrust::make_zip_iterator(
           thrust::make_tuple(
               tempNodeForceXUnreduced.begin(),
               tempNodeForceYUnreduced.begin(),
               tempNodeForceZUnreduced.begin())),
               tempMembraneIdReduced.begin(),
       thrust::make_zip_iterator(
           thrust::make_tuple(
               tempNodeForceXReduced.begin(),
               tempNodeForceYReduced.begin(),
               tempNodeForceZReduced.begin())),
       thrust::equal_to<int>(), CVec3Add())) - tempMembraneIdReduced.begin();//binary_pred, binary_op 


    
    
   
    thrust::for_each(
        thrust::make_zip_iterator(//1st begin
            thrust::make_tuple(
                tempMembraneIdReduced.begin(),
                tempNodeForceXReduced.begin(),
                tempNodeForceYReduced.begin(),
                tempNodeForceZReduced.begin())),
        thrust::make_zip_iterator(//1st end
            thrust::make_tuple(
                tempMembraneIdReduced.begin(),
                tempNodeForceXReduced.begin(),
                tempNodeForceYReduced.begin(),
                tempNodeForceZReduced.begin())) + endKey,
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
        UCVec3 init(0, 0.0,0.0,0.0);
        UCVec3 temp = thrust::reduce(        
            thrust::make_zip_iterator(//1st begin
                thrust::make_tuple(
                    tempMembraneIdReduced.begin(),
                    tempNodeForceXReduced.begin(),
                    tempNodeForceYReduced.begin(),
                    tempNodeForceZReduced.begin())),
            thrust::make_zip_iterator(//1st end
                thrust::make_tuple(
                    tempMembraneIdReduced.begin(),
                    tempNodeForceXReduced.begin(),
                    tempNodeForceYReduced.begin(),
                    tempNodeForceZReduced.begin())) + endKey,
            init, UCVec3Add());

    //if (generalParams.iteration %10000 == 0) {
       // std::cout<<"linear applied to capsid: "<< thrust::get<0>(temp) << " "<<thrust::get<1>(temp) << " "<<thrust::get<2>(temp) << " "<<" "<< thrust::get<3>(temp) << std::endl;
            
      //  std::cout<<"pos: "<< capsidInfoVecs.nodeLocX[0] << " "<< capsidInfoVecs.nodeLocY[0] << " "<< capsidInfoVecs.nodeLocZ[0] << std::endl;
    //}
        capsidInfoVecs.forceX += -thrust::get<1>(temp);
        capsidInfoVecs.forceY += -thrust::get<2>(temp);
        capsidInfoVecs.forceZ += -thrust::get<3>(temp);
    
};

