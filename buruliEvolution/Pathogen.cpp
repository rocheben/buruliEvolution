//
//  Pathogen.cpp
//  buruliEvolution
//
//  Created by Benjamin ROCHE on 13/04/2018.
//  Copyright © 2018 Benjamin ROCHE. All rights reserved.
//

//
//  Cell.cpp
//  clonalEvolution
//
//  Created by Benjamin ROCHE on 14/03/2018.
//  Copyright © 2018 Benjamin ROCHE. All rights reserved.
//


#include "Pathogen.hpp"
using namespace std;
/** Crate new pathogen from another pathogen*/
Pathogen::Pathogen(Pathogen* pPathogenMother){
    viralLoadNeeded=pPathogenMother->getviralLoadNeeded();
    probInfection=pPathogenMother->getprobInfection();
    recoveryPeriod=pPathogenMother->getrecoveryPeriod();
    lifespan=pPathogenMother->getlifespan();
    excretionVolume=pPathogenMother->getexcretionVolume();
    mutationRate=pPathogenMother->getmutationRate();
    model=pPathogenMother->getmodel();
    virulence=pPathogenMother->getvirulence();
    waterInfect=false;
    ID=model->getNewPathogenId();
    sum=pPathogenMother->getSum();
    
    long lGenomeLength=model->getGenomeLength();
    genome=new int[lGenomeLength];
    parent=pPathogenMother;
    int* lGenomeParent=parent->getGenome();
    for(int lIndex=0;lIndex<lGenomeLength;lIndex++){
        genome[lIndex]=lGenomeParent[lIndex];
    }
}

/** Crate new pathogen with parameters*/
Pathogen::Pathogen(int pSum,float pviralLoadNeeded,float pprobInfection,float precoveryPeriod,float plifespan,float pexcretionVolume,float pvirulence, float pmutationRate,Model* pmodel){
    sum=pSum;
    viralLoadNeeded=pviralLoadNeeded;
    probInfection=pprobInfection;
    recoveryPeriod=precoveryPeriod;
    lifespan=plifespan;
    excretionVolume=pexcretionVolume;
    mutationRate=pmutationRate;
    waterInfect=false;
    virulence=pvirulence;
    model=pmodel;
    long lGenomeLength=model->getGenomeLength();
    genome=new int[lGenomeLength];
    for(int lIndex=0;lIndex<lGenomeLength;lIndex++){
        genome[lIndex]=ceil(model->getRandom()->getUnif()*4)-1;
    }
    ID=model->getNewPathogenId();
}


int* Pathogen::getGenome(){
    return genome;
}

/** mutation process for pathogen*/
void Pathogen::mutation(){
    //DEBUG
    float lProba=(1-exp(-mutationRate*model->gettimeStep()));
    if(model->getRandom()->getUnif()< lProba ){
        int baseMutation=ceil(model->getRandom()->getUnif()*(model->getGenomeLength()-4));
        genome[baseMutation+4]=genome[baseMutation]+1;
        if(genome[baseMutation+4]==4){
            genome[baseMutation+4]=0;
        }
    }

}

long Pathogen::getId(){
    return ID;
}

/*int Pathogen::getHammingDistance(double pSum1,double pSum2){
 unsigned lVal=(long)pSum1^(long)pSum2;
 unsigned dist = 0;
 while(lVal)
 {
 ++dist;
 lVal &= lVal - 1;
 }
 return dist;
 }*/

double Pathogen::getSum(){
    return sum;
}

/** Comparing 2 sequences*/
bool Pathogen::samePathogenSequence(Pathogen* pPathogen){
    bool lReturn=false;
    if(sum==pPathogen->getSum()){lReturn=true;}
    return lReturn;
}

/** Gettors*/
float Pathogen::getviralLoadNeeded(){return viralLoadNeeded;}
float Pathogen::getprobInfection(){return probInfection;}
float Pathogen::getrecoveryPeriod(){return recoveryPeriod;}
float Pathogen::getlifespan(){return lifespan;}
float Pathogen::getexcretionVolume(){return excretionVolume;}
float Pathogen::getmutationRate(){return mutationRate;}
bool Pathogen::getwaterInfect(){return waterInfect;}
float Pathogen::getvirulence(){return virulence;}
void Pathogen::setwaterInfect(bool pParam){waterInfect=pParam;}
Model* Pathogen::getmodel(){return model;}
Pathogen* Pathogen::getParent(){return parent;}
