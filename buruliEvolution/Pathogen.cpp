//
//  Pathogen.cpp
//  buruliEvolution
//
//  Created by Benjamin ROCHE on 13/04/2018.
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
    sum=pPathogenMother->getSum();
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
}

/** mutation process for pathogen*/
void Pathogen::mutation(){
    float lProba=(1-exp(-mutationRate*model->gettimeStep()));
    if(model->getRandom()->getUnif()< lProba ){
        if(sum==0){sum++;}
        else{
            if(model->getRandom()->getUnif()<0.5){
                sum++;
            }
            else{
                sum--;
            }
        }
    }
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
