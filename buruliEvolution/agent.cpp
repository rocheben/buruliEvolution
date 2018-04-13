//
//  agent.cpp
//  buruliEvolution
//
//  Created by Benjamin ROCHE on 13/04/2018.
//  Copyright © 2018 Benjamin ROCHE. All rights reserved.
//

#include "agent.hpp"

using namespace std;
/** Constructor from parameters*/
Agent::Agent(Model* pModel,float pdrinkingRate, float pdrinkingVolume, float pcontactRate,float pbirthRate,float pdeathRate,float pSeasonality,int pindexSpecies,int pindexAgent)
{
    model=pModel;
    drinkingRate=pdrinkingRate;
    drinkingVolume=pdrinkingVolume;
    contactRate=pcontactRate;
    lake=model->getLake();
    birthRate=pbirthRate;
    deathRate=pdeathRate;
    indexSpecies=pindexSpecies;
    seasonality=pSeasonality;
}

/** Constructor from another agent*/
Agent::Agent(Agent* pAgent,int pindexAgent)
{
    model=pAgent->getmodel();
    drinkingRate=pAgent->getdrinkingRate();
    drinkingVolume=pAgent->getdrinkingVolume();
    contactRate=pAgent->getcontactRate();
    lake=model->getLake();
    birthRate=pAgent->getbirthRate(true);
    deathRate=pAgent->getdeathRate(true);
    indexSpecies=pAgent->getindexSpecies();
    seasonality=pAgent->getseasonality();
}

Agent::~Agent(){
    
}
double Agent::absValue(double pValue){
    double lReturn=pValue;
    if(pValue<0){
        lReturn=0-pValue;
    }
    return lReturn;
}

int Agent::getHammingDistance(int pSum1,int pSum2){
    unsigned lVal=(long)pSum1^(long)pSum2;
    unsigned dist = 0;
    while(lVal)
    {
        ++dist;
        lVal &= lVal - 1;
    }
    return dist;
}

/** Compute cross-immunity*/
float Agent::getCrossImmunity(Pathogen* pPathogen){
    float lReturn=0;
    if(nextPathogen.length()==0 && currentPathogen.length()==0 && oldPathogen.length()==0){
        lReturn=1;
    }
    else{
        double lSum=pPathogen->getSum();
        double lSumMax=32768;
        double lMinDistance=lSumMax;
        for(int i=0;( (i<nextPathogen.length()) && (lMinDistance>0));i++){
            double lSumNew=((Pathogen*)(nextPathogen.getElement(i)))->getSum();
            if(lMinDistance>absValue(lSum-lSumNew)){
                lMinDistance=absValue(lSum-lSumNew);
            }
        }
        for(int i=0;( (i<currentPathogen.length()) && (lMinDistance>0));i++){
            double lSumNew=((Pathogen*)(currentPathogen.getElement(i)))->getSum();
            if(lMinDistance>absValue(lSum-lSumNew)){
                lMinDistance=absValue(lSum-lSumNew);
            }
        }
        //Looking for the closest pathogen in Recovered state
        for(int i=0;( (i<oldPathogen.length()) && (lMinDistance>0));i++){
            double lSumNew=((Pathogen*)(oldPathogen.getElement(i)))->getSum();
            if(lMinDistance>absValue(lSum-lSumNew)){
                lMinDistance=absValue(lSum-lSumNew);
            }
        }
        //double lMean=model->getmaxCrossImmunity()*(1-exp(-pow(lMinDistance,2)/5));
        double lMean=model->getmaxCrossImmunity()*(1-exp(-pow( ((lMinDistance)/2),2)));
        lReturn=lMean;
    }
    return lReturn;
}


/** Immunize host*/
void Agent::addOldPathogen(Pathogen* pPathogen){
    oldPathogen.add(pPathogen);
}

/** Add a pathogen which succesfully infect bird individual*/
Pathogen* Agent::addPathogen(Pathogen* pPathogen,float pCI,bool pWater){
    Pathogen* lPathogenTemp=new Pathogen(pPathogen);
    Pathogen* lPathogen=model->getPathogenAdd(lPathogenTemp);
    delete lPathogenTemp;
    nextPathogen.add(lPathogen);
    //float* lTempNew=new float;
    //*lTempNew=model->getIntercept()+(lPathogen->getrecoveryPeriod()-model->getIntercept())*pCI;
    /* *lTempNew=lPathogen->getrecoveryPeriod()*pCI;
     nextPathogenInfPeriod.add(lTempNew);*/
    float* lTempFlag=new float;
    *lTempFlag=0;
    if(pWater){*lTempFlag=1;}
    nextPathogenWb.add(lTempFlag);
    return lPathogen;
}

/** Function for drinking water and gets infection*/
void Agent::infectionWater(ListPerso* lPathog,ListPerso* lProbaInf){
    if(drinkingRate!=0){
        ListPerso lIndexList;
        
        for(int lIndex=0;lIndex<lPathog->length();lIndex++){float *lIndexTemp=new float;*lIndexTemp=lIndex;lIndexList.add(lIndexTemp);}
        
        //float* lProbaSum=new float[lPathog.length()];
        float lProba=0;
        //for(int lIndex=0;((lIndex<lPathog.length())&&(nextPathogen.length()+currentPathogen.length()<3));lIndex++){
        for(int lIndex=0;lIndex<lPathog->length();lIndex++){
            
            float lIndexRand=round(model->getRandom()->getUnif()*(lIndexList.length()-1));
            float* lIndexRand1=(float*)lIndexList.getElement((int)lIndexRand);
            int lIndexTemp=(int)(*lIndexRand1);
            lIndexList.remove(lIndexRand1);
            delete lIndexRand1;
            //int lIndexTemp=lIndex;
            
            lProba=(*(float*)lProbaInf->getElement(lIndexTemp));
            if(lProba>0){
                float lCI=getCrossImmunity((Pathogen*)lPathog->getElement(lIndexTemp));
                if(lCI>0){
                    float lProbaTemp=lProba*lCI;
                    if(model->getRandom()->getUnif()<lProbaTemp){
                        Pathogen* lTemp=(Pathogen*)lPathog->getElement((int)lIndexTemp);
                        Pathogen* lPathogen=addPathogen(lTemp,lCI,true);
                        if(lPathogen!=0){
                            //lPathogen->setwaterInfect(true);
                            //model->dataWrite(indexSpecies,lPathogen);
                            lake->addPathogen(lPathogen,0);
                        }
                    }
                }
            }
        }
        /*for(int lIndex=0;lIndex<lProbaInf->length();lIndex++){
         delete (float*)lProbaInf.getElement(lIndex);
         }*/
    }
}

/** Function for inter-individual transmission*/
void Agent::infectionDirect(double pRatio,ListPerso* lPathog,ListPerso* lViralLoad){
    if(contactRate!=0){
        float lProba;
        float* lIndexTempnew=new float[lPathog->length()];
        
        ListPerso lIndexList;
        
        
        for(int lIndex=0;lIndex<lPathog->length();lIndex++){
            lIndexTempnew[lIndex]=lIndex;
            lIndexList.add(&lIndexTempnew[lIndex]);
        }
        
        for(int lIndex=0;lIndex<lPathog->length();lIndex++){
            
            double lUnif=model->getRandom()->getUnif();
            
            double lLength=(lIndexList.length()-1);
            float lIndexRand=round(lUnif*lLength);
            
            float* lIndexRand1=(float*)lIndexList.getElement((int)lIndexRand);
            int lIndexTemp=(int)(*lIndexRand1);
            lIndexList.remove((void*)lIndexRand1);
            
            Pathogen* lTemp=(Pathogen*)lPathog->getElement(lIndexTemp);
            float lTempNbInf=(*(float*)lViralLoad->getElement(lIndexTemp));
            
            float lTempProbInf=lTemp->getprobInfection()*(1+model->getampTransm()*sin((2*PI* ( ((float)model->getCurrentTime())/((float)365)) )));
            float lProbaTemp=1-exp(-((lTempProbInf/365)*lTempNbInf)*model->gettimeStep());
            lProbaTemp=lProbaTemp/pRatio;
            
            if(lProbaTemp>0){
                float lCI=getCrossImmunity(lTemp);
                if(lCI>0){
                    lProba=lCI*lProbaTemp;
                    if(model->getRandom()->getUnif()<lProba){
                        if(currentPathogen.length()>1){
                            lCI=getCrossImmunity(lTemp);
                        }
                        Pathogen* lPathogen=addPathogen(lTemp,lCI,false);
                        if(lPathogen!=0){
                            //lPathogen->setwaterInfect(false);
                            lake->addPathogen(lPathogen,0);
                        }
                    }
                }
            }
        }
        delete [] lIndexTempnew;
    }
}

/** This function is the first called on object by the model scheduler: It sets 'time' of change (if it occurs) and implements all events*/
void Agent::agentStep(double pRatio,ListPerso* pList1Water,ListPerso* pList2Water,ListPerso* pList1Direct,ListPerso* pList2Direct)
{
    // Environmental transmission
    infectionWater(pList1Water,pList2Water);
    // Inter-individuals transmission
    infectionDirect(pRatio,pList1Direct,pList2Direct);
}

void Agent::reassortment(){
    /** Mutation process*/
    for(int lIndex=0;lIndex<currentPathogen.length();lIndex++){
        Pathogen* pPathogen=(Pathogen*)currentPathogen.getElement(lIndex);
        Pathogen* lPathogenTemp=new Pathogen(pPathogen);
        double lGenome=lPathogenTemp->getSum();
        lPathogenTemp->mutation();
        Pathogen* lPathogen=model->getPathogenAdd(lPathogenTemp);
        currentPathogen.setElement(lIndex,lPathogen);
        if(lPathogen->getSum()!=lGenome){
            model->writeEvolEvent(pPathogen,lPathogen,0);
        }
        delete lPathogenTemp;
    }
}

/** This function is called when function Update is called on all objects: It changes the state of agent initialize by function agentStep*/
void Agent::agentUpdate()
{
    ListPerso toEliminate;
    //ListPerso toEliminatePeriod;
    ListPerso toEliminateFlag;
    
    //From R to S
    for(int i=0;i<oldPathogen.length();i++){
        float lLatencyPeriod=model->getimmmunityPeriod();
        if(lLatencyPeriod>0){
            float lRate=1/(lLatencyPeriod);
            float lTimeStep=model->gettimeStep();
            float lProbaRecovery=1-exp(-lRate*lTimeStep);
            if(lProbaRecovery>0){
                if(model->getRandom()->getUnif()< lProbaRecovery){
                    toEliminate.add(oldPathogen.getElement(i));
                    //toEliminatePeriod.add(currentPathogenInfPeriod.getElement(i));
                }
            }
        }
    }
    if(toEliminate.length()>0){
        oldPathogen.removeAll(&toEliminate);
        toEliminate.removeAll();
    }
    //From I to R
    for(int i=0;i<currentPathogen.length();i++){
        float lLatencyPeriod=((Pathogen*)currentPathogen.getElement(i))->getrecoveryPeriod();
        //lLatencyPeriod=*((float*)currentPathogenInfPeriod.getElement(i));
        if(lLatencyPeriod>0){
            float lRate=1/(lLatencyPeriod);
            float lTimeStep=model->gettimeStep();
            float lProbaRecovery=1-exp(-lRate*lTimeStep);
            if(lProbaRecovery>0){
                if(model->getRandom()->getUnif()< lProbaRecovery){
                    toEliminate.add(currentPathogen.getElement(i));
                    //toEliminatePeriod.add(currentPathogenInfPeriod.getElement(i));
                    toEliminateFlag.add(currentPathogenWb.getElement(i));
                }
            }
        }
    }
    if(toEliminate.length()>0){
        oldPathogen.addAll(&toEliminate);
        currentPathogen.removeAll(&toEliminate);
        toEliminate.removeAll();
        
        //currentPathogenInfPeriod.removeAll(&toEliminatePeriod);
        /*for(int lIndex=0;lIndex<toEliminatePeriod.length();lIndex++){
         float* lTemp=(float*)toEliminatePeriod.getElement(lIndex);
         delete lTemp;
         }*/
        //toEliminatePeriod.removeAll();
        currentPathogenWb.removeAll(&toEliminateFlag);
        for(int lIndex=0;lIndex<toEliminateFlag.length();lIndex++){
            float* lTemp=(float*)toEliminateFlag.getElement(lIndex);
            delete lTemp;
        }
        toEliminateFlag.removeAll();
    }
    bool lNewCoInf=false;
    if(nextPathogen.length()+currentPathogen.length()>1 && nextPathogen.length()>0){
        lNewCoInf=true;
    }
    //From S to I
    for(int i=0;i<nextPathogen.length();i++){
        toEliminate.add(nextPathogen.getElement(i));
        //toEliminatePeriod.add(nextPathogenInfPeriod.getElement(i));
        toEliminateFlag.add(nextPathogenWb.getElement(i));
    }
    if(toEliminate.length()>0){
        currentPathogen.addAll(&toEliminate);
        //currentPathogenInfPeriod.addAll(&toEliminatePeriod);
        currentPathogenWb.addAll(&toEliminateFlag);
        nextPathogen.removeAll(&toEliminate);
        //nextPathogenInfPeriod.removeAll(&toEliminatePeriod);
        nextPathogenWb.removeAll(&toEliminateFlag);
        toEliminate.removeAll();
        //toEliminatePeriod.removeAll();
        toEliminateFlag.removeAll();
    }
    if(lNewCoInf){
        int lLength=currentPathogen.length();
        double* lPathogens=new double[lLength];
        for(int lIndex=0;lIndex<lLength;lIndex++){lPathogens[lIndex]=((Pathogen*)currentPathogen.getElement(lIndex))->getSum();}
        model->writeCoinf(lPathogens,lLength);
        delete[] lPathogens;
    }
    reassortment();
}

/** Gettors*/
Model* Agent::getmodel(){return model;}
float Agent::getdrinkingRate(){return drinkingRate;}
float Agent::getdrinkingVolume(){return drinkingVolume;}
float Agent::getcontactRate(){return contactRate;}
float Agent::getbirthRate(bool pOriginal){
    float lTemp=birthRate;
    if(!pOriginal){
        lTemp=birthRate*(1-seasonality*sin((2*PI* ( ((float)model->getCurrentTime())/((float)365)) )));
    }
    //lTemp=birthRate;
    return lTemp;
}
float Agent::getdeathRate(bool pOriginal){
    float lTemp=deathRate;
    if(!pOriginal){
        lTemp=deathRate*(1+seasonality*sin((2*PI* ( ((float)model->getCurrentTime())/((float)365)) )));
    }
    //float lReturn=deathRate;
    float lReturn=lTemp;
    
    /*for(int lIndex=0;lIndex<currentPathogen.length();lIndex++){
     float lVirulence=((Pathogen*)currentPathogen.getElement(lIndex))->getvirulence();
     lReturn=lReturn-lVirulence;
     }*/
    return lReturn;
}
int Agent::getindexSpecies(){return indexSpecies;}
float Agent::getseasonality(){return seasonality;}
ListPerso* Agent::getnextPathogens(){return &nextPathogen;}
ListPerso* Agent::getcurrentPathogens(){return &currentPathogen;}
ListPerso* Agent::getcurrentPathogensFlag(){return &currentPathogenWb;}
ListPerso* Agent::getoldPathogens(){return &oldPathogen;}
int Agent::getnbNextPathogen(){return nextPathogen.length();}
