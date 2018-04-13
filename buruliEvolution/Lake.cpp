//
//  Lake.cpp
//  buruliEvolution
//
//  Created by Benjamin ROCHE on 13/04/2018.
//  Copyright © 2018 Benjamin ROCHE. All rights reserved.
//
#include "Lake.hpp"
#include <math.h>
using namespace std;
/** Constructor*/
Lake::Lake(double pLakeVolume,float pSeasonality,Model* pModel){
    pathogens=new ListPerso();
    viralLoad=new ListPerso();
    seasonality=pSeasonality;
    lakeVolume=pLakeVolume;
    model=pModel;
}

/** Destructor*/
Lake::~Lake(){
    delete pathogens;
    delete viralLoad;
}

/** Return pathogens drunk in the lake regarding the drinking volume and lake volume*/
ListPerso* Lake::getPathogens(){
    return pathogens;
}

/** Return pathogens drunk in the lake regarding the drinking volume and lake volume*/
void Lake::getPathogens(ListPerso* pPathogens,ListPerso* pViralLoad){
    for(int i=0;i<pathogens->length();i++){
        float* lTemp=(float*)viralLoad->getElement(i);
        //if(*lTemp>0){
        pPathogens->add(pathogens->getElement(i));
        pViralLoad->add(lTemp);
        //}
    }
}

/** Return pathogens drunk in the lake regarding the drinking volume and lake volume*/
void Lake::getProbabilityInfection(float drinkingVolume,ListPerso* pPathogens,ListPerso* pProba){
    double lSum=0.0;
    for(int i=0;i<pathogens->length();i++){lSum+=*(float*)viralLoad->getElement(i);}
    for(int i=0;i<pathogens->length();i++){
        Pathogen* lTempPath=(Pathogen*)pathogens->getElement(i);
        
        float* lTemp=new float;
        float lV=(*(float*)viralLoad->getElement(i));
        float lKappa=lTempPath->getviralLoadNeeded();
        float lRate=( ( (float)drinkingVolume/365)/lakeVolume)*(lV/lSum)*(lV/(lV+lKappa));
        *lTemp=1-exp(-lRate*model->gettimeStep());
        if(*lTemp>0){
            pPathogens->add(pathogens->getElement(i));
            pProba->add(lTemp);
        }
    }
}

/** Return pathogens drunk in the lake regarding the drinking volume and lake volume*/
ListPerso* Lake::getViralLoad(){return viralLoad;}

/** Apply viral demography*/
void Lake::viralDemography(){
    ListPerso lPathog;
    ListPerso lListNbInf;
    model->getNbInfect(&lPathog,&lListNbInf);
    //For each pathogens in lake
    for(int i=0;i<pathogens->length();i++){
        //We get the pathogen and its viral load
        Pathogen* lPathogen=(Pathogen*)pathogens->getElement(i);
        float lViralLoad=(*(float*)viralLoad->getElement(i));
        float lTempInf=0;
        //We are looking for every infectious individuals
        for(int j=0;(j<(lPathog.length())&&(lTempInf==0));j++){
            Pathogen* lPathogen2=(Pathogen*)lPathog.getElement(j);
            if(lPathogen2->samePathogenSequence(lPathogen)){
                //If sequence is found, get the infectious population size
                lTempInf=(*(float*)lListNbInf.getElement(j));
            }
        }
        float epsilon=(lPathogen->getexcretionVolume()*365);
        float I1=lTempInf;
        float gamma1=(365/ApplyLakeCharateristics(lPathogen->getlifespan()));
        float V1=lViralLoad;
        float pT=model->gettimeStep()/(float)365;
        float lTemp=epsilon*I1/gamma1+exp(-gamma1*pT)*(V1-(epsilon*(I1)/gamma1));
        //float lTemp=V1/ka  epsilon*I1/gamma1+exp(-gamma1*pT)*(V1-(epsilon*(I1)/gamma1));
        viralLoad->setElement(i,lTemp);
    }
}

/** Apply lake characteristics on pathogen lifespan*/
float Lake::ApplyLakeCharateristics(float pLifespan){
    float lReturn=pLifespan*(1+seasonality*cos((2*PI*( (float)(model->getCurrentTime()) /(float)365))));
    //Apply lake characteristics
    return lReturn;
}

/** Update viral load*/
void Lake::update(){
    viralDemography();
}

void Lake::removeAllPathogens(){
    pathogens->removeAll();
    for(int lIndex=0;lIndex<viralLoad->length();lIndex++){
        delete (float*)viralLoad->getElement(lIndex);
    }
    viralLoad->removeAll();
}


/** Add a pathogen to the lake*/
void Lake::addPathogen(Pathogen* pPathogen,float pViralLoad){
    bool lFound=false;
    for(int i=0;i<pathogens->length();i++){
        if(pPathogen->samePathogenSequence( ((Pathogen*)pathogens->getElement(i))))
        {
            float lTempPointeur=(*(float*)viralLoad->getElement(i))+pViralLoad;
            viralLoad->setElement(i,lTempPointeur);
            lFound=true;
            break;
        }
    }
    if(!lFound){
        float* lTemp=new float;
        *lTemp=pViralLoad;
        pathogens->add(pPathogen);
        viralLoad->add(lTemp);
    }
}

int Lake::getNbPathogens(){
    return pathogens->length();
}
