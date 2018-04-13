//
//  Pathogen.hpp
//  buruliEvolution
//
//  Created by Benjamin ROCHE on 13/04/2018.
//  Copyright © 2018 Benjamin ROCHE. All rights reserved.
//
#ifndef PATHOGEN
#define PATHOGEN

class Pathogen;
#include "model.hpp"
#include <math.h>
#include <iostream>
#include <stdio.h>

class Pathogen{
private:
    float viralLoadNeeded;
    float probInfection;
    float recoveryPeriod;
    float lifespan;
    float mutationRate;
    float excretionVolume;
    bool waterInfect;
    int sequenceLength;
    float virulence;
    Model* model;
    double sum;
    double subtype;
public:
    /** Crate new pathogen*/
    Pathogen(Pathogen* pPathogenMother);
    Pathogen(int pSum,float pviralLoadNeeded,float pprobInfection,float precoveryPeriod,float plifespan,float pexcretionVolume,float pvirulence,float pmutationRate,Model* pmodel);
    //int getHammingDistance(double pSum1,double pSum2);
    void mutation();
    float getviralLoadNeeded();
    float getprobInfection();
    float getrecoveryPeriod();
    float getincubationPeriod();
    float getlifespan();
    float getexcretionVolume();
    float getmutationRate();
    bool getwaterInfect();
    void setwaterInfect(bool pParam);
    float getvirulence();
    Model* getmodel();
    bool samePathogenSequence(Pathogen* pPathogen);
    double getSum();
};
#endif
