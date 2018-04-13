//
//  agent.hpp
//  buruliEvolution
//
//  Created by Benjamin ROCHE on 13/04/2018.
//  Copyright © 2018 Benjamin ROCHE. All rights reserved.
//


#ifndef AGENT
#define AGENT

class Agent;

#include "list.hpp"
#include "model.hpp"
#include "Lake.hpp"
#include "Pathogen.hpp"
#include <iostream>
#include <math.h>
#define PI 3.14159265

/* */
class Agent{
private:
    /** Pathogens in Infectious state*/
    ListPerso currentPathogen;
    ListPerso currentPathogenWb;
    /** Pathogens in exposed state*/
    ListPerso nextPathogen;
    ListPerso nextPathogenWb;
    //ListPerso nextPathogenInfPeriod;
    //ListPerso currentPathogenInfPeriod;
    /** Pathogens in Recovered state*/
    ListPerso oldPathogen;
    /** Drinking rate of birds*/
    float drinkingRate;
    /** Drinking volume of birds*/
    float drinkingVolume;
    /** Contact rate between birds*/
    float contactRate;
    /** Birth rate*/
    float birthRate;
    /** Death rate*/
    float deathRate;
    /** Species index*/
    int indexSpecies;
    /** Computing variables*/
    Model* model;
    Lake* lake;
    float seasonality;
public:
    
    /** Constructor (first function called at the object's creation)*/
    Agent(Model* pModel,float pdrinkingRate, float pdrinkingVolume, float pcontactRate,float pbirthRate,float pdeathRate,float pSeasonality,int pindexSpecies,int indexAgent);
    Agent(Agent* pAgent,int indexAgent);
    ~Agent();
    double absValue(double pValue);
    int getHammingDistance(int lSum1,int lSum2);
    float getCrossImmunity(Pathogen* pPathogen);
    Pathogen* addPathogen(Pathogen* pPathogen,float pCI,bool pWater);
    void infectionWater(ListPerso* pList1,ListPerso* pList2);
    void infectionDirect(double pRatio,ListPerso* lPathog,ListPerso* lViralLoad);
    void agentStep(double pRatio,ListPerso* pList1Water,ListPerso* pList2Water,ListPerso* pList1Direct,ListPerso* pList2Direct);
    void agentUpdate();
    Model* getmodel();
    float getdrinkingRate();
    float getdrinkingVolume();
    float getcontactRate();
    float getbirthRate(bool pOriginal);
    float getdeathRate(bool pOriginal);
    int getindexSpecies();
    float getseasonality();
    ListPerso* getnextPathogens();
    ListPerso* getcurrentPathogens();
    ListPerso* getoldPathogens();
    ListPerso* getcurrentPathogensFlag();
    void reassortment();
    int getnbNextPathogen();
    void addOldPathogen(Pathogen* pPathogen);
};
#endif
