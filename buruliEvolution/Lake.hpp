//
//  Lake.hpp
//  buruliEvolution
//
//  Created by Benjamin ROCHE on 13/04/2018.
//  Copyright © 2018 Benjamin ROCHE. All rights reserved.
//
#ifndef LAKE
#define LAKE

class Lake;

#include "list.hpp"
#include "model.hpp"
#include "Pathogen.hpp"

class Lake{
private:
    ListPerso* pathogens;
    ListPerso* viralLoad;
    double lakeVolume;
    int verbose;
    float seasonality;
    Model* model;
public:
    Lake(double pLakeVolume,float pSeasonality,Model* pModel);
    ~Lake();
    void getPathogens(ListPerso* pPathogens,ListPerso* pViralLoad);
    void getProbabilityInfection(float drinkingVolume,ListPerso* pPathogens,ListPerso* pViralLoad);
    void viralDemography();
    float ApplyLakeCharateristics(float pLifespan);
    void update();
    void setSeasonality(float pSeasonality);
    void addPathogen(Pathogen* pPathogen,float pViralLoad);
    ListPerso* getPathogens();
    ListPerso* getViralLoad();
    void removeAllPathogens();
    int getNbPathogens();
};
#endif
