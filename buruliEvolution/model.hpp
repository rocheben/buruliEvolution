//
//  model.hpp
//  buruliEvolution
//
//  Created by Benjamin ROCHE on 13/04/2018.
//  Copyright © 2018 Benjamin ROCHE. All rights reserved.
//
#ifndef MODEL
#define MODEL

class Model;

#include "list.hpp"
#include "agent.hpp"
#include "Lake.hpp"
#include "Pathogen.hpp"
#include "random.h"
#define MAX_LENGTH 1000
class Model{
private:
    /** Parameters*/
    char birdsFile[256];
    char pathogensFile[256];
    Random* randomObject;
    /** Object Lake*/
    Lake* theLake;
    Lake* infectiousInd;
    Lake* infectiousIndWb;
    Lake* migrationLoads;
    /** ListPerso of all agents*/
    ListPerso* agentListPerso;
    ListPerso* pathogenListPerso;
    ListPerso* lIntrod;
    ListPerso* lIntrodRecov;
    ListPerso* pathogenIntrod;
    double tMax;
    char output[1024];
    char outputFolder[4096];
    double lakeVolume;
    double randomSeed;
    double currentTime;
    float tradeOffValues[8];
    int outputAll;
    int crossImmunityFlag;
    float timeStep;
    int year;
    float maxCrossImmunity;
    float ampLake;
    float ampTransm;
    float reassortment;
    int N;
    char dataOutput[MAX_LENGTH+1];
    char dataOutputAllAnti[MAX_LENGTH+1];
    char dataOutputEvents[MAX_LENGTH+1];
    char dataOutputCoinf[MAX_LENGTH+1];
    char dataOutputPop[MAX_LENGTH+1];
    int intervalMigration;
    int indexAgent;
    int immunity;
    float intervalWrite;
    double interceptInfPeriod;
    int modifyInfPeriod;
    float immmunityPeriod;
public:
    void convertSequence(int* pSequence,char* pOut,int pLength);
    Model(char* pSuffixe,float pTimeStep);
    ~Model();
    void loadBirdsSpecies();
    void loadPathogensSpecies();
    //void loadTradeOff(char* pFileName);
    void buildObjects();
    Lake* getLake();
    Random* getRandom();
    void addAgent(Agent* pAgent);
    void delAgent(Agent* pAgent);
    void goNormal();
    void goOptimized();
    bool loadParams(char*);
    void Init();
    bool allDataWrite(bool pWrite);
    bool getCrossImmun();
    void getNbInfect(ListPerso* pPathogens,ListPerso* pViralLoad);
    float getmaxCrossImmunity();
    double getCurrentTime();
    float getreassortment();
    Pathogen* getPathogenAdd(Pathogen* pPathogen);
    void closeBuffer();
    float gettimeStep();
    bool checkPreviousSimulations();
    int getPopSize();
    void updateNbInfectious(Agent* pAgent);
    int getimmunity();
    double getIntercept();
    int getmodifyInfPeriod();
    void immigration();
    void writeEvolEvent(Pathogen*,Pathogen*,Pathogen*);
    void writeCoinf(double* pPathogen,int pLength);
    float getampTransm();
    float getimmmunityPeriod();
    void writeBuffer(char* pFileName,char* pBuffer,char* pLine);
};
#endif
