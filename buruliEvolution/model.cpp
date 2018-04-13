//
//  model.cpp
//  buruliEvolution
//
//  Created by Benjamin ROCHE on 13/04/2018.
//  Copyright © 2018 Benjamin ROCHE. All rights reserved.
//
#include "model.hpp"
#include <string.h>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

using namespace std;

/** Constructor*/
Model::Model(char* pSuffixe,float pTimeStep){
    strcpy(output,pSuffixe);
    //agentListPerso=new ListPerso();
    strcpy(birdsFile,"/home/roche/Documents/Modeles/AvianFluStrains/avianFluStrains/species.dat");
    strcpy(pathogensFile,"/home/roche/Documents/Modeles/AvianFluStrains/avianFluStrains/pathogens.dat");
    strcpy(outputFolder,".");
    lakeVolume=10000;
    randomSeed=time(0);
    tMax=730;
    outputAll=1;
    crossImmunityFlag=1;
    lIntrod=new ListPerso();
    lIntrodRecov=new ListPerso();
    pathogenIntrod=new ListPerso();
    pathogenListPerso=new ListPerso();
    maxCrossImmunity=1;
    ampLake=0;
    ampTransm=0;
    reassortment=0;
    intervalMigration=150;
    timeStep=pTimeStep;
    indexAgent=0;
    intervalWrite=7;
    interceptInfPeriod=0.0;
}

bool Model::checkPreviousSimulations(){
    bool lReturn=true;
    char lFileName[2048];
    sprintf(lFileName,"%s/checkend_%s.dat",outputFolder,output);
    FILE* lFile=fopen(lFileName,"rt");
    if(lFile!=NULL){
        fgets(lFileName,2047,lFile);
        fclose(lFile);
        if(strstr(lFileName,"OK")!=NULL){lReturn=false;}
    }
    return lReturn;
}
/** Initialization function*/
void Model::Init(){
    //Object initialization
    theLake=new Lake(lakeVolume,ampLake,this);
    infectiousInd=new Lake(1,0,this);
    infectiousIndWb=new Lake(1,0,this);
    migrationLoads=new Lake(1,0,this);
    randomObject=new Random(randomSeed,output);
    
    //Files creation
    char lFileName[2048];
    sprintf(lFileName,"%s/seed_%s.dat",outputFolder,output);
    FILE* lFile=fopen(lFileName,"wt");
    sprintf(lFileName,"%u",randomObject->getSeed());
    fwrite(lFileName,1,strlen(lFileName),lFile);
    fclose(lFile);
    if(outputAll!=0){
        sprintf(lFileName,"%s/output_All_Anti_%s.dat",outputFolder,output);
        lFile=fopen(lFileName,"wt");
        fclose(lFile);
        sprintf(lFileName,"%s/output_All_Events_%s.dat",outputFolder,output);
        lFile=fopen(lFileName,"wt");
        fclose(lFile);
        sprintf(lFileName,"%s/output_All_Coinf_%s.dat",outputFolder,output);
        lFile=fopen(lFileName,"wt");
        fclose(lFile);
        sprintf(lFileName,"%s/output_All_Pop_%s.dat",outputFolder,output);
        lFile=fopen(lFileName,"wt");
        fclose(lFile);
    }
    year=0;
    strcpy(dataOutputAllAnti,"");
    strcpy(dataOutputEvents,"");
    strcpy(dataOutputCoinf,"");
    strcpy(dataOutputPop,"");
}

void Model::closeBuffer(){
    if(outputAll!=0){
        char lFileName[256];
        sprintf(lFileName,"%s/output_All_Anti_%s.dat",outputFolder,output);
        FILE* lFile=fopen(lFileName,"at");
        fwrite(dataOutputAllAnti,1,strlen(dataOutputAllAnti),lFile);
        fclose(lFile);
        sprintf(lFileName,"%s/output_All_Events_%s.dat",outputFolder,output);
        lFile=fopen(lFileName,"at");
        fwrite(dataOutputAllAnti,1,strlen(dataOutputAllAnti),lFile);
        fclose(lFile);
        sprintf(lFileName,"%s/output_All_Coinf_%s.dat",outputFolder,output);
        lFile=fopen(lFileName,"at");
        fwrite(dataOutputAllAnti,1,strlen(dataOutputAllAnti),lFile);
        fclose(lFile);
        sprintf(lFileName,"%s/output_All_Pop_%s.dat",outputFolder,output);
        lFile=fopen(lFileName,"at");
        fwrite(dataOutputAllAnti,1,strlen(dataOutputAllAnti),lFile);
        fclose(lFile);
    }
}

/** Destructor*/
Model::~Model(){
    closeBuffer();
    char lFileName[2048];
    sprintf(lFileName,"%s/checkend_%s.dat",outputFolder,output);
    FILE* lFile=fopen(lFileName,"wt");
    fwrite("OK",strlen("OK"),1,lFile);
    fclose(lFile);
    theLake->removeAllPathogens();
    delete theLake;
    infectiousInd->removeAllPathogens();
    infectiousIndWb->removeAllPathogens();
    migrationLoads->removeAllPathogens();
    delete infectiousInd;
    delete infectiousIndWb;
    delete migrationLoads;
    for(int i=0;i<agentListPerso->length();i++){delete (Agent*)agentListPerso->getElement(i);}
    delete agentListPerso;
    delete randomObject;
    delete lIntrod;
    delete lIntrodRecov;
    for(int i=0;i<pathogenListPerso->length();i++){delete (Pathogen*)pathogenListPerso->getElement(i);}
    delete pathogenListPerso;
}

/** Load bird species*/
void Model::loadBirdsSpecies(){
    char lLine[2048];
    FILE* lFile=fopen(birdsFile,"rt");
    int lIndexSpecies=0;
    //For each line within species file
    while(fgets(lLine,2047,lFile)){
        //If this line is not commented
        if( (strstr(lLine,"//")==NULL) && (strstr(lLine,"#")==NULL)){
            int lPopSize;
            float lParams[6];
            //Read parameter line
            sscanf(lLine,"%d\t%f\t%f\t%f\t%f\t%f\t%f",&lPopSize,&lParams[0],&lParams[1],&lParams[2],&lParams[3],&lParams[4],&lParams[5]);
            N=lPopSize;
            agentListPerso=new ListPerso(N);
            //For each individual of the population
            for(int j=0;j<lPopSize;j++){
                //cout<<j<<endl;
                //Agent creation
                Agent* lTemp=new Agent(this,lParams[0],lParams[1],lParams[2],lParams[3],lParams[4],lParams[5],lIndexSpecies,indexAgent);
                indexAgent++;
                //Adding to the list
                agentListPerso->add(lTemp);
                
            }
            lIndexSpecies++;
        }
    }
    fclose(lFile);
    //Pathogen introduction
    int previousTampon=0;
    for(int k=0;k<pathogenIntrod->length();k++){
        int lTampon=(int)(*(float*)lIntrod->getElement(k));
        if(lTampon>=1){
            for(int lIndex=previousTampon;lIndex<lTampon+previousTampon;lIndex++){
                Agent* lTemp=(Agent*)agentListPerso->getElement(lIndex);
                Pathogen* lPathogen=(Pathogen*)pathogenIntrod->getElement(k);
                lTemp->addPathogen(lPathogen,1,false);
                lTemp->agentUpdate();
                updateNbInfectious(lTemp);
            }
            previousTampon+=lTampon;
        }
        else{
            for(int lIndex=0;lIndex<agentListPerso->length();lIndex++){
                if(randomObject->getUnif()<(*(float*)lIntrod->getElement(k))){
                    Agent* lTemp=(Agent*)agentListPerso->getElement(lIndex);
                    lTemp->addPathogen((Pathogen*)pathogenIntrod->getElement(k),1,false);
                    lTemp->agentUpdate();
                    updateNbInfectious(lTemp);
                }
            }
        }
        int lTamponNew=(int)(*(float*)lIntrodRecov->getElement(k));
        for(int lIndex=previousTampon;lIndex<previousTampon+lTamponNew;lIndex++){
            Agent* lTemp=(Agent*)agentListPerso->getElement(lIndex);
            lTemp->addOldPathogen((Pathogen*)pathogenIntrod->getElement(k));
        }
    }
}

/** Load pathogens species*/
void Model::loadPathogensSpecies(){
    char lLine[4096];
    FILE* lFile=fopen(pathogensFile,"rt");
    //int lNbBits=-1;
    //Read the line
    while(fgets(lLine,4095,lFile)){
        if( (strstr(lLine,"//")==NULL) && (strstr(lLine,"#")==NULL)){
            //If sequence length is not yet determined
            /*if(lNbBits==-1){
             //Calculating sequence length
             char* lTemp=strstr(lLine,"\t");
             lNbBits=strlen(lLine)-strlen(lTemp);
             }
             
             //Get the sequence
             char* lSequenceTemp=new char[lNbBits+1];
             int* lSequence=new int[lNbBits];
             strncpy(lSequenceTemp,lLine,lNbBits);
             for(int k=0;k<lNbBits;k++){
             if(lSequenceTemp[k]=='0'){
             lSequence[k]=0;
             }
             else{
             lSequence[k]=1;
             }
             }
             char* lTempChar=lLine+lNbBits+1;*/
            float lParams[11];
            //Get the params
            sscanf(lLine,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",&lParams[0],&lParams[1],&lParams[2],&lParams[3],&lParams[4],&lParams[5],&lParams[6],&lParams[7],&lParams[8],&lParams[9],&lParams[10]);
            //Pathogen creation
            Pathogen* lTemp=new Pathogen(lParams[0],lParams[1],lParams[2],lParams[3],lParams[4],lParams[5],lParams[6],lParams[7],this);
            //Add to the lake with the given viral load
            Pathogen* lTemp1=getPathogenAdd(lTemp);
            theLake->addPathogen(lTemp1,lParams[8]);
            //Add to the list of pathogen which could be added initially to individuals
            pathogenIntrod->add(lTemp1);
            
            float* lTempFloat=new float;
            *lTempFloat=lParams[9];
            lIntrod->add(lTempFloat);
            
            float* lTempFloatRecov=new float;
            *lTempFloatRecov=lParams[10];
            lIntrodRecov->add(lTempFloatRecov);
            
            delete lTemp;
            //delete lSequenceTemp;
            //delete lSequence;
        }
    }
    fclose(lFile);
}

/** Load trade-off parameters*/
/*    void Model::loadTradeOff(char* pFileName){
 char lLine[2048];
 FILE* lFile=fopen(pFileName,"rt");
 if(lFile!=NULL){
 while(fgets(lLine,2047,lFile)){
 if( (strstr(lLine,"//")==NULL) && (strstr(lLine,"#")==NULL)){
 char lName[128];
 char lValue[128];
 sscanf(lLine,"%s\t%s",lName,lValue);
 if(strstr(lName,"b")!=0){
 tradeOffValues[0]=atof(lValue);
 }
 if(strstr(lName,"c")!=0){tradeOffValues[1]=atof(lValue);}
 if(strstr(lName,"d")!=0){tradeOffValues[2]=atof(lValue);}
 if(strstr(lName,"e")!=0){tradeOffValues[3]=atof(lValue);}
 if(strstr(lName,"f")!=0){tradeOffValues[4]=atof(lValue);}
 if(strstr(lName,"g")!=0){tradeOffValues[5]=atof(lValue);}
 if(strstr(lName,"h")!=0){tradeOffValues[6]=atof(lValue);}
 if(strstr(lName,"i")!=0){tradeOffValues[7]=atof(lValue);}
 }
 }
 }
 fclose(lFile);
 }*/

/** Building objects: World and agents creation */
void Model::buildObjects ()
{
    loadPathogensSpecies();
    loadBirdsSpecies();
}

/** Add a new agent*/
void Model::addAgent(Agent* pAgent){
    Agent* lTemp=new Agent(pAgent,indexAgent);
    indexAgent++;
    agentListPerso->add(lTemp);
}

/** Remove an existing agent*/
void Model::delAgent(Agent* pAgent){
    delete pAgent;
    agentListPerso->remove((void*)pAgent);
}

void Model::immigration(){
    ListPerso lPathogens;
    ListPerso lViralLoad;
    getNbInfect(&lPathogens,&lViralLoad);
    
    if(lPathogens.length()>0){
        migrationLoads->removeAllPathogens();
        for(int lIndex=0;lIndex<lPathogens.length();lIndex++){
            migrationLoads->addPathogen((Pathogen*)lPathogens.getElement(lIndex),(*(float*)lViralLoad.getElement(lIndex)));
        }
    }
    if(intervalMigration>0){
        if(randomObject->getUnif()<(1/(intervalMigration/timeStep))){
            
            migrationLoads->getPathogens(&lPathogens,&lViralLoad);
            float* lProbaSum=new float[lPathogens.length()];
            float* lProba=new float[lPathogens.length()];
            float lSum=0;
            for(int lIndex=0;lIndex<lPathogens.length();lIndex++){
                lProba[lIndex]=(*(float*)lViralLoad.getElement(lIndex));
                lSum+=lProba[lIndex];
            }
            for(int lIndex=0;lIndex<lPathogens.length();lIndex++){
                if(lIndex==0){lProbaSum[lIndex]=lProba[lIndex]/lSum;}
                else{lProbaSum[lIndex]=lProbaSum[lIndex-1]+lProba[lIndex]/lSum;}
            }
            float lChoose=randomObject->getUnif();
            Pathogen* lPathogenChoose=0;
            for(int lIndex=0;lIndex<lPathogens.length();lIndex++){
                if(lChoose<lProbaSum[lIndex]){
                    lPathogenChoose=(Pathogen*)lPathogens.getElement(lIndex);
                    break;
                }
            }
            
            for(int lIndex=0;lIndex<agentListPerso->length();lIndex++){
                Agent* lAgent=(Agent*)agentListPerso->getElement(lIndex);
                if(lAgent->getCrossImmunity(lPathogenChoose)==1){
                    lAgent->addPathogen(lPathogenChoose,1,false);
                    lAgent->agentUpdate();
                    break;
                }
            }
            delete[] lProba;
            delete[] lProbaSum;
        }
    }
}

void Model::writeBuffer(char* pExt,char* pBuffer,char* pLine){
    if(strlen(pBuffer)+strlen(pLine)<MAX_LENGTH){
        strcat(pBuffer,pLine);
    }
    else{
        char lFileName[2048];
        sprintf(lFileName,"%s/output_All_%s_%s.dat",outputFolder,pExt,output);
        FILE* lFile=fopen(lFileName,"at");
        fwrite(pBuffer,1,strlen(pBuffer),lFile);
        fwrite(pLine,1,strlen(pLine),lFile);
        fclose(lFile);
        strcpy(pBuffer,"");
    }
}

/** Simulation behavior for one thread version*/
void Model::goNormal(){
    //First filling arrays for immigration
    ListPerso lPathogens;
    ListPerso lViralLoad;
    getNbInfect(&lPathogens,&lViralLoad);
    for(int lIndex=0;lIndex<lPathogens.length();lIndex++){
        migrationLoads->addPathogen((Pathogen*)lPathogens.getElement(lIndex),(*(float*)lViralLoad.getElement(lIndex)));
    }
    //Start the simulation
    bool writeAll=true;
    for(double t=0;t<tMax;t+=timeStep){
        currentTime=t;
        cout<<currentTime<<endl;
        
        immigration();
        
        writeAll=allDataWrite(writeAll);
        if((int)currentTime%(int)intervalWrite!=0){
            writeAll=true;
        }
        
        char lLine[2048];
        sprintf(lLine,"%e\t%d\n",t,agentListPerso->length());
        writeBuffer("Pop",dataOutputPop,lLine);
        
        for(int lIndex=0;lIndex<agentListPerso->length();lIndex++){
            ((Agent*)(agentListPerso->getElement(lIndex)))->agentStep(1,NULL,NULL,NULL,NULL);
        }
        infectiousInd->removeAllPathogens();
        infectiousIndWb->removeAllPathogens();
        
        //Birth and death events
        ListPerso lAdd;
        ListPerso lRemove;
        
        //Call update function
        for(int lIndex=0;lIndex<agentListPerso->length();lIndex++){
            Agent* lAgent=((Agent*)(agentListPerso->getElement(lIndex)));
            lAgent->agentUpdate();
            updateNbInfectious(((Agent*)(agentListPerso->getElement(lIndex))));
            
            if(lAgent->getbirthRate(false)>0 && lAgent->getdeathRate(false)>0){
                float lRate=1/lAgent->getbirthRate(false);
                float lTimeStep=gettimeStep();
                float lProba=1-exp(-lRate*lTimeStep);
                float seasonality=lAgent->getseasonality();
                lProba=lProba*(1+seasonality*cos((2*PI* ( ((float)getCurrentTime())/((float)365)) )));
                if(randomObject->getUnif()<lProba){
                    lAdd.add(lAgent);
                }
                lRate=1/(lAgent->getdeathRate(false));
                float lProbaDC=1-exp(-lRate*lTimeStep);
                if(randomObject->getUnif()<lProbaDC){
                    lRemove.add(lAgent);
                }
            }
        }
        //Update and write lake status
        theLake->update();
        for(int lIndex=0;lIndex<lAdd.length();lIndex++){
            addAgent((Agent*)lAdd.getElement(lIndex));
        }
        for(int lIndex=0;lIndex<lRemove.length();lIndex++){
            delAgent((Agent*)lRemove.getElement(lIndex));
        }
    }
}

/** Simulation behavior for one thread version*/
void Model::goOptimized(){
    float lRatioOptimisation=2;
    //First filling arrays for immigration
    ListPerso lPathogens;
    ListPerso lViralLoad;
    getNbInfect(&lPathogens,&lViralLoad);
    for(int lIndex=0;lIndex<lPathogens.length();lIndex++){
        migrationLoads->addPathogen((Pathogen*)lPathogens.getElement(lIndex),(*(float*)lViralLoad.getElement(lIndex)));
    }
    //Start the simulation
    bool writeAll=true;
    for(double t=0;t<tMax;t+=timeStep){
        currentTime=t;
        cout<<currentTime<<endl;
        
        immigration();
        
        writeAll=allDataWrite(writeAll);
        if((int)currentTime%(int)intervalWrite!=0){
            writeAll=true;
        }
        
        char lLine[2048];
        sprintf(lLine,"%e\t%d\n",t,agentListPerso->length());
        writeBuffer("Pop",dataOutputPop,lLine);
        
        //Getting list of number of infectious individuals for each pathogen
        ListPerso lPathogDirect;
        ListPerso lViralLoad;
        getNbInfect(&lPathogDirect,&lViralLoad);
        
        // Computing the average number of infectious individuals by direct transmission
        float lProba;
        long lNbExpected=0;
        long lPopSize=agentListPerso->length();
        if(((Agent*)agentListPerso->getElement(0))->getcontactRate()>0){
            for(int lIndex=0;lIndex<pathogenListPerso->length();lIndex++){
                Pathogen* lTemp=(Pathogen*)pathogenListPerso->getElement(lIndex);
                lProba=lTemp->getprobInfection()*(1+getampTransm()*sin((2*PI* ( ((float)getCurrentTime())/((float)365)) )));
                lNbExpected+=lProba*agentListPerso->length();
            }
        }
        
        // Computing the average number of infectious individuals by environmental transmission
        ListPerso lPathog;
        ListPerso lProbaInf;
        theLake->getProbabilityInfection(((Agent*)agentListPerso->getElement(0))->getdrinkingVolume(),&lPathog,&lProbaInf);
        ListPerso lIndexList;
        
        if(((Agent*)agentListPerso->getElement(0))->getdrinkingRate()>0){
            for(int lIndex=0;lIndex<lPathog.length();lIndex++){
                lProba=(*(float*)lProbaInf.getElement(lIndex));
                lNbExpected+=lProba*agentListPerso->length();
            }
        }
        
        // Computing the ratio of individuals with transmission
        lNbExpected=lNbExpected*lRatioOptimisation;
        double lRatio=(double)lNbExpected/(double)lPopSize;
        if(lRatio>1){
            lRatio=1;
        }
        //Calling the step function for each of these individuals
        double lDone=0;
        double* lIndiceToStep=new double[lNbExpected];
        for(int lIndex=0;((lIndex<agentListPerso->length())&&(lDone<lNbExpected-1));lIndex++){
            if(getRandom()->getUnif()<lRatio){
                lIndiceToStep[(int)lDone]=lIndex;
                lDone++;
            }
        }
        //BOUCLE A PARALLELISER
        for(int lIndex=0;lIndex<lDone;lIndex++){
            ((Agent*)(agentListPerso->getElement(lIndiceToStep[lIndex])))->agentStep(lRatio,&lPathog,&lProbaInf,&lPathogDirect,&lViralLoad);
        }
        delete lIndiceToStep;
        infectiousInd->removeAllPathogens();
        infectiousIndWb->removeAllPathogens();
        
        //Birth and death events
        ListPerso lAdd;
        ListPerso lRemove;
        
        Agent* lAgent=((Agent*)(agentListPerso->getElement(0)));
        float seasonality=lAgent->getseasonality();
        float lRate=1/lAgent->getbirthRate(false);
        float lTimeStep=gettimeStep();
        lProba=1-exp(-lRate*lTimeStep);
        float lProbaB=lProba*(1+seasonality*cos((2*PI* ( ((float)getCurrentTime())/((float)365)) )));
        double lNbExpectedBirths=lProbaB*lPopSize*lRatioOptimisation;
        lRate=1/lAgent->getdeathRate(false);
        lProba=1-exp(-lRate*lTimeStep);
        float lProbaD=lProba*(1+seasonality*cos((2*PI* ( ((float)getCurrentTime())/((float)365)) )));
        double lNbExpectedDeaths=lProbaD*lPopSize*lRatioOptimisation;
        
        float lRatioDC=(double)lNbExpectedDeaths/(double)lPopSize;
        float lRatioB=(double)lNbExpectedBirths/(double)lPopSize;
        double lDoneB=0;
        double lDoneDC=0;
        //Call update function
        for(int lIndex=0;lIndex<agentListPerso->length();lIndex++){
            Agent* lAgent=((Agent*)(agentListPerso->getElement(lIndex)));
            lAgent->agentUpdate();
            updateNbInfectious(((Agent*)(agentListPerso->getElement(lIndex))));
            
            if(lDoneB<lNbExpectedBirths-1){
                if(getRandom()->getUnif()<lRatioB){
                    if(randomObject->getUnif()<(lProbaB/lRatioB)){
                        lAdd.add(lAgent);
                        lDoneB++;
                    }
                }
            }
            if(lDoneDC<lNbExpectedDeaths-1){
                if(getRandom()->getUnif()<lRatioDC){
                    if(randomObject->getUnif()<(lProbaD/lRatioDC)){
                        lRemove.add(lAgent);
                        lDoneDC++;
                    }
                }
            }
        }
        //Update and write lake status
        theLake->update();
        for(int lIndex=0;lIndex<lAdd.length();lIndex++){
            addAgent((Agent*)lAdd.getElement(lIndex));
        }
        for(int lIndex=0;lIndex<lRemove.length();lIndex++){
            delAgent((Agent*)lRemove.getElement(lIndex));
        }
        for(int lIndex=0;lIndex<lProbaInf.length();lIndex++){
            delete (float*)lProbaInf.getElement(lIndex);
        }
    }
}

void Model::updateNbInfectious(Agent* pAgent){
    ListPerso* lTempList=pAgent->getcurrentPathogens();
    ListPerso* lTempListFlag=pAgent->getcurrentPathogensFlag();
    for(int lIndex=0;lIndex<lTempList->length();lIndex++){
        Pathogen* lTemp=(Pathogen*)lTempList->getElement(lIndex);
        infectiousInd->addPathogen(lTemp,1);
        float lTempFlag=(*(float*)lTempListFlag->getElement(lIndex));
        if(lTempFlag==1){
            infectiousIndWb->addPathogen(lTemp,1);
        }
        else{
            infectiousIndWb->addPathogen(lTemp,0);
        }
    }
}

/** Write all data about each agent*/
bool Model::allDataWrite(bool pWriteAll){
    bool lWrite=pWriteAll;
    //Updating number of infectious individuals
    if(outputAll!=0){
        if( (int)currentTime%(int)intervalWrite==0 && currentTime>1 && lWrite){
            char lLine[2048];
            ListPerso lListPathogens;
            ListPerso lListInf;
            ListPerso lListPathogensWb;
            ListPerso lListInfWb;
            infectiousInd->getPathogens(&lListPathogens,&lListInf);
            infectiousIndWb->getPathogens(&lListPathogensWb,&lListInfWb);
            Pathogen* lPathogen=(Pathogen*)lListPathogens.getElement(0);
            for(int lIndex=0;lIndex<lListPathogens.length();lIndex++){
                lPathogen=(Pathogen*)lListPathogens.getElement(lIndex);
                double lNbInf=*((float*)lListInf.getElement(lIndex));
                double lNbInfWb=*((float*)lListInfWb.getElement(lIndex));
                if(lNbInf>0){
                    double lSum=lPathogen->getSum();
                    sprintf(lLine,"%f\t%f\t%f\t%f\n",currentTime,lSum,lNbInf,lNbInfWb);
                    writeBuffer("Anti",dataOutputAllAnti,lLine);
                }
            }
        }
        lWrite=false;
    }
    return lWrite;
}

/** Convert sequence from an integer array to a string*/
void Model::convertSequence(int* pSequence,char* pOut,int pLength){
    for(int lIndex=0;lIndex<pLength;lIndex++){if(pSequence[lIndex]==0){pOut[lIndex]='0';}else{pOut[lIndex]='1';}}
    pOut[pLength]='\0';
}

/** Load simulation parameters*/
bool Model::loadParams(char* pFileName){
    bool lReturn=true;
    FILE* lFile=fopen(pFileName,"rt");
    if(lFile!=0){
        char lLine[2048];
        while(fgets(lLine,2047,lFile)){
            if( (strstr(lLine,"//")==NULL) && (strstr(lLine,"#")==NULL)){
                char lName[128];
                char lValue[128];
                sscanf(lLine,"%s\t%s",lName,lValue);
                if(strstr(lName,"tMax")!=0){tMax=atoi(lValue);}
                if(strcmp(lName,"outputFolder")==0){
                    int lLength=strlen(lLine)-strlen("outputFolder")-2;
                    strncpy(outputFolder,lLine+strlen("outputFolder")+1,strlen(lLine)-strlen("outputFolder")-2);
                    outputFolder[lLength]='\0';
                }
                if(strstr(lName,"birdsFile")!=0){
                    int lLength=strlen(lLine)-strlen("birdsFile")-2;
                    strncpy(birdsFile,lLine+strlen("birdsFile")+1,strlen(lLine)-strlen("birdsFile")-2);
                    birdsFile[lLength]='\0';
                }
                if(strstr(lName,"pathogensFile")!=0){
                    int lLength=strlen(lLine)-strlen("pathogensFile")-2;
                    strncpy(pathogensFile,lLine+strlen("pathogensFile")+1,strlen(lLine)-strlen("pathogensFile")-2);
                    pathogensFile[lLength]='\0';
                }
                if(strstr(lName,"randomSeed")!=0){randomSeed=atoi(lValue);}
                if(strstr(lName,"lakeVolume")!=0){lakeVolume=atoi(lValue);}
                if(strstr(lName,"outputAll")!=0){outputAll=atoi(lValue);}
                if(strstr(lName,"crossImmunity")!=0){crossImmunityFlag=atoi(lValue);}
                if(strstr(lName,"maxCrossImmunity")!=0){maxCrossImmunity=atof(lValue);}
                if(strstr(lName,"ampLake")!=0){ampLake=atof(lValue);}
                if(strstr(lName,"ampTransm")!=0){ampTransm=atof(lValue);}
                if(strstr(lName,"reassortment")!=0){reassortment=atof(lValue);}
                if(strstr(lName,"intervalMigration")!=0){intervalMigration=atoi(lValue);}
                if(strstr(lName,"immunity")!=0){immunity=atoi(lValue);}
                if(strstr(lName,"intervalWrite")!=0){intervalWrite=atoi(lValue);}
                if(strstr(lName,"interceptInfPeriod")!=0){interceptInfPeriod=atof(lValue);}
                if(strstr(lName,"modifyInfPeriod")!=0){modifyInfPeriod=atoi(lValue);}
                if(strstr(lName,"immunityPeriod")!=0){immmunityPeriod=atof(lValue);}
            }
        }
        fclose(lFile);
    }
    if(checkPreviousSimulations()){
        Init();
    }
    else{
        lReturn=false;
    }
    return lReturn;
}

void Model::getNbInfect(ListPerso* pPathogens,ListPerso* pViralLoad){
    infectiousInd->getPathogens(pPathogens,pViralLoad);
}

Pathogen* Model::getPathogenAdd(Pathogen* pPathogen){
    Pathogen* lAddress=0;
    for(int lIndex=0;((lIndex<pathogenListPerso->length())&&(lAddress==0));lIndex++){
        Pathogen* lPathogen=(Pathogen*)pathogenListPerso->getElement(lIndex);
        if(lPathogen->samePathogenSequence(pPathogen)){
            lAddress=lPathogen;
        }
    }
    if(lAddress==0){
        Pathogen* lPathogen=new Pathogen(pPathogen);
        lAddress=lPathogen;
        pathogenListPerso->add(lAddress);
    }
    return lAddress;
}

void Model::writeEvolEvent(Pathogen* pPathogen1,Pathogen* pPathogen2,Pathogen* pPathogen3){
    char lLine[2048];
    if(pPathogen3==0){
        //Mutation
        sprintf(lLine,"%f\t%f\t%f\t%d\n",currentTime,pPathogen1->getSum(),pPathogen2->getSum(),pPathogen1->getwaterInfect());
    }
    else{
        //Reassortment
        sprintf(lLine,"%f\t%f\t%f\t%f\t%d\t%d\n",currentTime,pPathogen1->getSum(),pPathogen2->getSum(),pPathogen3->getSum(),pPathogen1->getwaterInfect(),pPathogen2->getwaterInfect());
    }
    writeBuffer("Events",dataOutputEvents,lLine);
}

void Model::writeCoinf(double* pPathogen,int pLength){
    char lLine[2048];
    sprintf(lLine,"%f\t",currentTime);
    for(int lIndex=0;lIndex<pLength;lIndex++){
        sprintf(lLine,"%s\t%f",lLine,pPathogen[lIndex]);
    }
    sprintf(lLine,"%s\n",lLine);
    writeBuffer("Coinf",dataOutputCoinf,lLine);
}

float Model::getimmmunityPeriod(){
    return immmunityPeriod;
}

/** Gettors*/
Lake* Model::getLake(){return theLake;}
Random* Model::getRandom(){return randomObject;}
bool Model::getCrossImmun(){return (bool)crossImmunityFlag;}
float Model::getmaxCrossImmunity(){return maxCrossImmunity;}
double Model::getCurrentTime(){return currentTime;}
float Model::getreassortment(){return reassortment;}
float Model::gettimeStep(){return timeStep;}
int Model::getPopSize(){return agentListPerso->length();}
int Model::getimmunity(){return immunity;}
double Model::getIntercept(){return interceptInfPeriod;}
int Model::getmodifyInfPeriod(){return modifyInfPeriod;}
float Model::getampTransm(){return ampTransm;}
