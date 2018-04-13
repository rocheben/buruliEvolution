//
//  main.cpp
//  buruliEvolution
//
//  Created by Benjamin ROCHE on 13/04/2018.
//  Copyright © 2018 Benjamin ROCHE. All rights reserved.
//
#include <iostream>
#include <stdio.h>
#include <string.h>

using namespace std;

#include "model.hpp"

int main(int argc,char* argv[])
{
    cout << "Individual-based model for avian flu multi strains" << endl;
    Model* lModel;
    char lParamFile[256];
    char lSuffixe[256];
    float lTimeStep=1;
    //If paramfile, timestep and outputfile are filled
    if(argc>=4){
        strcpy(lSuffixe,argv[3]);
        lTimeStep=atof(argv[2]);
        strcpy(lParamFile,argv[1]);
        
        //Build model object
        lModel=new Model(lSuffixe,lTimeStep);
        //Load parameters
        if(lModel->loadParams(lParamFile)){
            //Build hosts
            lModel->buildObjects();
            if(argc==5){
                if(atof(argv[4])==1){
                    cout<<"OPTIMIZED: Loading parameter file : "<<lParamFile<<" goes to "<<lSuffixe<<endl;
                    //Launch the simulation
                    lModel->goOptimized();
                }
                else{
                    cout<<"Loading parameter file : "<<lParamFile<<" goes to "<<lSuffixe<<endl;
                    //Launch the simulation
                    lModel->goOptimized();
                }
            }
            else{
                cout<<"Loading parameter file : "<<lParamFile<<" goes to "<<lSuffixe<<endl;
                //Launch the simulation
                lModel->goOptimized();
            }
            delete lModel;
        }
        else{
            cout<<"Simulations already done... Skipping"<<endl;
        }
    }
    else{
        cout<<"Syntax: "<<argv[0]<<" paramFile timeStep suffix"<<endl;
    }
    return 0;
}
