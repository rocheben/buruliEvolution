//
//  list.hpp
//  buruliEvolution
//
//  Created by Benjamin ROCHE on 13/04/2018.
//  Copyright © 2018 Benjamin ROCHE. All rights reserved.
//
#ifndef LIST
#define LIST
#define PRE_ALLOC 10
class ListPerso{
private:
    long preAlloc;
    void** vector;
    int lengthAllocated;
    int lengthUsed;
    
    void extendVector(){
        void** lNewVect=new void*[lengthUsed+preAlloc];
        for(int lIndex=0;lIndex<lengthUsed;lIndex++){lNewVect[lIndex]=vector[lIndex];}
        delete vector;
        vector=lNewVect;
        lengthAllocated=lengthUsed+preAlloc;
    }
    
public:
    
    ListPerso(){
        preAlloc=PRE_ALLOC;
        vector=new void*[preAlloc];
        lengthAllocated=preAlloc;
        lengthUsed=0;
    }
    
    ListPerso(double pNumber){
        preAlloc=pNumber;
        vector=new void*[preAlloc];
        lengthAllocated=preAlloc;
        lengthUsed=0;
    }
    
    ~ListPerso(){
        delete[] vector;
    }
    
    /*void flush(){
     void** lNewVector=new void*[lengthAllocated];
     int lIndexEcriture=0;
     int lNewLengthUsed=0;
     for(int lIndex=0;lIndex<lengthUsed;lIndex++){
     if(vector[lIndex]!=0){
     lNewVector[lIndexEcriture]=vector[lIndex];
     lIndexEcriture++;
     lNewLengthUsed++;
     }
     }
     void** oldVector=vector;
     vector=lNewVector;
     delete[] oldVector;
     vector=lNewVector;
     if(lNewLengthUsed==0){
     int tata=1;
     }
     lengthUsed=lNewLengthUsed;
     }*/
    
    int length(){return lengthUsed;}
    
    void* getElement(int pElement){return vector[pElement];}
    
    void add(void* pParam){
        if(lengthUsed==lengthAllocated){extendVector();}
        vector[lengthUsed]=pParam;
        lengthUsed++;
    }
    
    void addAll(ListPerso* pParam){
        while(lengthUsed+pParam->length()>lengthAllocated){extendVector();}
        for(int lIndex=0;lIndex<pParam->length();lIndex++){
            vector[lengthUsed]=pParam->getElement(lIndex);
            lengthUsed++;
        }
    }
    
    void removeAll(ListPerso* pParam){
        int lElemFound=0;
        int lOffset=0;
        for(int lIndexNew=0;lIndexNew<lengthUsed;lIndexNew++){
            bool lFound=false;
            for(int lIndex=lElemFound;((lIndex<pParam->length())&&(!lFound));lIndex++){
                if(vector[lIndexNew]==pParam->getElement(lIndex)){
                    lFound=true;
                    lElemFound++;
                }
            }
            if(lFound){lOffset++;}
            void* prevVal=vector[lIndexNew];
            vector[lIndexNew]=vector[lIndexNew+lOffset];
            if(lFound){lIndexNew--;lengthUsed--;}
            if(vector[lIndexNew]!=prevVal && !lFound){lIndexNew--;}
        }
    }
    
    void removeAll(){
        lengthUsed=0;
        lengthAllocated=preAlloc;
    }
    
    void* remove(void* pParam){
        void* lReturn=0;
        int lFound=-1;
        for(int lIndex=0;((lIndex<lengthUsed)&&(lFound==-1));lIndex++){
            if(vector[lIndex]==pParam){
                lReturn=vector[lIndex];
                lFound=lIndex;
                break;
            }
        }
        lengthUsed--;
        for(int lIndex=lFound;lIndex<lengthUsed;lIndex++){
            vector[lIndex]=vector[lIndex+1];
        }
        return lReturn;
    }
    
    /*void remove(float pParam){
     int lFound=-1;
     for(int lIndex=0;((lIndex<lengthUsed)&&(lFound==-1));lIndex++){
     float lValue=*(float*)vector[lIndex];
     if(lValue==pParam){
     lFound=lIndex;
     break;
     }
     }
     lengthUsed--;
     for(int lIndex=lFound;lIndex<lengthUsed;lIndex++){
     vector[lIndex]=vector[lIndex+1];
     }
     }*/
    
    void setElement(int pPlace,void* pParam){
        vector[pPlace]=pParam;
    }
    
    void setElement(int pPlace,float pParam){
        /*float* lTemp=new float;
         *lTemp=pParam;*/
        *(float*)vector[pPlace]=pParam;
    }
};
#endif
