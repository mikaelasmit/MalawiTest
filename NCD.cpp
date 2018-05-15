/////////////////////////////////////////////////////////////////
//    NCD.cpp                                                  //
//    Created by Mikaela Smit on 25/01/2018.				   //
//    Copyright (c) 2014 Mikaela Smit. All rights reserved.    //
//    These event relate to NCDs.                              //
/////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <ctime>
#include <math.h>
#include <vector>
#include <random>
#include <cmath>

#include "eventfunctions.h"
#include "event.h"
#include "eventQ.h"
#include "person.h"
#include "errorcoutmacro.h"
#include "CParamReader.hpp"
#include "LoadParams.h"
#include "CountryParams.hpp"
#include "Intervention.hpp"
#include "HPVInfection.hpp"
#include "NCD.hpp"

using namespace std;

// hard code to remove
// 0. Move NCD risk array snad mortality array to country param
// 1. ncd_index
// 2. CauseOfDeath
// 3. p_mort



//// --- OUTSIDE INFORMATION --- ////
extern double *p_GT;                            // Time related
extern int EndYear;

extern person** MyArrayOfPointersToPeople;      // Event related
extern priority_queue<event*, vector<event*>, timeComparison> *p_PQ;
extern vector<event*> Events;

extern double**  NCDArray;                      // Arrays
extern double** DeathArray_Women;
extern double** DeathArray_Men;

extern double MortRisk[6];                      // Mortality related
extern double MortRisk_Cancer[5];
extern double MortAdj;

extern double Risk_DiabHC;                      // NCD common causal pathways
extern double Risk_DiabHT;
extern double Risk_DiabCKD;
extern double Risk_DiabCVD;

extern double Risk_HCHT;
extern double Risk_HCCVD;

extern double Risk_HTCKD;
extern double Risk_HTCVD;

extern double Risk_NCD_Diabetes[5];
extern int relatedNCDs_Diab[5];
extern int nr_NCD_Diab;

extern double Risk_NCD_HT[3];
extern int relatedNCDs_HT[3];
extern int nr_NCD_HT;

extern double Risk_NCD_HC[3];
extern int relatedNCDs_HC[3];
extern int nr_NCD_HC;



//// --- Important Functions --- ////
int RandomMinMax_2(int min, int max){							// Random number generator between min and max number
    return rand()%(max-min+1)+min;							    // Note: if min=0 and max=4 it will generate 0,1,2,3,4
}

double randfrom_3(double min, double max){
    double range = (max - min);
    double div = RAND_MAX / range;
    return min + (rand() / div);
}




//// --- NCD EVENTS --- ////
void EventMyDepressionDate(person *MyPointerToPerson)
{
    E(cout << endl << endl << "This patient is going to develop depression!" << endl;)
    
    if(MyPointerToPerson->Alive==1 && MyPointerToPerson->Depression_status==0) {
        
        // 1. Update status
        MyPointerToPerson->Depression_status=1;
        
        // 2. Lets see if we need to update death date
        int ncd_index=1;
        int k=(MyPointerToPerson->DoB-1800);
        int j=0;
        double	d = ((double) rand() / (RAND_MAX)) ;
        double TestDeathDate=-998;
        
        if (MyPointerToPerson->Sex==1){
            while(d>DeathArray_Men[k][j]*MortRisk[ncd_index] && j<121){j++;}
            TestDeathDate=(MyPointerToPerson->DoB+j);
        }
        
        if (MyPointerToPerson->Sex==2) {
            while(d>DeathArray_Women[k][j]*MortRisk[ncd_index] && j<121){j++;}
            TestDeathDate=(MyPointerToPerson->DoB+j);
        }
        
        // 3. Feed into EventQ if applicable
        if (TestDeathDate<MyPointerToPerson->DateOfDeath && TestDeathDate>*p_GT){
            
            MyPointerToPerson->DateOfDeath=TestDeathDate;
            
            if (MyPointerToPerson->DateOfDeath<EndYear){
                int p=MyPointerToPerson->PersonID-1;
                event * DeathEvent = new event;
                Events.push_back(DeathEvent);
                DeathEvent->time = MyPointerToPerson->DateOfDeath;
                DeathEvent->p_fun = &EventMyDeathDate;
                DeathEvent->person_ID = MyArrayOfPointersToPeople[p];
                p_PQ->push(DeathEvent);
                
                MyPointerToPerson->CauseOfDeath=4;                          // Update cause of death
            }
        }
    }
    E(cout << endl << endl << "This patient has developed depression!" << endl;)
}

void EventMyAsthmaDate(person *MyPointerToPerson)
{
    E(cout << endl << endl << "This patient is going to develop asthma!" << endl;)
    
    if(MyPointerToPerson->Alive==1 && MyPointerToPerson->Asthma_status==0) {
        
        // 1. Update status
        MyPointerToPerson->Asthma_status=1;
        
        // 2. Lets see if they die earlier
        int ncd_index=2;
        int k=(MyPointerToPerson->DoB-1800);
        int j=0;
        double	d = ((double) rand() / (RAND_MAX)) ;
        double TestDeathDate=-998;
        
        if (MyPointerToPerson->Sex==1){
            while(d>DeathArray_Men[k][j]*MortRisk[ncd_index] && j<121){j++;}
            TestDeathDate=(MyPointerToPerson->DoB+j);
        }
        
        if (MyPointerToPerson->Sex==2) {
            while(d>DeathArray_Women[k][j]*MortRisk[ncd_index] && j<121){j++;}
            TestDeathDate=(MyPointerToPerson->DoB+j);
        }
        
        // 3. Feed into EventQ if applicable
        if (TestDeathDate<MyPointerToPerson->DateOfDeath && TestDeathDate>*p_GT){
            
            MyPointerToPerson->DateOfDeath=TestDeathDate;
            
            if (MyPointerToPerson->DateOfDeath<EndYear){
                int p=MyPointerToPerson->PersonID-1;
                event * DeathEvent = new event;
                Events.push_back(DeathEvent);
                DeathEvent->time = MyPointerToPerson->DateOfDeath;
                DeathEvent->p_fun = &EventMyDeathDate;
                DeathEvent->person_ID = MyArrayOfPointersToPeople[p];
                p_PQ->push(DeathEvent);
             
                MyPointerToPerson->CauseOfDeath=5;                          // Update cause of death
            }
        }
    }
    E(cout << endl << endl << "This patient has developed asthma!" << endl;)
}

void EventMyStrokeDate(person *MyPointerToPerson)
{
    E(cout << endl << endl << "This patient is going to hava a stroke!" << endl;)
    
    if(MyPointerToPerson->Alive==1 && MyPointerToPerson->Stroke_status==0) {
        
        // 1. Update status
        MyPointerToPerson->Stroke_status=1;
        
        // 2. Lets see if they die earlier
        int ncd_index=3;
        int k=(MyPointerToPerson->DoB-1800);
        int j=0;
        double	d = ((double) rand() / (RAND_MAX)) ;
        double TestDeathDate=-998;
        double p_mort=((double) rand () / (RAND_MAX));

        if (p_mort<0.6) {                                                    // See if they die within the next year
            double YearFraction=(RandomMinMax_2(1,12))/12.1;                //REFERENCE: http://stroke.ahajournals.org/content/32/9/2131.full
            TestDeathDate=*p_GT + YearFraction;
        }
        
        if (p_mort>=0.6){
            if (MyPointerToPerson->Sex==1){
                while(d>DeathArray_Men[k][j]*MortRisk[ncd_index] && j<121){j++;}
                TestDeathDate=(MyPointerToPerson->DoB+j);
            }
            
            if (MyPointerToPerson->Sex==2) {
                while(d>DeathArray_Women[k][j]*MortRisk[ncd_index] && j<121){j++;}
                TestDeathDate=(MyPointerToPerson->DoB+j);
            }
        }
        
        // 3. Feed into EventQ if applicable
        if (TestDeathDate<MyPointerToPerson->DateOfDeath && TestDeathDate>*p_GT){
            
            MyPointerToPerson->DateOfDeath=TestDeathDate;
            
            if (MyPointerToPerson->DateOfDeath<EndYear){
                int p=MyPointerToPerson->PersonID-1;
                event * DeathEvent = new event;
                Events.push_back(DeathEvent);
                DeathEvent->time = MyPointerToPerson->DateOfDeath;
                DeathEvent->p_fun = &EventMyDeathDate;
                DeathEvent->person_ID = MyArrayOfPointersToPeople[p];
                p_PQ->push(DeathEvent);
                
                MyPointerToPerson->CauseOfDeath=6;                              // Update cause of death
            }
        }
    }
    
    E(cout << endl << endl << "This patient has had a stroke!" << endl;)
}

void EventMyMIDate(person *MyPointerToPerson)
{
    E(cout << endl << endl << "This patient is going to have an MI!" << endl;)
    
    if(MyPointerToPerson->Alive==1 && MyPointerToPerson->MI_status==0) {
        
        // 1. Update stats
        MyPointerToPerson->MI_status=1;
        
        // 2. Lets see if they die earlier
        int ncd_index=6;
        int k=(MyPointerToPerson->DoB-1800);
        int j=0;
        double	d = ((double) rand() / (RAND_MAX)) ;
        double TestDeathDate=-998;
        double p_mort=((double) rand () / (RAND_MAX));
        
        if (p_mort<0.6){                                                    // REFERENCE: http://stroke.ahajournals.org/content/32/9/2131.full
            double YearFraction=(RandomMinMax_2(1,12))/12.1;                // See when during this year this person will die
            TestDeathDate=*p_GT + YearFraction;
        }
        
        if (p_mort>=0.6){
            if (MyPointerToPerson->Sex==1){
                while(d>DeathArray_Men[k][j]*MortRisk[ncd_index] && j<121){j++;}
                TestDeathDate=(MyPointerToPerson->DoB+j);
            }
            
            if (MyPointerToPerson->Sex==2) {
                while(d>DeathArray_Women[k][j]*MortRisk[ncd_index] && j<121){j++;}
                TestDeathDate=(MyPointerToPerson->DoB+j);
            }
        }
        
        // 3. Feed into EventQ if applicable
        if (TestDeathDate<MyPointerToPerson->DateOfDeath && TestDeathDate>*p_GT){
            
            MyPointerToPerson->DateOfDeath=TestDeathDate;
            
            if (MyPointerToPerson->DateOfDeath<EndYear){
                int p=MyPointerToPerson->PersonID-1;
                event * DeathEvent = new event;
                Events.push_back(DeathEvent);
                DeathEvent->time = MyPointerToPerson->DateOfDeath;
                DeathEvent->p_fun = &EventMyDeathDate;
                DeathEvent->person_ID = MyArrayOfPointersToPeople[p];
                p_PQ->push(DeathEvent);
                
                MyPointerToPerson->CauseOfDeath=6;                           // Update cause of death
            }
        }
    }
    E(cout << endl << endl << "This patient has had an MI!" << endl;)
}

void EventMyDiabetesDate(person *MyPointerToPerson){
    
    E(cout << "This patient is going to develop diabetes, let's see if they are at increased risk for other NCDs!" << endl;)
    
    
    if (MyPointerToPerson->Alive==1 && MyPointerToPerson->Diabetes_status==0)
    {
        // 1. Update status
        MyPointerToPerson->Diabetes_status=1;
        
        // 2. Lets see if they dies ealier
        int ncd_index=4;
        int k=(MyPointerToPerson->DoB-1800);
        int j=0;
        double	d = ((double) rand() / (RAND_MAX)) ;
        double TestDeathDate=-998;
        
        if (MyPointerToPerson->Sex==1){
            while(d>DeathArray_Men[k][j]*MortRisk[ncd_index] && j<121){j++;}
            TestDeathDate=(MyPointerToPerson->DoB+j);
        }
        
        if (MyPointerToPerson->Sex==2) {
            while(d>DeathArray_Women[k][j]*MortRisk[ncd_index] && j<121){j++;}
            TestDeathDate=(MyPointerToPerson->DoB+j);
        }
        
        // 3. Feed into EventQ if applicable
        if (TestDeathDate<MyPointerToPerson->DateOfDeath && TestDeathDate>*p_GT){
            
            MyPointerToPerson->DateOfDeath=TestDeathDate;
            
            if (MyPointerToPerson->DateOfDeath<EndYear){
                int p=MyPointerToPerson->PersonID-1;
                event * DeathEvent = new event;
                Events.push_back(DeathEvent);
                DeathEvent->time = MyPointerToPerson->DateOfDeath;
                DeathEvent->p_fun = &EventMyDeathDate;
                DeathEvent->person_ID = MyArrayOfPointersToPeople[p];
                p_PQ->push(DeathEvent);
                
                MyPointerToPerson->CauseOfDeath=7;                          // Update cause of death
            }
        }
        
        // 4. Re-evaluate risk of related NCDs
        int ncd_nr=0;
        int age_index=(*p_GT-MyPointerToPerson->DoB);
        
        while (ncd_nr<nr_NCD_Diab)
        {
            double DateNCD=-997;
            double r = randfrom_3(NCDArray[relatedNCDs_Diab[ncd_nr]][age_index]*Risk_NCD_Diabetes[ncd_nr] ,1*Risk_NCD_Diabetes[ncd_nr] );
            
            // 4.1 // If we are getting an NCD lets get the age and date of NCD
            if (r<NCDArray[relatedNCDs_Diab[ncd_nr]][120]*Risk_NCD_Diabetes[ncd_nr])
            {
                int i=0;
                while (r>NCDArray[relatedNCDs_Diab[ncd_nr]][i]*Risk_NCD_Diabetes[ncd_nr]){i++;}
                double YearFraction=(RandomMinMax_2(1,12))/12.1;
                DateNCD=MyPointerToPerson->DoB+i+YearFraction;
            }
            
            // 4.2 Lets see if this pushed forward the existing NCD date
            if (DateNCD>=*p_GT && DateNCD<MyPointerToPerson->NCD_DatesVector.at(relatedNCDs_Diab[ncd_nr]))                  {
                
                MyPointerToPerson->NCD_DatesVector.at(relatedNCDs_Diab[ncd_nr])=DateNCD;
                
                if (ncd_nr==0)
                {
                    MyPointerToPerson->HT=DateNCD;
                    
                    int p=MyPointerToPerson->PersonID-1;
                    event * HTEvent = new event;
                    Events.push_back(HTEvent);
                    HTEvent->time = MyPointerToPerson->HT;
                    HTEvent->p_fun = &EventMyHyptenDate;
                    HTEvent->person_ID = MyArrayOfPointersToPeople[p];
                    p_PQ->push(HTEvent);
                }
                
                if (ncd_nr==1)
                {
                    MyPointerToPerson->Stroke=DateNCD;
                    
                    int p=MyPointerToPerson->PersonID-1;
                    event * StrokeEvent = new event;
                    Events.push_back(StrokeEvent);
                    StrokeEvent->time = MyPointerToPerson->Stroke;
                    StrokeEvent->p_fun = &EventMyStrokeDate;
                    StrokeEvent->person_ID = MyArrayOfPointersToPeople[p];
                    p_PQ->push(StrokeEvent);
                }
                
                if (ncd_nr==2)
                {
                    MyPointerToPerson->CKD=DateNCD;
                    
                    int p=MyPointerToPerson->PersonID-1;
                    event * CKDEvent = new event;
                    Events.push_back(CKDEvent);
                    CKDEvent->time = MyPointerToPerson->CKD;
                    CKDEvent->p_fun = &EventMyCKDDate;
                    CKDEvent->person_ID = MyArrayOfPointersToPeople[p];
                    p_PQ->push(CKDEvent);
                }
                
                if (ncd_nr==3)
                {
                    MyPointerToPerson->MI=DateNCD;
                    
                    int p=MyPointerToPerson->PersonID-1;
                    event * MIEvent = new event;
                    Events.push_back(MIEvent);
                    MIEvent->time = MyPointerToPerson->MI;
                    MIEvent->p_fun = &EventMyMIDate;
                    MIEvent->person_ID = MyArrayOfPointersToPeople[p];
                    p_PQ->push(MIEvent);
                }
                
                if (ncd_nr==4)
                {
                    MyPointerToPerson->HC=DateNCD;
                    
                    int p=MyPointerToPerson->PersonID-1;
                    event * HCEvent = new event;
                    Events.push_back(HCEvent);
                    HCEvent->time = MyPointerToPerson->HC;
                    HCEvent->p_fun = &EventMyHypcholDate;
                    HCEvent->person_ID = MyArrayOfPointersToPeople[p];
                    p_PQ->push(HCEvent);
                }
            }
            ncd_nr++;
        }
    }
    E(cout << endl << endl << "This patient developed diabetes and related NCD risk has been evaluated!" << endl;)
}

void EventMyHyptenDate(person *MyPointerToPerson)
{
    E(cout << "This patient is going to develop HT, let's see if they are at increased risk for other NCDs" << endl;)
    
    if (MyPointerToPerson->Alive==1 && MyPointerToPerson->HT_status==0)
    {
        // 1. Update status
        MyPointerToPerson->HT_status=1;
        
        // 2. Lets see if they die ealier
        int ncd_index=0;
        int k=(MyPointerToPerson->DoB-1800);
        int j=0;
        double	d = ((double) rand() / (RAND_MAX)) ;
        double TestDeathDate=-998;
        
        if (MyPointerToPerson->Sex==1){
            while(d>DeathArray_Men[k][j]*MortRisk[ncd_index] && j<121){j++;}
            TestDeathDate=(MyPointerToPerson->DoB+j);
        }
        
        if (MyPointerToPerson->Sex==2) {
            while(d>DeathArray_Women[k][j]*MortRisk[ncd_index] && j<121){j++;}
            TestDeathDate=(MyPointerToPerson->DoB+j);
        }
        
        // 3. Feed into EventQ if applicable
        if (TestDeathDate<MyPointerToPerson->DateOfDeath && TestDeathDate>*p_GT){
            
            MyPointerToPerson->DateOfDeath=TestDeathDate;
            
            if (MyPointerToPerson->DateOfDeath<EndYear){
                int p=MyPointerToPerson->PersonID-1;
                event * DeathEvent = new event;
                Events.push_back(DeathEvent);
                DeathEvent->time = MyPointerToPerson->DateOfDeath;
                DeathEvent->p_fun = &EventMyDeathDate;
                DeathEvent->person_ID = MyArrayOfPointersToPeople[p];
                p_PQ->push(DeathEvent);
                
                MyPointerToPerson->CauseOfDeath=3;                  // Update cause of death
            }
        }
        
        // 4. Re-evaluate risk of related NCDs
        int ncd_nr=0;
        double DateNCD=-997;
       
        while (ncd_nr<nr_NCD_HT)
        {
            double r = ((double) rand() / (RAND_MAX));
            
            // 4.1 // If we are getting an NCD lets get the age and date of NCD
            if (r<NCDArray[relatedNCDs_HT[ncd_nr]][120]*Risk_NCD_HT[ncd_nr])
            {
                int i=0;
                while (r>NCDArray[relatedNCDs_HT[ncd_nr]][i]*Risk_NCD_HT[ncd_nr]){i++;}
                double YearFraction=(RandomMinMax_2(1,12))/12.1;
                DateNCD=MyPointerToPerson->DoB+i+YearFraction;
            }
            
            // 4.2 Lets see if this pushed forward the existing NCD date
            if (DateNCD>=*p_GT && DateNCD<MyPointerToPerson->NCD_DatesVector.at(relatedNCDs_HT[ncd_nr]))
            {
                
                MyPointerToPerson->NCD_DatesVector.at(relatedNCDs_HT[ncd_nr])=DateNCD;
                
                if (ncd_nr==0)
                {
                    MyPointerToPerson->Stroke=DateNCD;
                    
                    int p=MyPointerToPerson->PersonID-1;
                    event * StrokeEvent = new event;
                    Events.push_back(StrokeEvent);
                    StrokeEvent->time = MyPointerToPerson->Stroke;
                    StrokeEvent->p_fun = &EventMyStrokeDate;
                    StrokeEvent->person_ID = MyArrayOfPointersToPeople[p];
                    p_PQ->push(StrokeEvent);
                }
                
                if (ncd_nr==1)
                {
                    MyPointerToPerson->CKD=DateNCD;
                    
                    int p=MyPointerToPerson->PersonID-1;
                    event * CKDEvent = new event;
                    Events.push_back(CKDEvent);
                    CKDEvent->time = MyPointerToPerson->CKD;
                    CKDEvent->p_fun = &EventMyCKDDate;
                    CKDEvent->person_ID = MyArrayOfPointersToPeople[p];
                    p_PQ->push(CKDEvent);
                }
                
                if (ncd_nr==2)
                {
                    MyPointerToPerson->MI=DateNCD;
                    
                    int p=MyPointerToPerson->PersonID-1;
                    event * MIEvent = new event;
                    Events.push_back(MIEvent);
                    MIEvent->time = MyPointerToPerson->MI;
                    MIEvent->p_fun = &EventMyMIDate;
                    MIEvent->person_ID = MyArrayOfPointersToPeople[p];
                    p_PQ->push(MIEvent);
                }
            }
            ncd_nr++;
        }
    }
    E(cout << endl << endl << "This patient developed HT and related NCD risk has been evaluated" << endl;)
}

void EventMyHypcholDate(person *MyPointerToPerson)
{
    E(cout << "This patient is going to develop HC, let's see if they are at increased risk for other NCD" << endl;)
    
    if (MyPointerToPerson->Alive==1 && MyPointerToPerson->HC_status==0)
    {
        // 1. Update status
        MyPointerToPerson->HC_status=1;
        
        // 2. Lets see if they die ealier
        int ncd_index=0;
        int k=(MyPointerToPerson->DoB-1800);
        int j=0;
        double	d = ((double) rand() / (RAND_MAX)) ;
        double TestDeathDate=-998;
        
        if (MyPointerToPerson->Sex==1){
            while(d>DeathArray_Men[k][j]*MortRisk[ncd_index] && j<121){j++;}
            TestDeathDate=(MyPointerToPerson->DoB+j);
        }
        
        if (MyPointerToPerson->Sex==2) {
            while(d>DeathArray_Women[k][j]*MortRisk[ncd_index] && j<121){j++;}
            TestDeathDate=(MyPointerToPerson->DoB+j);
        }
        
        // 3. Feed into the EventQ if applicable
        if (TestDeathDate<MyPointerToPerson->DateOfDeath && TestDeathDate>*p_GT){
            
            MyPointerToPerson->DateOfDeath=TestDeathDate;
            
            if (MyPointerToPerson->DateOfDeath<EndYear){
                int p=MyPointerToPerson->PersonID-1;
                event * DeathEvent = new event;
                Events.push_back(DeathEvent);
                DeathEvent->time = MyPointerToPerson->DateOfDeath;
                DeathEvent->p_fun = &EventMyDeathDate;
                DeathEvent->person_ID = MyArrayOfPointersToPeople[p];
                p_PQ->push(DeathEvent);
                
                MyPointerToPerson->CauseOfDeath=3;                  // Update cause of death
            }
        }
        
        // 4. Re-evaluate risk of related NCDs
        int ncd_nr=0;
        double DateNCD=-997;
        
        while (ncd_nr<nr_NCD_HC)
        {
            double r = ((double) rand() / (RAND_MAX));
            
            // 4.1 If we are getting an NCD lets get the age and date of NCD
            if (r<NCDArray[relatedNCDs_HC[ncd_nr]][120]*Risk_NCD_HC[ncd_nr]){
                int i=0;
                while (r>NCDArray[relatedNCDs_HC[ncd_nr]][i]*Risk_NCD_HC[ncd_nr]){i++;}
                double YearFraction=(RandomMinMax_2(1,12))/12.1;
                DateNCD=MyPointerToPerson->DoB+i+YearFraction;
            }
            
            // 4.2 Lets see if this pushed forward the existing NCD date
            if (DateNCD>=*p_GT && DateNCD<MyPointerToPerson->NCD_DatesVector.at(relatedNCDs_HC[ncd_nr]))
            {
                MyPointerToPerson->NCD_DatesVector.at(relatedNCDs_HC[ncd_nr])=DateNCD;
                
                if (ncd_nr==0)
                {
                    MyPointerToPerson->HT=DateNCD;
                    
                    int p=MyPointerToPerson->PersonID-1;
                    event * HTEvent = new event;
                    Events.push_back(HTEvent);
                    HTEvent->time = MyPointerToPerson->HT;
                    HTEvent->p_fun = &EventMyHyptenDate;
                    HTEvent->person_ID = MyArrayOfPointersToPeople[p];
                    p_PQ->push(HTEvent);
                }
                
                if (ncd_nr==1)
                {
                    MyPointerToPerson->Stroke=DateNCD;
                    
                    int p=MyPointerToPerson->PersonID-1;
                    event * StrokeEvent = new event;
                    Events.push_back(StrokeEvent);
                    StrokeEvent->time = MyPointerToPerson->Stroke;
                    StrokeEvent->p_fun = &EventMyStrokeDate;
                    StrokeEvent->person_ID = MyArrayOfPointersToPeople[p];
                    p_PQ->push(StrokeEvent);
                }
                
                if (ncd_nr==2)
                {
                    MyPointerToPerson->MI=DateNCD;
                    
                    int p=MyPointerToPerson->PersonID-1;
                    event * MIEvent = new event;
                    Events.push_back(MIEvent);
                    MIEvent->time = MyPointerToPerson->MI;
                    MIEvent->p_fun = &EventMyMIDate;
                    MIEvent->person_ID = MyArrayOfPointersToPeople[p];
                    p_PQ->push(MIEvent);
                }
            }
            
            ncd_nr++;
        }
    }
    E(cout << endl << endl << "This patient developed HC and related NCD risk has been evaluated" << endl;)
}


void EventMyCKDDate (person *MyPointerToPerson)
{
    
    E(cout << endl << endl << "This patient is going to develop CKD!" << endl;)
    
    if (MyPointerToPerson->CKD_status==0 && MyPointerToPerson->Alive == 1){
        
       // 1. Update status
        MyPointerToPerson->CKD_status=1;
        
        // L2. ets see if death happens earlier
        int ncd_index=5;
        int k=(MyPointerToPerson->DoB-1800);
        int j=0;
        double	d = ((double) rand() / (RAND_MAX)) ;
        double TestDeathDate=-998;
        
        if (MyPointerToPerson->Sex==1){
            while(d>DeathArray_Men[k][j]*MortRisk[ncd_index] && j<121){j++;}
            TestDeathDate=(MyPointerToPerson->DoB+j);
        }
        
        if (MyPointerToPerson->Sex==2) {
            while(d>DeathArray_Women[k][j]*MortRisk[ncd_index] && j<121){j++;}
            TestDeathDate=(MyPointerToPerson->DoB+j);
        }
        
        // 3. Feed into the EventQ if applicable
        if (TestDeathDate<MyPointerToPerson->DateOfDeath && TestDeathDate>*p_GT) {
            
            MyPointerToPerson->DateOfDeath=TestDeathDate;
            
            if (MyPointerToPerson->DateOfDeath<EndYear){
                int p=MyPointerToPerson->PersonID-1;
                event * DeathEvent = new event;
                Events.push_back(DeathEvent);
                DeathEvent->time = MyPointerToPerson->DateOfDeath;
                DeathEvent->p_fun = &EventMyDeathDate;
                DeathEvent->person_ID = MyArrayOfPointersToPeople[p];
                p_PQ->push(DeathEvent);
                
                MyPointerToPerson->CauseOfDeath=8;          // Update cause of death
            }
        }
    }
    E(cout << endl << endl << "This patient has developed CKD" << endl;)
}


void EventMyBreastDate (person *MyPointerToPerson)			// Function executed when person develops hypertension
{
    
    E(cout << endl << endl << "This patient just developed breast cancer!" << endl;)
    
    if (MyPointerToPerson->Breast_status==0 && MyPointerToPerson->Alive == 1){
        
        MyPointerToPerson->Breast_status=1;
        
        // Lets see if we need to update death date
        int cancer_index=0;
        
        // Lets see if they die earlier
        int k=(MyPointerToPerson->DoB-1800);					// To find corresponding year of birth from mortality array
        int j=0;												// This will be matched to probability taken from random number generator
        double	d = ((double) rand() / (RAND_MAX)) ;			// get a random number to determine Life Expectancy
        double TestDeathDate;
        
        if (MyPointerToPerson->Sex==1){
            while(d>DeathArray_Men[k][j]*MortRisk_Cancer[cancer_index] && j<121){j++;}
            TestDeathDate=(MyPointerToPerson->DoB+j);
        }
        
        if (MyPointerToPerson->Sex==2) {
            while(d>DeathArray_Women[k][j]*MortRisk_Cancer[cancer_index] && j<121){j++;}
            TestDeathDate=(MyPointerToPerson->DoB+j);
        }
        
        if (TestDeathDate<MyPointerToPerson->DateOfDeath && TestDeathDate>*p_GT) {
            
            MyPointerToPerson->DateOfDeath=TestDeathDate;
            
            // 2. Lets feed death into the eventQ
            if (MyPointerToPerson->DateOfDeath<EndYear){
                int p=MyPointerToPerson->PersonID-1;
                event * DeathEvent = new event;
                Events.push_back(DeathEvent);
                DeathEvent->time = MyPointerToPerson->DateOfDeath;
                DeathEvent->p_fun = &EventMyDeathDate;
                DeathEvent->person_ID = MyArrayOfPointersToPeople[p];
                p_PQ->push(DeathEvent);
                
                // Update cause of death
                MyPointerToPerson->CauseOfDeath=9;
            }
        }
    }
}

void EventMyCervicalDate (person *MyPointerToPerson)			// Function executed when person develops hypertension
{
    
    E(cout << endl << endl << "This patient just developed cervical cancer!" << endl;)
    
    if (MyPointerToPerson->Cervical_status==0 && MyPointerToPerson->Alive == 1){
        
        MyPointerToPerson->Cervical_status=1;
        
        // Lets see if we need to update death date
        int cancer_index=1;
        
        // Lets see if they die earlier
        int k=(MyPointerToPerson->DoB-1800);					// To find corresponding year of birth from mortality array
        int j=0;												// This will be matched to probability taken from random number generator
        double	d = ((double) rand() / (RAND_MAX)) ;			// get a random number to determine Life Expectancy
        double TestDeathDate;
        
        if (MyPointerToPerson->Sex==1){
            while(d>DeathArray_Men[k][j]*MortRisk_Cancer[cancer_index] && j<121){j++;}
            TestDeathDate=(MyPointerToPerson->DoB+j);
        }
        
        if (MyPointerToPerson->Sex==2) {
            while(d>DeathArray_Women[k][j]*MortRisk_Cancer[cancer_index] && j<121){j++;}
            TestDeathDate=(MyPointerToPerson->DoB+j);
        }
        
        if (TestDeathDate<MyPointerToPerson->DateOfDeath && TestDeathDate>*p_GT) {
            
            MyPointerToPerson->DateOfDeath=TestDeathDate;
            
            // 2. Lets feed death into the eventQ
            if (MyPointerToPerson->DateOfDeath<EndYear){
                int p=MyPointerToPerson->PersonID-1;
                event * DeathEvent = new event;
                Events.push_back(DeathEvent);
                DeathEvent->time = MyPointerToPerson->DateOfDeath;
                DeathEvent->p_fun = &EventMyDeathDate;
                DeathEvent->person_ID = MyArrayOfPointersToPeople[p];
                p_PQ->push(DeathEvent);
                
                // Update cause of death
                MyPointerToPerson->CauseOfDeath=10;
            }
        }
    }
}

void EventMyColoDate (person *MyPointerToPerson)			// Function executed when person develops hypertension
{
    
    E(cout << endl << endl << "This patient just developed breast cancer!" << endl;)
    
    if (MyPointerToPerson->Colo_status==0 && MyPointerToPerson->Alive == 1){
        
        MyPointerToPerson->Colo_status=1;
        
        // Lets see if we need to update death date
        int cancer_index=2;
        
        // Lets see if they die earlier
        int k=(MyPointerToPerson->DoB-1800);					// To find corresponding year of birth from mortality array
        int j=0;												// This will be matched to probability taken from random number generator
        double	d = ((double) rand() / (RAND_MAX)) ;			// get a random number to determine Life Expectancy
        double TestDeathDate;
        
        if (MyPointerToPerson->Sex==1){
            while(d>DeathArray_Men[k][j]*MortRisk_Cancer[cancer_index] && j<121){j++;}
            TestDeathDate=(MyPointerToPerson->DoB+j);
        }
        
        if (MyPointerToPerson->Sex==2) {
            while(d>DeathArray_Women[k][j]*MortRisk_Cancer[cancer_index] && j<121){j++;}
            TestDeathDate=(MyPointerToPerson->DoB+j);
        }
        
        if (TestDeathDate<MyPointerToPerson->DateOfDeath && TestDeathDate>*p_GT) {
            
            MyPointerToPerson->DateOfDeath=TestDeathDate;
            
            // 2. Lets feed death into the eventQ
            if (MyPointerToPerson->DateOfDeath<EndYear){
                int p=MyPointerToPerson->PersonID-1;
                event * DeathEvent = new event;
                Events.push_back(DeathEvent);
                DeathEvent->time = MyPointerToPerson->DateOfDeath;
                DeathEvent->p_fun = &EventMyDeathDate;
                DeathEvent->person_ID = MyArrayOfPointersToPeople[p];
                p_PQ->push(DeathEvent);
                
                // Update cause of death
                MyPointerToPerson->CauseOfDeath=11;
            }
        }
    }
}



void EventMyLiverDate (person *MyPointerToPerson)			// Function executed when person develops hypertension
{
    
    E(cout << endl << endl << "This patient just developed Liver cancer!" << endl;)
    
    if (MyPointerToPerson->Liver_status==0 && MyPointerToPerson->Alive == 1){
        
        MyPointerToPerson->Liver_status=1;
        
        //cout << "We have Liver " << MyPointerToPerson->Liver << endl;
        
        // Lets see if we need to update death date
        int cancer_index=3;
        
        // Lets see if they die earlier
        int k=(MyPointerToPerson->DoB-1800);					// To find corresponding year of birth from mortality array
        int j=0;												// This will be matched to probability taken from random number generator
        double	d = ((double) rand() / (RAND_MAX)) ;			// get a random number to determine Life Expectancy
        double TestDeathDate;
        
        if (MyPointerToPerson->Sex==1){
            while(d>DeathArray_Men[k][j]*MortRisk_Cancer[cancer_index] && j<121){j++;}
            TestDeathDate=(MyPointerToPerson->DoB+j);
        }
        
        if (MyPointerToPerson->Sex==2) {
            while(d>DeathArray_Women[k][j]*MortRisk_Cancer[cancer_index] && j<121){j++;}
            TestDeathDate=(MyPointerToPerson->DoB+j);
        }
        
        
        if (TestDeathDate<MyPointerToPerson->DateOfDeath && TestDeathDate>*p_GT) {
            
            MyPointerToPerson->DateOfDeath=TestDeathDate;
            
            // 2. Lets feed death into the eventQ
            if (MyPointerToPerson->DateOfDeath<EndYear){
                int p=MyPointerToPerson->PersonID-1;
                event * DeathEvent = new event;
                Events.push_back(DeathEvent);
                DeathEvent->time = MyPointerToPerson->DateOfDeath;
                DeathEvent->p_fun = &EventMyDeathDate;
                DeathEvent->person_ID = MyArrayOfPointersToPeople[p];
                p_PQ->push(DeathEvent);
                
                // Update cause of death
                MyPointerToPerson->CauseOfDeath=12;
            }
        }
    }
}


void EventMyOesoDate (person *MyPointerToPerson)			// Function executed when person develops hypertension
{
    
    E(cout << endl << endl << "This patient just developed Oesophageal cancer!" << endl;)
    
    if (MyPointerToPerson->Oeso_status==0 && MyPointerToPerson->Alive == 1){
        
        MyPointerToPerson->Oeso_status=1;
        
        //cout << "We have Oeso " << MyPointerToPerson->Oeso << endl;
        
        // Lets see if we need to update death date
        int cancer_index=4;
        
        // Lets see if they die earlier
        int k=(MyPointerToPerson->DoB-1800);					// To find corresponding year of birth from mortality array
        int j=0;												// This will be matched to probability taken from random number generator
        double	d = ((double) rand() / (RAND_MAX)) ;			// get a random number to determine Life Expectancy
        double TestDeathDate;
        
        if (MyPointerToPerson->Sex==1){
            while(d>DeathArray_Men[k][j]*MortRisk_Cancer[cancer_index] && j<121){j++;}
            TestDeathDate=(MyPointerToPerson->DoB+j);
        }
        
        if (MyPointerToPerson->Sex==2) {
            while(d>DeathArray_Women[k][j]*MortRisk_Cancer[cancer_index] && j<121){j++;}
            TestDeathDate=(MyPointerToPerson->DoB+j);
        }
        
        if (TestDeathDate<MyPointerToPerson->DateOfDeath && TestDeathDate>*p_GT) {
            
            MyPointerToPerson->DateOfDeath=TestDeathDate;
            
            // 2. Lets feed death into the eventQ
            if (MyPointerToPerson->DateOfDeath<EndYear){
                int p=MyPointerToPerson->PersonID-1;
                event * DeathEvent = new event;
                Events.push_back(DeathEvent);
                DeathEvent->time = MyPointerToPerson->DateOfDeath;
                DeathEvent->p_fun = &EventMyDeathDate;
                DeathEvent->person_ID = MyArrayOfPointersToPeople[p];
                p_PQ->push(DeathEvent);
                
                // Update cause of death
                MyPointerToPerson->CauseOfDeath=13;
            }
        }
    }
}


void EventMyProstateDate (person *MyPointerToPerson)			// Function executed when person develops hypertension
{
    
    E(cout << endl << endl << "This patient just developed prostate cancer!" << endl;)
    
    if (MyPointerToPerson->Prostate==0 && MyPointerToPerson->Alive == 1){
        
        MyPointerToPerson->Prostate_status=1;
        
        // Lets see if we need to update death date
        int cancer_index=5;
        
        // Lets see if they die earlier
        int k=(MyPointerToPerson->DoB-1800);					// To find corresponding year of birth from mortality array
        int j=0;												// This will be matched to probability taken from random number generator
        double	d = ((double) rand() / (RAND_MAX)) ;			// get a random number to determine Life Expectancy
        double TestDeathDate;
        
        if (MyPointerToPerson->Sex==1){
            while(d>DeathArray_Men[k][j]*MortRisk_Cancer[cancer_index] && j<121){j++;}
            TestDeathDate=(MyPointerToPerson->DoB+j);
        }
        
        if (MyPointerToPerson->Sex==2) {
            while(d>DeathArray_Women[k][j]*MortRisk_Cancer[cancer_index] && j<121){j++;}
            TestDeathDate=(MyPointerToPerson->DoB+j);
        }
        
        if (TestDeathDate<MyPointerToPerson->DateOfDeath && TestDeathDate>*p_GT) {
            
            MyPointerToPerson->DateOfDeath=TestDeathDate;
            
            // 2. Lets feed death into the eventQ
            if (MyPointerToPerson->DateOfDeath<EndYear){
                int p=MyPointerToPerson->PersonID-1;
                event * DeathEvent = new event;
                Events.push_back(DeathEvent);
                DeathEvent->time = MyPointerToPerson->DateOfDeath;
                DeathEvent->p_fun = &EventMyDeathDate;
                DeathEvent->person_ID = MyArrayOfPointersToPeople[p];
                p_PQ->push(DeathEvent);
                
                // Update cause of death
                MyPointerToPerson->CauseOfDeath=14;
            }
        }
    }
}


void EventMyOtherCanDate (person *MyPointerToPerson)			// Function executed when person develops hypertension
{
    
    E(cout << endl << endl << "This patient just developed Other cancer!" << endl;)
    
    if (MyPointerToPerson->OtherCan_status==0 && MyPointerToPerson->Alive == 1){
        
        MyPointerToPerson->OtherCan_status=1;
        
        // Lets see if we need to update death date
        int cancer_index=6;
        
        // Lets see if they die earlier
        int k=(MyPointerToPerson->DoB-1800);					// To find corresponding year of birth from mortality array
        int j=0;												// This will be matched to probability taken from random number generator
        double	d = ((double) rand() / (RAND_MAX)) ;			// get a random number to determine Life Expectancy
        double TestDeathDate;
        
        if (MyPointerToPerson->Sex==1){
            while(d>DeathArray_Men[k][j]*MortRisk_Cancer[cancer_index] && j<121){j++;}
            TestDeathDate=(MyPointerToPerson->DoB+j);
        }
        
        if (MyPointerToPerson->Sex==2) {
            while(d>DeathArray_Women[k][j]*MortRisk_Cancer[cancer_index] && j<121){j++;}
            TestDeathDate=(MyPointerToPerson->DoB+j);
        }
        
        if (TestDeathDate<MyPointerToPerson->DateOfDeath && TestDeathDate>*p_GT) {
            
            MyPointerToPerson->DateOfDeath=TestDeathDate;
            
            // 2. Lets feed death into the eventQ
            if (MyPointerToPerson->DateOfDeath<EndYear){
                int p=MyPointerToPerson->PersonID-1;
                event * DeathEvent = new event;
                Events.push_back(DeathEvent);
                DeathEvent->time = MyPointerToPerson->DateOfDeath;
                DeathEvent->p_fun = &EventMyDeathDate;
                DeathEvent->person_ID = MyArrayOfPointersToPeople[p];
                p_PQ->push(DeathEvent);
                
                // Update cause of death
                MyPointerToPerson->CauseOfDeath=15;
            }
        }
    }
}




