/////////////////////////////////////////////////////////////////
//    eventfunction.cpp                                        //
//    Created by Mikaela Smit on 07/11/2014.				   //
//    Copyright (c) 2014 Mikaela Smit. All rights reserved.    //
//    These event are added to the EventQ.					   //
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
#include "event.h"
#include "eventQ.h"
#include "person.h"
#include "NCD.hpp"
#include "LoadParams.h"
#include "eventfunctions.h"
#include "errorcoutmacro.h"
#include "CParamReader.hpp"
#include "CountryParams.hpp"
#include "Intervention.hpp"
#include "HPVInfection.hpp"
#include "HIVInfection.hpp"

using namespace std;

// Clean this code!
// Remove all unecessary headers and variables


//// --- OUTSIDE INFORMATION --- ////
extern double *p_GT;
extern person** MyArrayOfPointersToPeople;
extern priority_queue<event*, vector<event*>, timeComparison> *p_PQ;
extern vector<event*> Events;


extern double*** CD4_startarray;
extern double**  CD4_prog_rates;
extern double*** Death_CD4_rates;
extern double*** ART_CD4_rates;
extern double**  NCDArray;
extern double**  HPVarray;
extern int** ARTMen;
extern int** ARTWomen;
extern int*  ARTKids;
extern int*  ARTMen_sum;
extern int*  ARTWomen_sum;
extern int   age_atrisk_hpv;
extern int   max_age_atrisk_hpv;
extern double HPV_Array_buffer;
extern double HPV_hiv_Array_buffer;
extern int    HPV_Status_Recovered;
extern int int_CVDIntervention;


//// --- Important Internal informtaion --- ////
extern int RandomMinMax_2(int min, int max);
double randfrom_1(double min, double max){
    double range = (max - min);
    double div = RAND_MAX / range;
    return min + (rand() / div);
}


/// HIV Increased risk for NCD Parameter ///
extern double Risk_HIVHT;                                 // Increased risk of HT given HIV (and CKD given HIV below) - from Schouten et al CID 2016
extern double Risk_HIVCKD;
extern double Risk_NCDHIV[2];
extern int relatedNCD_HIV[2];
extern int nr_NCD_HIV;
extern int ageAdult;


/// HIV related
int count_ARTMen[7]={0, 0, 0, 0, 0, 0, 0};              // Count by CD4 count category
int count_ARTAdult_Women[7]={0, 0, 0, 0, 0, 0, 0};      // Count by CD4 count category
int count_ARTMen_sum=0;                                 // Count the sum for the year - men
int count_ARTWomen_sum=0;                               // Count the sum for the year - women
int count_ARTKids=0;
int HIV_Ref_PersonID[500000]={};
int countHIVRef=0;
int KIDS_HIV_Ref_PersonID[70000]={};
int countKIDSHIVRef=0;
int Elig_Men=0;
int Elig_Women=0;
int Elig_Kids=0;


//// --- HIV EVENT --- ////
void EventMyHIVInfection(person *MyPointerToPerson){
    
    E(cout << "Somebody is about to get infected with HIV: " << endl;)			// Error message - can be switched on/off
    
    if(MyPointerToPerson->Alive == 1) {											// Only execute this is patient is still alove
        
        if (MyPointerToPerson->DateOfDeath<MyPointerToPerson->HIV){
            cout << "Error! GT " << *p_GT << " PersonID: " << MyPointerToPerson->PersonID << " Death: " << MyPointerToPerson->DateOfDeath << " HIV: " << MyPointerToPerson->HIV<< endl;
        }
        
        MyPointerToPerson->Age= (*p_GT - MyPointerToPerson->DoB);				// Update age to get correct parameter below
        
        //// Update key parameters for ART
        if (MyPointerToPerson->Age>ageAdult){
            HIV_Ref_PersonID[countHIVRef]=MyPointerToPerson->PersonID;
            countHIVRef++;
            if (MyPointerToPerson->Sex==1){Elig_Men++;}
            if (MyPointerToPerson->Sex==2){Elig_Women++;}
        }
        
        if (MyPointerToPerson->Age<ageAdult) {
            KIDS_HIV_Ref_PersonID[countKIDSHIVRef]=MyPointerToPerson->PersonID;
            countKIDSHIVRef++;
            Elig_Kids++;
            
            // Now we need to schedule an event where we move them to adult category
            double DateCD4CatSwitch=(MyPointerToPerson->DoB+ageAdult);
            
            // Schedule event for moving from kids ART to adult ART
            event * ARTCatMove = new event;
            Events.push_back(ARTCatMove);
            ARTCatMove->time = DateCD4CatSwitch;
            ARTCatMove->p_fun = &EventARTCatSwitch;
            ARTCatMove->person_ID = MyPointerToPerson;
            p_PQ->push(ARTCatMove);
        }
        
        //// --- Get my CD4 count at start --- ////
        double	h = ((double)rand() / (RAND_MAX));								// Gets a random number between 0 and 1.
        int i=0;
        int j=0;
        int a=25;
        
        while (MyPointerToPerson->Age>=a && a<46){a=a+10; i++;};				// To get the right age-specific row in the above sex-specific arrays
        while (h>CD4_startarray[MyPointerToPerson->Sex-1][i][j] && j<3){j++;}	// To get the corresponding CD4 count for correct age and sex from random 'h' generated above
        
        MyPointerToPerson->CD4_cat=0+j;											// CD4 count cat (variable over time)
        MyPointerToPerson->CD4_cat_start=0+j;									// CD4 count cat at start (to keep in records)
        MyPointerToPerson->CD4_change.at(MyPointerToPerson->CD4_cat)=*p_GT;
        
        //// --- Let's see what will happen next (Death, CD4 count progression or ART initiation) ---- ////
        double FindART_CD4_rate = ART_CD4_rates[MyPointerToPerson->Sex-1][i][MyPointerToPerson->CD4_cat];
        double FindCD4_rate = CD4_prog_rates[MyPointerToPerson->Sex-1][MyPointerToPerson->CD4_cat];
        double FindDeath_CD4_rate = Death_CD4_rates[MyPointerToPerson->Sex-1][i][MyPointerToPerson->CD4_cat];
        
        // Lets see when CD4 count progression would start
        double CD4_test = 0;
        double cd4 = ((double)rand() / (RAND_MAX));
        CD4_test = (-1/FindCD4_rate) * log(cd4);
        
        // Lets see when death would happen
        double death_test = 0;
        double dd = ((double)rand() / (RAND_MAX));
        death_test = (-1/FindDeath_CD4_rate) * log(dd);
        
        
        //  If they do not start ART then we need to evaluate if they die or decrease in CD4 count
        if (MyPointerToPerson->ART<0){
            
            if (CD4_test<death_test){
                event * CD4change = new event;
                Events.push_back(CD4change);
                CD4change->time = *p_GT+CD4_test;
                CD4change->p_fun = &EventCD4change;
                CD4change->person_ID = MyPointerToPerson;
                p_PQ->push(CD4change);
                
            }
            
            if (death_test<=CD4_test){
                
                double death_test_date = *p_GT +death_test;                  // Get the actual date, not just time until death
                
                if (death_test_date<MyPointerToPerson->DateOfDeath && death_test_date>*p_GT){   // Check HIV deaths happens before natural death
                    
                    MyPointerToPerson->DateOfDeath=death_test_date;
                    
                    event * DeathEvent = new event;						    // In that case we need to add the new death date to the EventQ
                    Events.push_back(DeathEvent);
                    DeathEvent->time = MyPointerToPerson->DateOfDeath;
                    DeathEvent->p_fun = &EventMyDeathDate;
                    DeathEvent->person_ID = MyPointerToPerson;
                    p_PQ->push(DeathEvent);
                    
                    MyPointerToPerson->CauseOfDeath=2;                      // and we will need to update Cause of Death
                }
            }
        }
        
        // Now we need to update NCD risks for HT and CKD
        E(cout << endl << endl << "We are assigning NCDs for HIV+ upon infection!" << endl;)
        
        double age_at_death=MyPointerToPerson->DateOfDeath-MyPointerToPerson->DoB;
        
        if (age_at_death>18)
        {
            
            // Some basic code and finding index for not getting NCDs
            int ncd_nr=0;                                                  // Assisgn all the possible NCDs in this code
            double DateNCD=-997;
            
            
            // Re-evaluate NCD related to HIV
            while (ncd_nr<nr_NCD_HIV){
                double r = ((double) rand() / (RAND_MAX));                 // Get a random number for each NCD
                
                // If we are getting an NCD lets get age and date
                if (r<=NCDArray[relatedNCD_HIV[ncd_nr]][120]*Risk_NCDHIV[ncd_nr])             // If they will get and NCD lets get the age and date and also update mortality
                {
                    
                    // Lets get the index for age at NCD
                    int i=0;
                    while (r>NCDArray[relatedNCD_HIV[ncd_nr]][i]*Risk_NCDHIV[ncd_nr]){i++;}
                    
                    // Lets get the age and date they will have the NCD
                    double YearFraction=(RandomMinMax_2(1,12))/12.1;                          // This gets month of birth as a fraction of a year
                    DateNCD=MyPointerToPerson->DoB+i+YearFraction;
                    
                    // Lets see if this pushed NCD date forward
                    if (DateNCD>=*p_GT && DateNCD<MyPointerToPerson->NCD_DatesVector.at(relatedNCD_HIV[ncd_nr]))
                    {
                        // Lets update the Date everywhere and add to the queue
                        MyPointerToPerson->NCD_DatesVector.at(relatedNCD_HIV[ncd_nr])=DateNCD;
                        
                        
                        if (ncd_nr==0)
                        {
                            MyPointerToPerson->HT=DateNCD;
                            //// --- Lets feed Hypertension into the eventQ --- ////
                            int p=MyPointerToPerson->PersonID-1;
                            event * HTEvent = new event;
                            Events.push_back(HTEvent);
                            HTEvent->time = MyPointerToPerson->HT;
                            HTEvent->p_fun = &EventMyHyptenDate;
                            HTEvent->person_ID = MyArrayOfPointersToPeople[p];
                            p_PQ->push(HTEvent);
                        }
                        
                        
                        else if (ncd_nr==1)
                        {
                            MyPointerToPerson->CKD=DateNCD;
                            //// --- Lets feed MI into the eventQ --- ////
                            int p=MyPointerToPerson->PersonID-1;
                            event * CKDEvent = new event;
                            Events.push_back(CKDEvent);
                            CKDEvent->time = MyPointerToPerson->CKD;
                            CKDEvent->p_fun = &EventMyCKDDate;
                            CKDEvent->person_ID = MyArrayOfPointersToPeople[p];
                            p_PQ->push(CKDEvent);
                        }
                        
                        
                    }
                    
                    ncd_nr++;                                                    // Lets do the next NCD
                    
                }
            }
        }
        E(cout << "We finished assigning NCDs for HIV+!" << endl;)
        
        // Now we need to re-evaluate the risk of HPV infection for women with HIV
        
        E(cout << endl << endl << "We're re-evaluating the probability of HPV for HIV+ women upon infection!" << endl);
        
        if (MyPointerToPerson->Sex==2 && MyPointerToPerson->Age<=max_age_atrisk_hpv){
            
            double ReTestHPVDate=-999;
            double age_now=*p_GT-MyPointerToPerson->DoB;
            int f = floor(age_now);
            double h = randfrom_1(HPVarray[0][f]*HPV_hiv_Array_buffer,HPVarray[0][65]*HPV_hiv_Array_buffer);
            
            double r = ((double) rand() / (RAND_MAX));
            if (r<=HPVarray[0][65]*HPV_hiv_Array_buffer){                           // --- Lifetime probability = 90.7% for HIV-pos
                int i = 0;                                                          // --- (array from Yamada et al. 2008)
                while (h>(HPVarray[0][i]*HPV_hiv_Array_buffer)){i++;}
                double YearFraction=(RandomMinMax_2(1,12))/ 12.1;
                ReTestHPVDate = MyPointerToPerson->DoB + i + YearFraction;
            }
            
            if(ReTestHPVDate>=*p_GT && MyPointerToPerson->HPVcount==0){
                // --- She'll get it!
                MyPointerToPerson->HPV[MyPointerToPerson->HPVcount]=ReTestHPVDate;
                int p=MyPointerToPerson->PersonID-1;
                event * HPV_DateofInfectionEvent = new event;
                Events.push_back(HPV_DateofInfectionEvent);
                HPV_DateofInfectionEvent->time = MyPointerToPerson->HPV[MyPointerToPerson->HPVcount];
                HPV_DateofInfectionEvent->p_fun = &EventMyHPVInfection;
                HPV_DateofInfectionEvent->person_ID = MyArrayOfPointersToPeople[p];
                p_PQ->push(HPV_DateofInfectionEvent);
            }
            
            if(ReTestHPVDate>=*p_GT && MyPointerToPerson->HPVcount>0 && MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1]==HPV_Status_Recovered){
                // --- She'll get it!
                MyPointerToPerson->HPV[MyPointerToPerson->HPVcount]=ReTestHPVDate;
                int p=MyPointerToPerson->PersonID-1;
                event * HPV_DateofInfectionEvent = new event;
                Events.push_back(HPV_DateofInfectionEvent);
                HPV_DateofInfectionEvent->time = MyPointerToPerson->HPV[MyPointerToPerson->HPVcount];
                HPV_DateofInfectionEvent->p_fun = &EventMyHPVInfection;
                HPV_DateofInfectionEvent->person_ID = MyArrayOfPointersToPeople[p];
                p_PQ->push(HPV_DateofInfectionEvent);
            }
        }
    }
    
    E(cout << "Somebody has just been infected with HIV!" << endl;)				// Error message - can be switched on/off
}

void EventCD4change(person *MyPointerToPerson){
    
    E(cout << "Somebody is about to experience a drop in CD4 count: " << endl;)	// Error message - can be switched on/off
    
    if(MyPointerToPerson->Alive == 1) {											// Add any additional things to progression of CD4 count
        
        //// --- Let's get the right index for all relevant arrays used here and update important info ---- ////
        int i=0;
        int a=25;
        
        while (MyPointerToPerson->Age>=a && a<46){a=a+10; i++;};				// To get the right age-specific row in the above sex-specific arrays
        MyPointerToPerson->Age= (*p_GT - MyPointerToPerson->DoB);				// Update age to get correct parameter below
        MyPointerToPerson->CD4_cat=MyPointerToPerson->CD4_cat+1;                // Update CD4 count
        MyPointerToPerson->CD4_change.at(MyPointerToPerson->CD4_cat)=*p_GT;
        
        
        //// --- When CD4 count hits the lowest possible value and no ART is started --- ////
        if (MyPointerToPerson->CD4_cat==6 && MyPointerToPerson->ART<0){
            
            double FindDeath_CD4_rate = Death_CD4_rates[MyPointerToPerson->Sex-1][i][MyPointerToPerson->CD4_cat];
            
            // Lets see when death would happen
            double death_test = 0;
            double dd = ((double)rand() / (RAND_MAX));
            death_test = (-1/FindDeath_CD4_rate) * log(dd);
            //cout <<  death_test << endl;
            
            double death_test_date = *p_GT +death_test;                      // Get the actual date, not just time until death
            
            if (death_test_date<MyPointerToPerson->DateOfDeath && death_test_date>*p_GT){            // Check HIV deaths happens before natural death
                
                MyPointerToPerson->DateOfDeath=death_test_date;
                
                event * DeathEvent = new event;								// In that case we need to add the new death date to the EventQ
                Events.push_back(DeathEvent);
                DeathEvent->time = MyPointerToPerson->DateOfDeath;
                DeathEvent->p_fun = &EventMyDeathDate;
                DeathEvent->person_ID = MyPointerToPerson;
                p_PQ->push(DeathEvent);
                
                MyPointerToPerson->CauseOfDeath=2;                          // and we will need to update Cause of Death
            }
        }
        
        //// --- In case CD4 count is higher than minimum possible category.  They can either die OR Progress AND they still have a 1 year risk of death if on ART --- ///
        if (MyPointerToPerson->CD4_cat<6){
            
            //// --- Let's see what will happen next (Death, CD4 count progression or ART initiation) ---- ////
            double FindART_CD4_rate = ART_CD4_rates[MyPointerToPerson->Sex-1][i][MyPointerToPerson->CD4_cat];
            double FindCD4_rate = CD4_prog_rates[MyPointerToPerson->Sex-1][MyPointerToPerson->CD4_cat];
            double FindDeath_CD4_rate = Death_CD4_rates[MyPointerToPerson->Sex-1][i][MyPointerToPerson->CD4_cat];
            
            // Lets see when CD4 count progression would start
            double CD4_test = 0;
            double cd4 = ((double)rand() / (RAND_MAX));
            CD4_test = (-1/FindCD4_rate) * log(cd4);
            
            
            // Lets see when death would happen
            double death_test = 0;
            double dd = ((double)rand() / (RAND_MAX));
            death_test = (-1/FindDeath_CD4_rate) * log(dd);
            
            
            // If patient hasn't started ART yet
            if (MyPointerToPerson->ART<0){
                
                if (CD4_test<death_test){
                    event * CD4change = new event;
                    Events.push_back(CD4change);
                    CD4change->time = *p_GT+CD4_test;
                    CD4change->p_fun = &EventCD4change;
                    CD4change->person_ID = MyPointerToPerson;
                    p_PQ->push(CD4change);
                }
                
                if (death_test<=CD4_test){
                    
                    double death_test_date = *p_GT +death_test;                      // Get the actual date, not just time until death
                    
                    if (death_test_date<MyPointerToPerson->DateOfDeath && death_test_date>*p_GT){            // Check HIV deaths happens before natural death
                        MyPointerToPerson->DateOfDeath=death_test_date;
                        event * DeathEvent = new event;								// In that case we need to add the new death date to the EventQ
                        Events.push_back(DeathEvent);
                        DeathEvent->time = MyPointerToPerson->DateOfDeath;
                        DeathEvent->p_fun = &EventMyDeathDate;
                        DeathEvent->person_ID = MyPointerToPerson;
                        p_PQ->push(DeathEvent);
                        
                        MyPointerToPerson->CauseOfDeath=2;                          // and we will need to update Cause of Death
                    }
                    
                }
                
            }
            
        }
        
    }
    
    E(cout << "Somebody has just experiences a drop in CD4 count!" << endl;)	// Error message - can be switched on/off
}



//// --- ART Category switch from Kids to Adult --- ////
void EventARTCatSwitch(person *MyPointerToPerson){
    
    E(cout << "We are switching our kids to adult ART count: " << endl;)       // Error message - can be switched on/off
    
    if(MyPointerToPerson->Alive == 1) {
        
        
        // If they are not on ART
        if (MyPointerToPerson->ART==-999){
            
            // First we need to add them to eligible adults for selection
            HIV_Ref_PersonID[countHIVRef]=MyPointerToPerson->PersonID;         // And the patientID reference
            countHIVRef++;                                                     // Then we need to updated the adult counter
            Elig_Kids=Elig_Kids-1;
            if (MyPointerToPerson->Sex==1){Elig_Men++;}
            if (MyPointerToPerson->Sex==2){Elig_Women++;}
        }
        
        // If they are on ART
        if (MyPointerToPerson->ART>0){
            
            if (count_ARTKids>0) {count_ARTKids=count_ARTKids-1;}              // This count how many children are still on ART
            
            if (MyPointerToPerson->Sex==1){
                count_ARTMen_sum++;                                             // Lets add a count to the adult category
                count_ARTMen[4]++;                                        // Lets update the ART CD4 array counter
            }
            
            if (MyPointerToPerson->Sex==2){
                count_ARTWomen_sum++;                                           // Lets add a count to the adult category
                count_ARTAdult_Women[4]++;                                      // Lets update the ART CD4 array counter
            }
            
            if (int_CVDIntervention==1){
                double TimeOfScreen = 0;
                TimeOfScreen = MyPointerToPerson->DoB+18;
                //Push a CVD intervention date
                event * CVD_Screen_Date = new event;
                Events.push_back(CVD_Screen_Date);
                CVD_Screen_Date->time = TimeOfScreen;
                CVD_Screen_Date->p_fun = &EventCVDPrevIntervention;
                CVD_Screen_Date->person_ID = MyPointerToPerson;
                p_PQ->push(CVD_Screen_Date);
                
            }
            
        }
    }
}
