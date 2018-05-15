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
#include "NCD.hpp"
#include "person.h"
#include "LoadParams.h"
#include "eventfunctions.h"
#include "errorcoutmacro.h"
#include "CParamReader.hpp"
#include "CountryParams.hpp"
#include "Intervention.hpp"
#include "HPVInfection.hpp"

using namespace std;

// remove all unecessary varibales and arrays and references
// move ART start numbers to HIVInfection .hpp


//// --- OUTSIDE INFORMATION --- ////
extern int *p_PY;
extern double *p_GT;
extern int ART_start_yr;
extern int total_population;
extern person** MyArrayOfPointersToPeople;
extern priority_queue<event*, vector<event*>, timeComparison> *p_PQ;
extern vector<event*> Events;

extern int* ARTKids;
extern int** ARTMen;
extern int* ARTMen_sum;
extern int** ARTWomen;
extern int* ARTWomen_sum;


//// --- Important Internal informtaion --- ////
extern int RandomMinMax_2(int min, int max);


/// HIV Increased risk for NCD Parameter ///
extern int ageAdult;
extern double ARTbuffer;
extern int count_ARTMen[7];              // Count by CD4 count category
extern int count_ARTAdult_Women[7];      // Count by CD4 count category
extern int count_ARTMen_sum;             // Count the sum for the year - men
extern int count_ARTWomen_sum;           // Count the sum for the year - women
extern int count_ARTKids;
int count_AdultsART;
int ART_index;


/// HIV related
extern int HIV_Ref_PersonID[500000];
extern int countHIVRef;
extern int KIDS_HIV_Ref_PersonID[70000];
extern int countKIDSHIVRef;
extern int Elig_Men;
extern int Elig_Women;
extern int Elig_Kids;

// HPV variables
extern int int_CCscreening;
extern int age_HPVvaccination;
extern int int_HPVvaccination;
extern int yearintervention_start;

//CVD Intervention variables
extern int int_CVDIntervention;

//// --- NEW YEAR FUNCTION --- ////                                     /
void EventTellNewYear(person *MyPointerToPerson){
    cout << "A new year has started, it is now " << *p_GT << endl;      // Tells us a new year has started
    
    // Lets set the first people on ART
    if (*p_GT>=ART_start_yr){
        
        // Add people on ART until we hit our aim
        while (count_ARTKids<ARTKids[ART_index] || ARTbuffer*count_ARTMen_sum<ARTMen_sum[ART_index] || ARTbuffer*count_ARTWomen_sum<ARTWomen_sum[ART_index]){
            
           if (count_ARTKids<ARTKids[ART_index] ){
                
                int a=(RandomMinMax_2(1,countKIDSHIVRef-1));    // Get a random person and update age
                int i=KIDS_HIV_Ref_PersonID[a]-1;
                MyArrayOfPointersToPeople[i]->Age=(*p_GT - MyArrayOfPointersToPeople[i]->DoB);
                
                
                // Lets get some kids onto ART
                if (((MyArrayOfPointersToPeople[i]->Age<ageAdult && *p_GT<2012)|| (MyArrayOfPointersToPeople[i]->Age<ageAdult && *p_GT>=2012)) && MyArrayOfPointersToPeople[i]->HIV>0 && MyArrayOfPointersToPeople[i]->HIV<=*p_GT && MyArrayOfPointersToPeople[i]->ART==-999 && MyArrayOfPointersToPeople[i]->Alive==1 && count_ARTKids<ARTKids[ART_index] ){
                    
                    MyArrayOfPointersToPeople[i]->ART=*p_GT;            // Lets set ART date
                    MyArrayOfPointersToPeople[i]->CD4_cat_ARTstart=MyArrayOfPointersToPeople[i]->CD4_cat;       // Lets set CD4 cat at ART start
                    count_ARTKids++;                                    // Update our counter
                    Elig_Kids=Elig_Kids-1;
                    
                    
                    if (int_CCscreening==1 && *p_GT<=yearintervention_start && MyArrayOfPointersToPeople[i]->Sex==2){
                        double TimeOfScreen = 0;
                        if (MyArrayOfPointersToPeople[i]->DoB<=2001) {TimeOfScreen = 2018;}
                        if (MyArrayOfPointersToPeople[i]->DoB>2001) {TimeOfScreen = MyArrayOfPointersToPeople[i]->DoB+17;}
                        event * CC_First_screen = new event;
                        Events.push_back(CC_First_screen);
                        CC_First_screen->person_ID = MyArrayOfPointersToPeople[i];
                        CC_First_screen->time = TimeOfScreen;
                        CC_First_screen->p_fun = &EventVIA_Screening;
                        p_PQ->push(CC_First_screen);
                        p_PQ->push(CC_First_screen);
                    }
                }
            }
            
            
            // Lets see if adults start ART
            if (ARTbuffer*count_ARTMen_sum<ARTMen_sum[ART_index] || ARTbuffer*count_ARTWomen_sum<ARTWomen_sum[ART_index]){
                
                int a=(RandomMinMax_2(1,countHIVRef-1));
                int i=HIV_Ref_PersonID[a]-1;
                MyArrayOfPointersToPeople[i]->Age=(*p_GT - MyArrayOfPointersToPeople[i]->DoB);
                
                // Lets let men start ART
                if (MyArrayOfPointersToPeople[i]->Sex==1 && MyArrayOfPointersToPeople[i]->Age>=ageAdult && MyArrayOfPointersToPeople[i]->HIV>0 && MyArrayOfPointersToPeople[i]->HIV<*p_GT && MyArrayOfPointersToPeople[i]->ART==-999 &&MyArrayOfPointersToPeople[i]->Alive==1 && count_ARTMen_sum<ARTMen_sum[ART_index] ){
                    
                    if (ARTbuffer*ARTMen[ART_index][MyArrayOfPointersToPeople[i]->CD4_cat] > count_ARTMen[MyArrayOfPointersToPeople[i]->CD4_cat]){
                        
                        MyArrayOfPointersToPeople[i]->ART=*p_GT;                        // Lets set ART date
                        MyArrayOfPointersToPeople[i]->CD4_cat_ARTstart=MyArrayOfPointersToPeople[i]->CD4_cat;       // Lets set CD4 cat at ART start
                        count_ARTMen[MyArrayOfPointersToPeople[i]->CD4_cat]++;    // Update our counter CD4/ART array
                        count_ARTMen_sum++;                                             // Update the sum counter
                        count_AdultsART++;
                        Elig_Men=Elig_Men-1;
                        
                        if (int_CVDIntervention==1){
                            double TimeOfScreen = 0;
                            if (MyArrayOfPointersToPeople[i]->DoB<=2000) {TimeOfScreen = 2018;}
                            if (MyArrayOfPointersToPeople[i]->DoB>2000) {TimeOfScreen = MyArrayOfPointersToPeople[i]->DoB+18;}
                            //Push a CVD intervention date
                            event * CVD_Screen_Date = new event;
                            Events.push_back(CVD_Screen_Date);
                            CVD_Screen_Date->time = TimeOfScreen;
                            CVD_Screen_Date->p_fun = &EventCVDPrevIntervention;
                            CVD_Screen_Date->person_ID = MyArrayOfPointersToPeople[i];
                            p_PQ->push(CVD_Screen_Date);
                            
                        }
                        
                    }
                }
            
                
                
                // Lets start women on ART
                if (MyArrayOfPointersToPeople[i]->Sex==2 && MyArrayOfPointersToPeople[i]->Age>=ageAdult && MyArrayOfPointersToPeople[i]->HIV>0 && MyArrayOfPointersToPeople[i]->HIV<*p_GT && MyArrayOfPointersToPeople[i]->ART==-999 &&MyArrayOfPointersToPeople[i]->Alive==1 && count_ARTWomen_sum<ARTWomen_sum[ART_index]){
                    
                    
                    if (ARTbuffer*ARTWomen[ART_index][MyArrayOfPointersToPeople[i]->CD4_cat] > count_ARTAdult_Women[MyArrayOfPointersToPeople[i]->CD4_cat]){
                        
                        
                        MyArrayOfPointersToPeople[i]->ART=*p_GT;            // Lets set ART date
                        MyArrayOfPointersToPeople[i]->CD4_cat_ARTstart=MyArrayOfPointersToPeople[i]->CD4_cat;       // Lets set CD4 cat at ART start
                        count_ARTAdult_Women[MyArrayOfPointersToPeople[i]->CD4_cat]++;  // Update our counter CD4/ART array
                        count_ARTWomen_sum++;                               // Update the sum counter
                        count_AdultsART++;
                        Elig_Women=Elig_Women-1;
                        
                        if (int_CVDIntervention==1){
                            double TimeOfScreen = 0;
                            if (MyArrayOfPointersToPeople[i]->DoB<=2000) {TimeOfScreen = 2018;}
                            if (MyArrayOfPointersToPeople[i]->DoB>2000) {TimeOfScreen = MyArrayOfPointersToPeople[i]->DoB+18;}
                            
                            //Push a CVD intervention date
                            event * CVD_Screen_Date = new event;
                            Events.push_back(CVD_Screen_Date);
                            CVD_Screen_Date->time = TimeOfScreen;
                            CVD_Screen_Date->p_fun = &EventCVDPrevIntervention;
                            CVD_Screen_Date->person_ID = MyArrayOfPointersToPeople[i];
                            p_PQ->push(CVD_Screen_Date);
                                
                        }
                    }
                        
                    if (int_CCscreening==1){
                        double TimeOfScreen = 0;
                        if (*p_GT>=2018) {TimeOfScreen = *p_GT;}
                        if (*p_GT<2018) {TimeOfScreen = 2018;}
                        event * CC_First_screen = new event;
                        Events.push_back(CC_First_screen);
                        CC_First_screen->person_ID = MyArrayOfPointersToPeople[i];
                        CC_First_screen->time = TimeOfScreen;
                        CC_First_screen->p_fun = &EventVIA_Screening;
                        p_PQ->push(CC_First_screen);
                        p_PQ->push(CC_First_screen);
                    }
                }
            }
        }
        
        // Lets update the ART index
        ART_index++;
        if (ART_index>12){ART_index=12;}
    }
    
    // Lets get the pointer to the right year range
    if (*p_GT<1955){*p_PY = 0; };
    if (*p_GT >= 1955 && *p_GT<1960){*p_PY = 1; };
    if (*p_GT >= 1960 && *p_GT<1965){*p_PY = 2; };
    if (*p_GT >= 1965 && *p_GT<1970){*p_PY = 3; };
    if (*p_GT >= 1970 && *p_GT<1975){*p_PY = 4; };
    if (*p_GT >= 1975 && *p_GT<1980){*p_PY = 5; };
    if (*p_GT >= 1980 && *p_GT<1985){*p_PY = 6; };
    if (*p_GT >= 1985 && *p_GT<1990){*p_PY = 7; };
    if (*p_GT >= 1990 && *p_GT<1995){*p_PY = 8; };
    if (*p_GT >= 1995 && *p_GT<2000){*p_PY = 9; };
    if (*p_GT >= 2000 && *p_GT<2005){*p_PY = 10;};
    if (*p_GT >= 2005){*p_PY = 11; };
    E(cout << "Let's check relevant things have been updated... *p_PY: " << *p_PY << " and Global Time: " << *p_GT << endl;)
    
    // Schedule event for next year
    event * RecurrentTellNewYear = new event;
    Events.push_back(RecurrentTellNewYear);
    RecurrentTellNewYear->time = *p_GT + 1;
    RecurrentTellNewYear->p_fun = &EventTellNewYear;
    p_PQ->push(RecurrentTellNewYear);
        
        cout << "Next year is been scheduled " << *p_GT+1 << endl;
    
    E(cout << "We have finished telling you the new year and setting fertility variables for the year." << endl;)	// Error message - can be switched on/off
    
}


//// --- BIRTH EVENT AND MAKING NEW PERSON --- ////
void EventBirth(person *MyPointerToPerson){
    
    E(cout << "A birth is about to happen" << endl;)
    
    if(MyPointerToPerson->Alive == 1) {
        
        total_population=total_population+1;									// Update total population for output and for next new entry
        
        // Creating a new person
        MyArrayOfPointersToPeople[total_population-1]=new person();
        (MyArrayOfPointersToPeople[total_population-1])->PersonIDAssign(total_population-1);
        (MyArrayOfPointersToPeople[total_population-1])->Alive=1;
        (MyArrayOfPointersToPeople[total_population-1])->GenderDistribution();
        (MyArrayOfPointersToPeople[total_population-1])->GetMyDoBNewEntry();
        (MyArrayOfPointersToPeople[total_population-1])->GetDateOfDeath();
        (MyArrayOfPointersToPeople[total_population-1])->GetDateOfBaby();
        (MyArrayOfPointersToPeople[total_population-1])->GetMyDateNCD();
        (MyArrayOfPointersToPeople[total_population-1])->GetMyDateCancers();
        (MyArrayOfPointersToPeople[total_population-1])->GetMyDateOfHPVInfection();
        (MyArrayOfPointersToPeople[total_population-1])->GetMyDateOfHIVInfection();
        (MyArrayOfPointersToPeople[total_population-1])->MotherID=MyPointerToPerson->PersonID;			// Give child their mothers ID
        
        MyPointerToPerson->ChildIDVector.push_back((MyArrayOfPointersToPeople[total_population-1]));	// Give mothers their child's ID
        
        
        // This code is for the HPV vaccination Intervention
        if (*p_GT>=(yearintervention_start-age_HPVvaccination) && int_HPVvaccination==1 && MyArrayOfPointersToPeople[total_population-1]->Sex==2)
        {
            event * HPV_VaccinationEvent = new event;
            Events.push_back(HPV_VaccinationEvent);
            HPV_VaccinationEvent->time = MyArrayOfPointersToPeople[total_population-1]->DoB + age_HPVvaccination;
            HPV_VaccinationEvent->p_fun = &EventMyHPVVaccination;
            HPV_VaccinationEvent->person_ID = MyArrayOfPointersToPeople[total_population-1];
            p_PQ->push(HPV_VaccinationEvent);
        }
        
        E(cout << "We have finished giving birth " << endl;)
    }
}

//// --- DEATH EVENT --- ////
void EventMyDeathDate(person *MyPointerToPerson){
    
    E(cout << "A death is about to happen " << endl;)
    
    if (MyPointerToPerson->Alive==1){
        
        MyPointerToPerson->Alive=0;
        
        if (MyPointerToPerson->CauseOfDeath==-999)      // Updates cause of death to "other" if not assigned already
        {
            MyPointerToPerson->CauseOfDeath=1;
        }
        
        
        // Lets also update ART information as people die
        // 1. If they are not on ART
        if (MyPointerToPerson->HIV>0 && MyPointerToPerson->HIV<*p_GT && MyPointerToPerson->ART==-999){
            
            // First we need to update age
            MyPointerToPerson->Age=(*p_GT - MyPointerToPerson->DoB);  // First we update age
            
            if (MyPointerToPerson->Sex==1 && MyPointerToPerson->Age>=ageAdult){Elig_Men=Elig_Men-1;}
            if (MyPointerToPerson->Sex==2 && MyPointerToPerson->Age>=ageAdult){Elig_Women=Elig_Women-1;}
            if (MyPointerToPerson->Age<ageAdult){Elig_Kids=Elig_Kids-1;}
        }
        
        
        // 2. If they are on ART
        if (MyPointerToPerson->HIV>0 && MyPointerToPerson->HIV<*p_GT && MyPointerToPerson->ART>0){
            
            // First we need to update age
            MyPointerToPerson->Age=(*p_GT - MyPointerToPerson->DoB);  // First we update age
            
            // If it is a child
            if (MyPointerToPerson->Age<ageAdult){
                if (count_ARTKids>0) {count_ARTKids=count_ARTKids-1;}
            }
            
            // If it is a man
            if (MyPointerToPerson->Sex==1 && MyPointerToPerson->Age>=ageAdult){
                count_ARTMen[MyPointerToPerson->CD4_cat_start]=count_ARTMen[MyPointerToPerson->CD4_cat_start]-1;
                count_ARTMen_sum=count_ARTMen_sum-1;
            }
            
            // If it is a woman
            if (MyPointerToPerson->Sex==2 && MyPointerToPerson->Age>=ageAdult){
                count_ARTAdult_Women[MyPointerToPerson->CD4_cat_start]=count_ARTAdult_Women[MyPointerToPerson->CD4_cat_start]-1;
                count_ARTWomen_sum=count_ARTWomen_sum-1;
            }
        }
    }
    
    E(cout << "Person " << MyPointerToPerson->PersonID << " just died. Their life status now is: " << MyPointerToPerson->Alive << endl;)
}
