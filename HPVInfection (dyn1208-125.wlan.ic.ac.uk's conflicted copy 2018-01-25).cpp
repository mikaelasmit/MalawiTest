//
//  HPVInfection.cpp
//  HIVModelZimbabwe
//
//  Created by Mikaela Smit on 25/01/2018.
//  Copyright © 2018 Mikaela Smit. All rights reserved.
//


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

#include "HPVInfection.hpp"

using namespace std;

//// --- OUTSIDE INFORMATION --- ////
extern double *p_GT;
extern int EndYear;
extern priority_queue<event*, vector<event*>, timeComparison> *p_PQ;
extern vector<event*> Events;

extern int HPV_Status_HPV;
extern int HPV_Status_CIN1;
extern int HPV_Status_CIN2_3;
extern int HPV_Status_CIS;
extern int HPV_Status_ICC;
extern int HPV_Status_Recovered;


//// --- Parameters TO MOVE!!!!
double          CIN1_Rates[2]= {0.2,0.8};
double          CIN2_3_Rates[2]= {0.4,0.6};
double          CIS_Rates[2]= {0.65,0.35};
double          ICC_Rates[2]= {1.0,0.0};
extern double   hpv_date_after_death;
extern double   no_hpv_infection;
extern double** HPVarray;


//// --- Important functions --- ////
int RandomMinMax_2(int min, int max){                            // Provide function for random number generator between min and max number
    return rand()%(max-min+1)+min;                                // Note: if min=0 and max=4 it will generate 0,1,2,3,4
}


//// --- HPV EVENT --- ////

void EventMyHPVInfection(person *MyPointerToPerson){                    // Function executed when somebody develops HPV infection
    
    E(cout << "Somebody just got infected with HPV and will either progress to CIN1 or recover: " << endl;)
    
    if(MyPointerToPerson->Alive == 1 && MyPointerToPerson->HPV_Status<1 && MyPointerToPerson->HPVvaccination_status==0){                                  // Only execute this function if the person is still alive
        MyPointerToPerson->Age= (*p_GT - MyPointerToPerson->DoB);       // Update age to get correct parameter belowhp
        
        if (MyPointerToPerson->HPV_DateofInfection>0){
            MyPointerToPerson->HPV_Status=HPV_Status_HPV;
            
            int year = floor(*p_GT);
            double months = floor(((1-(*p_GT-year+0.01))*12));
            
            int j=0;                                                    // Find out how many years it will take to progress/recover (uniform dist from 0 to 3)
            float TestCIN1Date=0;
            std::random_device rd;
            std::mt19937 gen{rd()};
            std::uniform_int_distribution<> dis{0, 3};
            j = dis(gen);
            
            double YearFraction=-999;                                   // Get the date when the progression/recovery will take place
            if(months>=1){YearFraction=(RandomMinMax_2(0,months))/12.1;}
            if(months<1){YearFraction=0;}
            TestCIN1Date=(MyPointerToPerson->HPV_DateofInfection+j)+YearFraction;
            
            double    h = ((double)rand() / (RAND_MAX));                // Get a random number for the probability of progression vs recovery
            
            //// In case they recover from HPV
            if (h>CIN1_Rates[1]){
                MyPointerToPerson->CIN1_DateofProgression=no_hpv_infection;
                MyPointerToPerson->HPV_DateofRecovery=TestCIN1Date;
                
                //// Let's push recovery into the events Q
                event * HPV_DateofRecoveryEvent = new event;
                Events.push_back(HPV_DateofRecoveryEvent);
                HPV_DateofRecoveryEvent->time = MyPointerToPerson->HPV_DateofRecovery;
                HPV_DateofRecoveryEvent->p_fun = &EventMyHPVRecovery;
                HPV_DateofRecoveryEvent->person_ID = MyPointerToPerson;
                p_PQ->push(HPV_DateofRecoveryEvent);
                
                //cout << MyPointerToPerson << " got infected with HPV on " << MyPointerToPerson->HPV_DateofInfection << ". Given H = " << h << ", she recovered on " << MyPointerToPerson->HPV_DateofRecovery << " (" << MyPointerToPerson->CIN1_DateofProgression << ")" << endl;
            }
            
            //// In case they progress to CIN1
            if (h<=CIN1_Rates[1]){
                if (TestCIN1Date<MyPointerToPerson->DateOfDeath){MyPointerToPerson->CIN1_DateofProgression=TestCIN1Date;}{MyPointerToPerson->HPV_DateofRecovery=-977;}
                if (TestCIN1Date>=MyPointerToPerson->DateOfDeath) {MyPointerToPerson->CIN1_DateofProgression=hpv_date_after_death;}
                
                //// Let's push progression into the events Q
                event * CIN1_DateofProgressionEvent = new event;
                Events.push_back(CIN1_DateofProgressionEvent);
                CIN1_DateofProgressionEvent->time = MyPointerToPerson->CIN1_DateofProgression;
                CIN1_DateofProgressionEvent->p_fun = &EventMyCIN1Status;
                CIN1_DateofProgressionEvent->person_ID = MyPointerToPerson;
                p_PQ->push(CIN1_DateofProgressionEvent);
                
                //cout << MyPointerToPerson << " got infected with HPV on " << MyPointerToPerson->HPV_DateofInfection << ". Given H = " << h << ", she progressed to CIN1 on " << MyPointerToPerson->CIN1_DateofProgression << endl;
            }
        }
    }
    E(cout << "Somebody with HPV just progressed to CIN1 or recovered!" << endl;)
}

void EventMyCIN1Status(person *MyPointerToPerson){
    
    E(cout << "Somebody with CIN1 is about to progress to CIN2/3 or recover: " << endl;)
    
    if(MyPointerToPerson->Alive == 1){
        MyPointerToPerson->Age= (*p_GT - MyPointerToPerson->DoB);
        
        ///Decide if person recovers or moves on to CIN2/3
        if (MyPointerToPerson->CIN1_DateofProgression>0){
            MyPointerToPerson->HPV_Status=HPV_Status_CIN1;
            int year = floor(*p_GT);
            double months = floor(((1-(*p_GT-year+0.01))*12));
            
            int j=0;
            float TestCIN2_3Date=0;
            std::random_device rd;
            std::mt19937 gen{rd()};
            std::uniform_int_distribution<> dis{0, 3};
            j = dis(gen);
            
            double YearFraction=-999;
            if(months>=1){YearFraction=(RandomMinMax_2(0,months))/12.1;}
            if(months<1){YearFraction=0;}
            double    h = ((double)rand() / (RAND_MAX));
            TestCIN2_3Date=(MyPointerToPerson->CIN1_DateofProgression+j)+YearFraction;
            
            //// In case they recover from CIN1
            if (h>CIN2_3_Rates[1]){
                MyPointerToPerson->CIN2_3_DateofProgression=no_hpv_infection;
                MyPointerToPerson->CIN1_DateofRecovery=TestCIN2_3Date;
                
                //// Let's push recovery into the events Q
                event * CIN1_DateofRecoveryEvent = new event;
                Events.push_back(CIN1_DateofRecoveryEvent);
                CIN1_DateofRecoveryEvent->time = MyPointerToPerson->CIN1_DateofRecovery;
                CIN1_DateofRecoveryEvent->p_fun = &EventMyHPVRecovery;
                CIN1_DateofRecoveryEvent->person_ID = MyPointerToPerson;
                p_PQ->push(CIN1_DateofRecoveryEvent);
                
                //cout << MyPointerToPerson << " got infected with HPV on " << MyPointerToPerson->HPV_DateofInfection << " and progressed to CIN1 on " << MyPointerToPerson->CIN1_DateofProgression << ". Given H = " << h << ", she recovered from CIN1 on " << MyPointerToPerson->CIN1_DateofRecovery << " (" << MyPointerToPerson->CIN2_3_DateofProgression << ")" << endl;
            }
            
            //// In case they progress to CIN2/3
            if (h<=CIN2_3_Rates[1]){
                if (TestCIN2_3Date<MyPointerToPerson->DateOfDeath){MyPointerToPerson->CIN2_3_DateofProgression=TestCIN2_3Date;}{MyPointerToPerson->CIN1_DateofRecovery=-977;}
                if (TestCIN2_3Date>=MyPointerToPerson->DateOfDeath) {MyPointerToPerson->CIN2_3_DateofProgression=hpv_date_after_death;}
                
                //// Let's push progression into the events Q
                event * CIN2_3_DateofProgressionEvent = new event;
                Events.push_back(CIN2_3_DateofProgressionEvent);
                CIN2_3_DateofProgressionEvent->time = MyPointerToPerson->CIN2_3_DateofProgression;
                CIN2_3_DateofProgressionEvent->p_fun = &EventMyCIN2_3Status;
                CIN2_3_DateofProgressionEvent->person_ID = MyPointerToPerson;
                p_PQ->push(CIN2_3_DateofProgressionEvent);
                
                //cout << MyPointerToPerson << " got infected with HPV on " << MyPointerToPerson->HPV_DateofInfection << " and progressed to CIN1 on " << MyPointerToPerson->CIN1_DateofProgression << ". Given H = " << h << ", she will now progress from CIN1 to CIN 2/3 on " << MyPointerToPerson->CIN2_3_DateofProgression << endl;
            }
        }
    }
    E(cout << "Somebody with CIN1 just progressed to CIN2/3 or recovered!" << endl;)
}

void EventMyCIN2_3Status(person *MyPointerToPerson){
    
    E(cout << "Somebody with CIN2_3 is about to progress to CIS or recover: " << endl;)
    
    if(MyPointerToPerson->Alive == 1){
        MyPointerToPerson->Age= (*p_GT - MyPointerToPerson->DoB);
        
        ///Decide if person recovers or moves on to CIS
        if (MyPointerToPerson->CIN2_3_DateofProgression>0){
            MyPointerToPerson->HPV_Status=HPV_Status_CIN2_3;
            int year = floor(*p_GT);
            double months = floor(((1-(*p_GT-year+0.01))*12));
            
            int j=0;
            float TestCISDate=0;
            std::random_device rd;
            std::mt19937 gen{rd()};
            std::uniform_int_distribution<> dis{0, 3};
            j = dis(gen);
            
            double YearFraction=-999;
            if(months>=1){YearFraction=(RandomMinMax_2(0,months))/12.1;}
            if(months<1){YearFraction=0;}
            double    h = ((double)rand() / (RAND_MAX));
            TestCISDate=(MyPointerToPerson->CIN2_3_DateofProgression+j)+YearFraction;
            
            //// In case they recover from CIN2/3
            if (h>CIS_Rates[1]){
                MyPointerToPerson->CIS_DateofProgression=no_hpv_infection;
                MyPointerToPerson->CIN2_3_DateofRecovery=TestCISDate;
                
                //// Let's push recovery into the events Q
                event * CIN2_3_DateofRecoveryEvent = new event;
                Events.push_back(CIN2_3_DateofRecoveryEvent);
                CIN2_3_DateofRecoveryEvent->time = MyPointerToPerson->CIN2_3_DateofRecovery;
                CIN2_3_DateofRecoveryEvent->p_fun = &EventMyHPVRecovery;
                CIN2_3_DateofRecoveryEvent->person_ID = MyPointerToPerson;
                p_PQ->push(CIN2_3_DateofRecoveryEvent);
                
                //cout << MyPointerToPerson << " got infected with HPV on " << MyPointerToPerson->HPV_DateofInfection << ", progressed to CIN1 on " << MyPointerToPerson->CIN1_DateofProgression << " and to CIN2/3 on " << MyPointerToPerson->CIN2_3_DateofProgression << ". Given H = " << h << ", she recovered from CIN2/3 on " << MyPointerToPerson->CIN2_3_DateofRecovery << " (" << MyPointerToPerson->CIS_DateofProgression << ")" << endl;
            }
            
            //// In case they progress to CIS
            if (h<=CIS_Rates[1]){
                if (TestCISDate<MyPointerToPerson->DateOfDeath){MyPointerToPerson->CIS_DateofProgression=TestCISDate;}{MyPointerToPerson->CIN2_3_DateofRecovery=-977;}
                if (TestCISDate>=MyPointerToPerson->DateOfDeath) {MyPointerToPerson->CIS_DateofProgression=hpv_date_after_death;}
                
                //// Let's push progression into the events Q
                event * CIS_DateofProgressionEvent = new event;
                Events.push_back(CIS_DateofProgressionEvent);
                CIS_DateofProgressionEvent->time = MyPointerToPerson->CIS_DateofProgression;
                CIS_DateofProgressionEvent->p_fun = &EventMyCISStatus;
                CIS_DateofProgressionEvent->person_ID = MyPointerToPerson;
                p_PQ->push(CIS_DateofProgressionEvent);
                
                //cout << MyPointerToPerson << " got infected with HPV on " << MyPointerToPerson->HPV_DateofInfection << ", progressed to CIN1 on " << MyPointerToPerson->CIN1_DateofProgression << " and to CIN2/3 on " << MyPointerToPerson->CIN2_3_DateofProgression << ". Given H = " << h << ", she will now progress from CIN2/3 to CIS on " << MyPointerToPerson->CIS_DateofProgression << endl;
            }
        }
    }
    E(cout << "Somebody with CIN2_3 just progressed to CIS or recovered!" << endl;)
}

void EventMyCISStatus(person *MyPointerToPerson){
    
    E(cout << "Somebody with CIS is about to progress to ICC or recover(?): " << endl;)
    
    if(MyPointerToPerson->Alive == 1){
        MyPointerToPerson->Age= (*p_GT - MyPointerToPerson->DoB);
        
        ///Decide if person recovers or moves on to ICC
        if (MyPointerToPerson->CIS_DateofProgression>0){
            MyPointerToPerson->HPV_Status=HPV_Status_CIS;
            int year = floor(*p_GT);
            double months = floor(((1-(*p_GT-year+0.01))*12));
            
            int j=0;
            float TestICCDate=0;
            std::random_device rd;
            std::mt19937 gen{rd()};
            std::uniform_int_distribution<> dis{0, 3};
            j = dis(gen);
            
            double YearFraction=-999;
            if(months>=1){YearFraction=(RandomMinMax_2(0,months))/12.1;}
            if(months<1){YearFraction=0;}
            double    h = ((double)rand() / (RAND_MAX));
            TestICCDate=(MyPointerToPerson->CIS_DateofProgression+j)+YearFraction;
            
            //// In case they recover from CIS there's a bug! Check with cout below
            if (h>ICC_Rates[1]){
                MyPointerToPerson->ICC_DateofInfection=no_hpv_infection;
                MyPointerToPerson->CIS_DateofRecovery=TestICCDate;
                
                //// Feed recovery into events Q
                event * CIS_DateofRecoveryEvent = new event;
                Events.push_back(CIS_DateofRecoveryEvent);
                CIS_DateofRecoveryEvent->time = MyPointerToPerson->CIS_DateofRecovery;
                CIS_DateofRecoveryEvent->p_fun = &EventMyHPVRecovery;
                CIS_DateofRecoveryEvent->person_ID = MyPointerToPerson;
                p_PQ->push(CIS_DateofRecoveryEvent);
                
                //cout << MyPointerToPerson << " got infected with HPV on " << MyPointerToPerson->HPV_DateofInfection << ", progressed to CIN1 on " << MyPointerToPerson->CIN1_DateofProgression << ", to CIN2/3 on " << MyPointerToPerson->CIN2_3_DateofProgression << " and to CIS on" << MyPointerToPerson->CIS_DateofProgression << ". Given H = " << h << " there is a bug, fix it! (" << MyPointerToPerson->CIS_DateofRecovery << ")" << endl;
            }
            
            //// Get date of progression to ICC
            if (h<=ICC_Rates[1]){
                if (TestICCDate<MyPointerToPerson->DateOfDeath){MyPointerToPerson->ICC_DateofInfection=TestICCDate;}{MyPointerToPerson->CIS_DateofRecovery=-977;}
                if (TestICCDate>=MyPointerToPerson->DateOfDeath){MyPointerToPerson->ICC_DateofInfection=hpv_date_after_death;}
                
                
                //// Feed progression into events Q
                event * ICC_DateofInfectionEvent = new event;
                Events.push_back(ICC_DateofInfectionEvent);
                ICC_DateofInfectionEvent->time = MyPointerToPerson->ICC_DateofInfection;
                ICC_DateofInfectionEvent->p_fun = &EventMyICCStatus;
                ICC_DateofInfectionEvent->person_ID = MyPointerToPerson;
                p_PQ->push(ICC_DateofInfectionEvent);
                
                //cout << MyPointerToPerson << " got infected with HPV on " << MyPointerToPerson->HPV_DateofInfection << ", progressed to CIN1 on " << MyPointerToPerson->CIN1_DateofProgression << ", to CIN2/3 on " << MyPointerToPerson->CIN2_3_DateofProgression << " and to CIS on" << MyPointerToPerson->CIS_DateofProgression << ". Given H = " << h << " she will now progress to ICC on " << MyPointerToPerson->ICC_DateofInfection << ". She spent a total of " << MyPointerToPerson->ICC_DateofInfection - MyPointerToPerson->HPV_DateofInfection << " years with a curable disease!" << " Her HPV status is " << MyPointerToPerson->HPV_Status << endl;
            }
        }
    }
    E(cout << "Somebody with CIS just progressed to ICC!" << endl;)
}


void EventMyICCStatus(person *MyPointerToPerson){
    E(cout << "Somebody with ICC is about to get its HPV status updated: " << endl;)
    if(MyPointerToPerson->Alive == 1){MyPointerToPerson->HPV_Status=HPV_Status_ICC;}
    E(cout << "Somebody with ICC just got its status updated!" << endl;)
}


void EventMyHPVRecovery(person *MyPointerToPerson){
    E(cout << "Somebody is about to recover from a stage of HPV infection: " << endl;)
    if(MyPointerToPerson->Alive == 1){
        if (MyPointerToPerson->HPV_DateofRecovery>0){MyPointerToPerson->HPV_Status=HPV_Status_Recovered;}
        if (MyPointerToPerson->CIN1_DateofRecovery>0){MyPointerToPerson->HPV_Status=HPV_Status_Recovered;}
        if (MyPointerToPerson->CIN2_3_DateofRecovery>0){MyPointerToPerson->HPV_Status=HPV_Status_Recovered;}
        if (MyPointerToPerson->CIS_DateofRecovery>0){MyPointerToPerson->HPV_Status=HPV_Status_Recovered;}}
}
