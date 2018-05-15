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
#include  "HIVInfection.hpp"
#include "HPVInfection.hpp"

using namespace std;

//// --- OUTSIDE INFORMATION --- ////
extern double *p_GT;
extern int EndYear;
extern person** MyArrayOfPointersToPeople;
extern priority_queue<event*, vector<event*>, timeComparison> *p_PQ;
extern vector<event*> Events;

extern int HPV_Status_HPV;
extern int HPV_Status_CIN1;
extern int HPV_Status_CIN2_3;
extern int HPV_Status_CIS;
extern int HPV_Status_ICC;
extern int HPV_Status_Recovered;
extern int HPV_Status_ReInfected;


//// --- Parameters TO MOVE!!!!
double          CIN1_Rates[2]= {0.1,0.6}; // 10% progress to CIN1, 50% clear HPV by 0.5y and an 40% within 1-4y (Rodrigurez et al. 2010; Winer et al., 2011; Gravitt et al., 2017)
double          CIN2_3_Rates[2]= {0.12,0.4}; // 40% progress to CIN2/3 - 12% will progress fast, since the HPV16 prevalence is 11.9% in women with CIN1 in Kenya (HPV centre - 2017)
double          CIS_Rates[2]= {0.75,0.25};
double          ICC_Rates[2]= {1.0,0.0};
extern double   hpv_date_after_death;
extern double   no_hpv_infection;
extern double** HPVarray;
extern double   HPV_Array_buffer;
extern double   HPV_hiv_Array_buffer;


//// --- Important functions --- ////
extern int RandomMinMax_2(int min, int max);
double randfrom_2(double min, double max){
    double range = (max - min);
    double div = RAND_MAX / range;
    return min + (rand() / div);
}




//// --- HPV EVENT --- ////

void EventMyHPVInfection(person *MyPointerToPerson){                    // Function executed when somebody develops HPV infection
    
    E(cout << "Somebody just got infected with HPV and will either progress to CIN1 or recover: " << endl;)
    
    if(MyPointerToPerson->Alive == 1 && (MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount]<1 || MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount]==HPV_Status_ReInfected) && MyPointerToPerson->HPVvaccination_status==0){
        MyPointerToPerson->HPVcount++;
        MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1]=HPV_Status_HPV;
        //double age_now = *p_GT - MyPointerToPerson->DoB;
        
        // Let's get a random test date and a random probability for progression/recovery
        int     j = RandomMinMax_2(1,3);
        float   TestDate=0;
        double YearFraction=-999;
        YearFraction=(RandomMinMax_2(1,12))/12.1;
        TestDate=*p_GT+j+YearFraction;
        double    h = ((double)rand() / (RAND_MAX));
        
        // --- In case they recover from HPV
        if (h>CIN1_Rates[0] && h<CIN1_Rates[1]){  // 50% of clearance by 6 months (Rdz et al., 2010; Winer et al., 2011; Gravitt et al., 2017)
            MyPointerToPerson->HPV_DateofRecovery[MyPointerToPerson->HPVcount-1]=*p_GT+0.5;
            MyPointerToPerson->HPV_StatusAtRecovery[MyPointerToPerson->HPVcount-1]=1; // Clean later; 1 = recovered from HPV
            // Push recovery into the events Q
            event * HPV_DateofRecoveryEvent = new event;
            Events.push_back(HPV_DateofRecoveryEvent);
            HPV_DateofRecoveryEvent->time = MyPointerToPerson->HPV_DateofRecovery[MyPointerToPerson->HPVcount-1];
            HPV_DateofRecoveryEvent->p_fun = &EventMyHPVRecovery;
            HPV_DateofRecoveryEvent->person_ID = MyPointerToPerson;
            p_PQ->push(HPV_DateofRecoveryEvent);
        }
        
        if (h>CIN1_Rates[1]){ // The remaining 40% of clearance will occur within 4 years, although 7% will have an LTP (Castle et al., 2011)
            int r = randfrom_2(0,1);
            if (r<=0.07){MyPointerToPerson->HPV_DateofRecovery[MyPointerToPerson->HPVcount-1]=*p_GT+7;}
            if(r>0.07){MyPointerToPerson->HPV_DateofRecovery[MyPointerToPerson->HPVcount-1]=TestDate;}
            MyPointerToPerson->HPV_StatusAtRecovery[MyPointerToPerson->HPVcount-1]=1; // Clean later; 1 = recovered from HPV
            // Push recovery into the events Q
            event * HPV_DateofRecoveryEvent = new event;
            Events.push_back(HPV_DateofRecoveryEvent);
            HPV_DateofRecoveryEvent->time = MyPointerToPerson->HPV_DateofRecovery[MyPointerToPerson->HPVcount-1];
            HPV_DateofRecoveryEvent->p_fun = &EventMyHPVRecovery;
            HPV_DateofRecoveryEvent->person_ID = MyPointerToPerson;
            p_PQ->push(HPV_DateofRecoveryEvent);
        }
        
        // --- In case they progress to CIN1
        if (h<=CIN1_Rates[0]){
            double  j = randfrom_2(6.5,7.3);      // --- Progression will ocur between 6.5 and 7.3 years (Schlecht et al., 2003)
            MyPointerToPerson->HPV_DateofProgression[MyPointerToPerson->HPVcount-1]=*p_GT+j;
            // Push progression into the events Q
            event * CIN1_DateofProgressionEvent = new event;
            Events.push_back(CIN1_DateofProgressionEvent);
            CIN1_DateofProgressionEvent->time = MyPointerToPerson->HPV_DateofProgression[MyPointerToPerson->HPVcount-1];
            CIN1_DateofProgressionEvent->p_fun = &EventMyCIN1Status;
            CIN1_DateofProgressionEvent->person_ID = MyPointerToPerson;
            p_PQ->push(CIN1_DateofProgressionEvent);
        }
    }
    E(cout << "Somebody with HPV just progressed to CIN1 or recovered!" << endl;)
}

void EventMyCIN1Status(person *MyPointerToPerson){
    
    E(cout << "Somebody with CIN1 is about to progress to CIN2/3 or recover: " << endl;)
    
    if(MyPointerToPerson->Alive == 1 && MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1]==HPV_Status_HPV){
        
        MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1]=HPV_Status_CIN1;
        
        int j=0;
        float TestDate=0;
        random_device rd;
        mt19937 gen{rd()};
        uniform_int_distribution<> dis{1, 2};
        j = dis(gen);
        double YearFraction=-999;
        YearFraction=(RandomMinMax_2(1,12))/12.1;
        double    h = ((double)rand() / (RAND_MAX));
        TestDate=*p_GT+j+YearFraction;
        
        // --- In case they recover from CIN1
        if (h>CIN2_3_Rates[1]){
            MyPointerToPerson->HPV_DateofRecovery[MyPointerToPerson->HPVcount-1]=TestDate;
            MyPointerToPerson->HPV_StatusAtRecovery[MyPointerToPerson->HPVcount-1]=2;      // Clean later; 2 = recovered from CIN1
            // Push recovery into the events Q
            event * CIN1_DateofRecoveryEvent = new event;
            Events.push_back(CIN1_DateofRecoveryEvent);
            CIN1_DateofRecoveryEvent->time = MyPointerToPerson->HPV_DateofRecovery[MyPointerToPerson->HPVcount-1];
            CIN1_DateofRecoveryEvent->p_fun = &EventMyHPVRecovery;
            CIN1_DateofRecoveryEvent->person_ID = MyPointerToPerson;
            p_PQ->push(CIN1_DateofRecoveryEvent);
        }
        
        // --- In case they progress to CIN2/3
        if (h>CIN2_3_Rates[0] && h<=CIN2_3_Rates[1]){
            double  j=randfrom_2(7,7.85);      // --- Progression will ocur between 7.0 and 7.85 years (Schlecht et al., 2003)
            MyPointerToPerson->HPV_DateofProgression[MyPointerToPerson->HPVcount-1]=*p_GT+j;
            // Push progression into the events Q
            event * CIN2_3_DateofProgressionEvent = new event;
            Events.push_back(CIN2_3_DateofProgressionEvent);
            CIN2_3_DateofProgressionEvent->time = MyPointerToPerson->HPV_DateofProgression[MyPointerToPerson->HPVcount-1];
            CIN2_3_DateofProgressionEvent->p_fun = &EventMyCIN2_3Status;
            CIN2_3_DateofProgressionEvent->person_ID = MyPointerToPerson;
            p_PQ->push(CIN2_3_DateofProgressionEvent);
        }
        
        if (h<=CIN2_3_Rates[0]){
            MyPointerToPerson->HPV_DateofProgression[MyPointerToPerson->HPVcount-1]=*p_GT;
            // Push progression into the events Q
            event * CIN2_3_DateofProgressionEvent = new event;
            Events.push_back(CIN2_3_DateofProgressionEvent);
            CIN2_3_DateofProgressionEvent->time = MyPointerToPerson->HPV_DateofProgression[MyPointerToPerson->HPVcount-1];
            CIN2_3_DateofProgressionEvent->p_fun = &EventMyCIN2_3Status;
            CIN2_3_DateofProgressionEvent->person_ID = MyPointerToPerson;
            p_PQ->push(CIN2_3_DateofProgressionEvent);
        }
    }
    E(cout << "Somebody with CIN1 just progressed to CIN2/3 or recovered!" << endl;)
}

void EventMyCIN2_3Status(person *MyPointerToPerson){
    
    E(cout << "Somebody with CIN2_3 is about to progress to CIS or recover: " << endl;)
    
    if(MyPointerToPerson->Alive == 1 && MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1]==HPV_Status_CIN1){
        MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1]=HPV_Status_CIN2_3;
        
        // Let's get a random test date and a random probability for progression/recovery
        float TestDate=0;
        double YearFraction=-999;
        YearFraction=(RandomMinMax_2(1,12))/12.1;
        double    h = ((double)rand() / (RAND_MAX));
        TestDate=*p_GT+YearFraction;
        
        // --- In case they recover from CIN2/3
        if (h>CIS_Rates[0]){
            MyPointerToPerson->HPV_DateofRecovery[MyPointerToPerson->HPVcount-1]=TestDate;
            MyPointerToPerson->HPV_StatusAtRecovery[MyPointerToPerson->HPVcount-1]=3;      // Clean later; 3 = recovered from CIN2/3
            // --- Push recovery into the events Q
            event * CIN2_3_DateofRecoveryEvent = new event;
            Events.push_back(CIN2_3_DateofRecoveryEvent);
            CIN2_3_DateofRecoveryEvent->time = MyPointerToPerson->HPV_DateofRecovery[MyPointerToPerson->HPVcount-1];
            CIN2_3_DateofRecoveryEvent->p_fun = &EventMyHPVRecovery;
            CIN2_3_DateofRecoveryEvent->person_ID = MyPointerToPerson;
            p_PQ->push(CIN2_3_DateofRecoveryEvent);
        }
        
        // --- In case they progress to CIS
        if (h<=CIS_Rates[0]){
            MyPointerToPerson->HPV_DateofProgression[MyPointerToPerson->HPVcount-1]=TestDate;
            // Push progression into the events Q
            event * CIS_DateofProgressionEvent = new event;
            Events.push_back(CIS_DateofProgressionEvent);
            CIS_DateofProgressionEvent->time = MyPointerToPerson->HPV_DateofProgression[MyPointerToPerson->HPVcount-1];
            CIS_DateofProgressionEvent->p_fun = &EventMyCISStatus;
            CIS_DateofProgressionEvent->person_ID = MyPointerToPerson;
            p_PQ->push(CIS_DateofProgressionEvent);
        }
    }
    E(cout << "Somebody with CIN2_3 just progressed to CIS or recovered!" << endl;)
}

void EventMyCISStatus(person *MyPointerToPerson){
    
    E(cout << "Somebody with CIS is about to progress to ICC or recover(?): " << endl;)
    
    if(MyPointerToPerson->Alive == 1 && MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1]==HPV_Status_CIN2_3){
        MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1]=HPV_Status_CIS;
        MyPointerToPerson->HPV_DateofRecovery[MyPointerToPerson->HPVcount-1]=2040;
        MyPointerToPerson->HPV_DateOfCIS=*p_GT;
        
        // Let's get a random test date and a random probability for progression/recovery
        float TestDate=0;
        double YearFraction=-999;
        YearFraction=(RandomMinMax_2(1,12))/12.1;
        double    h = ((double)rand() / (RAND_MAX));
        TestDate=*p_GT+YearFraction;
        
        //// In case they recover from CIS there's a bug! Check with cout below
        if (h>ICC_Rates[0]){
            MyPointerToPerson->HPV_DateofRecovery[MyPointerToPerson->HPVcount-1]=TestDate;
            MyPointerToPerson->HPV_StatusAtRecovery[MyPointerToPerson->HPVcount-1]=4;  // Clean later; 4 = recovered from CIS
            //// Feed recovery into events Q
            event * CIS_DateofRecoveryEvent = new event;
            Events.push_back(CIS_DateofRecoveryEvent);
            CIS_DateofRecoveryEvent->time = MyPointerToPerson->HPV_DateofRecovery[MyPointerToPerson->HPVcount-1];
            CIS_DateofRecoveryEvent->p_fun = &EventMyHPVRecovery;
            CIS_DateofRecoveryEvent->person_ID = MyPointerToPerson;
            p_PQ->push(CIS_DateofRecoveryEvent);
            cout << "There's a bug; nobody should recover from CIS" << endl;
        }
        
        //// Get date of progression to ICC
        if (h<=ICC_Rates[0]){
            MyPointerToPerson->HPV_DateofProgression[MyPointerToPerson->HPVcount-1]=TestDate;
            //// Feed progression into events Q
            event * ICC_DateofProgressionEvent = new event;
            Events.push_back(ICC_DateofProgressionEvent);
            ICC_DateofProgressionEvent->time = MyPointerToPerson->HPV_DateofProgression[MyPointerToPerson->HPVcount-1];
            ICC_DateofProgressionEvent->p_fun = &EventMyICCStatus;
            ICC_DateofProgressionEvent->person_ID = MyPointerToPerson;
            p_PQ->push(ICC_DateofProgressionEvent);
        }
    }
    E(cout << "Somebody with CIS just progressed to ICC!" << endl;)
}


void EventMyICCStatus(person *MyPointerToPerson){
    E(cout << "Somebody with ICC is about to get its HPV status updated: " << endl;)
    if(MyPointerToPerson->Alive == 1 && MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1]==HPV_Status_CIS){
        MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1]=HPV_Status_ICC;
        
        double  j=randfrom_2(0, 5);
        double  TestDeathDate=*p_GT+j;
        if(TestDeathDate<MyPointerToPerson->DateOfDeath){
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
                MyPointerToPerson->CauseOfDeath=16;
            }
        }
    }
    E(cout << "Somebody with ICC just got its status updated!" << endl;)
}


void EventMyHPVRecovery(person *MyPointerToPerson){
    E(cout << "Somebody is about to recover from a stage of HPV infection: " << endl;)
    //cout << "Status " << MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1] << endl;
    
    if(MyPointerToPerson->Alive == 1 && MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1]<HPV_Status_Recovered){
        MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1]=HPV_Status_Recovered;
        
        double age_now = *p_GT-MyPointerToPerson->DoB;
        
        // --- An HIV negative woman will be re-infected/re-activate with HPV
        if (MyPointerToPerson->HPVcount<=4 && (MyPointerToPerson->HIV<0 || MyPointerToPerson->HIV>*p_GT)){
            double TestDate=-997;
            int f = floor(age_now);
            double h = randfrom_2(HPVarray[1][f]*HPV_Array_buffer,HPVarray[1][65]*HPV_Array_buffer);
            double r = ((double) rand() / (RAND_MAX));
            if (r<=0.4274){
                int i=0;
                while (h>(HPVarray[1][i])*HPV_Array_buffer){i++;}
                double YearFraction=(RandomMinMax_2(1,12))/ 12.1;
                TestDate = MyPointerToPerson->DoB + i + YearFraction;
                if(TestDate >= *p_GT){
                    MyPointerToPerson->HPV[MyPointerToPerson->HPVcount] = TestDate;
                    MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount] = HPV_Status_ReInfected;
                    event * HPV_ReInfection = new event;
                    Events.push_back(HPV_ReInfection);
                    HPV_ReInfection->time = MyPointerToPerson->HPV[MyPointerToPerson->HPVcount];
                    HPV_ReInfection->p_fun = &EventMyHPVInfection;
                    HPV_ReInfection->person_ID = MyPointerToPerson;
                    p_PQ->push(HPV_ReInfection);
                }
            }
        }
        // --- An HIV positive woman will be re-infected/re-activate with HPV
        MyPointerToPerson->Age = (*p_GT-MyPointerToPerson->DoB);
        if (MyPointerToPerson->HIV<=*p_GT && MyPointerToPerson->HIV>0 && MyPointerToPerson->HPVcount<=4){
            
            double TestDate=-999;
            int f = floor(age_now);
            double h = randfrom_2(HPVarray[0][f]*HPV_hiv_Array_buffer,HPVarray[0][f+4]*HPV_hiv_Array_buffer);
            
            double r = ((double) rand() / (RAND_MAX));
            if (r<=0.8){   // HPVarray[2][65]*HPV_hiv_Array_buffer
                int i=0;
                while (h>(HPVarray[0][i]*HPV_hiv_Array_buffer)){i++;}
                double YearFraction=(RandomMinMax_2(1,12))/ 12.1;
                TestDate = MyPointerToPerson->DoB + i + YearFraction;
                if(TestDate >= *p_GT){
                    MyPointerToPerson->HPV[MyPointerToPerson->HPVcount] = TestDate;
                    MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount] = HPV_Status_ReInfected;
                    event * HPV_ReInfection = new event;
                    Events.push_back(HPV_ReInfection);
                    HPV_ReInfection->time = MyPointerToPerson->HPV[MyPointerToPerson->HPVcount];
                    HPV_ReInfection->p_fun = &EventMyHPVInfection;
                    HPV_ReInfection->person_ID = MyPointerToPerson;
                    p_PQ->push(HPV_ReInfection);
                }
            }
        }
    }
}
