////////////////////////////////////////////////////////////////
//    Mikaela Epidemiology Model							  //
//    Created by Mikaela Smit on 7/11/2014.				      //
//    Copyright (c) 2014 Mikaela Smit. All rights reserved.   //
////////////////////////////////////////////////////////////////

#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <ctime>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <cmath>
#include <math.h>
#include "person.h"
#include "event.h"
#include "eventQ.h"
#include "eventfunctions.h"
#include "errorcoutmacro.h"
#include "LoadParams.h"
#include "CParamReader.hpp"
#include "CountryParams.hpp"
#include "Intervention.hpp"
#include "NCD.hpp"

using namespace std;

// TIDY UP
// 1. Clean up NCD risk from NCD.hpp
// 2. Clean up HIV.cpp
// Clean up person.cpp
// Add ART events for HPV vacciantion for ART
// Clean eventfunction
// clean intervemtiom
// make sure that we DO book in events even if they happen after death!!

//// --- 1. Choose country --- ////
int country=1;                  // 1=KENYA      2=ZIMBABWE      3=MALAWI      4=KENYA - UG

//// --- 2. Update file paths --- ////
string InputFileDirectory="/Users/Monkeyface/Dropbox/Ageing in Kenya and Zimbabwe - project/Model Clean_Miki/HIVModelZimbabwe";
string OutputFileDirectory="/Users/Monkeyface/Dropbox/Ageing in Kenya and Zimbabwe - project/Model Clean_Miki/Kenya_File.csv";

//// --- 3. Interventions Switch and Details --- ////
int yearintervention_start=2018;
int int_HPVvaccination=0;
int age_HPVvaccination=9;
int int_CVDIntervention=0;
int int_CCscreening=0;


//// --- 4. Central Model Parameters --- ////
double StartYear=1950;
int    EndYear=2035;
int    factor=100;                  // Ffraction of population run in the mode
const long long int final_number_people=100000000;


double Risk_DiabHC=1.12;            // NCD parameters
double Risk_DiabHT=1.4;
double Risk_DiabCKD=1.5;
double Risk_DiabCVD=2.31;
double Risk_HCHT=1.277;
double Risk_HCCVD=1.41;
double Risk_HTCKD=1.69;
double Risk_HTCVD=1.26;

double Risk_NCD_Diabetes[5]={Risk_DiabHT, Risk_DiabCVD, Risk_DiabCKD, Risk_DiabCVD, Risk_DiabHC};
int relatedNCDs_Diab[5]={0, 3, 5, 6, 7};
int nr_NCD_Diab=sizeof(relatedNCDs_Diab)/sizeof(relatedNCDs_Diab[0]);

double Risk_NCD_HT[3]={Risk_HTCVD, Risk_HTCKD, Risk_HTCVD};
int relatedNCDs_HT[3]={3, 5, 6};
int nr_NCD_HT=sizeof(relatedNCDs_HT)/sizeof(relatedNCDs_HT[0]);

double Risk_NCD_HC[3]={Risk_HCHT, Risk_HCCVD, Risk_HCCVD};
int relatedNCDs_HC[3]={0, 3, 6};
int nr_NCD_HC=sizeof(relatedNCDs_HC)/sizeof(relatedNCDs_HC[0]);


double Risk_HIVHT=1.49;                       // HIV increases NCD risk  - from Schouten et al CID 2016
double Risk_HIVCKD=2;
double Risk_NCDHIV[2]={Risk_HIVHT, Risk_HIVCKD};
int relatedNCD_HIV[2]={0, 5};
int nr_NCD_HIV=sizeof(relatedNCD_HIV)/sizeof(relatedNCD_HIV[0]);
extern int ageAdult;


double MortRisk[6]= {0, 0, 0.85, 1.3, 1.1, 0.8}; // Mortality parameters//{0.087, 0, 1.4, 670.87, 12.23, 5};
double MortRisk_Cancer[5]= {1, 1, 1, 1, 1.05};
// Original values from Smith et al Factors associated with : 1.52 (HT), 1.77 (diabetes)
//{0.087, 0, 1.4, 670.87, 12.23};   // Both this and above needs to be fitted


//// --- Pointers and external information
double      *p_GT;
int         *p_PY;
int         PY=0;
int         age_tostart_CCscreening;
extern int  total_population;
priority_queue<event*, vector<event*>, timeComparison> *p_PQ;
person** MyArrayOfPointersToPeople = new person*[final_number_people];
vector<event *> Events;

extern double background_d;        // Mortality parameters
extern double HIV_d;
extern double IHD_d;
extern double Depression_d;
extern double Asthma_d;
extern double Stroke_d;
extern double Diabetes_d;
extern double CKD_d;
extern double Colo_d;
extern double Liver_d;
extern double Oeso_d;
extern double Prostate_d;
extern double OtherCan_d;
extern double MortAdj;







int main(){
    
    srand(time(NULL));														 // Random Number generator using PC time
    
    cout << endl << "Jambo / Hello / Hola!" << endl << endl ;                // Check if model is running
    
    
    //// --- Load parameters --- ////
    cout << "Section 1 - Loading country related items.  " << endl << endl;
    getParamsString(country);                                                // Gets the ParamDirectory string to pull country-specific files
    loadCountryParams(country);                                              // Call the function that loads country-specific parameters
    
    
    cout << "Section 2 - We are loading the arrays.  " << endl;
    // Load HIV Arrays
    loadCD4StartArray();
    loadCD4ProgArray();
    loadCD4DeathArray();
    loadCD4ARTArray();
    
    // Load ART Arrays
    loadARTKidsArray();
    loadARTWomenArray();
    loadARTMenArray();
    loadARTMen_sumArray();
    loadARTWomen_sumArray();
    
    // Load Demographic Arrays
    loadAgeDistribution();
    loadAgeMin();
    loadAgeMax();
    loadNrChildren();
    loadNrChildrenProb();
    
    // Load Large Arrays
    loadBirthArray();
    loadDeathArray_Women();
    loadDeathArray_Men();
    loadHIVArray_Women();
    loadHIVArray_Men();
    loadNCDArray();
    loadCancerArray();
    loadHPVarray();
    
    cout << "Section 3 - All arrays and country-specific items were loaded successfully.  " << endl;
    
    
    //// ---- Warning Code --- ////
    E(cout << "NOTE: The Macro for error code is working" << endl << endl;)
    if (StartYear > 1954){cout << "CAREFUL!! Fertility at the beginning is based on a TRF_init of 0, based on GlobalTime being between 1950 and 1954" << endl << endl;}
    if (EndYear > 2010){cout << "CAREFUL!! The model runs after 2010 and output things need to be updated to avoid errors" << endl << endl;}
    
    
    //// --- Some pointers ---- ////
    double GlobalTime=StartYear;											// Define Global Time and set it to 0 at the beginning of the model
    p_GT=&GlobalTime;														// Define the location the pointer to Global time is pointing to
    
    priority_queue<event*, vector<event*>, timeComparison> iQ;				// Define the Priority Q
    p_PQ=&iQ;																// Define pointer to event Q
    p_PY=&PY;
    
    
    //// --- MAKING POPULATION--- ////
    
    cout << "Section 4 - We're going to create a population.  " << endl;
    
    for(int i=0; i<total_population; i++){									// REMEMBER: this needs to stay "final_number_people" or it will give error with CSV FILES!!!!
        MyArrayOfPointersToPeople[i]=new person();							// The 'new person' the actual new person
        int a=i;
        (MyArrayOfPointersToPeople[i])->PersonIDAssign(a);					// --- Assign PersonID ---
    }
    // cout << "We have assigned a new ID to " << endl;
    
    for(int i=0; i<total_population; i++){
        (MyArrayOfPointersToPeople[i])->Alive=1;							// --- Assign Life Status ---
        (MyArrayOfPointersToPeople[i])->GenderDistribution();				// --- Assign Sex- ---
        (MyArrayOfPointersToPeople[i])->GetMyDoB();							// --- Assign DoB/Age ---
        (MyArrayOfPointersToPeople[i])->GetDateOfDeath();					// --- Assign date of death ---
        (MyArrayOfPointersToPeople[i])->GetMyDateOfHPVInfection();            // --- Assign date of HPV infection ---
        if (MyArrayOfPointersToPeople[i]->Sex == 2 && MyArrayOfPointersToPeople[i]->Age<50 && MyArrayOfPointersToPeople[i]->AgeAtDeath>=15) {(MyArrayOfPointersToPeople[i])->GetDateOfBaby();}              // --- Assign Birth of all Children- ---
        (MyArrayOfPointersToPeople[i])->GetMyDateNCD();                     // --- Get date of NCDs ---
        (MyArrayOfPointersToPeople[i])->GetMyDateCancers();                     // --- Get date of NCDs ---
        (MyArrayOfPointersToPeople[i])->GetMyDateOfHIVInfection();           // ---Get date of HIV infection ---
    }
    
    cout << "Section 5 - We've finished creating a population.  " << endl;
    
    
    //// --- EVENTQ --- ////
    cout << "Section 6 - We are going to create key events.  " << endl;
    
    event * TellNewYear = new event;										//// --- Tell me every time  a new year start --- ////
    Events.push_back(TellNewYear);
    TellNewYear->time = StartYear;
    TellNewYear->p_fun = &EventTellNewYear;
    iQ.push(TellNewYear);
    
    /*event * InterventionEvent = new event;                                  //// --- Push in interventions ---////
    Events.push_back(InterventionEvent);
    InterventionEvent->time = yearintervention_start;
    InterventionEvent->p_fun = &EventStartIntervention;
    iQ.push(InterventionEvent);
    
    event * CC_First_screen = new event;                                    //// --- Cervical Cancer screening --- ///
    Events.push_back(CC_First_screen);
    CC_First_screen->time = 2018;
    CC_First_screen->p_fun = &EventMyFirst_VIA_Screening;
    iQ.push(CC_First_screen);
    
    
    event * PreCVDInterventionEvent = new event;                                  //// --- Push in CVD intervention ---////
    Events.push_back(PreCVDInterventionEvent);
    PreCVDInterventionEvent->time = yearintervention_start;
    PreCVDInterventionEvent->p_fun = &EventCVDPrevIntervention;
    iQ.push(PreCVDInterventionEvent);
    */
    
    //// --- LETS RUN THE EVENTQ --- ////
    cout << endl << endl << "The characteristics of the event queue:" << endl;
    cout << "The first event will ocurr in " << iQ.top()->time << ".  " << endl;
    cout << "The size of the event queue is " << iQ.size() << endl;
    
    while(iQ.top()->time< EndYear +1){                                      // this loop throws up error because no recurrent birthday pushing gt over 5 yrs and iq.pop means gt cannot be updated after pop
        GlobalTime=iQ.top()->time;											// careful with order of global time update - do not touch or touch and check!!
        iQ.top()-> p_fun(iQ.top()->person_ID);
        iQ.pop();
    }
    
    
    //// --- Output the results in a csv file --- ////
    cout << "Section 7 - We are going to save the output. " << endl;
    FILE* ProjectZim;
    ProjectZim = fopen(OutputFileDirectory.c_str(),"w");
    

     for (int i=0; i<total_population; i++) {								// Note: If adding more variables to be output, need to adapt the %x
         fprintf(ProjectZim,"%d, %d, %f, %f, %d, %d, %f, %d, %f, %d, %d, %f, %f, %f, %f, %f, %d, %f, %f, %f, %f, %f, %f, %f, %d, %f, %f, %f, %f, %f, %f, %f, %f, %d, %f, %f, %d, %f, %f, %f, %f, %f, %d, %d, %d, %d, %d, %f, %f, %f, %f, %f, %d, %d, %d, %d, %d, %f, %f, %f, %f, %f, %f, %d, %d, %f, %f, %f, %f \n",
                 MyArrayOfPointersToPeople[i]->PersonID,
                 MyArrayOfPointersToPeople[i]->Sex,
                 MyArrayOfPointersToPeople[i]->DoB,
                 MyArrayOfPointersToPeople[i]->Age,
                 MyArrayOfPointersToPeople[i]->MotherID,                         // Number 5
                 MyArrayOfPointersToPeople[i]->DatesBirth.size(),
                 MyArrayOfPointersToPeople[i]->DateOfDeath,
                 MyArrayOfPointersToPeople[i]->AgeAtDeath,
                 MyArrayOfPointersToPeople[i]->HIV,
                 MyArrayOfPointersToPeople[i]->CD4_cat,                          // Number 10
                 MyArrayOfPointersToPeople[i]->ART,
                 MyArrayOfPointersToPeople[i]->HT,
                 MyArrayOfPointersToPeople[i]->Depression,
                 MyArrayOfPointersToPeople[i]->Asthma,
                 MyArrayOfPointersToPeople[i]->Stroke,                           // Number 15
                 MyArrayOfPointersToPeople[i]->Diabetes,
                 MyArrayOfPointersToPeople[i]->CD4_cat_ARTstart,
                 MyArrayOfPointersToPeople[i]->CD4_change.at(0),
                 MyArrayOfPointersToPeople[i]->CD4_change.at(1),
                 MyArrayOfPointersToPeople[i]->CD4_change.at(2),                 // Number 20
                 MyArrayOfPointersToPeople[i]->CD4_change.at(3),
                 MyArrayOfPointersToPeople[i]->CD4_change.at(4),
                 MyArrayOfPointersToPeople[i]->CD4_change.at(5),
                 MyArrayOfPointersToPeople[i]->CD4_change.at(6),
                 MyArrayOfPointersToPeople[i]->CauseOfDeath,                     // Number 25
                 MyArrayOfPointersToPeople[i]->CKD,
                 MyArrayOfPointersToPeople[i]->Breast,
                 MyArrayOfPointersToPeople[i]->Cervical,
                 MyArrayOfPointersToPeople[i]->Colo,
                 MyArrayOfPointersToPeople[i]->Liver,                            // Number 30
                 MyArrayOfPointersToPeople[i]->Oeso,
                 MyArrayOfPointersToPeople[i]->Prostate,
                 MyArrayOfPointersToPeople[i]->OtherCan,
                 MyArrayOfPointersToPeople[i]->Stroke_status,             // Check if used and, if not, remove
                 MyArrayOfPointersToPeople[i]->MI,                               // Number 35
                 MyArrayOfPointersToPeople[i]->HC,
                 MyArrayOfPointersToPeople[i]->HPVcount,
                 MyArrayOfPointersToPeople[i]->HPV[0],
                 MyArrayOfPointersToPeople[i]->HPV[1],
                 MyArrayOfPointersToPeople[i]->HPV[2],                               // Number 40
                 MyArrayOfPointersToPeople[i]->HPV[3],
                 MyArrayOfPointersToPeople[i]->HPV[4],
                 MyArrayOfPointersToPeople[i]->HPV_Status[0],
                 MyArrayOfPointersToPeople[i]->HPV_Status[1],
                 MyArrayOfPointersToPeople[i]->HPV_Status[2],                        // Number 45
                 MyArrayOfPointersToPeople[i]->HPV_Status[3],
                 MyArrayOfPointersToPeople[i]->HPV_Status[4],
                 MyArrayOfPointersToPeople[i]->HPV_DateofRecovery[0],
                 MyArrayOfPointersToPeople[i]->HPV_DateofRecovery[1],
                 MyArrayOfPointersToPeople[i]->HPV_DateofRecovery[2],                               // Number 50
                 MyArrayOfPointersToPeople[i]->HPV_DateofRecovery[3],
                 MyArrayOfPointersToPeople[i]->HPV_DateofRecovery[4],
                 MyArrayOfPointersToPeople[i]->HPV_StatusAtRecovery[0],
                 MyArrayOfPointersToPeople[i]->HPV_StatusAtRecovery[1],
                 MyArrayOfPointersToPeople[i]->HPV_StatusAtRecovery[2],                               // Number 55
                 MyArrayOfPointersToPeople[i]->HPV_StatusAtRecovery[3],
                 MyArrayOfPointersToPeople[i]->HPV_StatusAtRecovery[4],
                 MyArrayOfPointersToPeople[i]->HPV_DateofProgression[0],
                 MyArrayOfPointersToPeople[i]->HPV_DateofProgression[1],
                 MyArrayOfPointersToPeople[i]->HPV_DateofProgression[2],                              // Number 60
                 MyArrayOfPointersToPeople[i]->HPV_DateofProgression[3],
                 MyArrayOfPointersToPeople[i]->HPV_DateofProgression[4],
                 MyArrayOfPointersToPeople[i]->VIA,
                 MyArrayOfPointersToPeople[i]->VIAcount,
                 MyArrayOfPointersToPeople[i]->HPVvaccination_status,                        // Number 65
                 MyArrayOfPointersToPeople[i]->HPVvaccination_date,
                 MyArrayOfPointersToPeople[i]->HPV_DateOfCIS,
                 MyArrayOfPointersToPeople[i]->CVD_HT_Treat_Date,
                 MyArrayOfPointersToPeople[i]->CVD_HC_Treat_Date
                 );}
    fclose(ProjectZim);
    
    
    // COUNT OUTPUT FOR FITTING
    cout << "Section 7 - Mortality fitting output. " << endl;
    int count_2016deaths=0;
    int count_causeofdeath[14]={0};
    
    for (int i=0; i<total_population; i++) {
        if (MyArrayOfPointersToPeople[i]->DateOfDeath>=2013 && MyArrayOfPointersToPeople[i]->DateOfDeath<2014)
        {
            count_2016deaths++;
            count_causeofdeath[MyArrayOfPointersToPeople[i]->CauseOfDeath-1]++;
        }
    }
    
    // Calculate percentages
    cout << "FINAL: total: " << count_2016deaths << " and nr background " << count_causeofdeath[0] << endl;
    
    double background_m =(count_causeofdeath[0]/(double)count_2016deaths)*100;
    double HIV_m        =(count_causeofdeath[1]/(double)count_2016deaths)*100;
    double HT_m         =(count_causeofdeath[2]/(double)count_2016deaths)*100;
    double Depression_m =(count_causeofdeath[3]/(double)count_2016deaths)*100;
    double Asthma_m     =(count_causeofdeath[4]/(double)count_2016deaths)*100;
    double Stroke_m     =(count_causeofdeath[5]/(double)count_2016deaths)*100;
    double Diabetes_m   =(count_causeofdeath[6]/(double)count_2016deaths)*100;
    double CKD_m        =(count_causeofdeath[7]/(double)count_2016deaths)*100;
    double Colo_m       =(count_causeofdeath[8]/(double)count_2016deaths)*100;
    double Liver_m      =(count_causeofdeath[9]/(double)count_2016deaths)*100;
    double Oeso_m       =(count_causeofdeath[10]/(double)count_2016deaths)*100;
    double Prostate_m    =(count_causeofdeath[11]/(double)count_2016deaths)*100;
    double OtherCan_m   =(count_causeofdeath[12]/(double)count_2016deaths)*100;
    double HPV_Status_m        =(count_causeofdeath[13]/(double)count_2016deaths)*100;
    double HPV_DateofInfection_m        =(count_causeofdeath[13]/(double)count_2016deaths)*100;
    double CIN1_DateofProgression_m        =(count_causeofdeath[13]/(double)count_2016deaths)*100;
    double CIN2_3_DateofProgression_m        =(count_causeofdeath[13]/(double)count_2016deaths)*100;
    double CIS_DateofProgression_m        =(count_causeofdeath[13]/(double)count_2016deaths)*100;
    double ICC_DateofInfection_m        =(count_causeofdeath[13]/(double)count_2016deaths)*100;
    double HPV_DateofRecovery_m        =(count_causeofdeath[13]/(double)count_2016deaths)*100;
    double CIN1_DateofRecovery_m        =(count_causeofdeath[13]/(double)count_2016deaths)*100;
    double CIN2_3_DateofRecovery_m        =(count_causeofdeath[13]/(double)count_2016deaths)*100;
    double CIS_DateofRecovery_m        =(count_causeofdeath[13]/(double)count_2016deaths)*100;
    
    // Output the model percentages
    cout << "Background " << background_m << endl;
    cout << "HIV "        << HIV_m << endl;
    cout << "HT "         << HT_m << endl;
    cout << "Depression " << Depression_m << endl;
    cout << "Asthma "     << Asthma_m << endl;
    cout << "Stroke "     << Stroke_m << endl;
    cout << "Diabetes "   << Diabetes_m << endl;
    cout << "CKD "        << CKD_m << endl;
    cout << "Colo "       << Colo_m << endl;
    cout << "Liver "      << Liver_m << endl;
    cout << "Oeso "       << Oeso_m << endl;
    cout << "Prostate "    << Prostate_m << endl;
    cout << "OtherCan "   << OtherCan_m << endl;
    cout << "HPV_Status_m "        << HPV_Status_m << endl;
    cout << "HPV_DateofInfection "        << HPV_DateofInfection_m << endl;
    cout << "CIN1_DateofProgression "        << CIN1_DateofProgression_m << endl;
    cout << "CIN2_3_DateofProgression "        << CIN2_3_DateofProgression_m << endl;
    cout << "CIS_DateofProgression "        << CIS_DateofProgression_m << endl;
    cout << "ICC_DateofInfection "        << ICC_DateofInfection_m << endl;
    cout << "HPV_DateofRecovery "        << HPV_DateofRecovery_m << endl;
    cout << "CIN1_DateofRecovery "        << CIN1_DateofRecovery_m << endl;
    cout << "CIN2_3_DateofRecovery "        << CIN2_3_DateofRecovery_m << endl;
    cout << "CIS_DateofRecovery "        << CIS_DateofRecovery_m << endl;
    
    
    // Least square calculation
    double sum_MinLik=  pow ((background_m  - background_d),2) +
    pow ((HIV_m         - HIV_d),2) +
    pow ((HT_m          - IHD_d),2) +
    pow ((Depression_m  - Depression_d),2) +
    pow ((Asthma_m      - Asthma_d),2) +
    pow ((Stroke_m      - Stroke_d),2) +
    pow ((Diabetes_m    - Diabetes_d),2) +
    pow ((CKD_m         - CKD_d),2) +
    pow ((Colo_m        - Colo_d),2) +
    pow ((Liver_m       - Liver_d),2) +
    pow ((Oeso_m        - Oeso_d),2) +
    pow ((Prostate_m     - Prostate_d),2) +
    pow ((OtherCan_m    - OtherCan_d),2) ;
    
    
    cout << "Least Square " << sum_MinLik << endl;
    

    //// --- LETS AVOID MEMORY LEAKS AND DELETE ALL NEW EVENTS --- ////
    cout << "Lets delete the heap! " << endl;
    for(int i=0; i<Events.size()-1; i++){
        E(cout << "Event " << Events.size() << " is " << Events.at(i)->time << endl;)
        delete Events.at(i);
        E(cout << "Event " << Events.size() << " is " << Events.at(i)->time << endl;)
    }
    
    for(int i=0; i<total_population; i++){									// REMEMBER: this needs to stay "final_number_people" or it will give error with CSV FILES!!!!
        delete MyArrayOfPointersToPeople[i];								// The 'new person' the actual new person
    }
    
    // --- End of code ---
    cout << endl << "Kwaheri / Goodbye / Adios!\n";
    system("pause");
    return 0;
    
    
}
