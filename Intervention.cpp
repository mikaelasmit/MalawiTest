

///////////////////////////////////////////////////////////
//  Intervention.cpp                                     //
//  hivmodelzimbabwe                                     //
//                                                       //
//  Created by Mikaela Smit on 12/01/2018.               //
//  Copyright © 2018 Mikaela Smit. All rights reserved.  //
//  File for executing interventions                     //
//                                                       //
///////////////////////////////////////////////////////////


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
#include "NCD.hpp"

using namespace std;

// clean HPV intervention

//// --- Important Internal informtaion --- ////
extern int RandomMinMax_2(int min, int max);

//// --- OUTSIDE INFORMATION --- ////
extern double *p_GT;
extern int int_HPVvaccination;
extern int int_CCscreening;
extern person** MyArrayOfPointersToPeople;
extern priority_queue<event*, vector<event*>, timeComparison> *p_PQ;
extern vector<event*> Events;

extern double** HPVarray;
extern int HPV_Status_HPV;
extern int HPV_Status_CIN1;
extern int HPV_Status_CIN2_3;
extern int HPV_Status_CIS;
extern int HPV_Status_ICC;
extern int HPV_Status_Recovered;
extern int HPV_Status_ReInfected;
extern int CC_Screening_Count;
extern int CC_ScreenOutcome;
extern int CC_CryoOutcome;
extern int age_atrisk_hpv;
extern int max_age_atrisk_hpv;
extern double VIAsens;

extern int count_AdultsART;
extern int HIV_Ref_PersonID[500000];

extern int int_CVDIntervention;
extern int CVD_First_Screening;
extern double CVD_First_Screening_Date;
extern int CVD_Screening_Count;
extern int CVD_Treat_Outcome;
extern double CVD_HT_Treat_Date;
extern double CVD_HC_Treat_Date;

extern double**  NCDArray;
extern double Risk_DiabCVD;
extern double Risk_NCD_HT[3];
extern int relatedNCDs_HT[3];
extern double Risk_NCD_HC[3];
extern int relatedNCDs_HC[3];



/// --- Turning on all interventions to be rolled out --- ///

void EventStartIntervention(person *MyPointerToPerson){
    
    if (int_HPVvaccination==1)          // --- Rolling out HPV vaccination
    {
        cout << "We are rolling out vaccination for HPV and it is the year: " << *p_GT << endl;
    }
    
    if (int_CVDIntervention==1)          // --- Rolling out CVD intervention for just HIV+
    {
        cout << "We are rolling out CVD intervention for HIV+" << *p_GT << endl;
    }
    
    if (int_CVDIntervention==2)          // --- Rolling out CVD intervention for just HIV+
    {
        cout << "We are rolling out CVD intervention for everyone" << *p_GT << endl;
    }
    
    if (int_CCscreening==1)
    {
        cout << "We are rolling out CC screen-and-treat with VIA and Cryo " << *p_GT << endl;
    }
}



////////CVD Intervention //////////////////

void EventCVDPrevIntervention(person *MyPointerToPerson){
    
    if (MyPointerToPerson->Alive==1){
        
        if (MyPointerToPerson->CVD_First_Screening==0){
            MyPointerToPerson->CVD_First_Screening=1;
            MyPointerToPerson->CVD_First_Screening=*p_GT;
            
        }

        MyPointerToPerson->Age=(*p_GT - MyPointerToPerson->DoB);
        MyPointerToPerson->CVD_Screening_Count++;
        double additionalrisk=1;
        
        if (MyPointerToPerson->Diabetes_status==1){
            additionalrisk=Risk_DiabCVD;
            
        }
        
    //    cout << "ID " << MyPointerToPerson->PersonID << "No Screens " << MyPointerToPerson->CVD_Screening_Count << endl;
//        cout << "personid " << MyPointerToPerson->PersonID << "HT " << MyPointerToPerson->HT << "HC " << MyPointerToPerson->HC << "MI " << MyPointerToPerson->MI << "Stroke " << MyPointerToPerson->Stroke << "alive " <<  MyPointerToPerson->Alive << endl;
        
        /////////////////////////////////////////////////////////// --- A person without any CVD risk factors is about to get screened ---
        
        if(MyPointerToPerson->HT_status==0 && MyPointerToPerson->HC_status==0){
            
        }
        
        /////////////////////////////////////////////////////////// --- A person with HT is about to get screened ---
        if(MyPointerToPerson->HT_status==1 && MyPointerToPerson->HC_status==0 && MyPointerToPerson->CVD_Treat_Outcome!=1){
            MyPointerToPerson->CVD_Treat_Outcome=1;
            MyPointerToPerson->CVD_HT_Treat_Date=*p_GT;
            
            double h = ((double)rand()/(RAND_MAX));
            if (h<=0.40){
                int ncd_nr=0;
                double DateNCD=-997;
                
                while (ncd_nr<=2){
                    double r = ((double)rand() / (RAND_MAX));
                    
                    if (r<NCDArray[relatedNCDs_HT[ncd_nr]][120]*additionalrisk){
                        int i=1;
                        
                        while (r>NCDArray[relatedNCDs_HT[ncd_nr]][i]*additionalrisk){i++;}
                        double YearFraction=(RandomMinMax_2(1,12))/12.1;
                        DateNCD=MyPointerToPerson->DoB+i+YearFraction;
                    }
                    
                    if ((DateNCD>=*p_GT && DateNCD>MyPointerToPerson->NCD_DatesVector.at(relatedNCDs_HT[ncd_nr]) && MyPointerToPerson->NCD_DatesVector.at(relatedNCDs_HT[ncd_nr])>0)|| (DateNCD<0)){
                        MyPointerToPerson->NCD_DatesVector.at(relatedNCDs_HT[ncd_nr])=DateNCD;
                        
                        if (ncd_nr==0){
                            MyPointerToPerson->Stroke=DateNCD;
                            
                            if (DateNCD>0){
                                int p=MyPointerToPerson->PersonID-1;
                                event * StrokeEvent = new event;
                                Events.push_back(StrokeEvent);
                                StrokeEvent->time = MyPointerToPerson->Stroke;
                                StrokeEvent->p_fun = &EventMyStrokeDate;
                                StrokeEvent->person_ID = MyArrayOfPointersToPeople[p];
                                p_PQ->push(StrokeEvent);
                            }
                        }
                        
                        if (ncd_nr==1){
                            MyPointerToPerson->CKD=DateNCD;
                            
                            if (DateNCD>0){
                                int p=MyPointerToPerson->PersonID-1;
                                event * CKDEvent = new event;
                                Events.push_back(CKDEvent);
                                CKDEvent->time = MyPointerToPerson->CKD;
                                CKDEvent->p_fun = &EventMyCKDDate;
                                CKDEvent->person_ID = MyArrayOfPointersToPeople[p];
                                p_PQ->push(CKDEvent);
                            }
                        }
                        
                        if (ncd_nr==2){
                            MyPointerToPerson->MI=DateNCD;
                            
                            if (DateNCD>0){
                                int p=MyPointerToPerson->PersonID-1;
                                event * MIEvent = new event;
                                Events.push_back(MIEvent);
                                MIEvent->time = MyPointerToPerson->MI;
                                MIEvent->p_fun = &EventMyMIDate;
                                MIEvent->person_ID = MyArrayOfPointersToPeople[p];
                                p_PQ->push(MIEvent);
                            }
                        }
                    }
                    ncd_nr++;
                }
            }
        }
        
        /////////////////////////////////////////////////////////// --- A person with HC is about to get screened ---
        if(MyPointerToPerson->HT_status==0 && MyPointerToPerson->HC_status==1 && MyPointerToPerson->CVD_Treat_Outcome!=2){
            MyPointerToPerson->CVD_Treat_Outcome=2;
            MyPointerToPerson->CVD_HC_Treat_Date=*p_GT;
            
            double h = ((double)rand()/(RAND_MAX));
            if (h<=0.40){
                int ncd_nr=0;
                double DateNCD=-997;
                
                while (ncd_nr<=2){
                    double r = ((double)rand() / (RAND_MAX));
                    
                    if (r<NCDArray[relatedNCDs_HC[ncd_nr]][120]*additionalrisk){
                        int i=1;
                        
                        while (r>NCDArray[relatedNCDs_HC[ncd_nr]][i]*additionalrisk){i++;}
                        double YearFraction=(RandomMinMax_2(1,12))/12.1;
                        DateNCD=MyPointerToPerson->DoB+i+YearFraction;
                    }
                    
                    if ((DateNCD>=*p_GT && DateNCD>MyPointerToPerson->NCD_DatesVector.at(relatedNCDs_HC[ncd_nr]) && MyPointerToPerson->NCD_DatesVector.at(relatedNCDs_HC[ncd_nr])>0)||(DateNCD<0)){
                        MyPointerToPerson->NCD_DatesVector.at(relatedNCDs_HC[ncd_nr])=DateNCD;
                        
                        if (ncd_nr==0){
                            MyPointerToPerson->HT=DateNCD;
                            
                            if (DateNCD>0){
                                int p=MyPointerToPerson->PersonID-1;
                                event * HTEvent = new event;
                                Events.push_back(HTEvent);
                                HTEvent->time = MyPointerToPerson->HT;
                                HTEvent->p_fun = &EventMyHyptenDate;
                                HTEvent->person_ID = MyArrayOfPointersToPeople[p];
                                p_PQ->push(HTEvent);
                            }
                        }
                        
                        if (ncd_nr==1){
                            MyPointerToPerson->Stroke=DateNCD;
                            
                            if (DateNCD>0){
                                int p=MyPointerToPerson->PersonID-1;
                                event * StrokeEvent = new event;
                                Events.push_back(StrokeEvent);
                                StrokeEvent->time = MyPointerToPerson->Stroke;
                                StrokeEvent->p_fun = &EventMyStrokeDate;
                                StrokeEvent->person_ID = MyArrayOfPointersToPeople[p];
                                p_PQ->push(StrokeEvent);
                            }
                        }
                        
                        if (ncd_nr==2){
                            MyPointerToPerson->MI=DateNCD;
                            
                            if (DateNCD>0){
                                int p=MyPointerToPerson->PersonID-1;
                                event * MIEvent = new event;
                                Events.push_back(MIEvent);
                                MIEvent->time = MyPointerToPerson->MI;
                                MIEvent->p_fun = &EventMyMIDate;
                                MIEvent->person_ID = MyArrayOfPointersToPeople[p];
                                p_PQ->push(MIEvent);
                            }
                        }
                    }
                    ncd_nr++;
                }
            }
        }
        
        /////////////////////////////////////////////////////////// --- A person with HT and HC is about to get screened ---
        if(MyPointerToPerson->HT_status==1 && MyPointerToPerson->HC_status==1 && MyPointerToPerson->CVD_Treat_Outcome!=3){
            if(MyPointerToPerson->CVD_Treat_Outcome==1||MyPointerToPerson->CVD_Treat_Outcome==0){
                MyPointerToPerson->CVD_HC_Treat_Date=*p_GT;
                
            }
            if(MyPointerToPerson->CVD_Treat_Outcome==2||MyPointerToPerson->CVD_Treat_Outcome==0){
                MyPointerToPerson->CVD_HT_Treat_Date=*p_GT;
            }
            
            double h = ((double)rand()/(RAND_MAX));
            if (h<=0.40){
                int ncd_nr=0;
                double DateNCD=-997;
                
                while (ncd_nr<=2){
                    double r = ((double)rand() / (RAND_MAX));
                    
                    if (r<NCDArray[relatedNCDs_HT[ncd_nr]][120]*additionalrisk){
                        int i=1;
                        
                        while (r>NCDArray[relatedNCDs_HT[ncd_nr]][i]*additionalrisk){i++;}
                        double YearFraction=(RandomMinMax_2(1,12))/12.1;
                        DateNCD=MyPointerToPerson->DoB+i+YearFraction;
                    }
                    
                    if ((DateNCD>=*p_GT && DateNCD>MyPointerToPerson->NCD_DatesVector.at(relatedNCDs_HT[ncd_nr]) && MyPointerToPerson->NCD_DatesVector.at(relatedNCDs_HT[ncd_nr])>0)|| (DateNCD<0)){
                        MyPointerToPerson->NCD_DatesVector.at(relatedNCDs_HT[ncd_nr])=DateNCD;
                        
                        if (ncd_nr==0){
                            MyPointerToPerson->Stroke=DateNCD;
                            
                            if (DateNCD>0){
                                int p=MyPointerToPerson->PersonID-1;
                                event * StrokeEvent = new event;
                                Events.push_back(StrokeEvent);
                                StrokeEvent->time = MyPointerToPerson->Stroke;
                                StrokeEvent->p_fun = &EventMyStrokeDate;
                                StrokeEvent->person_ID = MyArrayOfPointersToPeople[p];
                                p_PQ->push(StrokeEvent);
                            }
                        }
                        
                        if (ncd_nr==1){
                            if (MyPointerToPerson->CVD_Treat_Outcome==2|| MyPointerToPerson->CVD_Treat_Outcome==0){
                                MyPointerToPerson->CKD=DateNCD;
                                
                                if (DateNCD>0){
                                    int p=MyPointerToPerson->PersonID-1;
                                    event * CKDEvent = new event;
                                    Events.push_back(CKDEvent);
                                    CKDEvent->time = MyPointerToPerson->CKD;
                                    CKDEvent->p_fun = &EventMyCKDDate;
                                    CKDEvent->person_ID = MyArrayOfPointersToPeople[p];
                                    p_PQ->push(CKDEvent);
                                }
                            }
                        }
                        
                        if (ncd_nr==2){
                            MyPointerToPerson->MI=DateNCD;
                            
                            if (DateNCD>0){
                                int p=MyPointerToPerson->PersonID-1;
                                event * MIEvent = new event;
                                Events.push_back(MIEvent);
                                MIEvent->time = MyPointerToPerson->MI;
                                MIEvent->p_fun = &EventMyMIDate;
                                MIEvent->person_ID = MyArrayOfPointersToPeople[p];
                                p_PQ->push(MIEvent);
                            }
                        }
                    }
                    ncd_nr++;
                }
            }
            MyPointerToPerson->CVD_Treat_Outcome=3;
        }
    }
    
  //          cout << "personid " << MyPointerToPerson->PersonID << "HT " << MyPointerToPerson->HT << "HC " << MyPointerToPerson->HC << "MI " << MyPointerToPerson->MI << "Stroke " << MyPointerToPerson->Stroke << endl;
    
    int p=MyPointerToPerson->PersonID-1;
    event * RecurrentCVD_Screen = new event;
    Events.push_back(RecurrentCVD_Screen);
    RecurrentCVD_Screen->time = *p_GT + 3;
    RecurrentCVD_Screen->p_fun = &EventCVDPrevIntervention;
    RecurrentCVD_Screen->person_ID = MyArrayOfPointersToPeople[p];
    p_PQ->push(RecurrentCVD_Screen);
    
 //   cout << "personid " << MyPointerToPerson->PersonID << "next screen " << *p_GT + 3 << endl;
}





/// --- HPV vaccination event --- ///

void EventMyHPVVaccination(person *MyPointerToPerson){
    
    if (MyPointerToPerson->Alive == 1){
        
        MyPointerToPerson->HPVvaccination_status=1;
        MyPointerToPerson->HPVvaccination_date=*p_GT;
        
    }
}



/// --- CC screen-and-treat events --- ///

void EventVIA_Screening(person *MyPointerToPerson){
    
    if (int_CCscreening==1){
        
        // --- Set elegibility criteria for first screening with VIA
        if (MyPointerToPerson->Alive==1 && MyPointerToPerson->Sex==2 && MyPointerToPerson->VIAcount==0 && MyPointerToPerson->ART<=*p_GT && MyPointerToPerson->ART != -999){
            
            // --- Get a date of screening and updated counter
            MyPointerToPerson->VIA=*p_GT;
            MyPointerToPerson->VIAcount++;
            
            
            // --- A woman with no HPV infection or that has recovered from a previous infection is about to be screened
            if(MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount]<=0 || MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1]==HPV_Status_Recovered)
            {
                // Test will be normal, let's book her for a second screening in 6 months
                int p=MyPointerToPerson->PersonID-1;
                event * DateOfSecondScreening = new event;
                Events.push_back(DateOfSecondScreening);
                DateOfSecondScreening->time = *p_GT+0.5;
                DateOfSecondScreening->p_fun = &EventMy_VIA_FollowUp;
                DateOfSecondScreening->person_ID = MyArrayOfPointersToPeople[p];
                p_PQ->push(DateOfSecondScreening);
            }
            
            // --- A woman with HPV infection, but no cervical abnormalities is about to be screened
            if (MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1]==HPV_Status_HPV)
            {
                // Test will be normal, let's book her for a second screening in 6 months
                int p=MyPointerToPerson->PersonID-1;
                event * DateOfSecondScreening = new event;
                Events.push_back(DateOfSecondScreening);
                DateOfSecondScreening->time = *p_GT+0.5;
                DateOfSecondScreening->p_fun = &EventMy_VIA_FollowUp;
                DateOfSecondScreening->person_ID = MyArrayOfPointersToPeople[p];
                p_PQ->push(DateOfSecondScreening);
            }
            
            // --- A woman with CIN1 is about to be screened
            if (MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1]==2)
            {
                double h = ((double) rand() / (RAND_MAX));  // Let's see if she's diagnosed
                if(h>VIAsens){    // No luck =(
                    MyPointerToPerson->CC_ScreenOutcome=2;   // 2 = False negative
                    // Let's book her for a second screening in 6 months and hope she's lucky next time
                    int p=MyPointerToPerson->PersonID-1;
                    event * DateOfSecondScreening = new event;
                    Events.push_back(DateOfSecondScreening);
                    DateOfSecondScreening->time = *p_GT+0.5;
                    DateOfSecondScreening->p_fun = &EventMy_VIA_FollowUp;
                    DateOfSecondScreening->person_ID = MyArrayOfPointersToPeople[p];
                    p_PQ->push(DateOfSecondScreening);
                }
                
                if(h<=VIAsens) {  // She was lucky!
                    MyPointerToPerson->CC_ScreenOutcome=1;   // 1 = True positive
                    // Let's give her cryotherapy right now and see if she's cured
                    double j = ((double) rand() / (RAND_MAX));
                    
                    if(j>0.875){    // Not cured
                        MyPointerToPerson->CC_CryoOutcome=2;   // 2 = Cryo not successful
                        // Let's book her for follow-up in one year and hope she's lucky
                        int p=MyPointerToPerson->PersonID-1;
                        event * DateOfSecondScreening = new event;
                        Events.push_back(DateOfSecondScreening);
                        DateOfSecondScreening->time = *p_GT+1;
                        DateOfSecondScreening->p_fun = &EventMy_VIA_FollowUp;
                        DateOfSecondScreening->person_ID = MyArrayOfPointersToPeople[p];
                        p_PQ->push(DateOfSecondScreening);
                    }
                    
                    if(j<=0.875){   // She'll be cured with cryo!
                        MyPointerToPerson->CC_CryoOutcome=1;   // 1 = Cryo successful
                        MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1]=HPV_Status_Recovered;       // Change her status to recovered
                        MyPointerToPerson->HPV_DateofRecovery[MyPointerToPerson->HPVcount-1]=*p_GT;
                        // Let's book her for follow-up in one year, according to guidelines
                        int p=MyPointerToPerson->PersonID-1;
                        event * DateOfSecondScreening = new event;
                        Events.push_back(DateOfSecondScreening);
                        DateOfSecondScreening->time = *p_GT+1;
                        DateOfSecondScreening->p_fun = &EventMy_VIA_FollowUp;
                        DateOfSecondScreening->person_ID = MyArrayOfPointersToPeople[p];
                        p_PQ->push(DateOfSecondScreening);

                    }
                }
            }
            
            // --- A woman with CIN2/3 is about to be screened
            if (MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1]==3)
            {
                double h = ((double) rand() / (RAND_MAX));  // Let's see if she's diagnosed
                if(h>VIAsens){    // No luck =(
                    MyPointerToPerson->CC_ScreenOutcome=2;   // 2 = False negative
                    // Let's book her for a second screening in 6 months and hope she's lucky next time
                    int p=MyPointerToPerson->PersonID-1;
                    event * DateOfSecondScreening = new event;
                    Events.push_back(DateOfSecondScreening);
                    DateOfSecondScreening->time = *p_GT+0.5;
                    DateOfSecondScreening->p_fun = &EventMy_VIA_FollowUp;
                    DateOfSecondScreening->person_ID = MyArrayOfPointersToPeople[p];
                    p_PQ->push(DateOfSecondScreening);
                }
                
                if(h<=VIAsens) {  // She was lucky!
                    MyPointerToPerson->CC_ScreenOutcome=1;   // 1 = True positive
                    // Let's give her cryotherapy right now and see if she's cured
                    double j = ((double) rand() / (RAND_MAX));
                    
                    if(j>0.875){    // Not cured
                        MyPointerToPerson->CC_CryoOutcome=2;   // 2 = Cryo not successful
                        // Let's book her for follow-up in one year and hope she's lucky
                        int p=MyPointerToPerson->PersonID-1;
                        event * DateOfSecondScreening = new event;
                        Events.push_back(DateOfSecondScreening);
                        DateOfSecondScreening->time = *p_GT+1;
                        DateOfSecondScreening->p_fun = &EventMy_VIA_FollowUp;
                        DateOfSecondScreening->person_ID = MyArrayOfPointersToPeople[p];
                        p_PQ->push(DateOfSecondScreening);
                    }
                    
                    if(j<=0.875){   // She'll be cured with cryo!
                        MyPointerToPerson->CC_CryoOutcome=1;   // 1 = Cryo successful
                        MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1]=HPV_Status_Recovered;       // Change her status to recovered
                        MyPointerToPerson->HPV_DateofRecovery[MyPointerToPerson->HPVcount-1]=*p_GT;
                        // Let's book her for follow-up in one year, according to guidelines
                        int p=MyPointerToPerson->PersonID-1;
                        event * DateOfSecondScreening = new event;
                        Events.push_back(DateOfSecondScreening);
                        DateOfSecondScreening->time = *p_GT+1;
                        DateOfSecondScreening->p_fun = &EventMy_VIA_FollowUp;
                        DateOfSecondScreening->person_ID = MyArrayOfPointersToPeople[p];
                        p_PQ->push(DateOfSecondScreening);
                    }
                }
                
                cout << MyPointerToPerson->PersonID << " is HPV status " << MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1] << " with a counter of " << MyPointerToPerson->HPVcount << endl;
                
            }
            // --- A woman with CIS is about to be screened
            if (MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1]==4)
            {
                double h = ((double) rand() / (RAND_MAX));  // Let's see if she's diagnosed
                if(h>VIAsens){    // No luck =(
                    MyPointerToPerson->CC_ScreenOutcome=2;   // 2 = False negative
                    // Let's book her for a second screening in 6 months and hope she's lucky next time
                    int p=MyPointerToPerson->PersonID-1;
                    event * DateOfSecondScreening = new event;
                    Events.push_back(DateOfSecondScreening);
                    DateOfSecondScreening->time = *p_GT+0.5;
                    DateOfSecondScreening->p_fun = &EventMy_VIA_FollowUp;
                    DateOfSecondScreening->person_ID = MyArrayOfPointersToPeople[p];
                    p_PQ->push(DateOfSecondScreening);
                }
                
                if(h<=VIAsens) {  // She was diagnosed; needs referral!
                    MyPointerToPerson->CC_ScreenOutcome=1;   // 1 = True positive
                    // Let's see what her chances are of having access to an appropriate referral centre
                    double j = ((double) rand() / (RAND_MAX));
                    
                    if(j>0.175){    // Check this assumption... i.e. only 17.5% of women in Kenya will access colpo+biop+leep services?
                        MyPointerToPerson->CC_CryoOutcome=4;   // 4 = fell through the loop-holes of the system... she will progress with her disease
                    }
                    
                    if(j<=0.175){   // She'll access specialised care!
                        MyPointerToPerson->CC_CryoOutcome=3;   // 3 = referred to specialised care
                        // Let's book her for follow-up in one year, according to guidelines
                        int p=MyPointerToPerson->PersonID-1;
                        event * Referral = new event;
                        Events.push_back(Referral);
                        Referral->time = *p_GT;
                        Referral->p_fun = &EventMy_CIS_Referral;
                        Referral->person_ID = MyArrayOfPointersToPeople[p];
                        p_PQ->push(Referral);
                    }
                }
            }
            
            // --- A woman with ICC is about to be screened
            if (MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1]==5)
            {
                double h = ((double) rand() / (RAND_MAX));  // Let's see if she's diagnosed
                if(h>VIAsens){    // No luck =(
                    MyPointerToPerson->CC_ScreenOutcome=2;   // 2 = False negative
                    // Let's book her for a second screening in 6 months and hope she's lucky next time
                    int p=MyPointerToPerson->PersonID-1;
                    event * DateOfSecondScreening = new event;
                    Events.push_back(DateOfSecondScreening);
                    DateOfSecondScreening->time = *p_GT+0.5;
                    DateOfSecondScreening->p_fun = &EventMy_VIA_FollowUp;
                    DateOfSecondScreening->person_ID = MyArrayOfPointersToPeople[p];
                    p_PQ->push(DateOfSecondScreening);
                }
                
                if(h<=VIAsens) {  // She was diagnosed; needs referral!
                    MyPointerToPerson->CC_ScreenOutcome=1;   // 1 = True positive
                    // Let's see what her chances are of having access to an appropriate referral centre
                    double j = ((double) rand() / (RAND_MAX));
                    if(j>0.175){    // Check this assumption... i.e. only 17.5% of women in Kenya will access specialised services?
                        MyPointerToPerson->CC_CryoOutcome=4;   // 4 = fell through the loop-holes of the system... she will progress with her disease
                    }
                    
                    if(j<=0.175){   // She'll access specialised care!
                        MyPointerToPerson->CC_CryoOutcome=3;   // 3 = referred to specialised care
                        // Let's book her for follow-up in one year, according to guidelines
                        int p=MyPointerToPerson->PersonID-1;
                        event * Referral = new event;
                        Events.push_back(Referral);
                        Referral->time = *p_GT;
                        Referral->p_fun = &EventMy_ICC_Referral;
                        Referral->person_ID = MyArrayOfPointersToPeople[p];
                        p_PQ->push(Referral);
                    }
                }
            }
        }
    }
}



void EventMy_VIA_FollowUp(person *MyPointerToPerson){
    
    double age_now = *p_GT-MyPointerToPerson->DoB;
    
    if (MyPointerToPerson->Alive==1 && age_now<=49 && MyPointerToPerson->VIAcount>0) {          // Continue on executing this function only while the person is alive and are still <49 years old
        MyPointerToPerson->VIA=*p_GT;
        MyPointerToPerson->VIAcount++;
        
        if (MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount]<=0 || MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1]==6){
            // Let's book her for a subsequent screening in one year
            int p=MyPointerToPerson->PersonID-1;
            event * DateOfReScreening = new event;
            Events.push_back(DateOfReScreening);
            DateOfReScreening->time = *p_GT+1;
            DateOfReScreening->p_fun = &EventMy_VIA_FollowUp;
            DateOfReScreening->person_ID = MyArrayOfPointersToPeople[p];
            p_PQ->push(DateOfReScreening);
        }
        
        if (MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1]==1){   // HPV infection, but no abnormalities on VIA
            MyPointerToPerson->VIAcount++;
            // Let's book her for a subsequent screening in one year
            int p=MyPointerToPerson->PersonID-1;
            event * DateOfReScreening = new event;
            Events.push_back(DateOfReScreening);
            DateOfReScreening->time = *p_GT+1;
            DateOfReScreening->p_fun = &EventMy_VIA_FollowUp;
            DateOfReScreening->person_ID = MyArrayOfPointersToPeople[p];
            p_PQ->push(DateOfReScreening);
        }
        
        if (MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1]==2){   // CIN 1
            MyPointerToPerson->VIAcount++;
            // Let's see what hapens
            double h = ((double) rand() / (RAND_MAX));
            if(h>VIAsens){    // No luck =(
                MyPointerToPerson->CC_ScreenOutcome=2;   // 2 = False negative
                // Let's book her for re-screening in 1 year and hope she's lucky next time
                int p=MyPointerToPerson->PersonID-1;
                event * DateOfSecondScreening = new event;
                Events.push_back(DateOfSecondScreening);
                DateOfSecondScreening->time = *p_GT+1;
                DateOfSecondScreening->p_fun = &EventMy_VIA_FollowUp;
                DateOfSecondScreening->person_ID = MyArrayOfPointersToPeople[p];
                p_PQ->push(DateOfSecondScreening);
            }
            if(h<=VIAsens) {  // She was lucky!
                MyPointerToPerson->CC_ScreenOutcome=1;   // 1 = True positive
                // Let's give her cryotherapy right now and see if she's cured
                double j = ((double) rand() / (RAND_MAX));
                
                if(j>0.875){    // Not cured
                    MyPointerToPerson->CC_CryoOutcome=2;   // 2 = Cryo not successful
                    // Let's book her for follow-up in one year and hope she's lucky
                    int p=MyPointerToPerson->PersonID-1;
                    event * DateOfSecondScreening = new event;
                    Events.push_back(DateOfSecondScreening);
                    DateOfSecondScreening->time = *p_GT+1;
                    DateOfSecondScreening->p_fun = &EventMy_VIA_FollowUp;
                    DateOfSecondScreening->person_ID = MyArrayOfPointersToPeople[p];
                    p_PQ->push(DateOfSecondScreening);
                }
                if(j<=0.875){   // She'll be cured with cryo!
                    MyPointerToPerson->CC_CryoOutcome=1;   // 1 = Cryo successful
                    MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1]=HPV_Status_Recovered;
                    MyPointerToPerson->HPV_DateofRecovery[MyPointerToPerson->HPVcount-1]=*p_GT;
                    // Let's book her for follow-up in one year, according to guidelines
                    int p=MyPointerToPerson->PersonID-1;
                    event * DateOfSecondScreening = new event;
                    Events.push_back(DateOfSecondScreening);
                    DateOfSecondScreening->time = *p_GT+1;
                    DateOfSecondScreening->p_fun = &EventMy_VIA_FollowUp;
                    DateOfSecondScreening->person_ID = MyArrayOfPointersToPeople[p];
                    p_PQ->push(DateOfSecondScreening);
                }
            }
        }
        
        if(MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1]==3){   // CIN 2/3
            MyPointerToPerson->VIAcount++;
            // Let's see what hapens
            double h = ((double) rand() / (RAND_MAX));
            if(h>VIAsens){    // No luck =(
                MyPointerToPerson->CC_ScreenOutcome=2;   // 2 = False negative
                // Let's book her for re-screening in 1 year and hope she's lucky next time
                int p=MyPointerToPerson->PersonID-1;
                event * DateOfSecondScreening = new event;
                Events.push_back(DateOfSecondScreening);
                DateOfSecondScreening->time = *p_GT+1;
                DateOfSecondScreening->p_fun = &EventMy_VIA_FollowUp;
                DateOfSecondScreening->person_ID = MyArrayOfPointersToPeople[p];
                p_PQ->push(DateOfSecondScreening);
            }
            if(h<=VIAsens) {  // She was lucky!
                MyPointerToPerson->CC_ScreenOutcome=1;   // 1 = True positive
                // Let's give her cryotherapy right now and see if she's cured
                double j = ((double) rand() / (RAND_MAX));
                
                if(j>0.875){    // Not cured
                    MyPointerToPerson->CC_CryoOutcome=2;   // 2 = Cryo not successful
                    // Let's book her for follow-up in one year and hope she's lucky
                    int p=MyPointerToPerson->PersonID-1;
                    event * DateOfSecondScreening = new event;
                    Events.push_back(DateOfSecondScreening);
                    DateOfSecondScreening->time = *p_GT+1;
                    DateOfSecondScreening->p_fun = &EventMy_VIA_FollowUp;
                    DateOfSecondScreening->person_ID = MyArrayOfPointersToPeople[p];
                    p_PQ->push(DateOfSecondScreening);
                }
                if(j<=0.875){   // She'll be cured with cryo!
                    MyPointerToPerson->CC_CryoOutcome=1;   // 1 = Cryo successful
                    MyPointerToPerson->HPV_Status[MyPointerToPerson->HPVcount-1]=HPV_Status_Recovered;
                    MyPointerToPerson->HPV_DateofRecovery[MyPointerToPerson->HPVcount-1]=*p_GT;
                    // Let's book her for follow-up in one year, according to guidelines
                    int p=MyPointerToPerson->PersonID-1;
                    event * DateOfSecondScreening = new event;
                    Events.push_back(DateOfSecondScreening);
                    DateOfSecondScreening->time = *p_GT+1;
                    DateOfSecondScreening->p_fun = &EventMy_VIA_FollowUp;
                    DateOfSecondScreening->person_ID = MyArrayOfPointersToPeople[p];
                    p_PQ->push(DateOfSecondScreening);
                }
            }
        }
    }
}

void EventMy_CIS_Referral(person *MyPointerToPerson){
}

void EventMy_ICC_Referral(person *MyPointerToPerson){
}



