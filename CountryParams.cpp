//////////////////////////////////////////////////////////
//                                                      //
//                  CountryParams.cpp                   //
//            Created by Pablo on 20/10/2017.           //
//  Copyright © 2017 Mikaela Smit. All rights reserved. //
//                                                      //
//////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "errorcoutmacro.h"
#include "CountryParams.hpp"
#include "LoadParams.h"

using namespace std;

// clean up all hpv related items

// General parameter
int UN_Pop;
int init_pop;
int total_population;
double Sex_ratio;
extern int factor;
int ageAdult;
int minAgeBirth;
int maxAgeBirth;


// HIV parameters
int ART_start_yr;
double ARTbuffer;

// HPV parameters
int age_atrisk_hpv;
int max_age_atrisk_hpv;
extern int age_tostart_CCscreening;
double HPV_Array_buffer;
double HPV_hiv_Array_buffer;
int HPV_Status_HPV;
int HPV_Status_CIN1;
int HPV_Status_CIN2_3;
int HPV_Status_CIS;
int HPV_Status_ICC;
int HPV_Status_Recovered;
int HPV_Status_ReInfected;
double VIAsens;
double hpv_date_after_death;
double no_hpv_infection;

// Directory parameters
string ParamDirectory;
extern string InputFileDirectory;

// Mortality parameters
double background_d;
double HIV_d;
double IHD_d;
double Depression_d;
double Asthma_d;
double Stroke_d;
double Diabetes_d;
double CKD_d;
double Colo_d;
double Liver_d;
double Oeso_d;
double Prostate_d;
double OtherCan_d;
extern double MortRisk[6];
extern double MortRisk_Cancer[5];
double MortAdj;


/////////////////// FUNCTION: IF LOOP FOR COUNTRY-SPECIFIC PARAMETERS //////////////////////
void loadCountryParams(int x){
    if (x == 1){                            // KENYA
        
        cout << "You have selected Kenya" << endl;
        
        // Central parameters
        UN_Pop=5909800;
        init_pop=UN_Pop/factor;
        total_population=init_pop;
        Sex_ratio=0.495639296;
        ageAdult=15;
        minAgeBirth=15;
        maxAgeBirth=50;
        
        // HIV parameters
        ART_start_yr=2004;
        ARTbuffer=1.01;
        
        // HPV parameters
        age_atrisk_hpv=15;
        max_age_atrisk_hpv=65;
        age_tostart_CCscreening=25;
        HPV_Array_buffer=1.3;
        HPV_hiv_Array_buffer=1;
        VIAsens=0.76;               // Vary sensitivity from 48 to 76% (Viviano et al. 2017)
        no_hpv_infection=-988;
        hpv_date_after_death = -977;
        HPV_Status_HPV=1;
        HPV_Status_CIN1=2;
        HPV_Status_CIN2_3=3;
        HPV_Status_CIS=4;
        HPV_Status_ICC=5;
        HPV_Status_Recovered=6;
        HPV_Status_ReInfected=7;
        
        // Mortality parameters
        MortAdj=1;
        background_d =71.32;                // Mortality percentages from GBD 2016
        HIV_d        =15.56;
        IHD_d         =3.99;
        Depression_d =0.0;
        Asthma_d     =0.45;
        Stroke_d     =3.92;
        Diabetes_d   =1.27;
        CKD_d        =1.5;
        Colo_d       =0.39;
        Liver_d      =0.34;
        Oeso_d       =0.32;
        Prostate_d    =0.31;
        OtherCan_d   =0.71;
    }
    
    else if (x == 2){                                          // ZIMBABWE
        
        cout << "You have selected Zimbabwe" << endl;
        
        // Central parameters
        UN_Pop=2565000;
        init_pop=UN_Pop/factor;
        total_population=init_pop;
        Sex_ratio=0.4986;
        ageAdult=15;
        minAgeBirth=15;
        maxAgeBirth=50;
        
        // HIV parameters
        ART_start_yr=2004;
        ARTbuffer=1.01;
        
        // HPV parameters
        age_atrisk_hpv=15;
        age_tostart_CCscreening=25;
        HPV_Array_buffer=1.3;
        HPV_hiv_Array_buffer=1;
        VIAsens=0.76;
        no_hpv_infection=-988;
        hpv_date_after_death = -977;
        HPV_Status_HPV=1;
        HPV_Status_CIN1=2;
        HPV_Status_CIN2_3=3;
        HPV_Status_CIS=4;
        HPV_Status_ICC=5;
        HPV_Status_Recovered=6;
        
        // Mortality percentages from GBD 2013
        MortAdj=1;
        background_d =56.6;
        HIV_d        =29.6;
        IHD_d         =1.00;
        Depression_d =0.0;
        Asthma_d     =0.7;
        Stroke_d     =4.6;
        Diabetes_d   =1.7;
        CKD_d        =1.3;
        Colo_d       =0.3;
        Liver_d      =0.3;
        Oeso_d       =0.6;
        Prostate_d    =0.2;
        OtherCan_d   =3.1;
    }
    
    
    else if (x == 3){                                          // MALAWI need to replace params
        
        cout << "You have selected MALAWI" << endl;
        
        
        // Central Parameters
        UN_Pop=2565000;
        init_pop=UN_Pop/factor;
        total_population=init_pop;
        Sex_ratio=0.4986;
        ageAdult=15;
        minAgeBirth=15;
        maxAgeBirth=50;
        
        // HIV parameters
        ART_start_yr=2004;
        ARTbuffer=1.01;
        
        // HPV parameters
        age_atrisk_hpv=15;
        age_tostart_CCscreening=25;
        HPV_Array_buffer=1.3;
        HPV_hiv_Array_buffer=1;
        VIAsens=0.76;
        no_hpv_infection=-988;
        hpv_date_after_death = -977;
        HPV_Status_HPV=1;
        HPV_Status_CIN1=2;
        HPV_Status_CIN2_3=3;
        HPV_Status_CIS=4;
        HPV_Status_ICC=5;
        HPV_Status_Recovered=6;
        
        // Mortality percentages from GBD 2013
        MortAdj=1;
        background_d =56.6;
        HIV_d        =29.6;
        IHD_d         =1.00;
        Depression_d =0.0;
        Asthma_d     =0.7;
        Stroke_d     =4.6;
        Diabetes_d   =1.7;
        CKD_d        =1.3;
        Colo_d       =0.3;
        Liver_d      =0.3;
        Oeso_d       =0.6;
        Prostate_d    =0.2;
        OtherCan_d   =3.1;
    }
    
    else if (x == 4){                                          // Uasin Gishu - KENYA
        
        cout << "You have selected Kenya - Uasin Gishu" << endl;
        
        // Central Parameters
        UN_Pop=99075;
        init_pop=UN_Pop/factor;
        total_population=init_pop;
        Sex_ratio=0.50237783022;
        ageAdult=15;
        minAgeBirth=15;
        maxAgeBirth=50;
        
        // HIV parameters
        ART_start_yr=2004;
        ARTbuffer=1.01;
        
        // HPV parameters
        age_atrisk_hpv=15;
        age_tostart_CCscreening=25;
        HPV_Array_buffer=1.3;
        HPV_hiv_Array_buffer=1;
        VIAsens=0.76;
        no_hpv_infection=-988;
        hpv_date_after_death = -977;
        HPV_Status_HPV=1;
        HPV_Status_CIN1=2;
        HPV_Status_CIN2_3=3;
        HPV_Status_CIS=4;
        HPV_Status_ICC=5;
        HPV_Status_Recovered=6;
        
        // Mortality parameters
        MortAdj=1;
        background_d =71.32;
        HIV_d        =15.56;
        IHD_d         =3.99;
        Depression_d =0.0;
        Asthma_d     =0.45;
        Stroke_d     =3.92;
        Diabetes_d   =1.27;
        CKD_d        =1.5;
        Colo_d       =0.39;
        Liver_d      =0.34;
        Oeso_d       =0.32;
        Prostate_d    =0.31;
        OtherCan_d   =0.71;
    }
}


void getParamsString(int x){
    if (x == 1){
        ParamDirectory=InputFileDirectory + "/Kenya/";
    }
    
    else if (x == 2){
        ParamDirectory=InputFileDirectory + "/Zimbabwe/";
    }
    else if (x == 3){
        ParamDirectory=InputFileDirectory + "/Malawi/";
    }
    
    else if (x == 4){
        ParamDirectory=InputFileDirectory + "/Malawi/";
    }
    
}

