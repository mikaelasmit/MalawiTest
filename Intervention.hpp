///////////////////////////////////////////////////////////
//  Intervention.hpp                                     //
//  hivmodelzimbabwe                                     //
//                                                       //
//  Created by Mikaela Smit on 12/01/2018.               //
//  Copyright © 2018 Mikaela Smit. All rights reserved.  //
//  File for executing interventions                     //
//                                                       //
///////////////////////////////////////////////////////////


#include <stdio.h>
#include "person.h"

using namespace std;


//// --- This file contains a number of interventions --- ////
void EventStartIntervention(person *MyPointerToPerson);        // This function initialises all interventions we want to roll out

    /// --- HPV functions
    void EventMyHPVVaccination(person *MyPointerToPerson);          // Vaccinate girls for HPV
    void EventVIA_Screening(person *MyPointerToPerson);             // Initial screening for CC with VIA/Cry
    void EventMy_VIA_FollowUp(person *MyPointerToPerson);           // Follow up with VIA/Cry per national guidelines
    void EventMy_CIS_Referral(person *MyPointerToPerson);           // Refer women with CIS to specialised services
    void EventMy_ICC_Referral(person *MyPointerToPerson);           // Refer women with ICC to specialised services

    /// --- CVD functions
    void EventCVDPrevIntervention(person *MyPointerToPerson);       // Intervention for CVD (HT and HC)
