/////////////////////////////////////////////////////////////////
//    HPVInfection.hh                                          //
//    Created by Mikaela Smit on 25/01/2018.                   //
//    Copyright (c) 2014 Mikaela Smit. All rights reserved.    //
//    These event are added to the EventQ.                     //
/////////////////////////////////////////////////////////////////

#ifndef HPVInfection_hpp
#define HPVInfection_hpp

#include <stdio.h>
#include "person.h"

#endif /* HPVInfection_hpp */

//// --- HPV Events ---
//void EventMyHPVRe_EvaluationUponHIV_infection(person *MyPointerToPerson);   // Function executed to re-evaluate risk of HPV upon HIV infection
void EventMyHPVInfection(person *MyPointerToPerson);        // Function executed once person gets infected with HPV
void EventMyHPVRecovery(person *MyPointerToPerson);         // Function executed once person recovers from HPV
void EventMyCIN1Status(person *MyPointerToPerson);          // Function executed once person progresses from HPV to CIN1
void EventMyCIN2_3Status(person *MyPointerToPerson);        // Function executed once person progresses from CIN1 to CIN2_3
void EventMyCISStatus(person *MyPointerToPerson);           // Function executed once person progresses from CIN2_3 to CIS
void EventMyICCStatus(person *MyPointerToPerson);           // Function executed once person progresses from CIS to ICC
