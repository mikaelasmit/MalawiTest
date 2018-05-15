//////////////////////////////////////////////////////////////
//    eventfunction.h                                       //
//    Created by Mikaela Smit on 07/11/2014.                //
//    Copyright (c) 2014 Mikaela Smit. All rights reserved. //
//    These event are added to the EventQ.                  //
//////////////////////////////////////////////////////////////

#include <stdio.h>
#include "person.h"

using namespace std;

#ifndef HIVInfection_hpp
#define HIVInfection_hpp

#include <stdio.h>

#endif /* HIVInfection_hpp */

//// --- HIV Events ---
void EventMyHIVInfection(person *MyPointerToPerson);        // Function executed once person gets infecte with HIV
void EventCD4change(person *MyPointerToPerson);			   // Function executed when CD4 count is due to change (go down)
void EventARTCatSwitch(person *MyPointerToPerson);         // Function to switch people from Kids ART cat to Adult ART cat

