/////////////////////////////////////////////////////////////////
//    NCD.hpp                                                 //
//    Created by Mikaela Smit on 25/01/2018.				   //
//    Copyright (c) 2014 Mikaela Smit. All rights reserved.    //
//    These event relate to NCDs.                              //
/////////////////////////////////////////////////////////////////

#ifndef NCD_hpp
#define NCD_hpp

#endif /* NCD_hpp */

#include <stdio.h>
#include "person.h"

using namespace std;


//// --- NCD Events ---
void EventMyAsthmaDate(person *MyPointerToPerson);			// Function executed when person develops asthma
void EventMyCKDDate(person *MyPointerToPerson);             // Function executed when person develops CKD
void EventMyDepressionDate(person *MyPointerToPerson);		// Function executed when person develops depression
void EventMyDiabetesDate(person *MyPointerToPerson);        // Function executed when person develops diabetes
void EventMyHypcholDate(person *MyPointerToPerson);			// Function executed when person develops high-cholesterol
void EventMyHyptenDate(person *MyPointerToPerson);			// Function executed when person develops hypertension
void EventMyMIDate(person *MyPointerToPerson);              // Function executed when person develops MI
void EventMyStrokeDate(person *MyPointerToPerson);			// Function executed when person develops stroke


/// --- Cancer events ---
void EventMyBreastDate(person *MyPointerToPerson);          // Function executed when person develops breast cancer
void EventMyCervicalDate(person *MyPointerToPerson);        // Function executed when person develops Cervical cancer
void EventMyColoDate(person *MyPointerToPerson);            // Function executed when person develops Colorectum cancer
void EventMyLiverDate(person *MyPointerToPerson);			// Function executed when person develops Liver cancer
void EventMyOesoDate(person *MyPointerToPerson);            // Function executed when person develops Oesophageal cancer
void EventMyProstateDate(person *MyPointerToPerson);		// Function executed when person develops prostate cancer
void EventMyOtherCanDate(person *MyPointerToPerson);        // Function executed when person develops other cancers


