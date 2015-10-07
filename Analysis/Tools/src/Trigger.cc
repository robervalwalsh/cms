/**\class Trigger Trigger.cc Analysis/Tools/src/Trigger.cc

 Description: [Base Class for Trigger objects]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rostyslav Shevchenko
//         Created:  Fri, 25 Sep 2015 13:51:04 GMT
//
//

// system include files
//
// user include files

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Analysis/Tools/interface/Trigger.h"

//
// class decalaration
//

using namespace analysis;
using namespace analysis::tools;

Trigger::Trigger() : Candidate()
{

}

Trigger::~Trigger()
{

}

// ................methods ................
int Trigger::fired() {return fired_;}

void Trigger::fired(const int & fired) { fired_ = fired}
