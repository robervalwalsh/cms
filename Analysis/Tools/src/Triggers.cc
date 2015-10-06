/**\class Triggers Triggers.cc Analysis/Tools/src/Triggers.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Rostyslav Shevchenko
//         Created:  Fri, 25 Sep 2015 14:08:26 GMT
//
//

// system include files
#include <iostream>
//
// user include files
#include "FWCore/Framework/interface/Event.h"
//
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Analysis/Tools/interface/Triggers.h"


//
// class declaration
//

using namespace analysis;
using namespace analysis::tools;

//
// constructors and destructor
//
Triggers::Triggers(TChain * mainTree, TChain *firedTree, const std::string unique_name) : Candidates(mainTree)
{

  firedTree -> SetBranchAddress(unique_name,fired_);

}

Triggers::~Triggers()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------

Trigger Triggers::at(const int& i)
{
   using namespace edm;

   Trigger trig;

   if ( i >= n_ || pt_ < 0 )
   {
      return trig;
   }

   trig.set(pt_[i],eta_[i],phi_[i],e_[i],q_[i]);
   trig.fired(fired_[i]);

   return trig;

}
