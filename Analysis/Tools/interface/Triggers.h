#ifndef Analysis_Tools_Triggers_h
#define Analysis_Tools_Triggers_h 1

// -*- C++ -*-
//
// Package:    Analysis/Tools
// Class:      Triggers
//
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
#include <memory>
#include <vector>
// user include files

#include "TTree.h"
#include "TChain.h"
#include "Analysis/Tools/interface/Candidates.h"
#include "Analysis/Tools/interface/Trigger.h"

//
// class declaration
//

namespace analysis {
   namespace tools {

      class Triggers : public Candidates {
         public:
            Triggers(TChain *, TChain *, const std::string);
           ~Triggers();

            Trigger at(const int &);

         private:
            // ----------member data ---------------------------
            int fired_[max_];
      };
   }
}

#endif  // Analysis_Tools_Triggers_h
