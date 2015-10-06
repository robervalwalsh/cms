#ifndef Analysis_Tools_Trigger_h
#define Analysis_Tools_Trigger_h 1

// -*- C++ -*-
//
// Package:    Analysis/Tools
// Class:      Trigger
//
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
#include <memory>
//
// user include files
#include "Analysis/Tools/interface/Candidate.h"
//
// class declaration
//

namespace analysis {
   namespace tools {

      class Trigger : public Candidate {
         public:
            Trigger();
           ~Trigger();
            using Candidate::set; // in case needed to overload the function set

         private:
            // ----------member data ---------------------------

            //
      };
   }
}

#endif  // Analysis_Tools_Trigger_h
