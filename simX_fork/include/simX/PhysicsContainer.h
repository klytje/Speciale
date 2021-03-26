#ifndef PHYSICSCONTAINER_H
#define PHYSICSCONTAINER_H

#include "simX/PhysicsEvent.h"

// C++ and ROOT header files
#include <string>
#include <ausa/TFileWrapper.h>

class TTree;
class TBranch;
class TFile;

namespace simX {

    class PhysicsContainer {

        public:
            PhysicsContainer( std::shared_ptr<AUSA::TFileWrapper> wrapper );
            PhysicsContainer( std::string path );
            virtual ~PhysicsContainer();
            
            void fill( PhysicsEvent& physEvent );

        private:
        	TTree * t0; // TFile will delete this.
            std::shared_ptr<AUSA::TFileWrapper> f;
            TFile* file;


	        UChar_t _MU;

            std::vector<UChar_t> _Z, _A;
            std::vector<Float_t> _E, _TH, _PHI;

            TBranch *bZ, *bA, *bE, *bTH, *bPHI;

            void ensureCapacity(size_t cap);

            void init(TFile* f);
    };
}

#endif	/* PHYSICSCONTAINER_H */
