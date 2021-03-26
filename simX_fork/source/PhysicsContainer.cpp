
#include "simX/PhysicsContainer.h"
#include "simX/git.h"
#include <ausa/util/FileUtil.h>
#include <ausa/AUSA.h>
#include <TTree.h>
#include <TFile.h>
#include <TObjString.h>
#include <TVector3.h>
#include <TMath.h>
#include <ausa/util/memory>

using namespace std;
using namespace simX;

namespace {
	size_t START_CAP = 10;
}

PhysicsContainer::PhysicsContainer(std::shared_ptr<AUSA::TFileWrapper> wrapper)
		: _MU(0), _Z(START_CAP), _A(START_CAP), _E(START_CAP), _TH(START_CAP), _PHI(START_CAP)
{
	TDirectory::TContext cxt{nullptr};
	f = wrapper;
	init(&(wrapper->get()));
}

PhysicsContainer::PhysicsContainer( string path ) : PhysicsContainer(std::make_unique<AUSA::TFileWrapper>(path, "RECREATE"))
{
}

void PhysicsContainer::init(TFile* f) {
	file = f;

	// pointer to tree
	t0 = new TTree( "phys", "simX Physics" );
	t0->SetDirectory(f);

	// b : an 8 bit unsigned integer (UChar_t)  0-255
	// F : a 32 bit floating point (Float_t)

	t0 -> Branch ( "MULT" , &_MU , "_MU/b" );
	bZ = t0 -> Branch ( "Z" , _Z.data() , "_Z[_MU]/b" );
	bA = t0 -> Branch ( "A" , _A.data() , "_A[_MU]/b" );
	bE = t0 -> Branch ( "ENERGY" , _E.data() , "_E[_MU]/F" );
	bTH = t0 -> Branch ( "THETA" , _TH.data() , "_TH[_MU]/F" );
	bPHI = t0 -> Branch ( "PHI" , _PHI.data() , "_PHI[_MU]/F" );

	TObjString ausaSha(AUSA::GIT_HASH);
	TObjString ausaBranch(AUSA::GIT_BRANCH);
	TObjString simXSha(simX::GIT_HASH);
	TObjString simXBranch(simX::GIT_BRANCH);

	f -> WriteTObject(&ausaBranch, "AUSALIB_HASH");
	f -> WriteTObject(&ausaSha, "AUSALIB_BRANCH");
	f -> WriteTObject(&simXSha, "SIMX_HASH");
	f -> WriteTObject(&simXBranch, "SIMX_BRANCH");
}

PhysicsContainer::~PhysicsContainer() 
{
	if (file != nullptr) {
		file->WriteTObject(t0);
		file->Flush();
	}
}


void PhysicsContainer::fill( PhysicsEvent& physEvent ) {
	ensureCapacity(physEvent.size());

	_MU = 0;

    // loop over particles in event
    for (auto p : physEvent) {
        _Z[_MU]  = p -> getZ();
        _A[_MU]  = p -> getA();
        _E[_MU]  = p -> getKineticEnergyLab();
        _TH[_MU] = p -> getDirectionLab().Theta() * 180./TMath::Pi();
        _PHI[_MU] = p -> getDirectionLab().Phi() * 180./TMath::Pi();
        _MU++;
    }    
    // pass data to ROOT tree
	t0 -> Fill();
}

void PhysicsContainer::ensureCapacity(size_t required) {
	auto current = _Z.capacity();
	if (required < current) return;

	auto next = 2*current+1;
	_Z.resize(next);
	bZ->SetAddress(_Z.data());
	_A.resize(next);
	bA->SetAddress(_A.data());
	_E.resize(next);
	bE->SetAddress(_E.data());
	_TH.resize(next);
	bTH->SetAddress(_TH.data());
	_PHI.resize(next);
	bPHI->SetAddress(_PHI.data());
}
