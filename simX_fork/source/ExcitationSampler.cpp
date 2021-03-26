
#include "simX/ExcitationSampler.h"
#include "simX/TwoBodyDecay.h"
#include "simX/weight/SimpleBreitWigner.h"
#include "simX/weight/BreitWignerWithPenetrability.h"

#include <TVector3.h>
#include <vector>
#include <math.h>
#include <simX/CompoundFormation.h>
#include <simX/Logger.h>
#include "ausa/util/memory"

using namespace std;
using namespace simX;


ExcitationSampler::ExcitationSampler( ProcessChain& chain )
 :  chain( chain ),
    DIM(0),
    isInitialized(false),
    verbose(0) {

    static auto log = log::getLogger("ExcitationSampler");


    // Select sampler type
    ROOT::Math::DistSamplerOptions::SetDefaultSampler("Foam");   
   
    // Below, we loop over all the decays and count the number of 
    // daughters with finite widths that need sampling.
    // Pointers to the daugthers are kept in the vector
    // 'fatDaughters' for later use.
    for (auto& p : chain) { // loop over process chain

        auto& d = p.getDaughters();
        for (auto& daughter : d) {  // loop over daughters 
        
            if (daughter->hasFiniteWidth())
                fatDaughters.push_back(daughter);
            else
                slimDaughters.push_back(daughter);
        }
    }
    DIM = static_cast<int>(fatDaughters.size());
    log->debug("Dimension of excitation sampling problem {}", DIM);
        
    // Create sampler
    sampler = ROOT::Math::Factory::CreateDistSampler();
    if (sampler==0) 
        throw invalid_argument("Foam sampler could not be created - exit.");
    
    // Determine min and max values allowed by energy conservation
    determineLimits();    
}


ExcitationSampler::~ExcitationSampler() {
    delete sampler;
}


void ExcitationSampler::sampleExcitationEnergies() {

    for (auto d : slimDaughters)
        d->setExcitationEnergy(d->getNominalExcitationEnergy());

    // sampling only required if DIM>0
 	if (DIM==0) return;

    // Initialize sampler if needed
    if (!isInitialized) init();	

    // Do sampling
    double v[DIM];
    sampler -> Sample(v);	
    
    // Update excitation energies
    for (int i=0; i<DIM; i++)
        fatDaughters[i] -> setExcitationEnergy( v[i] );
        
        
}


double ExcitationSampler::bwprod(double *x, double *p) {

	// Return value
	double fval=1;

	// Variables
    for (int i=0; i<DIM; i++)
        fatDaughters[i] -> setExcitationEnergy( x[i] );
        
    // Evaluate weight functions
    int c=0;
    for (auto&process : chain) {
		double w = process.getWeight();
        fval *= w;
        c++;
 		if (w==0) break; // terminate loop if weight is zero
	}

    return fval;
}


void ExcitationSampler::init() {

	// Sampler only needed if DIM>0
	if (DIM==0) return;
	
    // Set function
    int NPAR = 0;
	TF1 * f = new TF1("weightfunc", this, &ExcitationSampler::bwprod, 0, 1, NPAR);
	sampler -> SetFunction(*f,DIM);
	
    // Set ranges
	int i=0;
	double xmin[DIM], xmax[DIM];		
	for (auto& daughter : fatDaughters) {
        xmin[i] = daughter -> getMinExcitationEnergy();
        xmax[i] = daughter -> getMaxExcitationEnergy();
        if (xmin[i]>=xmax[i])
            throw invalid_argument("Insufficient energy for reaction to proceed.");
        i++;
    }
	sampler -> SetRange(xmin,xmax);

    // Display some info on screen		
    cout << "Initializing excitation-energy sampler ... (this may take a few minutes)" << endl;
	if (verbose) {
	    cout << "\n DIM = " << DIM;
	    cout << "\n xmin = ( ";
		for (int i=0; i<DIM; i++) {
		    if (i>0) cout << ", ";
		    cout << xmin[i];
		}
		cout << ")";
		cout << "\n xmax = ( ";
		for (int i=0; i<DIM; i++) {
   		    if (i>0) cout << ", ";
		    cout << xmax[i];
        }
        cout << " )";
    }
		
    // initialize
    bool ret = sampler -> Init();
    	
    if (!ret)
        throw invalid_argument("Error: excitation-energy sampler could not be initialized");
    if (verbose) cout << "\nInitialization completed" << endl;
    
    isInitialized = true;
}


void ExcitationSampler::determineLimits() {
    
    // Access tree
    auto& t = chain.getTree();

    // Depth (0-1-2- ... -N)
    int N = t.max_depth();
    
    // First process of chain
    auto first = chain.getFirstProcess();
    bool beginsWithCompound = dynamic_cast<CompoundFormation*>(first) != nullptr;

    // Iteration limits (used further down)
    int nbu, ntd;    
    if (beginsWithCompound) {
        nbu = 2;
        ntd = 1;
    }
    else {
        nbu = 0;
        ntd = 0;    
    }

    // Determine min and max excitation of compound nucleus (if relevant)
    if (beginsWithCompound) {
        auto& d = first -> getDaughters();
        double ex0 = d[0] -> getNominalExcitationEnergy();   
        d[0] -> setMinExcitationEnergy( ex0 );   
        d[0] -> setMaxExcitationEnergy( ex0 );
    }
            
    // Leaf iterators
    auto b = t.begin_leaf();
    auto e = t.end_leaf();
    
    // Loop over leafs and check that daughters all have zero width
    // (using iterator syntax)
    for (auto it=b; it!=e; ++it) {
        ProcessChain::DecayPtr& proc = *it;
        auto& d = proc->getDaughters();
        for (auto& daughter : d) { // loop over daugthers
            if (daughter->hasFiniteWidth()) 
                throw invalid_argument("Final-state particles must all have zero width.");
        }
    }
    
    // Bottom-up iteration : set min excitation (of parents)
    for (auto n=N; n>=nbu; --n) {
        // Fixed depth iterators
        auto it = t.begin_fixed(n);
        // Loop from left to right
        while (t.is_valid(it)) {
            ProcessChain::DecayPtr& proc = *it;
            auto& d = proc -> getDaughters();
            auto parent = d[0] -> getParent();
            double mp = parent -> getMass();
            double exmin = -mp;
            for (auto& daughter : d) { // loop over daugthers
                exmin += daughter -> getMass();
                exmin += daughter -> getMinExcitationEnergy();
            }
            if (parent->hasFiniteWidth()) {
                parent -> setMinExcitationEnergy(exmin); // set min excitation of parent
            }
            ++it;
        }
    }

    // Top-down iteration : set max excitation (of daughters)
    for (auto n=ntd; n<=N; n++) {
        // Fixed depth iterators
        auto it = t.begin_fixed(n);
        // Loop from left to right
        while (t.is_valid(it) && !it.isLeaf()) {
            ProcessChain::DecayPtr& proc = *it;
            auto& d = proc -> getDaughters();
            auto parent = d[0] -> getParent();
            for (auto& di : d) {
                double exmax = parent->getMass() + parent->getMaxExcitationEnergy();
                for (auto& dj : d) {
                    exmax -= dj->getMass();
                    if (dj!=di)  exmax -= dj->getMinExcitationEnergy();
                }
                if (di->hasFiniteWidth()) {
                    di -> setMaxExcitationEnergy(exmax); // set max excitation of daughter
                }
            }
            ++it;
        }
    }
    
    // Check that energy-conserving process is possible
    for (auto n=ntd; n<=N; n++) {
        // Fixed depth iterators
        auto it = t.begin_fixed(n);
        // Loop from left to right
        while (t.is_valid(it)) {
            ProcessChain::DecayPtr& proc = *it;
            auto& d = proc -> getDaughters();
            auto parent = d[0] -> getParent();
            double eimax = parent->getMass() + parent->getMaxExcitationEnergy();  
            double efmin = 0;         
            for (auto& di : d) {
                efmin += di->getMass() + di->getMinExcitationEnergy();
            }
            if (efmin>eimax) {
                throw invalid_argument("Process violates energy conservation!");
            }
            ++it;
        }
    }
}




