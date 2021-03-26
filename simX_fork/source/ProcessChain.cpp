#include "simX/ProcessChain.h"
#include <iostream>
#include <unordered_set>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/random_access_index.hpp>

using namespace std;
using namespace simX;


ProcessChain::ProcessChain( std::unique_ptr<Particle> beam, Tree& chain )
 :  beam(move(beam)), t(move(chain))
{
    // We move everything to a vector, such that nodes lies sequentially in memory
    for (auto& item : t) {
        ProcessChain::chain.push_back(item.get());
    }

    validate();
}

ProcessChain::~ProcessChain() {/*Nothing*/}

ProcessChain::iterator_::iterator_() {/*Nothing*/}

ProcessChain::iterator_::iterator_(sub_iterator iter) : iter(iter)
{/*Nothing*/}

NuclearProcess &ProcessChain::iterator_::operator*() const {
    return **iter;
}

bool ProcessChain::iterator_::operator==(const ProcessChain::iterator_ &it) const {
    return iter == it.iter;
}

bool ProcessChain::iterator_::operator!=(const ProcessChain::iterator_& it) const {
    return !((*this) == it);
}

ProcessChain::iterator_ & ProcessChain::iterator_::operator++() {
    iter++;
    return *this;
}

NuclearProcess* const ProcessChain::findDecayOf(const std::string& parent) {
    for (auto& c : chain) {
        auto& d = c->getDaughters();

        for (auto& daughter : d) {
            auto p = daughter->getParent();
            if (p != nullptr && p->getName() == parent) return c;
        }
    }
    return nullptr;
}

NuclearProcess* ProcessChain::getFirstProcess() {
    return !chain.empty() ? chain[0] : nullptr;
}

std::vector<Particle*> ProcessChain::findFinalStates(){
    using namespace ::boost;
    using namespace ::boost::multi_index;

    // Here we define a container which allows for fast lookup + insertion order
    // ie. a merge of Set and vector.
    using Set = multi_index_container<
            Particle*,
            indexed_by<
                    random_access<>,  // keep insertion order
                    hashed_unique<identity<Particle*>> // Only accept 1 entry of each
    >>;

    Set parent, all;
    std::vector<Particle*> result;

    for (auto& p : chain) {
        for (auto& d : p->getDaughters()) {
            all.push_back(d);
            parent.push_back(d->getParent());
        }
    }

    for (auto p : all.get<0>()) {
        if (!std::count(std::begin(parent), std::end(parent), p)) result.push_back(p);
    }
    return result;
}

void ProcessChain::validate() const {
    for (auto& c : chain) {
        
        // initial A, Z for process 
        int Ai = c -> getNumberOfNucleons();
        int Zi = c -> getNumberOfProtonsInInitialState();
        
        // final A, Z for process
        int Af=0, Zf=0;
        auto& d = c->getDaughters();
        for (auto& daughter : d) {
            Af += daughter->getA();
            Zf += daughter->getZ();
        }
//        if (Ai!=Af || Zi!=Zf) return false; // ok for beta decay if electron is assigned Z=-1
        if (Ai != Af || Zi != Zf) {
            string message = "Decay of ";
            message += d[0]->getParent()->getName();
            message += " -> {";
            for (auto p : d) message += p->getName() + ", ";
            message.resize(message.size()-2);
            message += "} violates conservation of A and Z!";
            throw invalid_argument(message);
        }
    }
}


// Iterators
namespace simX {
    ProcessChain::iterator begin(const ProcessChain &c) {
        return ProcessChain::iterator(c.chain.begin());
    }

    ProcessChain::iterator end(const ProcessChain &c) {
        return ProcessChain::iterator(c.chain.end());
    }
}
