
#ifndef DECAYCHAIN_h
#define DECAYCHAIN_h 1

// simX headers
#include "simX/Particle.h"
#include "simX/NuclearProcess.h"

#include <tree.hh>
#include <memory>

namespace simX {

    class ProcessChain {

        public:
            using DecayPtr = std::unique_ptr<NuclearProcess>;
            using Tree = tree<DecayPtr>;

            /**
             * Create a ProcessChain with a specifec beam particle and reaction tree.
             * @param beam The beam particle.
             * @param chain The reaction chain this chain will be based on.
             *
             * @post This constructor will ensure that it is sound in terms of A and Z conservation.
             */
            ProcessChain( std::unique_ptr<Particle> beam, Tree& chain );
            ~ProcessChain();

            /**
            * Get first particle.
            */
            inline Particle& getBeam() { return *beam; }

            class iterator_ {
            public:
                using sub_iterator = std::vector<NuclearProcess*>::const_iterator;

                iterator_();
                iterator_(sub_iterator iter);

                NuclearProcess &     operator*() const;
                bool       operator==(const iterator_ &) const;
                bool       operator!=(const iterator_ &) const;
                iterator_ &  operator++();

            private:
                sub_iterator iter;
            };

            typedef iterator_ iterator;
            friend iterator begin(const ProcessChain & c);
            friend iterator end(const ProcessChain & c);

            const Tree& getTree() const { return t;}

            /**
             * Find first NuclearProcess where parent is named the input
             */
            NuclearProcess* const findDecayOf(const std::string& parent);
            
            /**
             * Returns pointer to first process
             */            
            NuclearProcess* getFirstProcess();

            std::vector<Particle*> findFinalStates();

        private:
            friend class iterator_;
            std::unique_ptr<Particle> beam;
            std::vector<NuclearProcess*> chain;
            Tree t;


            /**
             * Checks if A,Z add up for each nuclear process
             */
            void validate() const;
    };
}

#endif	/* DECAYCHAIN_H */
