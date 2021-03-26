
#ifndef WEIGHTCALCULATOR_H
#define	WEIGHTCALCULATOR_H


namespace simX {
    
    /**
    * Base class for calculating statistical of events/processes weights
    */
    class WeightCalculator {
        public:
            WeightCalculator() = default;
            virtual ~WeightCalculator() = default;

            /**
            * Returns weight
            */
            virtual double getWeight() const = 0;

        protected:

        private:
    };
}

#endif	/* WEIGHTCALCULATOR_H */

