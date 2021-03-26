#ifndef SIMX_DETECTION_W1_H
#define SIMX_DETECTION_W1_H 1

// simX header files
#include "simX/Detection/SegmentedDetector.h"

// C++ and ROOT header files
#include <TVector3.h>



namespace simX {

    /**
    * Class for W1 detectors
    */            
    class W1 : public SegmentedDetector {

        public:
            W1( std::string name,
                UInt_t nFrontStrip, UInt_t nBackStrip,
                TVector3 center, TVector3 normal, TVector3 up,
                bool sharing,
                double activeVolumeThickness
                    , double deadlayerFrontThickness
                    , double deadlayerBackThickness
                    , double gridThickness
                    , double contactBackThick
                    , double stripFrontWidth, double stripBackWidth
                    , double stripFrontGap, double stripBackGap
                    , double gridWidth
                    , bool reverseOrderingFront = false, bool reverseOrderingBack = false );


        UInt_t frontStripCount() const;

        UInt_t backStripCount() const;

        virtual std::string description() override;

        /**
         * Override simX defaults in order to comply with AUSAlib behaviour.
         * In essence AUSAlib calculates horizontal strip numbers top down.
         */
        virtual void reverseStripOrdering(size_t i, bool b ) override;

    protected:
        virtual bool strikesGrid(double u) override;
        virtual int gridIndex() override;

        struct Args;

    private:
        double sideLength, gridWidth;
        double stripWidthFront, stripWidthBack;
        double stripGapFront, stripGapBack;
        UInt_t nFrontStrip, nBackStrip;
        int gridI;
    };
}

#endif	/* SIMX_DETECTION_W1_H */
