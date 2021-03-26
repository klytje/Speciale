#ifndef SIMX_S3_H
#define SIMX_S3_H 1

// simX header files
#include <simX/Detection/SegmentedDetector.h>

// C++ and ROOT header files
#include <TVector3.h>


namespace simX {

    /**
    * Class for S3 detectors
    */
    class S3 : public SegmentedDetector {

    public:
        S3(std::string name, TVector3 center, TVector3 normal, TVector3 up,
           UInt_t nSpokes, UInt_t nRings, double innerRadius, double ringWidth,
           double ringGap, double spokeGap,
           bool sharing,
           double activeThickness, double deadlayerFrontThickness,
           double deadlayerBackThickness,
           bool reverseSpokesOrdering = false, bool reverseRingOrdering = false);

        virtual ~S3();


        UInt_t spokeCount();

        UInt_t ringCount();

        virtual std::string description() override;


        UInt_t getNSpokes() const {
            return nSpokes;
        }

        UInt_t getNRings() const {
            return nRings;
        }

        double getRingWidth() const {
            return ringWidth;
        }

        double getRingGap() const {
            return ringGap;
        }

        double getSpokeGap() const {
            return spokeGap;
        }

        double getInnerRadius() const {
            return innerRadius;
        }

    protected:
        virtual bool strikesGrid(double u) override;

        virtual int gridIndex() override;

        struct Args;

    private:
        UInt_t nSpokes, nRings;
        double ringWidth, ringGap, spokeGap, innerRadius;
    };
}

#endif	/* SIMX_S3_H */
