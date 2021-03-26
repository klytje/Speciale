
#include "simX/Detection/SegmentedDetector.h"
#include <ausa/util/memory>
#include <simX/Logger.h>

using namespace AUSA::EnergyLoss;
using namespace AUSA::Calibration;
using namespace std;
using namespace simX;

// all length dimensions in units of mm
namespace {
    double boringFluctuation(int chan, double E) {
        (void) chan;
        (void) E;
        return 0;
    }

    DetectorCalibration buildDummyCalibration(size_t n) {
        return DetectorCalibration{vector<LinearCalibration>{n, LinearCalibration{0, 1}}};
    }
}


SegmentedDetector::SegmentedDetector(std::string name, std::unique_ptr<Args> args)
        : center0(TVector3(0,0,0)), normal(args->buildNormal()), up(args->buildOrientation()),
          doSharing(args->doSharing()), doIonizing(args->doIonizing()),

          detVol(args->buildDetectionVolume()), layers(args->buildLayers()),
          name(name), fluctuationFunction(boringFluctuation),
          calibration(buildDummyCalibration(static_cast<size_t>(args->nChannels()))),
          segmentation(args->buildSegmentation()),
          numberOfChannels(args->nChannels())
{
    // left-to-right unit vector
    right = up.Cross(normal);
    
    // move layers and detVol to specified position
    setPosition( args->buildCenter() );
}

const Detector::Volume& SegmentedDetector::getDetectorVolume() const {
    return *detVol;
}

const simX::Detector::DetectorOutput& SegmentedDetector::detect ( Particle& part ) {

    static auto log = log::getLogger("SegmentedDetector");

    // Clear output buffer
    detOut.clear();
    
    // Direction of motion
    auto& dir = part.getDirectionLab();

    // Which layer (if any) is the electrode grid
    int g = gridIndex();

    // Assume normal ordering (particle hits detector from front)
    // s=start, e=end, j=ordering
    int s = 0;
    int e = (int) layers.size();
    int j = 1;

    if (dir.Dot(normal) > 0) {  // Oh, we were wrong... (inverse ordering)
        j = -1;
        s = (int) (layers.size() - 1);
        e = -1;
    }

    log->debug("'{}' enters {} with {}keV", part.getName(), getName(), part.getKineticEnergyLab());
    // Loop over layers
    for (int i=s; i!=e; i+=j) {
        auto& l = layers[i];

        // Determine point of intersection and effective thickness
        TVector3 intersection;
        double t, dist;
        bool hit = l -> getIntersection( part.getPosition(), part.getDirectionLab(), intersection, t, dist);
        log->debug("'{}' hits {}mm layer", part.getName(), t);
        
        // If particle does not hit layer, move on to next layer
        if (!hit) break;
        
        // If dist>0, propagate detector to surface of layer (through vacuum)
        if (dist>0) part.setPosition( part.getPosition() + dist*part.getDirectionLab() );
        
        // get projected coordinates (u,v)   (u=left-right, v=down-up)
        TVector3 displ = part.getPosition() - l->getCenter();
        TVector3 proj = displ-displ.Dot(normal)*normal;
        double u = proj.Dot(right); 
        double v = proj.Dot(up);
        TVector2 uv(u,v);
        
        // Perform energy-loss calculation?
        double eloss = part.getKineticEnergyLab();
        bool doEloss = (i!=g) || strikesGrid(u); // grid is special case

        // non-ionizing energy loss
        double elossNonIoni = 0;

        // If no, just update particle position
        if (!doEloss) {
            eloss = 0;
            part.setPosition( part.getPosition()+t*part.getDirectionLab() );        
        }
        else { // Otherwise, do energy-loss calculation
            part.propagate(*l, part, -1, &elossNonIoni);
            eloss -= part.getKineticEnergyLab();
        }

        // ionizing energy loss
        double elossIoni = eloss;
        if (doIonizing) elossIoni -= elossNonIoni;
        
        // For the active volume, determine which strips have been affected,
        // and, for each, the amount of charge collected (energy signal) and the timing
        if (l->isActive()) 
            chargeCollection(uv, elossIoni, detOut);

        log->debug("'{}' loses {}keV in {} layer {}", part.getName(), eloss, layers[i]->getMaterial().getName(), i);
        // if particle has lost all its energy, end loop
        if (part.getKineticEnergyLab()<=0) break;
    }

    if (part.getKineticEnergyLab()!=0) {
        log->debug("'{}' leaves {} with {}keV", part.getName(), getName(), part.getKineticEnergyLab());
    } else {
        log->debug("'{}' is stopped in {}", part.getName(), getName());
    }
    return detOut;
}


void SegmentedDetector::chargeCollection( TVector2& uv, double& eloss, simX::Detector::DetectorOutput& detOut ) {

    // RECTANGULAR  =>  vertical and horizontal
    // CIRCULAR     =>  spokes and rings

    UInt_t prev = 0;

    for (int j = 0; j < segmentation.size(); ++j) {
        auto& seg = segmentation[j];

        vector<int> c;
        c.reserve(2);


        auto x0 = seg.coordinateStart;
        auto w = seg.segmentWidth;
        auto rev = seg.reversed;
        int n = seg.nSegments;

        auto x = seg.mapToCoordinate(uv);

        // Which strip has been hit?
        int i = static_cast<int>(std::floor((x-x0)/w));
        if (i < 0) return;
        int ch = i;
        if (rev) ch = n-1-i; // reverse strip ordering
        c.push_back(ch);


        // *** Sharing is assumed to occur if, and only if,
        // *** particle hits interstrip region. This is the simplest
        // *** possible model (more sophisticated models could
        // *** probably be designed which additionally would depend
        // *** on particle ID, energy, angle of entry, etc ...

        double f = 0;
        if (doSharing) {
            auto g = seg.calculateGap(uv);

            double dx = x - (x0+i*w);
            auto halfg = g/2;

            if (dx<halfg) { // sharing with previous strip
                int iprev = seg.modulate(i-1);
                int ch = iprev;
                if (rev) ch = n-1-iprev; // reverse strip ordering
                c.push_back(ch);
                f = 1.-dx/halfg;
            }
            else if (dx>w-halfg) { // sharing with next strip
                int inext = seg.modulate(i+1);
                int ch = inext;
                if (rev) ch = n-1-inext; // reverse strip ordering
                c.push_back(ch);
                f = (dx-(w-halfg))/halfg;
            }
        }

        // time
        double t = 0;

        // in case of no sharing
        if (c.size()==1) {
            int ch = c[0] + prev;
            double e = eloss + fluctuation(ch, eloss);
            DataTuple dat{ch, e, t};
            detOut.push_back(dat);
        }

            // in case of sharing
        else if (c.size()==2) {
            int offset = prev;
            int ch1 = c[0] + offset;
            int ch2 = c[1] + offset;
            double e1 = (1.0 - sharing(f)) * eloss;
            double e2 = sharing(f) * eloss;
            e1 += fluctuation(ch1, e1);
            DataTuple dat1{ch1, e1, t};
            detOut.push_back(dat1);

            if (0<=c[1] && c[1]<n) {
                e2 += fluctuation(ch2, e2);
                DataTuple dat2{ch2, e2, t};
                detOut.push_back(dat2);
            }
        }

        prev += n;
    }
}


int SegmentedDetector::getNumberOfChannels() {
    return numberOfChannels;
}


// Sharing distribution (deduced from CMAM 2008 data, talk to O. Kirsebom)
// x gives the position  (0:edge of strip, 1:midway between strips)
double SegmentedDetector::sharing( double x ) {
  double y=0.5*x-0.0423*sin(TMath::Pi()*x);
  return y;
}

// fluctuation in number of charges collected
double SegmentedDetector::fluctuation( int chan, double& e ) {
    return fluctuationFunction(chan, e);
}


void SegmentedDetector::setPosition( TVector3 c ) {
    // displacement vector
    TVector3 dc = c-center0;
    // move layers
    for (auto& l : layers) {
        TVector3 c0 = l->getCenter();
        l->setCenter( c0+dc );
    }
    // move detVol
    TVector3 c0 = detVol -> getCenter();
    detVol -> setCenter( c0+dc );
    // new center position
    center0 = c;
}


TVector3 SegmentedDetector::getPosition() const {
    return center0;
}

TVector3 SegmentedDetector::getDirection() const {
    return normal;
}

TVector3 SegmentedDetector::getUp() const {
    return up;
}


const std::string& SegmentedDetector::getName() const {
    return name;
}

Detector::FluctuationFunction SegmentedDetector::getFluctuationFunction() const {
    return fluctuationFunction;
}

void SegmentedDetector::setFluctuationFunction(FluctuationFunction fluctuationFunction) {
    SegmentedDetector::fluctuationFunction = fluctuationFunction;
}

unsigned int SegmentedDetector::energyToChannel(int channel, double e) {
    return calibration.energy2Adc(e, channel+1);
}

void SegmentedDetector::setCalibration(AUSA::Calibration::DetectorCalibration calibration) {
    if (calibration.size() != numberOfChannels) throw invalid_argument("Calibration must have " + to_string(numberOfChannels) + " channels.");

    SegmentedDetector::calibration = calibration;
}

void SegmentedDetector::reverseStripOrdering(size_t i, bool b) {
    segmentation[i].reversed = b;
}
