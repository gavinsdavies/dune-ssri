#ifndef _TMS_RECO_H_SEEN_
#define _TMS_RECO_H_SEEN_

#include <vector>
#include <utility>
#include <cmath>

// For Hough line visualisation
#include "TF1.h"
#include "TH2D.h"

#include "TMS_Hit.h"
#include "TMS_Event.h"

// Maybe put this inside a separate namespace
class TMS_TrackFinder {

  public:

    static TMS_TrackFinder& GetFinder() {
      static TMS_TrackFinder Instance;
      return Instance;
    }

    void FindTracks(TMS_Event &event);
    const std::vector<TMS_Hit> & GetCandidates() { return Candidates; };

    TF1* GetHoughLine_zy() { return HoughLine_zy; };
    TF1* GetHoughLine_zx() { return HoughLine_zx; };

    int **GetAccumulator_zy() { return Accumulator_zy; };
    int **GetAccumulator_zx() { return Accumulator_zx; };

    TH2D *AccumulatorToTH2D(bool zy);

    void SetZMaxHough(double z) { zMaxHough = z;};


  private:
    TMS_TrackFinder();
    TMS_TrackFinder(TMS_TrackFinder const &) = delete;
    void operator=(TMS_TrackFinder const &) = delete;

    int FindBin(double Rho);
    // The candidates for each particle
    std::vector<TMS_Hit> Candidates;

    void Accumulate(double xvalue, double yvalue, double zvalue, TMS_Bar::BarType Type);

    // Number of theta bins
    int nTheta;
    int nRho;
    double RhoMin;
    double RhoMax;
    double ThetaMin;
    double ThetaMax;
    double ThetaWidth;
    double RhoWidth;
    double zMaxHough;

    int **Accumulator_zy;
    int **Accumulator_zx;

    TF1 *HoughLine_zy;
    TF1 *HoughLine_zx;
};

#endif
