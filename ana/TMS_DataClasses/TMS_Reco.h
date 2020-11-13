#ifndef _TMS_RECO_H_SEEN_
#define _TMS_RECO_H_SEEN_

#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>
#include <queue>
#include <list>
#include <unordered_map>

// For Hough line visualisation
#include "TF1.h"
#include "TH2D.h"

#include "TMS_Hit.h"
#include "TMS_Event.h"

// Utility class struct to store the node for track finding using A* or Best-First
class aNode {
  public:

    enum HeuristicType { kManhattan, kEuclidean, kUnkown };

    aNode(double xval, double yval, double ywval): 
      x(xval), y(yval), yw(ywval), 
      HeuristicCost(-999), NodeID(-999),
      Heuristic(kEuclidean) { // what calculator
    };

    aNode(double xval, double yval, double ywval, int ID): aNode(xval, yval, ywval) {
      NodeID = ID;
    };

    bool operator==(aNode const &other) {
      return (x == other.x && y == other.y);
    }

    double Calculate(aNode const &other) {
      // x is the plane number, y is in mm
      // jumping one plane incurs 10 ground, so reflect that here; jumping 2 planes (i.e. adjacent) should be 10 ground, jumping 4 planes (i.e. next to adjacent) is double that
      double deltax = (x-other.x)*5;
      // Moving 1 plane up is 10 ground cost, so reflect that here too
      double deltay = ((y-other.y)/yw)*10;
      if (Heuristic == kManhattan) return std::abs(deltax)+std::abs(deltay);
      else if (Heuristic == kEuclidean) return sqrt(deltax*deltax+deltay*deltay);
      else return 999;
      return 999;
    }

    void SetHeuristicCost(aNode const &other) {
      HeuristicCost = Calculate(other);
    }

    void SetHeuristicCost(double val) {
      HeuristicCost = val;
    }

    void SetHeuristic(HeuristicType a) {
      Heuristic = a;
    }

    // Position
    double x;
    double y;
    double yw;
    // Costs
    double HeuristicCost;
    // Neighbours
    std::unordered_map<aNode*, double> Neighbours;
    // self ID
    int NodeID;
    // What Heuristic do we use
    HeuristicType Heuristic;

    void Print() {
      std::cout << "NodeID: " << NodeID << std::endl;
      std::cout << "x, y = " << x << ", " << y << std::endl;
      std::cout << "Heuristic: " << HeuristicCost << std::endl; //" Ground: " << GroundCost << std::endl;
      std::cout << "Number of neighbours: " << Neighbours.size() << std::endl;
    }
};

inline bool operator<(aNode const &a, aNode const &b) {
  return a.HeuristicCost < b.HeuristicCost;
}

// Maybe put this inside a separate namespace
class TMS_TrackFinder {

  public:

    static TMS_TrackFinder& GetFinder() {
      static TMS_TrackFinder Instance;
      return Instance;
    }

    void FindTracks(TMS_Event &event);
    const std::vector<TMS_Hit> & GetCandidates() { return Candidates; };
    const std::vector<std::vector<TMS_Hit> >& GetTotalCandidates() { return TotalCandidates; };

    TF1* GetHoughLine_zy() { return HoughLine_zy; };
    TF1* GetHoughLine_zx() { return HoughLine_zx; };

    std::vector<TF1*> GetHoughLines_zy() { return HoughLines_zy; };
    std::vector<TF1*> GetHoughLines_zx() { return HoughLines_zx; };

    int **GetAccumulator_zy() { return Accumulator_zy; };
    int **GetAccumulator_zx() { return Accumulator_zx; };

    TH2D *AccumulatorToTH2D(bool zy);

    void SetZMaxHough(double z) { zMaxHough = z;};

    // Run a best first search
    void BestFirstSearch(const std::vector<TMS_Hit> &TMS_Hits);

    void HoughTransform(const std::vector<TMS_Hit> &TMS_Hits);

    // Clean up the hits, removing duplicates and zero entries
    std::vector<TMS_Hit> CleanHits(const std::vector<TMS_Hit> &TMS_Hits);
    // Get hits projected onto xz or yz
    std::vector<TMS_Hit> ProjectHits(const std::vector<TMS_Hit> &TMS_Hits, TMS_Bar::BarType bartype = TMS_Bar::kXBar);
    std::vector<TMS_Hit> RunAstar(const std::vector<TMS_Hit> &TMS_Hits);

    // Helper function to check if a hit is next to a gap
    bool NextToGap(double, double);

    void SpatialPrio(std::vector<TMS_Hit> &TMS_Hits);

  private:
    TMS_TrackFinder();
    TMS_TrackFinder(TMS_TrackFinder const &) = delete;
    void operator=(TMS_TrackFinder const &) = delete;

    int FindBin(double Rho);
    // The candidates for each particle
    std::vector<TMS_Hit> Candidates;
    std::vector<TMS_Hit> RawHits;

    std::vector<std::vector<TMS_Hit> > TotalCandidates;
    std::vector<TF1*> HoughLines_zy;
    std::vector<TF1*> HoughLines_zx;

    void Accumulate(double xvalue, double yvalue, double zvalue, TMS_Bar::BarType Type);

    // Number of theta bins
    /*
       int nTheta;
       int nRho;
       double RhoMin;
       double RhoMax;
       double ThetaMin;
       double ThetaMax;
       double ThetaWidth;
       double RhoWidth;
       */

    int nIntercept;
    int nSlope;
    double InterceptMin;
    double InterceptMax;
    double SlopeMin;
    double SlopeMax;
    double InterceptWidth;
    double SlopeWidth;

    int **Accumulator_zy;
    int **Accumulator_zx;

    TF1 *HoughLine_zy;
    TF1 *HoughLine_zx;
    //std::vector<TF1*> HoughLines_zy;
    //std::vector<TF1*> HoughLines_zx;

    double zMaxHough;

    unsigned int nMinHits;
    unsigned int nMaxMerges;

    double HighestCost;
    bool IsGreedy;
    // Which planes are next to the gaps (i.e. may cause discontinuities)?
    std::vector<int> PlanesNearGap;
};

#endif
