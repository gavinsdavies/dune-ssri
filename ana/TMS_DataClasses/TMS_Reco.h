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

    aNode(double xval, double yval): 
      x(xval), y(yval), 
      HeuristicCost(-999), GroundCost(-999),
      ParentNodeID(-999), NodeID(-999) {
    };

    aNode(double xval, double yval, int ID): aNode(xval, yval) {
      NodeID = ID;
    };

    bool operator==(aNode const &other) {
      return (x == other.x && y == other.y);
    }

    /*
    bool operator<(aNode const &other) {
      return other.HeuristicCost < HeuristicCost;
    }
    */

    double Calculate(aNode const &other) {
      double deltax = x-other.x;
      double deltay = y-other.y;
      // All these units are in mm
      return (std::abs(deltax)+std::abs(deltay))/100;
    }

    void SetHeuristic(aNode const &other) {
      HeuristicCost = Calculate(other);
    }

    void SetGround(aNode const &other) {
      GroundCost = Calculate(other);
    }

    void SetHeuristic(double val) {
      HeuristicCost = val;
    }

    void SetGround(double val) {
      GroundCost = val;
    }

    // Position
    double x;
    double y;
    // Costs
    double HeuristicCost;
    double GroundCost;
    // Neighbours
    std::list<aNode> Neighbours;
    // Parent node
    int ParentNodeID;
    // self ID
    int NodeID;

    void Print() {
      std::cout << "NodeID: " << NodeID << std::endl;
      std::cout << "x, y = " << x << ", " << y << std::endl;
      std::cout << "Heuristic: " << HeuristicCost << " Ground: " << GroundCost << std::endl;
      std::cout << "Number of neighbours: " << Neighbours.size() << std::endl;
      std::cout << "ParentNodeID: " << ParentNodeID << std::endl;
    }
};

inline bool operator<(aNode const &a, aNode const &b) {
  return a.HeuristicCost < b.HeuristicCost;
}

/*
   bool CompareHeuristic(aNode const &a, aNode const &b) {
   return a.HeuristicCost > b.HeuristicCost;
   }
   */

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

    void BestFirstSearch(std::vector<TMS_Hit> TMS_Hits);


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
};

#endif
