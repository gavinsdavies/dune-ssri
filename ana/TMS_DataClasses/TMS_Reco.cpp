#include "TMS_Reco.h"

TMS_TrackFinder::TMS_TrackFinder() :
  /*
  nTheta(2E3), // 1 degree accuracy
  nRho(2E3), // Somewhat arbitrary: comes from checking overlays with events
  RhoMin(-20E3),
  RhoMax(20E3),
  ThetaMin(0),
  ThetaMax(2*M_PI),
  ThetaWidth((ThetaMax-ThetaMin)/nTheta),
  RhoWidth((RhoMax-RhoMin)/nRho),
  */

  nIntercept(2E3),
  nSlope(2E3),
  InterceptMin(-25E3),
  InterceptMax(+25E3),
  SlopeMin(-1.5),
  SlopeMax(1.5),
  InterceptWidth((InterceptMax-InterceptMin)/nIntercept),
  SlopeWidth((SlopeMax-SlopeMin)/nSlope),
  zMaxHough((TMS_Const::TMS_Trans_Start+TMS_Const::TMS_Det_Offset[2])*10), // Max z for us to do Hough in, here choose transition layer
  nMinHits(20), // Minimum number of hits required to run track finding
  nMaxMerges(1), // Maximum number of merges for one hit
  // Initialise Highest cost to be very large
  HighestCost(999)
{
  // Apply the maximum Hough transform to the zx not zy: all bending happens in zx

  /*
  Accumulator_zy = new int*[nTheta];
  Accumulator_zx = new int*[nTheta];
  for (int i = 0; i < nTheta; ++i) {
    Accumulator_zy[i] = new int[nRho];
    Accumulator_zx[i] = new int[nRho];
  }
  */
  Accumulator_zy = new int*[nSlope];
  Accumulator_zx = new int*[nSlope];
  for (int i = 0; i < nSlope; ++i) {
    Accumulator_zy[i] = new int[nIntercept];
    Accumulator_zx[i] = new int[nIntercept];
  }

  // Make the TF1
  //HoughLine_zy = new TF1("LinearHough_zy", "[0]*x+[2]/[1]", (735+411)*10, (1450+411)*10);
  //HoughLine_zx = new TF1("LinearHough_zx", "[0]*x+[2]/[1]", (735+411)*10, zMaxHough);

  HoughLine_zy = new TF1("LinearHough_zy", "[0]+[1]*x", (TMS_Const::TMS_Thin_Start+TMS_Const::TMS_Det_Offset[2])*10, (1450+TMS_Const::TMS_Det_Offset[2])*10);
  HoughLine_zx = new TF1("LinearHough_zx", "[0]+[1]*x", (TMS_Const::TMS_Thin_Start+TMS_Const::TMS_Det_Offset[2])*10, zMaxHough);
  HoughLine_zy->SetLineStyle(kDashed);
  HoughLine_zx->SetLineStyle(kDashed);
  HoughLine_zy->SetLineColor(kRed);
  HoughLine_zx->SetLineColor(kRed);
}

// The generic track finder
void TMS_TrackFinder::FindTracks(TMS_Event &event) {
  std::cout << "Event: " << event.GetEventNumber() << std::endl;
  //

  // Reset the candidate vector
  Candidates.clear();
  RawHits.clear();
  TotalCandidates.clear();
  HoughLines_zy.clear();
  HoughLines_zx.clear();

  // Get the raw unmerged and untracked hits
  //RawHits = event.GetHits();

  std::vector<TMS_Hit> TMS_Hits = event.GetHits();
  // Require 20 hits
  if (TMS_Hits.size() < nMinHits) return;
  BestFirstSearch(TMS_Hits);
  return;

  // Maybe split up into zx and zy hits here

  // Multi-thread this?
  /*
  for (int i = 0; i < nTheta; ++i) {
    for (int j = 0; j < nRho; ++j) {
      Accumulator_zy[i][j] = 0;
      Accumulator_zx[i][j] = 0;
    }
  }
  */

  for (int i = 0; i < nSlope; ++i) {
    for (int j = 0; j < nIntercept; ++j) {
      Accumulator_zy[i][j] = 0;
      Accumulator_zx[i][j] = 0;
    }
  }

  // First run a simple Hough Transform
  for (std::vector<TMS_Hit>::iterator it = TMS_Hits.begin(); it!=TMS_Hits.end(); ++it) {
    TMS_Hit hit = (*it);
    TMS_Bar bar = hit.GetBar();
    //double EnergyDeposit = hit.GetE();
    //double Time = hit.GetT();
    double xhit = bar.GetX();
    double yhit = bar.GetY();
    double zhit = bar.GetZ();

    TMS_Bar::BarType Type = bar.GetBarType();
    // If z position is above region of interest, ignore hit
    if (Type == TMS_Bar::kYBar && zhit > zMaxHough) continue;

    Accumulate(xhit, yhit, zhit, Type);
  }

  // Get the maximum in the Accumulator
  /*
  double max_zy = 0;
  double max_zx = 0;
  int max_zy_theta_bin = 0;
  int max_zy_rho_bin = 0;
  int max_zx_theta_bin = 0;
  int max_zx_rho_bin = 0;
  //for (int i = 0; i < nTheta; ++i) {
    //for (int j = 0; j < nRho; ++j) {
  for (int i = 0; i < nSlope; ++i) {
    for (int j = 0; j < nIntercept; ++j) {
      if (Accumulator_zy[i][j] > max_zy) {
        max_zy = Accumulator_zy[i][j];
        max_zy_theta_bin = i;
        max_zy_rho_bin = j;
      }
      if (Accumulator_zx[i][j] > max_zx) {
        max_zx = Accumulator_zx[i][j];
        max_zx_theta_bin = i;
        max_zx_rho_bin = j;
      }
    }
  }
  */
  double max_zy = 0;
  double max_zx = 0;
  int max_zy_slope_bin = 0;
  int max_zy_inter_bin = 0;
  int max_zx_slope_bin = 0;
  int max_zx_inter_bin = 0;
  for (int i = 0; i < nSlope; ++i) {
    for (int j = 0; j < nIntercept; ++j) {
      if (Accumulator_zy[i][j] > max_zy) {
        max_zy = Accumulator_zy[i][j];
        max_zy_slope_bin = i;
        max_zy_inter_bin = j;
      }
      if (Accumulator_zx[i][j] > max_zx) {
        max_zx = Accumulator_zx[i][j];
        max_zx_slope_bin = i;
        max_zx_inter_bin = j;
      }
    }
  }

  // Transform back to a y=mx+c equivalent
  /*
  double ThetaOpt_zx = ThetaMin+max_zx_theta_bin*2*M_PI/nTheta;
  double RhoOpt_zx = RhoMin+max_zx_rho_bin*(RhoMax-RhoMin)/nRho;

  double ThetaOpt_zy = ThetaMin+max_zy_theta_bin*2*M_PI/nTheta;
  double RhoOpt_zy = RhoMin+max_zy_rho_bin*(RhoMax-RhoMin)/nRho;

  HoughLine_zy->SetParameter(0, tan(ThetaOpt_zy));
  HoughLine_zy->SetParameter(1, cos(ThetaOpt_zy));
  HoughLine_zy->SetParameter(2, RhoOpt_zy);

  HoughLine_zx->SetParameter(0, tan(ThetaOpt_zx));
  HoughLine_zx->SetParameter(1, cos(ThetaOpt_zx));
  HoughLine_zx->SetParameter(2, RhoOpt_zx);
  */

  double InterceptOpt_zx = InterceptMin+max_zx_inter_bin*(InterceptMax-InterceptMin)/nIntercept;
  double SlopeOpt_zx = SlopeMin+max_zx_slope_bin*(SlopeMax-SlopeMin)/nSlope;

  double InterceptOpt_zy = InterceptMin+max_zy_inter_bin*(InterceptMax-InterceptMin)/nIntercept;
  double SlopeOpt_zy = SlopeMin+max_zy_slope_bin*(SlopeMax-SlopeMin)/nSlope;
  HoughLine_zy->SetParameter(0, InterceptOpt_zy);
  HoughLine_zy->SetParameter(1, SlopeOpt_zy);

  HoughLine_zx->SetParameter(0, InterceptOpt_zx);
  HoughLine_zx->SetParameter(1, SlopeOpt_zx);

  HoughLines_zy.push_back(HoughLine_zy);
  HoughLines_zx.push_back(HoughLine_zx);

  // Then run a clustering on the Hough Transform
  // Hough transform is most likely to pick out straigh features at begining of track candidate, so start looking there
  TF1 *HoughLine = NULL;

  // Make a hard copy since we don't want to remove the original hits
  // Just make a pool of hits, and tracked hits
  // We'll pull hit out of the HitPool and put them into Candidates
  std::vector<TMS_Hit> HitPool = TMS_Hits;

  // Move hits from the pool into the candidates, and remove the candidates from the pool
  // HitPool is going to shrink or stay the same here
  for (std::vector<TMS_Hit>::iterator it = HitPool.begin(); it!=HitPool.end();) {
    TMS_Hit hit = (*it);
    TMS_Bar bar = hit.GetBar();
    double zhit = bar.GetZ();
    TMS_Bar::BarType Type = bar.GetBarType();
    // If z position is above region of interest, ignore hit
    if (Type == TMS_Bar::kYBar && zhit > zMaxHough) {
      ++it;
      continue;
    }

    if (Type == TMS_Bar::kXBar) {
      HoughLine = HoughLine_zy;
    } else if (Type == TMS_Bar::kYBar) {
      HoughLine = HoughLine_zx;
    }

    double HoughPoint = HoughLine->Eval(zhit);

    // Hough point is inside bar -> start clustering around bar
    if (bar.Contains(HoughPoint, zhit)) {
      Candidates.push_back(std::move(hit));
      // Remove from pool of hits
      it = HitPool.erase(it);
    } else {
      ++it;
    }
  }

  //MergeAdjacent(Candidates, HitPool);

  // Loop over the candidates, and add new adjacent candidates to the end
  size_t CandSize = Candidates.size();
  for (size_t i = 0; i < CandSize; ++i) {
    TMS_Hit Candidate = Candidates[i];
    TMS_Bar CandidateBar = Candidate.GetBar();
    int CandidatePlaneNumber = CandidateBar.GetPlaneNumber();
    TMS_Bar::BarType CandidateBarType = CandidateBar.GetBarType();

    // Count the number of times a candidate merges adjacent hits
    unsigned int nMerges = 0;
    bool InGap = false;
    double CandPos = CandidateBar.GetNotZ();
    double CandPosWidth = CandidateBar.GetNotZw();
    // First check if this hit may be affected by gap region issues
    // If so, allow search range to extend in z
    if (CandidateBarType == TMS_Bar::kYBar) {
      /*
      std::cout << "CandPos: " << CandidateBar.GetNotZ() << std::endl;
      std::cout << "Plane: " << CandidatePlaneNumber << std::endl;
      */
      double pos = CandPos + CandPosWidth;
      double neg = CandPos - CandPosWidth;
      // Check the top
      if ((pos > TMS_Const::TMS_Dead_Top_T[0] && pos < TMS_Const::TMS_Dead_Top_T[1]) ||
          (neg < TMS_Const::TMS_Dead_Top_T[1] && neg > TMS_Const::TMS_Dead_Top_T[0])) InGap = true;
      // Check the center
      else if ((pos > TMS_Const::TMS_Dead_Center_T[0] && pos < TMS_Const::TMS_Dead_Center_T[1]) ||
          (neg < TMS_Const::TMS_Dead_Center_T[1] && neg > TMS_Const::TMS_Dead_Center_T[0])) InGap = true;
      // Check the bottom
      else if ((pos > TMS_Const::TMS_Dead_Bottom_T[0] && pos < TMS_Const::TMS_Dead_Bottom_T[1]) ||
          (neg < TMS_Const::TMS_Dead_Bottom_T[1] && neg > TMS_Const::TMS_Dead_Bottom_T[0])) InGap = true;
      else InGap = false;
    }

    /*
    if (InGap) {
      std::cout << CandPos << " is close to gap" << std::endl;
    }
    */

    // Now loop over each hit
    for (std::vector<TMS_Hit>::iterator jt = HitPool.begin(); jt != HitPool.end();) {

      // If we've already merged more than we're allowed to
      if (nMerges >= nMaxMerges) break;
      TMS_Hit PoolHit = (*jt);
      TMS_Bar PoolBar = PoolHit.GetBar();
      int PoolPlaneNumber = PoolBar.GetPlaneNumber();
      TMS_Bar::BarType PoolBarType = PoolBar.GetBarType();

      // Only match the same type of bars (x bars with x bars, y bars with ybars)
      if (PoolBarType != CandidateBarType) {
        ++jt;
        continue;
      }

      // Now check the distance in x or y depending on bar
      double PoolPos = PoolBar.GetNotZ();
      double PoolPosWidth = PoolBar.GetNotZw();

      /*
      if (CandidateBarType == TMS_Bar::kYBar) {
        std::cout << PoolPos << "+/-" << PoolPosWidth << std::endl;
        std::cout << PoolPlaneNumber << std::endl;
      }
      */

      // Ensure adjacent or matching in z
      // Make an exception for the airgaps; allow any range in z
      if (!InGap && abs(CandidatePlaneNumber-PoolPlaneNumber) > 2) {
        ++jt;
        continue;
      }

      // Should we merge this hit
      bool Merge = false;
      // If this hasn't been merged already, attempt to do so
      if (nMerges < nMaxMerges) {
        // If the Candidate bar contains the pool bar + it's width, they're adjacent
        // Or even better, the adjacent hits have maxing not z

        // Preferentially merge in z (same NotZ position)
        if (CandidateBar.Contains(PoolPos, CandidateBar.GetZ())) {
          Merge = true;
          // Then check merge in NotZ
        } else if (CandidateBar.Contains(PoolPos + PoolPosWidth, CandidateBar.GetZ()) || 
            CandidateBar.Contains(PoolPos - PoolPosWidth, CandidateBar.GetZ())) {
          Merge = true;
          // Then check two bars away for the xz view
        } else if (CandidateBarType == TMS_Bar::kYBar && 
            (CandidateBar.Contains(PoolPos + 2*PoolPosWidth, CandidateBar.GetZ()) || 
             CandidateBar.Contains(PoolPos - 2*PoolPosWidth, CandidateBar.GetZ()))) {
          Merge = true;
        } /*else if (CandidateBarType == TMS_Bar::kYBar && 
            (CandidateBar.Contains(PoolPos + 3*PoolPosWidth, CandidateBar.GetZ()) || 
             CandidateBar.Contains(PoolPos - 3*PoolPosWidth, CandidateBar.GetZ()))) {
          Merge = true;
        }
        */

        // Make a special arrangement for if we're next to the gap
        else if (InGap && ( 
              PoolBar.Contains(TMS_Const::TMS_Dead_Top_T[1]+PoolPosWidth, PoolBar.GetZ()) ||
              PoolBar.Contains(TMS_Const::TMS_Dead_Top_T[0]-PoolPosWidth, PoolBar.GetZ()) ||
              PoolBar.Contains(TMS_Const::TMS_Dead_Center_T[1]+PoolPosWidth, PoolBar.GetZ()) ||
              PoolBar.Contains(TMS_Const::TMS_Dead_Center_T[0]-PoolPosWidth, PoolBar.GetZ()) ||
              PoolBar.Contains(TMS_Const::TMS_Dead_Bottom_T[1]+PoolPosWidth, PoolBar.GetZ()) ||
              PoolBar.Contains(TMS_Const::TMS_Dead_Bottom_T[0]-PoolPosWidth, PoolBar.GetZ()) )) {
          Merge = true;
        }

      }

      if (Merge) {
        //if (CandidateBarType == TMS_Bar::kYBar) std::cout << "  Merging " << PoolPos << "+/-" << PoolPosWidth << std::endl;
        Candidates.push_back(std::move(PoolHit));
        // Increment the number of candidates in the vector
        CandSize++;
        jt = HitPool.erase(jt);
        // Increment the merges for this candidate
        nMerges++;
      } else {
        ++jt;
      }
    }
  }

  // Now see if the yz and xz views roughly agree on stop and start positions of the main track
  // If not, it implies there may be two tracks and the yz view has selected one whereas xz selected theother, or a broken track in one of the views
  int StartPlane_xz = 99999999;
  int EndPlane_xz = -99999999;
  int StartPlane_yz = 99999999;
  int EndPlane_yz = -99999999;

  for (std::vector<TMS_Hit>::iterator it = Candidates.begin(); it != Candidates.end(); ++it) {
    TMS_Bar Bar = (*it).GetBar();
    // Compare the plane number
    int PlaneNumber = Bar.GetPlaneNumber();
    // if this bar is an x-bar, look at yz hits
    if (Bar.GetBarType() == TMS_Bar::kXBar) {
      if (PlaneNumber < StartPlane_yz) StartPlane_yz = PlaneNumber;
      if (PlaneNumber > EndPlane_yz) EndPlane_yz = PlaneNumber;

      // if this bar is an x-bar, look at yz hits
    } else if (Bar.GetBarType() == TMS_Bar::kYBar) {
      if (PlaneNumber < StartPlane_xz) StartPlane_xz = PlaneNumber;
      if (PlaneNumber > EndPlane_xz) EndPlane_xz = PlaneNumber;
    }

  }

  //bool Matched_yz_xz = true;
  if (abs(StartPlane_yz - StartPlane_xz) != 1 || abs(EndPlane_yz - EndPlane_xz) != 1) {
    //Matched_yz_xz = false;
    std::cout << "Not matching end/start plane" << std::endl;
    std::cout << "StartPlane yz: " << StartPlane_yz << ", EndPlane yz: " << EndPlane_yz << std::endl;
    std::cout << "StartPlane xz: " << StartPlane_xz << ", EndPlane xz: " << EndPlane_xz << std::endl;
  }

  /*
  // Do a Hough transform on the remaining hits
  // If this Hough transform matches within some number the original Hough transform parameters, merge the tracks
  // Reset the accumulator (or make a new one?)
  for (int i = 0; i < nSlope; ++i) {
    for (int j = 0; j < nIntercept; ++j) {
      Accumulator_zy[i][j] = 0;
      Accumulator_zx[i][j] = 0;
    }
  }

  // First run a simple Hough Transform
  for (std::vector<TMS_Hit>::iterator it = HitPool.begin(); it!=HitPool.end(); ++it) {
    TMS_Hit hit = (*it);
    TMS_Bar bar = hit.GetBar();
    //double EnergyDeposit = hit.GetE();
    //double Time = hit.GetT();
    double xhit = bar.GetX();
    double yhit = bar.GetY();
    double zhit = bar.GetZ();

    TMS_Bar::BarType Type = bar.GetBarType();
    // Scan over the entire range this time
    // Or maybe the z region in which the other view's track has been found?

    Accumulate(xhit, yhit, zhit, Type);
  }

  max_zy = 0;
  max_zx = 0;
  max_zy_slope_bin = 0;
  max_zy_inter_bin = 0;
  max_zx_slope_bin = 0;
  max_zx_inter_bin = 0;
  for (int i = 0; i < nSlope; ++i) {
    for (int j = 0; j < nIntercept; ++j) {
      if (Accumulator_zy[i][j] > max_zy) {
        max_zy = Accumulator_zy[i][j];
        max_zy_slope_bin = i;
        max_zy_inter_bin = j;
      }
      if (Accumulator_zx[i][j] > max_zx) {
        max_zx = Accumulator_zx[i][j];
        max_zx_slope_bin = i;
        max_zx_inter_bin = j;
      }
    }
  }

  double InterceptOpt_zx_2 = InterceptMin+max_zx_inter_bin*(InterceptMax-InterceptMin)/nIntercept;
  double SlopeOpt_zx_2 = SlopeMin+max_zx_slope_bin*(SlopeMax-SlopeMin)/nSlope;

  double InterceptOpt_zy_2 = InterceptMin+max_zy_inter_bin*(InterceptMax-InterceptMin)/nIntercept;
  double SlopeOpt_zy_2 = SlopeMin+max_zy_slope_bin*(SlopeMax-SlopeMin)/nSlope;

  TF1 *HoughLine_zy_2 = (TF1*)HoughLine_zy->Clone("HZY_2");
  HoughLine_zy_2->SetLineColor(kWhite);
  HoughLine_zy_2->SetParameter(0, InterceptOpt_zy_2);
  HoughLine_zy_2->SetParameter(1, SlopeOpt_zy_2);

  TF1 *HoughLine_zx_2 = (TF1*)HoughLine_zx->Clone("HZX_2");
  HoughLine_zx_2->SetLineColor(kWhite);
  HoughLine_zx_2->SetParameter(0, InterceptOpt_zx_2);
  HoughLine_zx_2->SetParameter(1, SlopeOpt_zx_2);
  */

  /*
     std::cout << InterceptOpt_zx << " " << InterceptOpt_zx_2 << std::endl;
     std::cout << SlopeOpt_zx << " " << SlopeOpt_zx_2 << std::endl;

     std::cout << InterceptOpt_zy << " " << InterceptOpt_zy_2 << std::endl;
     std::cout << SlopeOpt_zy << " " << SlopeOpt_zy_2 << std::endl;
     */

  /*
  // Get the hits that intersect the second hough line
  // First make a second candidate vector
  std::vector<TMS_Hit> Candidates2;
  for (std::vector<TMS_Hit>::iterator it = HitPool.begin(); it!=HitPool.end(); ) {
    TMS_Hit hit = (*it);
    TMS_Bar bar = hit.GetBar();
    double zhit = bar.GetZ();
    TMS_Bar::BarType Type = bar.GetBarType();
    if (Type == TMS_Bar::kXBar) {
      HoughLine = HoughLine_zy_2;
    } else if (Type == TMS_Bar::kYBar) {
      HoughLine = HoughLine_zx_2;
    }

    double HoughPoint = HoughLine->Eval(zhit);

    // Hough point is inside bar -> start clustering around bar
    if (bar.Contains(HoughPoint, zhit)) {
      Candidates2.push_back(std::move(hit));
      // Remove from pool of hits
      it = HitPool.erase(it);
    } else {
      ++it;
    }
  }

  // Loop over the candidates, and add new adjacent candidates to the end
  CandSize = Candidates2.size();
  for (size_t i = 0; i < CandSize; ++i) {
    TMS_Hit Candidate = Candidates2[i];
    TMS_Bar CandidateBar = Candidate.GetBar();
    int CandidatePlaneNumber = CandidateBar.GetPlaneNumber();
    TMS_Bar::BarType CandidateBarType = CandidateBar.GetBarType();

    // Count the number of times a candidate merges adjacent hits
    unsigned int nMerges = 0;
    // Now loop over each hit
    for (std::vector<TMS_Hit>::iterator jt = HitPool.begin(); jt != HitPool.end();) {
      TMS_Hit PoolHit = (*jt);
      TMS_Bar PoolBar = PoolHit.GetBar();
      int PoolPlaneNumber = PoolBar.GetPlaneNumber();
      TMS_Bar::BarType PoolBarType = PoolBar.GetBarType();

      // Only match the same type of bars (x bars with x bars, y bars with ybars)
      if (PoolBarType != CandidateBarType) {
        ++jt;
        continue;
      }

      // Ensure adjacent or matching in z
      // Need some exception here for the airgaps! FIGURE IT OUT
      if (abs(CandidatePlaneNumber-PoolPlaneNumber) > 2) {
        ++jt;
        continue;
      }

      // Now check the distance in x or y depending on bar
      double PoolPos = PoolBar.GetNotZ();
      double PoolPosWidth = PoolBar.GetNotZw();
      // If the Candidate bar contains the pool bar + it's width, they're adjacent
      // Or even better, the adjacent hits have maxing not z
      if (CandidateBar.Contains(PoolPos + PoolPosWidth, CandidateBar.GetZ()) || CandidateBar.Contains(PoolPos - PoolPosWidth, CandidateBar.GetZ()) || CandidateBar.Contains(PoolPos, CandidateBar.GetZ())) {
        Candidates2.push_back(std::move(PoolHit));
        // Increment the number of candidates in the vector
        CandSize++;
        jt = HitPool.erase(jt);
        // Increment the merges for this candidate
        nMerges++;
      } else {
        ++jt;
      }
    }
  }

  // Now we *should* have two candidates?
  std::cout << "Total hits: " << TMS_Hits.size() << std::endl;
  std::cout << "Candidates (Hough+Cluster): " << Candidates.size() << std::endl;
  std::cout << "Candidates2 (Leftover+Hough+cluster): " << Candidates2.size() << std::endl;
  std::cout << "Hitpool: " << HitPool.size() << std::endl;

  // Intersection in z (x-axis) of the two Hough lines
  double intersection_zx = (InterceptOpt_zx-InterceptOpt_zx_2)/(SlopeOpt_zx_2-SlopeOpt_zx);
  double intersection_zy = (InterceptOpt_zy-InterceptOpt_zy_2)/(SlopeOpt_zy_2-SlopeOpt_zy);
  std::cout << "Intersection zx: " << intersection_zx << std::endl;
  std::cout << "Intersection zy: " << intersection_zy << std::endl;

  TotalCandidates.push_back(Candidates);
  //TotalCandidates.push_back(Candidates2);

  HoughLines_zy.push_back(HoughLine_zy);
  //HoughLines_zy.push_back(HoughLine_zy_2);

  HoughLines_zx.push_back(HoughLine_zx);
  //HoughLines_zx.push_back(HoughLine_zx_2);

  // If this Hough transform matches the start and end position in z of the other view, choose the longest track as the muon and make 2 candidate tracks



  // Last scan tries looks at hits that occur around the gap regions
  // Then use this to inform the yz reconstruction (since broken tracks will be due to holes in x, not y)


  // Run a Hough transform on the remaining HitPool to find second tracks
  */


  /*
     std::cout << "After Hough+adjacent checks: " << std::endl;
     std::cout << "Of " << TMS_Hits.size()  << " hits " << Candidates.size() << " are candidates and " << HitPool.size() << " are left" << std::endl;
     */


}

/*
   void TMS_TrackFinder::MergeAdjacent(std::vector<TMS_Hit> Candidates, std::vector<TMS_Hit> Pool) {

// Loop over the candidates, and add new adjacent candidates to the end
size_t CandSize = Candidates.size();
for (size_t i = 0; i < CandSize; ++i) {
TMS_Hit Candidate = Candidates[i];
TMS_Bar CandidateBar = Candidate.GetBar();
int CandidatePlaneNumber = CandidateBar.GetPlaneNumber();
TMS_Bar::BarType CandidateBarType = CandidateBar.GetBarType();

// Count the number of times a candidate merges adjacent hits
int nMerges = 0;
// Now loop over each hit
for (std::vector<TMS_Hit>::iterator jt = HitPool.begin(); jt != HitPool.end();) {
TMS_Hit PoolHit = (*jt);
TMS_Bar PoolBar = PoolHit.GetBar();
int PoolPlaneNumber = PoolBar.GetPlaneNumber();
TMS_Bar::BarType PoolBarType = PoolBar.GetBarType();

// Only match the same type of bars (x bars with x bars, y bars with ybars)
if (PoolBarType != CandidateBarType) {
++jt;
continue;
}

// Ensure adjacent or matching in z
// Need some exception here for the airgaps! FIGURE IT OUT
if (abs(CandidatePlaneNumber-PoolPlaneNumber) > 2) {
++jt;
continue;
}

// Now check the distance in x or y depending on bar
double PoolPos = PoolBar.GetNotZ();
double PoolPosWidth = PoolBar.GetNotZw();
// If the Candidate bar contains the pool bar + it's width, they're adjacent
// Or even better, the adjacent hits have maxing not z
if (CandidateBar.Contains(PoolPos + PoolPosWidth, CandidateBar.GetZ()) || CandidateBar.Contains(PoolPos - PoolPosWidth, CandidateBar.GetZ()) || CandidateBar.Contains(PoolPos, CandidateBar.GetZ())) {
Candidates.push_back(std::move(PoolHit));
// Increment the number of candidates in the vector
CandSize++;
jt = HitPool.erase(jt);
// Increment the merges for this candidate
nMerges++;
} else {
++jt;
}
}
}
}
*/

// Find the bin for the accumulator
//int TMS_TrackFinder::FindBin(double rho) {
int TMS_TrackFinder::FindBin(double c) {
  // Since we're using uniform binning no need for binary search or similar
  //int bin = (rho-RhoMin)/RhoWidth;
  int bin = (c-InterceptMin)/InterceptWidth;
  //std::cout << bin << " " << rho << " " << InterceptMin << " " << InterceptWidth << std::endl;
  return bin;
}

// xvalue is x-axis, y value is y-axis
void TMS_TrackFinder::Accumulate(double xhit, double yhit, double zhit, TMS_Bar::BarType Type) {

  double yvar = -9999999999;
  int **Accumulator = NULL;
  if (Type == TMS_Bar::kXBar) {
    yvar = yhit;
    Accumulator = Accumulator_zy;
  } else if (Type == TMS_Bar::kYBar) {
    yvar = xhit;
    Accumulator = Accumulator_zx;
  }

  // Could probably multi-thread this operation
  // Now do the Hough
  //for (int i = 0; i < nTheta; ++i) {
  for (int i = 0; i < nSlope; ++i) {
    //float theta = ThetaMin+i*ThetaWidth;
    double m = SlopeMin+i*SlopeWidth;

    // Now calculate rho
    //float rho = (yvar-zhit*tan(theta))*cos(theta);
    double c = yvar-m*zhit;

    // Find which rho bin this corresponds to
    //int rho_bin = FindBin(rho);
    int c_bin = FindBin(c);
    //std::cout << c << " " << c_bin << std::endl;
    //if (c > InterceptMax || c < InterceptMin) {
    //std::cout << "Found c greater or less than maxium/minimum" << std::endl;
    //}

    // Fill the accumulator
    //Accumulator[i][rho_bin]++;
    Accumulator[i][c_bin]++;
  }
}

// Convert Accumulator to a TH2D
TH2D *TMS_TrackFinder::AccumulatorToTH2D(bool zy) {

  std::string Name;
  int **Accumulator;
  if (zy) {
    Name = "TMS_TrackFinder_Accumulator_zy";
    Accumulator = Accumulator_zy;
  } else {
    Name = "TMS_TrackFinder_Accumulator_zx";
    Accumulator = Accumulator_zx;
  }
  //TH2D *accumulator = new TH2D(Name.c_str(), (Name+";#theta (rad.);#rho (mm)").c_str(), nTheta, ThetaMin, ThetaMax, nRho, RhoMin, RhoMax);
  TH2D *accumulator = new TH2D(Name.c_str(), (Name+";m (slope);c (intercept) (mm)").c_str(), nSlope, SlopeMin, SlopeMax, nIntercept, InterceptMin, InterceptMax);

  /*
     for (int i = 0; i < nTheta; ++i) {
     for (int j = 0; j < nRho; ++j) {
     accumulator->SetBinContent(i+1, j+1, Accumulator[i][j]);
     }
     }
     */
  for (int i = 0; i < nSlope; ++i) {
    for (int j = 0; j < nIntercept; ++j) {
      accumulator->SetBinContent(i+1, j+1, Accumulator[i][j]);
    }
  }

  int maxx, maxy, maxz;
  accumulator->GetMaximumBin(maxx, maxy, maxz);
  double maxtheta = accumulator->GetXaxis()->GetBinCenter(maxx);
  double maxrho = accumulator->GetYaxis()->GetBinCenter(maxy);
  //accumulator->SetTitle(Form("#splitline{%s}{#theta_{max}=%.2f #rho_{max}=%.2f}", accumulator->GetTitle(), maxtheta, maxrho));
  accumulator->SetTitle(Form("#splitline{%s}{m_{max}=%.2f c_{max}=%.2f}", accumulator->GetTitle(), maxtheta, maxrho));

  // Set the minimum (easier to draw)
  double maxcounts = accumulator->GetMaximum();
  accumulator->SetMinimum(maxcounts*0.8);

  return accumulator;
}

// Implement A* algorithm for track finding, starting with most upstream to most downstream hit
void TMS_TrackFinder::BestFirstSearch(std::vector<TMS_Hit> TMS_Hits) {

  // First remove duplicate hits
  std::vector<TMS_Hit> TMS_Hits_Cleaned;

  for (std::vector<TMS_Hit>::iterator it = TMS_Hits.begin(); it != TMS_Hits.end(); ) {
    TMS_Hit hit = *it;
    TMS_Bar bar = hit.GetBar();
    // Maybe this hit has already been counted
    double z = bar.GetZ();
    double y = bar.GetNotZ();
    std::pair<double, double> temp = std::make_pair(z,y);
    bool Duplicate = false;

    // Skim out the duplicates
    for (std::vector<TMS_Hit>::iterator jt = TMS_Hits_Cleaned.begin(); jt != TMS_Hits_Cleaned.end(); ++jt) {
      TMS_Hit hit2 = *jt;
      double z2 = hit2.GetZ();
      double y2 = hit2.GetNotZ();
      std::pair<double, double> temp2 = std::make_pair(z2,y2);
      if (temp == temp2) {
        Duplicate = true;
        break;
      }
    }

    if (!Duplicate) {
      TMS_Hits_Cleaned.push_back(std::move(hit));
      it = TMS_Hits.erase(it);
    } else {
      ++it;
    }
  }

  // Now split in yz and xz hits
  std::vector<TMS_Hit> TMS_xz;
  std::vector<TMS_Hit> TMS_yz;

  for (std::vector<TMS_Hit>::iterator it = TMS_Hits_Cleaned.begin(); it != TMS_Hits_Cleaned.end(); ) {
    TMS_Hit hit = (*it);
    if (hit.GetBar().GetBarType() == TMS_Bar::kXBar) {
      TMS_yz.push_back(std::move(hit));
    } else if (hit.GetBar().GetBarType() == TMS_Bar::kYBar) {
      TMS_xz.push_back(std::move(hit));
    }
    it = TMS_Hits_Cleaned.erase(it);
  }

  // Now sort in z
  std::sort(TMS_yz.begin(), TMS_yz.end(), TMS_Hit::SortByZ);
  std::sort(TMS_xz.begin(), TMS_xz.end(), TMS_Hit::SortByZ);

  aNode First(TMS_xz.front().GetZ(), TMS_xz.front().GetNotZ(), 0);
  aNode Last(TMS_xz.back().GetZ(), TMS_xz.back().GetNotZ(), TMS_xz.size());

  // Calculate the Heuristic cost for each of the nodes to the last point
  // Use the Manhattan distance or Euclidean distance?
  // Also convert the TMS_Hit to our path-node
  std::vector<aNode> Nodes;
  // Give each node an ID
  int NodeID = 0;
  for (std::vector<TMS_Hit>::iterator it = TMS_xz.begin(); it != TMS_xz.end(); ++it, NodeID++) {

    // Make the node
    double x = (*it).GetZ();
    double y = (*it).GetNotZ();

    aNode TempNode(x, y, NodeID);
    TempNode.HeuristicCost = TempNode.Calculate(Last);
    TempNode.GroundCost = TempNode.Calculate(First);

    // Now need to find the neighbours of the node to decide on where to go
    // To find neighbours, look at the plane number and y
    int PlaneNumber = (*it).GetPlaneNumber();
    int NodeID_n = 0;
    for (std::vector<TMS_Hit>::iterator jt = TMS_xz.begin(); jt != TMS_xz.end(); ++jt, NodeID_n++) {
      int PlaneNumber_n = (*jt).GetPlaneNumber();
      int DeltaPlane = PlaneNumber_n-PlaneNumber;
      // Only look at adjacent bars in z (remember planes are alternating so adjacent xz bars are every two planes)
      if (abs(DeltaPlane) > 2) continue;
      double x_n = (*jt).GetZ();
      double y_n = (*jt).GetNotZ();
      double DeltaY = fabs(y_n-y);
      // Get the number of vertical bars (NotZw gives bar width)
      int DeltaY_Trunc = int(DeltaY/(*jt).GetNotZw());
      // Allow for three bars
      if (abs(DeltaY_Trunc) > 3) continue;

      aNode Neighbour(x_n, y_n, NodeID_n);
      Neighbour.SetHeuristic(Last);

      // And the ground cost will update as we go along?
      // The ground cost depends on if a diagonal or not is connected
      double GroundCost = -999;
      // left/right
      if (abs(DeltaY_Trunc) == 0 && abs(DeltaPlane) == 2) GroundCost = 10;
      // up/down
      else if (abs(DeltaY_Trunc) == 1 && abs(DeltaPlane) == 0) GroundCost = 10;
      // diagonal->prefer left/right followed by up/down
      else if (abs(DeltaY_Trunc) == 1 && abs(DeltaPlane) == 2) GroundCost = 30;
      // up/down by 2 (missing one bar)
      else if (abs(DeltaY_Trunc) == 2 && abs(DeltaPlane) == 0) GroundCost = 20;
      // up/down by 3 (missing three bars)
      else if (abs(DeltaY_Trunc) == 3 && abs(DeltaPlane) == 0) GroundCost = 30;
      // diagonal and missing 1 bar
      else if (abs(DeltaY_Trunc) == 2 && abs(DeltaPlane) == 2) GroundCost = 60;
      // diagonal and missing 2 bar
      else if (abs(DeltaY_Trunc) == 3 && abs(DeltaPlane) == 2) GroundCost = 90;
      // If we've truncated
      else if (abs(DeltaY_Trunc) == 0 && abs(DeltaPlane) == 0) GroundCost = 10;

      if (GroundCost == -999) {
        std::cout << "Missed exceptioN!" << std::endl;
        std::cout << "DeltaY_Trunc: " << DeltaY_Trunc << std::endl;
        std::cout << "DeltaY: " << DeltaY << std::endl;
        std::cout << "DeltaPlane: " << DeltaPlane << std::endl;
        std::cout << "getnotzW: " << (*jt).GetNotZw() << std::endl;
        throw;
      }

      Neighbour.SetGround(GroundCost);
      Neighbour.ParentNodeID = NodeID;
      Neighbour.NodeID = NodeID_n;
      TempNode.Neighbours.push_back(Neighbour);
    }
    //TempNode.Print();
    Nodes.push_back(TempNode);
  }

  //


  // Keep a map of the cost so far for reaching an aNode
  std::unordered_map<int, double> cost_so_far;
  // Keep a map of the where a given aNode was reached from
  std::unordered_map<int, int> came_from;
  // The priority queue of which node to next evaluate
  std::priority_queue<aNode, std::vector<aNode> > pq;

  Nodes.front().Print();
  std::cout << "to" << std::endl;
  Nodes.back().Print();

  std::cout << "Start loop..." << std::endl;
  pq.push(Nodes.front());
  cost_so_far[0] = 0;
  came_from[0] = 0;
  while (!pq.empty()) {
    //std::cout << "Next in queue..." << std::endl;
    aNode current = pq.top();
    //current.Print();
    pq.pop();
    if (current == Nodes.back()) break;

    // Loop over this node's neighbours
    for (aNode neighbour : current.Neighbours) {
      // Check the cost for getting to this node
      double new_cost = cost_so_far[current.NodeID]+neighbour.HeuristicCost+neighbour.GroundCost;
      // If we've never reached this node, or if this path to the node has less cost
      if (cost_so_far.find(neighbour.NodeID) == cost_so_far.end() ||
          new_cost < cost_so_far[neighbour.NodeID]) {
        /*
        std::cout << "updating" << std::endl;
        std::cout << current.NodeID << " to " << neighbour.NodeID << std::endl;
        std::cout << "old cost: " << cost_so_far[neighbour.NodeID] << " new cost: " << new_cost << std::endl;
        */
        cost_so_far[neighbour.NodeID] = new_cost;
        pq.push(Nodes.at(neighbour.NodeID));
        came_from[neighbour.NodeID] = current.NodeID;
      }
    }
  }

  for (auto pair: came_from) {
    std::cout << "Node " << pair.first << " came from " << pair.second << std::endl;
  }
  std::cout << " ********* " << std::endl;

  /*
  for (auto node : Nodes) {
    for (auto neighbour : node.Neighbours) {
      neighbour.Print();
    }
  }
  */


  // Can now make graphs to use

  // Starting point is the first vector entry in each (lowest z)

  // End point is the last vector entry in each (highest z)

  // First build the graph
  // Allow for moving up/down, left/right and all diagonals
  // Moving one plane to the right incurs no cost
  // Moving one plane to the left incurs one cost
  // Moving one bar up/down incurs one cost
  // Moving diagonally incurs three cost (to prefer up/down followed my left/right)
  //
  // Then also explore skipping one? i.e. connect over missing 
  // Probably don't allow skipping a z layer, but allow skipping one/two/three notZ layers

















  // First implementation of best-first search
  // Make connection costs for each of the hits
  /*
  std::vector<std::vector<ThreePair> > PositionCost;
  for (std::vector<TMS_Hit>::iterator it = TMS_yz.begin(); it != TMS_yz.end(); ++it) {
    TMS_Hit hit = (*it);
    //double z = hit.GetZ();
    double y = hit.GetNotZ();
    int planeno = hit.GetPlaneNumber();

    // Make the vector of costs for this hit
    std::vector<ThreePair> ThisCost;
    //std::cout << "Planenumber: " << planeno << std::endl;
    // Calculate the cost for all hits
    for (std::vector<TMS_Hit>::iterator jt = TMS_yz.begin(); jt != TMS_yz.end(); ++jt) {
      TMS_Hit comp = (*jt);
      // The difference in plane number
      int DeltaPlane = comp.GetPlaneNumber()-planeno;
      // Units of mm
      double DeltaY = fabs(comp.GetNotZ()-y);
      int DeltaY_Trunc = int(DeltaY/comp.GetNotZw());
      // Check only adjacent planes in z and y
      if (abs(DeltaPlane) > 2) continue;
      // Check distance in y isn't crazy (4 times the bar width of 4 cm)
      if (DeltaY_Trunc > 4) continue;

         std::cout << "y: " << y << std::endl;
         std::cout << "compy: " << comp.GetNotZ() << std::endl;
         std::cout << "DeltaY_Trunc: " << DeltaY_Trunc << std::endl;

      double cost = HighestCost;
      // Moving forwads
      if (DeltaPlane == 2 && DeltaY_Trunc == 0) cost = 1;
      // Moving upwards/downwards only
      else if (DeltaPlane == 0 && abs(DeltaY_Trunc) == 1) cost = 1;
      // Moving backwards only
      else if (DeltaPlane == -2 && DeltaY_Trunc == 0) cost = 1;
      // Moving diagonally (discourage this move over upward/downward + backward/forward
      else if (abs(DeltaPlane) == 2 && abs(DeltaY_Trunc) == 1) cost = 3;
      // Moving upwards/downwards by 2-4 steps: discourage vs taking 2-4 one-steps
      else if (DeltaPlane == 0 && abs(DeltaY_Trunc) > 1) cost = 1+abs(DeltaY_Trunc)*2;
      // Moving diagonally by N steps
      else if (abs(DeltaPlane) == 2 && abs(DeltaY_Trunc) > 1) cost = 1+abs(DeltaY_Trunc)*3;

      // If it's the same hit choose to include it, or merge these before? eek
      //else if (DeltaPlane == 0 && DeltaY_Trunc == 0) cost = 0;

      //std::cout << "   move to plane: " << comp.GetPlaneNumber() << " notz: " << comp.GetNotZ() << " cost: " << cost << std::endl;

      // Make the pair
      ThisCost.push_back(ThreePair(comp.GetPlaneNumber(), comp.GetNotZ(), cost));
    }

    PositionCost.push_back(ThisCost);
  }

  // Should now have our great vector of costs for each of the hits
  // Start the search
  //int StartPlane = TMS_yz.front().GetPlaneNumber();
  int FinishPlane = TMS_yz.back().GetPlaneNumber();
  double FinishNotZ = TMS_yz.back().GetNotZ();

  std::priority_queue<ThreePair, std::vector<ThreePair> > pq;
  pq.push(ThreePair(TMS_yz.front().GetPlaneNumber(), TMS_yz.front().GetNotZ(), 0));

  // Remember if we've visited each node
  std::vector<bool> Visited;
  for (size_t i = 0; i < PositionCost.size(); ++i) Visited.push_back(false);
  Visited[0] = true;

  while (!pq.empty()) {

    ThreePair bla = pq.top();
    int zplane = bla.xplane;
    double notz = bla.y;
    //std::cout << zplane << ", " << notz << ", " << bla.cost << std::endl;
    pq.pop();

    if (zplane == FinishPlane && notz == FinishNotZ) break;

    // Now loop over each of the connections
    int HitNumber = 0; 
    for (std::vector<std::vector<ThreePair> >::iterator it = PositionCost.begin(); it != PositionCost.end(); ++it, HitNumber++) {
      for (std::vector<ThreePair>::iterator jt = (*it).begin(); jt != (*it).end(); ++jt) {
        ThreePair curr = (*jt);
        //curr.Print();
        if (!Visited[HitNumber]) {
          pq.push(curr);
        }
      }
      // Remember that this has been visited already
      Visited[HitNumber] = true;
    }

  }
  */

  /*
     std::cout << "***********" << std::endl;
     for (std::vector<TMS_Hit>::iterator it = TMS_yz.begin(); it != TMS_yz.end(); ++it) {
     std::cout << (*it).GetBar().GetZ() << std::endl;
     }
     */

}
