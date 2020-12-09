#ifndef _TMS_KALMAN_H_SEEN_
#define _TMS_KALMAN_H_SEEN_

#include <iostream>

#include "TMatrixD.h"

// Include the physics for propgataion
#include "BetheBloch.h"
#include "MultipleScattering.h"
#include "Material.h"

// And the geometry for the materials
#include "TMS_Geom.h"
// Need to understand what a hit is
#include "TMS_Hit.h"


// State vector for Kalman filter
// Not actually ever truly observed (not a measurement!)
class TMS_State {
  public:
    double x;
    double y;
    double dxdz;
    double dydz;
    double qp;
    double z; // the dependent variable of the state vector

};

class TMS_Kalman {
  public:
    TMS_Kalman();
    TMS_Kalman(std::vector<TMS_Hit> &Candidates);
    void Update();
    void Predict(TMS_Hit &hit1, TMS_Hit &hit2);

  private:
    // Energy-loss calculator
    BetheBloch_Calculator Bethe;
    MultipleScatter_Calculator MSC;

    // State vector
    // x, y, dx/dz, dy/dz, q/p
    TMS_State State;

    // Propagator matrix
    // Takes us from detector k-1 to detector k
    TMatrixD PropMatrix;

    // Random variable w(k-1) includes random disturbances of track between z(k-1) and z(k) from multiple scattering
    // Noise matrix
    TMatrixD NoiseMatrix;

    // Measurement matrix
    TMatrixD MeasurementMatrix;

    // Remember if we're forward fitting
    bool ForwardFitting;

    double en;
    double mass;
    double mom;
    double total_en_var;
};

#endif
