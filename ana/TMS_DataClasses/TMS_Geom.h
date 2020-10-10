#ifndef _TMS_GEOM_H_SEEN_
#define _TMS_GEOM_H_SEEN_

#include <iostream>

#include "TGeoManager.h"

// Define the TMS geometry singleton
class TMS_Geom {
  public:

    // Get the instance
    static TMS_Geom& GetInstance() {
      static TMS_Geom Instance;
      return Instance;
    }

    // Get the geometry
    TGeoManager* GetGeometry() {
      // Shout if the geometry isn't set
      if (this->geom == NULL) {
        std::cerr << "Geometry not set!" << std::endl;
        std::cerr << "Call TMS_Geom.SetGeometry(TGeoManager *) first!" << std::endl;
        throw;
      }
      return this->geom;
    }

    // Set the geometry
    void SetGeometry(TGeoManager *geometry) {
      geom = geometry;
      std::cout << "Global geometry set to " << geometry->GetName() << std::endl;
    }


    TMS_Geom(TMS_Geom const &) = delete;
    void operator=(TMS_Geom const &) = delete;

  private:
    // The empty constructor
    TMS_Geom() {
      geom = NULL;
    };
    //TMS_Geom(TMS_Geom const &);
    //void operator=(TMS_Geom const &);

    // The actual geometry
    TGeoManager *geom;
};

#endif
