// For drawing the event
#include "TCanvas.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TBox.h"

#include "TMS_Event.h"

class TMS_EventViewer {
  public:
    static TMS_EventViewer GetViewer() {
      static TMS_EventViewer Instance;
      return Instance;
    }

    void Draw(TMS_Event &event);
  private:
    TMS_EventViewer();

    TH2D *xz_view;
    TH2D *yz_view;
    TCanvas *Canvas;

    TBox *xz_box_FV;
    TBox *xz_box_Full;

    TBox *yz_box_FV;
    TBox *yz_box_Full;
};
