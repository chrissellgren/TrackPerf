#pragma once

#include <TH1.h>

namespace EVENT
{
  class TrackerHit;
  class SimTrackerHit;
}

namespace IMPL
{
  class TrackerHitPlaneImpl;
}

namespace TrackPerf
{
  //! Histograms for tracker hit analysis
  class TrackerHitResoHists
  {
  public:
    TrackerHitResoHists(const TrackerHitResoHists &) = delete ;
    TrackerHitResoHists& operator =(const TrackerHitResoHists &) = delete ;

    //! Initialize empty histograms
    TrackerHitResoHists() ;

    // Fill histograms 
    void fill(const EVENT::TrackerHit* trkhit, const EVENT::SimTrackerHit* simtrkhit, IMPL::TrackerHitPlaneImpl* trkhitplane);

  private:
    TH1* h_x_pull;
    TH1* h_y_pull;
    TH1* h_dx_wide;
    TH1* h_dy_wide;
    TH1* h_dz_wide;
    TH1* h_dx;
    TH1* h_dy;
    TH1* h_dz;
    TH1* h_dr;
    TH1* h_cov_x;
    TH1* h_cov_y;
    TH1* h_cov_r;
  };
}
