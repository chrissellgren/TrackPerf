#pragma once

#include <TH2.h>
#include <TH1.h>

//#include <ACTSTracking/GeometryIdMappingTool.hxx>
#include "/opt/ilcsoft/muonc/ACTSTracking/v1.1.0/ACTSTracking/GeometryIdMappingTool.hxx"


namespace EVENT
{
  class TrackerHit;
}

namespace TrackPerf
{
  //! Histograms for cluster analysis
  class ClusterHists
  {
  public:
    ClusterHists(const ClusterHists &) = delete ;
    ClusterHists& operator =(const ClusterHists &) = delete ;

    //! Initialize empty histograms
    ClusterHists() ;

    // Fill histograms with a single track hit
    void fill(const EVENT::TrackerHit* trkhit);

  private:
    TH2* h_size_theta;
    TH2* h_cluster_pos;
    TH2* h_cluster_pos_0;
    TH2* h_cluster_pos_1;
    TH2* h_cluster_pos_2;
    TH2* h_cluster_pos_3;
    TH1* hits_by_layer;
    TH1* h_theta;
    TH1* h_edep_0deg;
    TH1* h_edep_90deg;
  };
}
