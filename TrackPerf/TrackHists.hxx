#pragma once

#include <TH1.h>
#include <TH2.h>

//#include <ACTSTracking/GeometryIdMappingTool.hxx>
#include "/opt/ilcsoft/muonc/ACTSTracking/v1.1.0/ACTSTracking/GeometryIdMappingTool.hxx"

namespace EVENT
{
  class Track;
}

namespace TrackPerf
{
  //! Histograms for reconstructed tracks
  class TrackHists
  {
  public:
    TrackHists(const TrackHists &) = delete ;
    TrackHists& operator =(const TrackHists &) = delete ;

    //! Initialize empty histograms
    TrackHists() ;

    // Fill histograms with a single track
    void fill(const EVENT::Track* track);


  private:
    //! magnetic field to use for curvature -> pT conversion
    float _Bz=3.57;

    //! Reconstructed track pT
    TH1* h_pt;
    TH1* h_lambda;
    TH1* h_phi;
    TH1* h_d0;
    TH1* h_z0;
    TH1* h_nhit;
    TH1* h_chi2_spatial;
    TH1* h_chi2_temp_avg;
    TH1* h_chi2_temp_vtx;
    TH1* h_chi2_temp_inner;
    TH1* h_chi2_temp_outer;
    TH1* h_chi2_temp_max;
    TH1* h_chi2_temp_std;
    TH1* h_cov_z0;
    TH2* h_lambda_nhit;
    TH2* h_pt_nhit;
    TH2* h_pt_lambda;
    TH2* h_tempchi2_nhit;
    TH1* h_nhit_vtx;
    TH1* h_nhit_inner;
    TH1* h_nhit_outer;
    TH1* h_nhit4;
    TH1* h_nhit5;
    TH1* h_nhit6;

    TH1* h_cov_d0_d0;         
    TH1* h_cov_d0_phi;       
    TH1* h_cov_d0_omega;      
    TH1* h_cov_d0_z0;         
    TH1* h_cov_d0_lambda;     
    TH1* h_cov_phi_phi;       
    TH1* h_cov_phi_omega;     
    TH1* h_cov_phi_z0;        
    TH1* h_cov_phi_lambda;    
    TH1* h_cov_omega_omega;   
    TH1* h_cov_omega_z0;      
    TH1* h_cov_omega_lambda;  
    TH1* h_cov_z0_z0;        
    TH1* h_cov_z0_lambda;     
    TH1* h_cov_lambda_lambda; 

    TH1* h_low_pt_lambda;

    TH1* h_nhit_24_2;
    TH1* h_nhit_24_4;
    TH1* h_nhit_24_6;
  };
}
