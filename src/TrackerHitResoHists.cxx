#include "TrackPerf/TrackerHitResoHists.hxx"
#include "marlin/VerbosityLevels.h"

#include <EVENT/TrackerHit.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitPlaneImpl.h>


using namespace TrackPerf;

TrackerHitResoHists::TrackerHitResoHists()
{
  h_x_pull = new TH1F("x_pull", ";(x_{reco} - x_{truth})/#sigma_{x}; Hits" , 50, -2, 2); 
  h_y_pull = new TH1F("y_pull", ";(y_{reco} - y_{truth})/#sigma_{y}; Hits" , 50, -2, 2); 
  // large range histograms for debugging bad reconstruction
  h_dx_wide = new TH1F("dx_wide", ";(x_{reco} - x_{truth}) (mm); Hits" , 100, -20, 20); 
  h_dy_wide = new TH1F("dy_wide", ";(y_{reco} - y_{truth}) (mm); Hits" , 100, -20, 20);
  h_dz_wide = new TH1F("dz_wide", ";(z_{reco} - z_{truth}) (mm); Hits" , 100, -20, 20);
  // absolute position difference
  h_dx = new TH1F("dx", ";(x_{reco} - x_{truth}) (mm); Hits" , 100, -0.03, 0.03); 
  h_dy = new TH1F("dy", ";(y_{reco} - y_{truth}) (mm); Hits" , 100, -0.03, 0.03);
  h_dz = new TH1F("dz", ";(z_{reco} - z_{truth}) (mm); Hits" , 100, -0.03, 0.03);
  h_dr = new TH1F("dr", ";(r_{reco} - r_{truth}) (mm); Hits" , 100, -20, 20); // r = sqrt(dx^2+dy^2+dz^2)
  // uncertainties
  h_cov_x  = new TH1F("cov_x" , ";X variance; Hits"                        , 50, 0, 0.01);
  h_cov_y  = new TH1F("cov_y" , ";Y variance; Hits"                        , 50, 0, 0.01);
  h_cov_r  = new TH1F("cov_r" , ";r variance; Hits"                        , 50, 0, 0.01);
}

void TrackerHitResoHists::fill(const EVENT::TrackerHit* trkhit, const EVENT::SimTrackerHit* simtrkhit, IMPL::TrackerHitPlaneImpl* trkhitplane)
{
  //Calculating theta
  float x = trkhit->getPosition()[0];
  float y = trkhit->getPosition()[1];
  float z = trkhit->getPosition()[2];
  float r = sqrt(pow(x,2)+pow(y,2));
  float incidentTheta = std::atan(r/z);
  if(incidentTheta<0)
    incidentTheta += M_PI;

  float x_reco = trkhit->getPosition()[0];
  float y_reco = trkhit->getPosition()[1];
  float z_reco = trkhit->getPosition()[2];
        
  float x_truth = simtrkhit->getPosition()[0];
  float y_truth = simtrkhit->getPosition()[1];
  float z_truth = simtrkhit->getPosition()[2];

  if(x_reco-x_truth==0){
    streamlog_out(DEBUG3) << "x truth = x_reco" << std::endl;
    streamlog_out(DEBUG3) << "x truth: " << x_truth << std::endl;
    streamlog_out(DEBUG3) << "y truth: " << y_truth << std::endl;
    streamlog_out(DEBUG3) << "yreco-ytruth: " << y_reco-y_truth << std::endl;
    streamlog_out(DEBUG3) << "incident theta: " << incidentTheta << std::endl << std::endl;
  }

  float sigma_x = trkhitplane->getdU();
  float sigma_y = trkhitplane->getdV();

  float delta_r = sqrt(pow(x_reco-x_truth,2)+pow(y_reco-y_truth,2)); // distance along xy plane
  float sigma_r = 1/delta_r * ((x_reco-x_truth)*sigma_x + (y_reco-y_truth)*sigma_y); // from error propagation

  h_x_pull->Fill((x_reco-x_truth)/sigma_x);
  h_y_pull->Fill((y_reco-y_truth)/sigma_y);

  h_dx_wide->Fill((x_reco-x_truth));
  h_dy_wide->Fill((y_reco-y_truth));
  h_dz_wide->Fill((z_reco-z_truth));

  h_dx->Fill((x_reco-x_truth));
  h_dy->Fill((y_reco-y_truth));
  h_dz->Fill((z_reco-z_truth));
  h_dr->Fill(delta_r);
  
  h_cov_x->Fill(sigma_x);
  h_cov_y->Fill(sigma_y);
  h_cov_r->Fill(sigma_r);
}

