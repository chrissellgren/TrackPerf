#include "TrackPerf/ClusterHists.hxx"
#include "marlin/VerbosityLevels.h"

#include <EVENT/TrackerHit.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <UTIL/LCTrackerConf.h>


using namespace TrackPerf;

ClusterHists::ClusterHists()
{
  h_size_theta    = new TH2F("cluster_size_vs_theta" , ";Cluster #theta; Cluster size" , 100, 0,  3.14,  10,  -0.5,  10.5  );
  h_cluster_pos   = new TH2F("cluster_position"      , ";z; r"                         , 100, -500, 500, 100, 0, 200);
  h_cluster_pos_0 = new TH2F("cluster_position_0"    , ";z; r"                         , 100, -500, 500, 100, 0, 200);
  h_cluster_pos_1 = new TH2F("cluster_position_1"    , ";z; r"                         , 100, -500, 500, 100, 0, 200);
  h_cluster_pos_2 = new TH2F("cluster_position_2"    , ";z; r"                         , 100, -500, 500, 100, 0, 200);
  h_cluster_pos_3 = new TH2F("cluster_position_3"    , ";z; r"                         , 100, -500, 500, 100, 0, 200);
  hits_by_layer   = new TH1F("numhits_by_layer"      , ";Layer Index; Number of Clusters",8,0,8);
  h_theta         = new TH1F("theta"                 , ";Theta;Number of Clusters"       ,100,0,3.15);
  h_edep_0deg     = new TH1F("edep_0to5deg"          , ";Energy Deposited (GeV);Clusters" ,100,0,0.0005);
  h_edep_90deg    = new TH1F("edep_89to91deg"        , ";Energy Deposited (GeV);Clusters",100,0,0.0005);
}

void ClusterHists::fill(const EVENT::TrackerHit* trkhit)
{
  //Calculate energy deposited
  float EDep = trkhit->getEDep();

  //Calculating theta
  float x = trkhit->getPosition()[0];
  float y = trkhit->getPosition()[1];
  float z = trkhit->getPosition()[2];
  float r = sqrt(pow(x,2)+pow(y,2));
  float incidentTheta = std::atan(r/z);
  streamlog_out(DEBUG9) << "theta before negative correction: " << incidentTheta << std::endl;

  if(incidentTheta<0)
    incidentTheta += M_PI;
  streamlog_out(DEBUG9) << "the value of theta is " << incidentTheta << std::endl;


  //Calculating cluster size
  const lcio::LCObjectVec &rawHits = trkhit->getRawHits(); 
  float max = -1000000;
  float min = 1000000; 

  float loopsize = rawHits.size();
  streamlog_out(DEBUG9) << "the size of rawhits is " << loopsize << std::endl;

  for (size_t j=0; j<loopsize; ++j) {
    streamlog_out(DEBUG9) << "the for loop is on iteration" << j << std::endl;
    lcio::SimTrackerHit *hitConstituent = dynamic_cast<lcio::SimTrackerHit*>( rawHits[j] );
    const double *localPos = hitConstituent->getPosition();
    float x_local = localPos[0];
    float y_local = localPos[1];
    streamlog_out(DEBUG9) << "y_local is: " << y_local << std::endl;
    if (y_local < min){
      min = y_local;
      }
    else if (y_local > max){
      max = y_local;          
      } 
    }
  streamlog_out(DEBUG9) << "the value of min and max are: " << min  << " and " << max << std::endl;
  float cluster_size = (max - min)+1;
  streamlog_out(DEBUG9) << "the value of cluster size is " << cluster_size << std::endl;

  //Get hit subdetector/layer 
  std::string _encoderString = lcio::LCTrackerCellID::encoding_string();
  UTIL::CellIDDecoder<lcio::TrackerHit> decoder(_encoderString);
  uint32_t systemID = decoder(trkhit)["system"];
  uint32_t layerID = decoder(trkhit)["layer"];
  // shift layer by 0.5 to resolve binning issue
  double layerID_adjusted = layerID + 0.5;

  // Fill for all hits
  h_size_theta->Fill(incidentTheta, cluster_size);
  h_theta->Fill(incidentTheta);
  h_cluster_pos->Fill(z,r);
  hits_by_layer->Fill(layerID_adjusted);

  // Fill energy deposition histograms based on angle
  float theta_deg = incidentTheta * (180/3.1416);
  if(theta_deg < 5 || theta_deg > 175){
    h_edep_0deg->Fill(EDep);
  }
  if(theta_deg > 89 && theta_deg < 91){
    h_edep_90deg->Fill(EDep);
  }
  
  // Fill based on which double layer region was hit
  if(layerID==0 or layerID==1){
    h_cluster_pos_0->Fill(z,r);
    }
  if(layerID==2 or layerID==3){
    h_cluster_pos_1->Fill(z,r);
    }
  if(layerID==4 or layerID==5){
    h_cluster_pos_2->Fill(z,r);
    }
  if(layerID==6 or layerID==7){
    h_cluster_pos_3->Fill(z,r);
    }   
}