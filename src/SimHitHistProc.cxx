#include "TrackPerf/SimHitHistProc.hxx"

#include "marlin/VerbosityLevels.h"

#include <set>

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <EVENT/SimTrackerHit.h>
#include <IMPL/TrackerHitPlaneImpl.h> // possibly unnecessary
#include <UTIL/LCTrackerConf.h> // possibly unnecessary


#include <marlin/AIDAProcessor.h>

#include <AIDA/ITree.h>

SimHitHistProc aSimHitHistProc ;

SimHitHistProc::SimHitHistProc() : Processor("SimHitHistProc")
{
  // modify processor description
  _description = "SimHitHistProc creates histograms mapping the location of Sim Tracker Hits." ;
  // register steering params: name, description, class-variable, default value
  // Sim Hits 
  registerInputCollection( LCIO::SIMTRACKERHIT,
			   "VBSimHitsCollection" , 
			   "Name of vertex barrel sim tracker hits collection",
			   _vbsimhitColName,
			   _vbsimhitColName
			   );     

  registerInputCollection( LCIO::SIMTRACKERHIT,
			   "IBSimHitsCollection" , 
			   "Name of inner barrel sim tracker hits collection",
			   _ibsimhitColName,
			   _ibsimhitColName
			   );   

  registerInputCollection( LCIO::SIMTRACKERHIT,
			   "OBSimHitsCollection" , 
			   "Name of outer barrel sim tracker hits collection",
			   _obsimhitColName,
			   _obsimhitColName
			   );  

  registerInputCollection( LCIO::SIMTRACKERHIT,
			   "VESimHitsCollection" , 
			   "Name of vertex endcap sim tracker hits collection",
			   _vesimhitColName,
			   _vesimhitColName
			   );  

  registerInputCollection( LCIO::SIMTRACKERHIT,
			   "IESimHitsCollection" , 
			   "Name of inner endcap sim tracker hits collection",
			   _iesimhitColName,
			   _iesimhitColName
			   );  

  registerInputCollection( LCIO::SIMTRACKERHIT,
			   "OESimHitsCollection" , 
			   "Name of outer endcap sim tracker hits collection",
			   _oesimhitColName,
			   _oesimhitColName
			   );  

  // Tracker Hits
  registerInputCollection( LCIO::TRACKERHIT,
			   "VBTrackerHitsCollection" , 
			   "Name of vertex barrel tracker hits collection",
			   _vbtrkhitColName,
			   _vbtrkhitColName
			   );     

  registerInputCollection( LCIO::TRACKERHIT,
			   "IBTrackerHitsCollection" , 
			   "Name of inner barrel tracker hits collection",
			   _ibtrkhitColName,
			   _ibtrkhitColName
			   );  

  registerInputCollection( LCIO::TRACKERHIT,
			   "OBTrackerHitsCollection" , 
			   "Name of outer barrel tracker hits collection",
			   _obtrkhitColName,
			   _obtrkhitColName
			   );  

  registerInputCollection( LCIO::TRACKERHIT,
			   "VETrackerHitsCollection" , 
			   "Name of vertex endcap tracker hits collection",
			   _vetrkhitColName,
			   _vetrkhitColName
			   );  

  registerInputCollection( LCIO::TRACKERHIT,
			   "IETrackerHitsCollection" , 
			   "Name of inner endcap tracker hits collection",
			   _ietrkhitColName,
			   _ietrkhitColName
			   );     

  registerInputCollection( LCIO::TRACKERHIT,
			   "OETrackerHitsCollection" , 
			   "Name of outer endcap tracker hits collection",
			   _oetrkhitColName,
			   _oetrkhitColName
			   );  
}

void SimHitHistProc::init()
{
  // Print the initial parameters
  printParameters() ;

  // Create histograms
  AIDA::ITree* tree=marlin::AIDAProcessor::tree(this);
  marlin::AIDAProcessor::histogramFactory(this);

  // Create ROOT histograms, with location setup by the above factory
  tree->mkdir("All" ); tree->cd("All" );
  int numbins_all = 1000;
  int rmin_all = 0;
  int rmax_all = 1600;
  int zmin_all = -2500;
  int zmax_all = 2500;
  h_x   = new TH1F("x  " , ";x   ; Num Hits" , numbins_all, -rmax_all, rmax_all);
  h_y   = new TH1F("y  " , ";y   ; Num Hits" , numbins_all, -rmax_all, rmax_all);
  h_z   = new TH1F("z  " , ";z   ; Num Hits" , numbins_all,  zmin_all, zmax_all);
  h_r   = new TH1F("r  " , ";r   ; Num Hits" , numbins_all,  rmin_all, rmax_all);
  h_t   = new TH1F("t  " , "Sim Hits from VXB and ITB;t   ; Num Hits" , numbins_all,  -2, 10); // ns
  h_z_r = new TH2F("z_r" , ";z_r ; r"        , numbins_all,  zmin_all, zmax_all, numbins_all, rmin_all, rmax_all);
  h_x_y = new TH2F("x_y" , ";x_y ; r"        , numbins_all, -rmax_all, rmax_all, numbins_all, -rmax_all, rmax_all);

  // vertex histograms
  tree->mkdir("../Vertex"); tree->cd("../Vertex");
  int numbins_vx = 500;
  int rmin_vx = 0;
  int rmax_vx = 120;
  int zmin_vx = -300;
  int zmax_vx = 300;
  h_x_vx   = new TH1F("x_vx  " , ";x   ; Num Hits" , numbins_vx, -rmax_vx, rmax_vx);
  h_y_vx   = new TH1F("y_vx  " , ";y   ; Num Hits" , numbins_vx, -rmax_vx, rmax_vx);
  h_z_vx   = new TH1F("z_vx  " , ";z   ; Num Hits" , numbins_vx,  zmin_vx, zmax_vx);
  h_r_vx   = new TH1F("r_vx  " , ";r   ; Num Hits" , numbins_vx,  rmin_vx, rmax_vx);
  h_t_vx   = new TH1F("t  " , ";t   ; Num Hits" , numbins_all,  -2, 10); // ns
  h_z_r_vx = new TH2F("z_r_vx" , ";z_r ; r"        , numbins_vx,  zmin_vx, zmax_vx, numbins_vx, rmin_vx, rmax_vx);
  h_x_y_vx = new TH2F("x_y_vx" , ";x_y ; r"        , numbins_vx, -rmax_vx, rmax_vx, numbins_vx, -rmax_vx, rmax_vx);
  h_t_tracker_vxb = new TH1F("t  " , ";t   ; Num Hits" , numbins_all,  -2, 10); // ns

  // IT histograms
  tree->mkdir("../IT"); tree->cd("../IT");
  int numbins_IT = 1000;
  int rmin_IT = 120;
  int rmax_IT = 600;
  int zmin_IT = 0;
  int zmax_IT = 1000;
  h_x_it   = new TH1F("x_vx  " , ";x   ; Num Hits" , numbins_IT, -rmax_IT, rmax_IT);
  h_y_it   = new TH1F("y_vx  " , ";y   ; Num Hits" , numbins_IT, -rmax_IT, rmax_IT);
  h_z_it   = new TH1F("z_vx  " , ";z   ; Num Hits" , numbins_IT,  zmin_IT, zmax_IT);
  h_r_it   = new TH1F("r_vx  " , ";r   ; Num Hits" , numbins_IT,  rmin_IT, rmax_IT);
  h_t_it   = new TH1F("t  " , ";t   ; Num Hits" , numbins_all,  -2, 10); // ns
  h_z_r_it = new TH2F("z_r_vx" , ";z_r ; r"        , numbins_IT,  zmin_IT, zmax_IT, numbins_IT, rmin_IT, rmax_IT);
  h_x_y_it = new TH2F("x_y_vx" , ";x_y ; r"        , numbins_IT, -rmax_IT, rmax_IT, numbins_IT, -rmax_IT, rmax_IT);

  // OT histograms
  tree->mkdir("../OT"); tree->cd("../OT");
  int numbins_OT = 500;
  int rmin_OT = 800;
  int rmax_OT = 1600;
  int zmin_OT = -2500;
  int zmax_OT = 2500;
  h_x_ot   = new TH1F("x_vx  " , ";x   ; Num Hits" , numbins_OT, -rmax_OT, rmax_OT);
  h_y_ot   = new TH1F("y_vx  " , ";y   ; Num Hits" , numbins_OT, -rmax_OT, rmax_OT);
  h_z_ot   = new TH1F("z_vx  " , ";z   ; Num Hits" , numbins_OT,  zmin_OT, zmax_OT);
  h_r_ot   = new TH1F("r_vx  " , ";r   ; Num Hits" , numbins_OT,  rmin_OT, rmax_OT);
  h_t_ot   = new TH1F("t  " , ";t   ; Num Hits" , numbins_all,  -2, 10); // ns
  h_z_r_ot = new TH2F("z_r_vx" , ";z_r ; r"        , numbins_OT,  zmin_OT, zmax_OT, numbins_OT, rmin_OT, rmax_OT);
  h_x_y_ot = new TH2F("x_y_vx" , ";x_y ; r"        , numbins_OT, -rmax_OT, rmax_OT, numbins_OT, -rmax_OT, rmax_OT);
  }

void SimHitHistProc::processRunHeader( LCRunHeader* /*run*/)
{ } 

void SimHitHistProc::processEvent (LCEvent * evt)
{
  // Get Sim Tracker Hits
  LCCollection* vbsimhitCol  =evt->getCollection(_vbsimhitColName);
  LCCollection* ibsimhitCol  =evt->getCollection(_ibsimhitColName);
  LCCollection* obsimhitCol  =evt->getCollection(_obsimhitColName);
  LCCollection* vesimhitCol  =evt->getCollection(_vesimhitColName);
  LCCollection* iesimhitCol  =evt->getCollection(_iesimhitColName);
  LCCollection* oesimhitCol  =evt->getCollection(_oesimhitColName);

  // Loop over hits from each collection, fill the histograms with the data from that collection
  for(int i=0; i<vbsimhitCol->getNumberOfElements(); ++i)
    {
      const EVENT::SimTrackerHit *simhit=static_cast<const EVENT::SimTrackerHit*>(vbsimhitCol->getElementAt(i));
      SimHitHistProc::fill(simhit, "vertexbarrel");}
  for(int i=0; i<ibsimhitCol->getNumberOfElements(); ++i)
    {
      const EVENT::SimTrackerHit *simhit=static_cast<const EVENT::SimTrackerHit*>(ibsimhitCol->getElementAt(i));
      SimHitHistProc::fill(simhit, "innerbarrel");}  
  for(int i=0; i<obsimhitCol->getNumberOfElements(); ++i)
    {
      const EVENT::SimTrackerHit *simhit=static_cast<const EVENT::SimTrackerHit*>(obsimhitCol->getElementAt(i));
      SimHitHistProc::fill(simhit, "outerbarrel");}
  for(int i=0; i<vesimhitCol->getNumberOfElements(); ++i)
    {
      const EVENT::SimTrackerHit *simhit=static_cast<const EVENT::SimTrackerHit*>(vesimhitCol->getElementAt(i));
      SimHitHistProc::fill(simhit, "vertexcap");}
  for(int i=0; i<iesimhitCol->getNumberOfElements(); ++i)
    {
      const EVENT::SimTrackerHit *simhit=static_cast<const EVENT::SimTrackerHit*>(iesimhitCol->getElementAt(i));
      SimHitHistProc::fill(simhit, "innercap");}
  for(int i=0; i<oesimhitCol->getNumberOfElements(); ++i)
    {
      const EVENT::SimTrackerHit *simhit=static_cast<const EVENT::SimTrackerHit*>(oesimhitCol->getElementAt(i));
      SimHitHistProc::fill(simhit, "outercap");}

  LCCollection* vbtrkhitCol  =evt->getCollection(_vbtrkhitColName);
  LCCollection* ibtrkhitCol  =evt->getCollection(_ibtrkhitColName);
  LCCollection* obtrkhitCol  =evt->getCollection(_obtrkhitColName);
  LCCollection* vetrkhitCol  =evt->getCollection(_vetrkhitColName);
  LCCollection* ietrkhitCol  =evt->getCollection(_ietrkhitColName);
  LCCollection* oetrkhitCol  =evt->getCollection(_oetrkhitColName);

  for(int i=0; i<vbtrkhitCol->getNumberOfElements(); ++i)
    {
      const EVENT::TrackerHit *trkhit=static_cast<const EVENT::TrackerHit*>(vbtrkhitCol->getElementAt(i));
      float t = trkhit->getTime();
      h_t_tracker_vxb->Fill(t);
      }
}

void SimHitHistProc::check( LCEvent * /*evt*/ )
{ }

void SimHitHistProc::end()
{ }


void SimHitHistProc::fill(const EVENT::SimTrackerHit* simhit, const std::string flag)
{
  //Get position
  float x = simhit->getPosition()[0];
  float y = simhit->getPosition()[1];
  float z = simhit->getPosition()[2];
  float r = sqrt(pow(x,2)+pow(y,2));
  float t = simhit->getTime();
  streamlog_out(DEBUG9) << "Sim Hit Parameters: x=" << x << ", y=" << y << ", z=" << z << ", r=" << r << ", t=" << t << std::endl;

  // Fill histograms with all hits
  h_x->Fill(x);  
  h_y->Fill(y);  
  h_z->Fill(z);  
  h_r->Fill(r);
  h_z_r->Fill(z,r);
  h_x_y->Fill(x,y);

  // Fill subdetector-specific histograms depending on input flag
  if (flag=="vertexbarrel" || flag=="vertexcap"){
    h_x_vx->Fill(x);  
    h_y_vx->Fill(y);  
    h_z_vx->Fill(z);  
    h_r_vx->Fill(r);  
    h_z_r_vx->Fill(z,r);
    h_x_y_vx->Fill(x,y);
  }

  if (flag=="innerbarrel" || flag=="innercap"){
    h_x_it->Fill(x);  
    h_y_it->Fill(y);  
    h_z_it->Fill(z);  
    h_r_it->Fill(r);
    h_z_r_it->Fill(z,r);
    h_x_y_it->Fill(x,y);
  }

  if (flag=="outerbarrel" || flag=="outercap"){
    h_x_ot->Fill(x);  
    h_y_ot->Fill(y);  
    h_z_ot->Fill(z);  
    h_r_ot->Fill(r);
    h_t_ot->Fill(t);  
    h_z_r_ot->Fill(z,r);
    h_x_y_ot->Fill(x,y);
  }
    // restrict vx histograms to only barrel
    if (flag=="vertexbarrel"){
        h_t_vx->Fill(t);
        h_t->Fill(t);
    }
    if (flag=="innerbarrel"){
        h_t_it->Fill(t); 
        h_t->Fill(t); 
    }
}

