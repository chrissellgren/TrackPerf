#include "TrackPerf/TrackPerfHistProc.hxx"

#include "TrackPerf/TrackHists.hxx"
#include "TrackPerf/TruthHists.hxx"
#include "TrackPerf/TrackResoHists.hxx"
#include "TrackPerf/TrackerHitResoHists.hxx"
#include "TrackPerf/ClusterHists.hxx"

#include <set>

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCRelation.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHitPlane.h>
#include <IMPL/TrackerHitPlaneImpl.h>

#include <marlin/AIDAProcessor.h>

#include <AIDA/ITree.h>

TrackPerfHistProc aTrackPerfHistProc ;

TrackPerfHistProc::TrackPerfHistProc()
  : Processor("TrackPerfHistProc")
{  
  // modify processor description
  _description = "TrackPerfHistProc creates a series of output histograms for tracking performance studies." ;

  // register steering parameters: name, description, class-variable, default value
  registerProcessorParameter("MatchProb",
                             "Minimum matching probabilty to be considered a good track-mc match.",
                             _matchProb,
                             _matchProb);

  registerInputCollection( LCIO::MCPARTICLE,
			   "MCParticleCollection" , 
			   "Name of the MCParticle collection"  ,
			   _mcpColName,
			   _mcpColName
			   );

  registerInputCollection( LCIO::TRACK,
			   "TrackCollection" ,
			   "Name of the Track collection" ,
			   _trkColName,
			   _trkColName
			   );

  registerInputCollection( LCIO::LCRELATION,
			   "MCTrackRelationCollection" , 
			   "Name of LCRelation collection with track to MC matching",
			   _trkMatchColName,
			   _trkMatchColName
			   );

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

  registerInputCollection( LCIO::LCRELATION,
		     "VBRelationCollection" ,
			   "Name of the input relation collection",
			   _VBRelationCollection,
		     _VBRelationCollection
		 	    );
}

void TrackPerfHistProc::init()
{
  // Print the initial parameters
  printParameters() ;

  // Create histograms
  AIDA::ITree* tree=marlin::AIDAProcessor::tree(this);
  marlin::AIDAProcessor::histogramFactory(this);

  // Create ROOT histograms, with location setup by the above factory
  tree->mkdir("all" ); tree->cd("all" );
  _allTracks =std::make_shared<TrackPerf::TrackHists>();
  _allTruths =std::make_shared<TrackPerf::TruthHists>();
  h_number_of_tracks = new TH1F("number_of_tracks", ";Number of seeds/tracks;Events", 100, 0, 300000);
  h_trackerhit_timing = new TH1F("time_of_flight", ";Time of flight for tracker hits [ns];Tracker hits;", 1000, -1, 1);
  tree->mkdir("../real"); tree->cd("../real");
  _realTracks=std::make_shared<TrackPerf::TrackHists>();
  _realTruths=std::make_shared<TrackPerf::TruthHists>();
  _realReso    =std::make_shared<TrackPerf::TrackResoHists>();
  h_relation_weight_real = new TH1F("relation_weight", ";Track-truth relation weight;", 100, 0, 2);
  tree->mkdir("../fake"); tree->cd("../fake");
  _fakeTracks=std::make_shared<TrackPerf::TrackHists>();
  h_number_of_fakes = new TH1F("number_of_fakes", ";Number of fake tracks;Events", 100,  0, 300000);
  h_relation_weight_fake = new TH1F("relation_weight", ";Track-truth relation weight;Tracks;", 100, 0, 2);
  tree->mkdir("../unmt"); tree->cd("../unmt");
  _unmtTruths=std::make_shared<TrackPerf::TruthHists>();

  tree->mkdir("../clusters_vb" ); tree->cd("../clusters_vb" );
  _uncertainties=std::make_shared<TrackPerf::TrackerHitResoHists>();
  _clusters_vb=std::make_shared<TrackPerf::ClusterHists>();

  tree->mkdir("../clusters_ib" ); tree->cd("../clusters_ib" );
  _uncertainties=std::make_shared<TrackPerf::TrackerHitResoHists>();
  _clusters_ib=std::make_shared<TrackPerf::ClusterHists>();

}

void TrackPerfHistProc::processRunHeader( LCRunHeader* /*run*/)
{ } 

void TrackPerfHistProc::processEventTracks( LCEvent * evt)
{
    // MCParticles

  LCCollection* mcpCol  =evt->getCollection(_mcpColName);

  if( mcpCol->getTypeName() != lcio::LCIO::MCPARTICLE )
    { throw EVENT::Exception( "Invalid collection type: " + mcpCol->getTypeName() ) ; }

  std::set<const EVENT::MCParticle*> mcpSet;
  for(uint32_t i=0;i<mcpCol->getNumberOfElements();i++)
    {
      const EVENT::MCParticle *mcp=static_cast<const EVENT::MCParticle*>(mcpCol->getElementAt(i));

      if(mcp->getGeneratorStatus()!=1)
	{ continue; }

      if(mcp->getCharge()==0)
	{ continue; }

      if(mcp->isDecayedInTracker())
	{ continue; }

      // Tracker acceptance
      const double* mom=mcp->getMomentum();
      double pt=std::sqrt(std::pow(mom[0],2)+std::pow(mom[1],2));
      double lambda=std::atan2(mom[2],pt);
      //if(fabs(lambda)>75./180*3.14)
      if(fabs(lambda)>0.8)
	{ continue; }

      mcpSet.insert(mcp);
      _allTruths->fill(mcp);
    }
    
    // Tracks

  LCCollection* trkCol  =evt->getCollection(_trkColName);
  if( trkCol == nullptr){
    std::cout << "Track collection pointer is a nullptr" << std::endl;
  }
  std::set<const EVENT::Track*> trkSet;

  if (trkCol != nullptr) {
    if( trkCol->getTypeName() != lcio::LCIO::TRACK )
      { throw EVENT::Exception( "Invalid collection type: " + trkCol->getTypeName() ) ; }
    
    for(uint32_t i=0;i<trkCol->getNumberOfElements();i++)
      {
        const EVENT::Track *trk=static_cast<const EVENT::Track*>(trkCol->getElementAt(i));
        
        float lambda=std::atan(trk->getTanLambda());
        if(fabs(lambda)>0.8)
        {continue;}

        trkSet.insert(trk);
        _allTracks->fill(trk);
      }
      h_number_of_tracks->Fill(trkSet.size());
  }

  // Relations

  LCCollection* relCol  =evt->getCollection(_trkMatchColName);
  std::set<const EVENT::LCRelation*> relSet;

  if( relCol != nullptr){
    if( relCol->getTypeName() != lcio::LCIO::LCRELATION )
      { throw EVENT::Exception( "Invalid collection type: " + relCol->getTypeName() ) ; }


    for(uint32_t i=0;i<relCol->getNumberOfElements();i++)
      {
        const EVENT::LCRelation *rel=static_cast<const EVENT::LCRelation*>(relCol->getElementAt(i));

        relSet.insert(rel);
      }
  }

  //
  // Loop over track to MC associations to save matched objects
  LCCollection* tr2mcCol=evt->getCollection(_trkMatchColName);
  if( trkCol != nullptr && relCol != nullptr && tr2mcCol != nullptr){
    if( tr2mcCol->getTypeName() != lcio::LCIO::LCRELATION )
      { throw EVENT::Exception( "Invalid collection type: "+ tr2mcCol->getTypeName() ); }

    for(int i=0; i<tr2mcCol->getNumberOfElements(); ++i)
      {
        const EVENT::LCRelation *rel=static_cast<const EVENT::LCRelation*>(tr2mcCol->getElementAt(i));
        const EVENT::MCParticle *mcp=static_cast<const EVENT::MCParticle*>(rel->getFrom());
        const EVENT::Track      *trk=static_cast<const EVENT::Track     *>(rel->getTo  ());

        if(mcpSet.count(mcp)==0)
    { continue; } // truth particle not selected

        if(rel->getWeight()>_matchProb)
    {  
      if(trkSet.find(trk) != trkSet.end())
      {
        _realTracks->fill(trk);
        _realTruths->fill(mcp);
        _realReso  ->fill(trk,mcp);
        h_relation_weight_real ->Fill(rel->getWeight());

        mcpSet.erase(mcp);
        trkSet.erase(trk);
        relSet.erase(rel);
      }
    }
      }  
    //
    // Save unmatched objects
    for(const EVENT::MCParticle* mcp : mcpSet)
      { _unmtTruths->fill(mcp); }
    for(const EVENT::Track* trk : trkSet)
      { _fakeTracks->fill(trk);} 
    for(const EVENT::LCRelation* rel : relSet)
      { h_relation_weight_fake->Fill(rel->getWeight());}
    h_number_of_fakes->Fill(trkSet.size());
}
}

void TrackPerfHistProc::processEventTrackerHits( LCEvent * evt)
{
  //Tracker hits

  LCCollection* vbtrkhitCol  =evt->getCollection(_vbtrkhitColName);
  LCCollection* ibtrkhitCol  =evt->getCollection(_ibtrkhitColName);
  LCCollection* obtrkhitCol  =evt->getCollection(_obtrkhitColName);
  LCCollection* vetrkhitCol  =evt->getCollection(_vetrkhitColName);
  LCCollection* ietrkhitCol  =evt->getCollection(_ietrkhitColName);
  LCCollection* oetrkhitCol  =evt->getCollection(_oetrkhitColName);
  LCCollection* VBRelationCollection =evt->getCollection(_VBRelationCollection);

  for(int i=0; i<VBRelationCollection->getNumberOfElements(); ++i)
    {
      streamlog_out(DEBUG3) << "Events in VB Relation Collection: " << VBRelationCollection->getNumberOfElements() << std::endl;
      EVENT::LCRelation *rel=static_cast<EVENT::LCRelation*>(VBRelationCollection->getElementAt(i));
      EVENT::TrackerHit *trkhit=dynamic_cast<EVENT::TrackerHit*>(rel->getFrom());
      EVENT::SimTrackerHit *simtrkhit=dynamic_cast<EVENT::SimTrackerHit*>(rel->getTo());

      IMPL::TrackerHitPlaneImpl *trkhitplane=dynamic_cast<IMPL::TrackerHitPlaneImpl*>(trkhit);

      if(trkhit==nullptr or simtrkhit==nullptr or trkhitplane==nullptr){
        std::cout << "Warning: Failed to dynamic cast to planar sensor" << std::endl;
        std::cout << "- Trackhit: " << trkhit << std::endl;
        std::cout << "- Simtrackhit: " << simtrkhit << std::endl;
        std::cout << "- Trackhitplane: " << trkhitplane << std::endl;
        continue;
      }

      _uncertainties->fill(trkhit,simtrkhit,trkhitplane);
    }

  streamlog_out(DEBUG3) << "Num Events in VB Hit Collection: " << vbtrkhitCol->getNumberOfElements() << std::endl;
  for(int i=0; i<vbtrkhitCol->getNumberOfElements(); ++i)
    {
      const EVENT::TrackerHit *trkhit=static_cast<const EVENT::TrackerHit*>(vbtrkhitCol->getElementAt(i));
      h_trackerhit_timing -> Fill(trkhit->getTime());
      _clusters_vb->fill(trkhit);}
  for(int i=0; i<ibtrkhitCol->getNumberOfElements(); ++i)
    {
      const EVENT::TrackerHit *trkhit=static_cast<const EVENT::TrackerHit*>(ibtrkhitCol->getElementAt(i));
      h_trackerhit_timing -> Fill(trkhit->getTime());
      _clusters_ib->fill(trkhit);}      
  for(int i=0; i<obtrkhitCol->getNumberOfElements(); ++i)
    {
      const EVENT::TrackerHit *trkhit=static_cast<const EVENT::TrackerHit*>(obtrkhitCol->getElementAt(i));
      h_trackerhit_timing -> Fill(trkhit->getTime());}
  for(int i=0; i<vetrkhitCol->getNumberOfElements(); ++i)
    {
      const EVENT::TrackerHit *trkhit=static_cast<const EVENT::TrackerHit*>(vetrkhitCol->getElementAt(i));
      h_trackerhit_timing -> Fill(trkhit->getTime());}
  for(int i=0; i<ietrkhitCol->getNumberOfElements(); ++i)
    {
      const EVENT::TrackerHit *trkhit=static_cast<const EVENT::TrackerHit*>(ietrkhitCol->getElementAt(i));
      h_trackerhit_timing -> Fill(trkhit->getTime());}
  for(int i=0; i<oetrkhitCol->getNumberOfElements(); ++i)
    {
      const EVENT::TrackerHit *trkhit=static_cast<const EVENT::TrackerHit*>(oetrkhitCol->getElementAt(i));
      h_trackerhit_timing -> Fill(trkhit->getTime());}
}

void TrackPerfHistProc::processEvent( LCEvent * evt )
{

  // if track collection and MCTrackRel collection are empty, disable functionality that relies on tracks
  bool trackingOn = true;
  if( _trkColName == "" && _trkMatchColName == "") {trackingOn = false;}
  //
  // Get object required collections and create lists
  // to keep track of unsaved objects.

  TrackPerfHistProc::processEventTrackerHits(evt);

  if (trackingOn){
    TrackPerfHistProc::processEventTracks(evt);
  }
  
} 

void TrackPerfHistProc::check( LCEvent * /*evt*/ )
{ }

void TrackPerfHistProc::end()
{ }
