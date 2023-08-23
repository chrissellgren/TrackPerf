#pragma once

#include <marlin/Processor.h>

#include <TH1.h>
#include <TH2.h>

namespace TrackPerf
{
  class TrackHists;
  class TruthHists;
  class TrackResoHists;
  class TrackerHitResoHists;
  class ClusterHists;
}

//! Creates a simple column wise ntuple in a HistProc from LCIO collections.
class TrackPerfHistProc : public marlin::Processor
{
public:
  virtual Processor*  newProcessor() { return new TrackPerfHistProc ; }

  TrackPerfHistProc(const TrackPerfHistProc &) = delete ;
  TrackPerfHistProc& operator =(const TrackPerfHistProc &) = delete ;
  TrackPerfHistProc() ;

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;

  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;

  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 

  /** Sub-functions of processEvent, seperated for readability if no tracking enabled*/
  virtual void processEventTrackerHits( LCEvent * evt);
  virtual void processEventTracks( LCEvent * evt );

  virtual void check( LCEvent * evt ) ; 

  /** Called after data processing for clean up.
   */
  virtual void end() ;  

protected:

  //virtual void FindLocalPosition( const EVENT::TrackerHit* hit, double *localPosition, double *localDirection);


private:
  //! Track Collection
  std::string _trkColName {};

  //! MC Particle Collection
  std::string _mcpColName {};

  //! Track to MC truth match collection
  std::string _trkMatchColName {};

  //! Tracker hit collections
  std::string _vbtrkhitColName {};
  std::string _ibtrkhitColName {};
  std::string _obtrkhitColName {};
  std::string _vetrkhitColName {};
  std::string _ietrkhitColName {};
  std::string _oetrkhitColName {};

  //! Tracker hit relation collection
  std::string _VBRelationCollection {};

  //! Determination of good vs bad match
  float _matchProb = 0.5;

  // Determine whether tracking output collections should be included
  //bool _trackingOn = true;
  
  //! Histograms
  std::shared_ptr<TrackPerf::TrackHists> _allTracks ;
  std::shared_ptr<TrackPerf::TrackHists> _realTracks;
  std::shared_ptr<TrackPerf::TrackHists> _fakeTracks;
  std::shared_ptr<TrackPerf::TruthHists> _allTruths ;
  std::shared_ptr<TrackPerf::TruthHists> _realTruths;
  std::shared_ptr<TrackPerf::TruthHists> _unmtTruths;
  std::shared_ptr<TrackPerf::TrackResoHists> _realReso;
  std::shared_ptr<TrackPerf::TrackerHitResoHists> _uncertainties;
  std::shared_ptr<TrackPerf::ClusterHists> _clusters_vb;
  std::shared_ptr<TrackPerf::ClusterHists> _clusters_ib;

  TH1 * h_number_of_fakes;
  TH1 * h_number_of_tracks;
  TH1 * h_relation_weight_real;
  TH1 * h_relation_weight_fake;
  TH1 * h_trackerhit_timing;
};
