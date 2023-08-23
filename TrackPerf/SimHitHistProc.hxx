#pragma once

#include <marlin/Processor.h>

#include <TH2.h>
#include <TH1.h>

//#include <ACTSTracking/GeometryIdMappingTool.hxx>
#include "/opt/ilcsoft/muonc/ACTSTracking/v1.1.0/ACTSTracking/GeometryIdMappingTool.hxx"


/* namespace SimHitPlots
{
  class SimHitHists;
} */

//! Creates a simple column wise ntuple in a HistProc from LCIO collections.
class SimHitHistProc : public marlin::Processor
{
  public:
    virtual Processor*  newProcessor() { return new SimHitHistProc ; }

    SimHitHistProc (const SimHitHistProc &) = delete ;
    SimHitHistProc & operator =(const SimHitHistProc &) = delete ;
    SimHitHistProc() ;

    // Fill histograms with a single sim hit
    void fill(const EVENT::SimTrackerHit* simhit, const std::string flag);

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

    virtual void check( LCEvent * evt ) ; 

    /** Called after data processing for clean up.
     */
    virtual void end() ;

  private:
    // SimHit collection
    std::string _vbsimhitColName {};
    std::string _vesimhitColName {};
    std::string _ibsimhitColName {};
    std::string _iesimhitColName {};
    std::string _obsimhitColName {};
    std::string _oesimhitColName {};

    //! Tracker hit collections
    std::string _vbtrkhitColName {};
    std::string _ibtrkhitColName {};
    std::string _obtrkhitColName {};
    std::string _vetrkhitColName {};
    std::string _ietrkhitColName {};
    std::string _oetrkhitColName {};

    // Histograms
    TH1* h_x;
    TH1* h_y;
    TH1* h_z;
    TH1* h_r;
    TH1* h_t;
    TH2* h_z_r;
    TH2* h_x_y;

    TH1* h_x_vx;  
    TH1* h_y_vx;  
    TH1* h_z_vx;  
    TH1* h_r_vx;
    TH1* h_t_vx;
    TH2* h_z_r_vx;
    TH2* h_x_y_vx;
    TH1* h_t_tracker_vxb;

    TH1* h_x_it;  
    TH1* h_y_it;  
    TH1* h_z_it;  
    TH1* h_r_it;
    TH1* h_t_it;  
    TH2* h_z_r_it;
    TH2* h_x_y_it;

    TH1* h_x_ot;  
    TH1* h_y_ot;  
    TH1* h_z_ot;  
    TH1* h_r_ot;
    TH1* h_t_ot;  
    TH2* h_z_r_ot;
    TH2* h_x_y_ot;
};