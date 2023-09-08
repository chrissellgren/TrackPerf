#pragma once

#include <marlin/Processor.h>
#include <vector>
#include <string>

//#include <ACTSTracking/GeometryIdMappingTool.hxx>
#include "/opt/ilcsoft/muonc/ACTSTracking/v1.1.0/ACTSTracking/GeometryIdMappingTool.hxx"

namespace TrackPerf
{
}

class FilterClusters : public marlin::Processor
{
public:
   virtual Processor* newProcessor() { return new FilterClusters ; }

   FilterClusters(const FilterClusters &) = delete ;
   FilterClusters& operator =(const FilterClusters &) = delete ;
   FilterClusters() ;

   /** Called at the begin of the job before anything is read.
    * Use to initialize the processor, e.g. book histograms.
    */
   virtual void init() ;

   /** Called for every run.
    */
   virtual void processRunHeader( LCRunHeader* run ) ;

   /** Called for every event - the working horse.
    */
   virtual void processEvent(LCEvent* evt) ;

   // Function to split long input vector for input ranges, cluster cuts into chunks
   std::vector<std::string>* splitVector(const std::vector<std::string>&, int, int);

   /** Called after data processing for clean up.
    */
   virtual void end() ;

private:
   //! Input track collection
   std::string _InTrackerHitCollection {};
   std::string _InRelationCollection {};

   //! Output track collection
   std::string _OutTrackerHitCollection {};
   std::string _OutRelationCollection {};
   std::string _DetectorType {};

   // Flags
   bool _FilterByLayer;
   bool isBarrel;
   bool isVertex;
   bool isInnerTracker;
   bool isOuterTracker;

   int numlayers;


   //! Ranges for theta or R
   std::vector<std::string> _InputRanges;

   //! Cut-offs for cluster size in various theta/R ranges
   std::vector<std::string> _ClusterSize;

   //! Layers to be filtered
   std::vector<std::string> _Layers;

   // Split vectors for different layers filtering
   std::vector<std::string>* splitInputRanges;
   std::vector<std::string>* splitClusterCuts;

   // sub-split vector for when looping over all layers
   std::vector<std::string> thisRanges;
   std::vector<std::string> thisClusterCuts;


};