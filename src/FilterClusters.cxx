#include "TrackPerf/FilterClusters.hxx"

#include <math.h>

#include <DD4hep/Detector.h>

#include <EVENT/Track.h>
#include <EVENT/TrackerHit.h>
#include <IMPL/LCCollectionVec.h>
#include <UTIL/LCTrackerConf.h>
#include <IMPL/LCRelationImpl.h>

#include <marlin/AIDAProcessor.h>

#include <AIDA/ITree.h>

#include "marlin/VerbosityLevels.h"


FilterClusters aFilterClusters ;

FilterClusters::FilterClusters()
  : Processor("FilterClusters")
{
  // modify processor description
  _description = "FilterClusters processor filters a collection of tracker hits based on cluster size and outputs a filtered collection";

  // register steering parameters: name, description, class-variable, default value
  registerProcessorParameter("InputRanges",
         "Divide theta or r into bins for different cluster size cuts",
         _InputRanges, 
         {}
          );
  registerProcessorParameter("DetectorType",
         "Indicate which geometry is being used",
         _DetectorType, 
         {}
          );
  registerProcessorParameter("FilterByLayer",
         "Indicate if filter should be uniform for each layer, or different ranges and cutoffs for each layer (double layer pair)",
         _FilterByLayer, 
         {}
          );
  registerProcessorParameter("ClusterSize",
		  	 "Maximum cluster size for each theta range",
			   _ClusterSize,
			   {}
			    );  
  registerProcessorParameter("Layers",
		  	 "Layers to be filtered",
			   _Layers,
			   {}
			    );  
  registerInputCollection( LCIO::TRACKERHIT,
		  	 "InTrackerHitCollection" ,
			   "Name of the input tracker hit collection",
			   _InTrackerHitCollection,
		     _InTrackerHitCollection
		 	    );

  registerInputCollection( LCIO::LCRELATION,
		     "InRelationCollection" ,
			   "Name of the input relation collection",
			   _InRelationCollection,
		     _InRelationCollection
		 	    );

  registerOutputCollection( LCIO::TRACKERHIT,
		  	 "OutTrackerHitCollection" ,
			   "Name of output tracker hit collection",
			   _OutTrackerHitCollection,
			   std::string("FilteredVBTrackerHits")
			    );

    registerOutputCollection( LCIO::LCRELATION,
		  	 "OutRelationCollection" ,
			   "Name of output relation collection",
			   _OutRelationCollection,
			   std::string("FilteredVBTrackerHitsRelations")
			    );

}

void FilterClusters::init()
{
  // Print the initial parameters
  printParameters() ;
}

void FilterClusters::processRunHeader( LCRunHeader* /*run*/)
{ }

// Function to split the input vector into chunks to specify cuts for each layer
/*std::vector<std::vector<std::string>> splitVector(const std::vector<std::string>& input, int numlayers, int originalLength) {
    std::vector<std::vector<std::string>> result;
    int startIndex = 0;
    int chunksize = originalLength / numlayers;

    // loop through all layers
    for (int i = 0; i < numlayers; i++) {
        std::vector<std::string> chunk;
        for (int j = startIndex; j < startIndex + chunksize && j < input.size(); j++) {
            chunk.push_back(input[j]);
        }
        result.push_back(chunk);
        startIndex += chunksize;
    }
    return result;
} */

void FilterClusters::processEvent( LCEvent * evt )
{
  // Determine if handling endcap or barrel
  if (_DetectorType.find("Barrel") != std::string::npos) {
      isBarrel=true;
    } else if (_DetectorType.find("Endcap") != std::string::npos) {
      isBarrel=false;
  } else {
      std::stringstream err  ; err << " Could not determine sub-detector type for: " << _DetectorType;
      throw Exception ( err.str() );
    }
  if (_DetectorType.find("Vertex") != std::string::npos) {
      isVertex=true;
      numlayers = 4;
  } else if (_DetectorType.find("InnerTracker") != std::string::npos) {
    isInnerTracker=true;
    if (isBarrel) numlayers = 3;
    // else numlayers = 8; // FIX: adjust for inner tracker endcap 
  } else if (_DetectorType.find("OuterTracker") != std::string::npos) {
    isOuterTracker=true;
    if (isBarrel) numlayers = 3;
    // else numlayers = 4; // FIX: adjust for outer tracker endcap 
  } else {
    std::stringstream err  ; err << " Could not determine sub-detector type for: " << _DetectorType;
    throw Exception ( err.str() );
  }
  streamlog_out(DEBUG8) << "Sub-detector element: Barrel: " << isBarrel << ", Vertex: " << isVertex << ", Inner Tracker: " << isInnerTracker << std::endl;
  streamlog_out(DEBUG8) << "Number of layers found for this sub-det: " << numlayers << std::endl;

  if (! _FilterByLayer) numlayers = 1;
  int rangesize = _InputRanges.size();
  int clustersize = _ClusterSize.size();

  // split the input vectors
  for (int layer = 0; layer < numlayers; layer++){
    std::vector<std::string> thisrangechunk;
    std::vector<std::string> thisclusterchunk;
    int chunksize_ranges = rangesize / numlayers;
    int chunksize_cluster = clustersize / numlayers;
    // split ranges vector
    for (int n = layer*chunksize_ranges; n < (layer+1)*chunksize_ranges; n++){
      thisrangechunk.push_back(_InputRanges[n]);
    }
    // cluster vector
    for (int n = layer*chunksize_cluster; n < (layer+1)*chunksize_cluster; n++){
      thisclusterchunk.push_back(_ClusterSize[n]);
    }
    splitInputRanges.push_back(thisrangechunk);
    splitClusterCuts.push_back(thisclusterchunk);
  }

  //splitInputRanges = splitVector(_InputRanges, numlayers, rangesize);
  //splitClusterCuts = splitVector(_ClusterSize, numlayers, clustersize);
    
  // Make the output track collection
  LCCollectionVec *OutTrackerHitCollection = new LCCollectionVec(LCIO::TRACKERHIT);
  OutTrackerHitCollection->setSubset(true);
  LCCollectionVec *OutRelationCollection   = new LCCollectionVec(LCIO::LCRELATION);
  OutRelationCollection->setSubset(true);

  // Get input collection
  LCCollection* InTrackerHitCollection  = evt->getCollection(_InTrackerHitCollection);
  LCCollection* InRelationCollection    = evt->getCollection(_InRelationCollection);

  if( InTrackerHitCollection->getTypeName() != lcio::LCIO::TRACKERHITPLANE )
    { throw EVENT::Exception( "Invalid collection type: " + InTrackerHitCollection->getTypeName() ) ; }
  if( InRelationCollection->getTypeName() != lcio::LCIO::LCRELATION )
    { throw EVENT::Exception( "Invalid collection type: " + InRelationCollection->getTypeName() ) ; }

  // Loop through each event in collection and filter
  for(int i=0; i<InTrackerHitCollection->getNumberOfElements(); ++i) //loop through all hits
    {
      streamlog_out(DEBUG0) << "Hit: " << i << std::endl;
      EVENT::TrackerHit *trkhit=static_cast<EVENT::TrackerHit*>(InTrackerHitCollection->getElementAt(i));
      EVENT::LCRelation *rel=static_cast<EVENT::LCRelation*>(InRelationCollection->getElementAt(i));
      
      //Calculating theta 
      float x = trkhit->getPosition()[0];
      float y = trkhit->getPosition()[1];
      float z = trkhit->getPosition()[2];
      float r = sqrt(pow(x,2)+pow(y,2));
      float incidentTheta = std::atan(r/z); 

      if(incidentTheta<0)
        incidentTheta += M_PI;

      //Calculating cluster size
      const lcio::LCObjectVec &rawHits = trkhit->getRawHits(); 
      float ymax = -1000000;
      float xmax = -1000000;
      float ymin =  1000000; 
      float xmin =  1000000; 

      float loopsize = rawHits.size();
      streamlog_out(DEBUG3) << "Number of raw hits: " << loopsize << std::endl;

      for (size_t j=0; j<loopsize; ++j) {
        lcio::SimTrackerHit *hitConstituent = dynamic_cast<lcio::SimTrackerHit*>( rawHits[j] );
        const double *localPos = hitConstituent->getPosition();
        float x_local = localPos[0];
        float y_local = localPos[1];
        streamlog_out(DEBUG6) << "Local y: " << y_local << ", local x: " << x_local << std::endl;
        if (y_local < ymin){
          ymin = y_local;
          }
        if (y_local > ymax){
          ymax = y_local;          
          } 

        if (x_local < xmin){
          xmin = x_local;
          }
        if (x_local > xmax){
          xmax = x_local;
          //streamlog_out(DEBUG2) << "Updated ymin: " << ymin << ", xmin: " << xmin << std::endl;
          //streamlog_out(DEBUG2) << "Updated ymax: " << ymax << ", xmax: " << xmax << std::endl;          
          }
        }
      streamlog_out(DEBUG2) << "Min y pos: " << ymin  << ", max y pos: " << ymax << std::endl;
      streamlog_out(DEBUG2) << "Min x pos: " << xmin  << ", max x pos: " << xmax << std::endl;
      float cluster_size_y = (ymax - ymin)+1;
      float cluster_size_x = (xmax - xmin)+1;
      float cluster_size_tot = loopsize;

      //Get hit subdetector/layer 
      std::string _encoderString = lcio::LCTrackerCellID::encoding_string();
      UTIL::CellIDDecoder<lcio::TrackerHit> decoder(_encoderString);
      uint32_t systemID = decoder(trkhit)["system"];
      uint32_t layerID = decoder(trkhit)["layer"];
      bool filter_layer = false;
      for (int i=0; i<_Layers.size(); ++i){
        if (layerID == std::stof(_Layers[i])) {
          filter_layer = true;
        }
      }

      // adjust for the double layer for vertex
      if (isVertex) layerID = std::floor(layerID/2);
      streamlog_out( DEBUG0 ) << "Layer ID (after dividing for DL): " << layerID << std::endl;
    
      std::vector<std::string> thisRanges;
      std::vector<std::string> thisClusterCuts;

      for (int layer = 0; layer < numlayers; ++layer){
        // choose ranges/cuts based on layer
        thisRanges = splitInputRanges[layer];
        thisClusterCuts = splitClusterCuts[layer];
        streamlog_out( DEBUG0 ) << "Layer: " << layer << std::endl;

        // only calculate filter if on the right layer
        if (layerID == layer){
          for (int i=0; i<thisRanges.size()-1; ++i) {
            float min = std::stof(thisRanges[i]);
            float max = std::stof(thisRanges[i+1]);
            streamlog_out( DEBUG0 ) << "Theta or R range: " << min << ", " << max << std::endl;
            
            if (isBarrel){
              streamlog_out( DEBUG0 ) << "Theta: " << incidentTheta << std::endl;
              if(incidentTheta > min && incidentTheta <= max && !filter_layer){
                streamlog_out( DEBUG0 ) << "Cluster size cut off: " << thisClusterCuts[i] << std::endl;
                streamlog_out( DEBUG0 ) << "Current y-cluster size (parallel to beam axis) : " << cluster_size_y << std::endl;
                if(cluster_size_y < std::stof(thisClusterCuts[i])) {
                  streamlog_out( DEBUG0 ) << "Cluster passes cut, hit added to output collection." << std::endl;
                  OutTrackerHitCollection->addElement(trkhit); 
                  OutRelationCollection->addElement(rel); }
                }
            }
            else {
              streamlog_out( DEBUG0 ) << "R: " << r << std::endl;
              if(r > min && r <= max && !filter_layer){
                streamlog_out( DEBUG0 ) << "Cluster size cut off: " << thisClusterCuts[i] << std::endl;
                streamlog_out( DEBUG0 ) << "Current total cluster size: " << cluster_size_tot << std::endl;
                if(cluster_size_tot < std::stof(thisClusterCuts[i])) {
                  streamlog_out( DEBUG0 ) << "Cluster passes cut, hit added to output collection." << std::endl;
                  OutTrackerHitCollection->addElement(trkhit); 
                  OutRelationCollection->addElement(rel); }
              }
            }
            }
        }
      }
      }
  // Save output track collection
  evt->addCollection(OutTrackerHitCollection, _OutTrackerHitCollection); 
  evt->addCollection(OutRelationCollection, _OutRelationCollection); 

  // Compare sizes of collections before/after filtering
  int nInTrackerHits = InTrackerHitCollection->getNumberOfElements();
  streamlog_out(MESSAGE) << "Number of Elements in Input Tracker Hits Collection: " << nInTrackerHits <<std::endl;
  int nOutTrackerHits = OutTrackerHitCollection->getNumberOfElements();
  streamlog_out(MESSAGE) << "Number of Elements in Output Tracker Hits Collection: " << nOutTrackerHits <<std::endl;
  if (nInTrackerHits != 0) {
    float retention = nOutTrackerHits/nInTrackerHits * 100;
    streamlog_out(MESSAGE) << "Percentage of hits retained: " << retention << std::endl;
  }
  
}


void FilterClusters::end()
{ }