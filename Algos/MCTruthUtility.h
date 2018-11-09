#ifndef MCTRUTHUTILITY_H
#define MCTRUTHUTILITY_H

//art includes
#include "canvas/Persistency/Common/FindManyP.h"
#include "art/Framework/Principal/Handle.h"


// larsoft includes
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Shower.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "larcore/Geometry/Geometry.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "uboone/AnalysisTree/MCTruth/IMCTruthMatching.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

// ROOT includes
#include "TLorentzVector.h"

// cc includes
#include <numeric>
#include <vector>
#include <functional>
#include <stdexcept>
#include <cmath>

namespace mctruthutil{

  class MCTruthUtility{
   
    public:


	    void FillPFParticleMCPVectors( std::unique_ptr<truth::IMCTruthMatching> &, art::FindManyP<recob::Cluster> &, art::FindManyP<recob::Hit> &,
							   std::vector< art::Ptr<recob::PFParticle> > &, art::Handle< std::vector<simb::MCParticle> > &, std::vector< art::Ptr<simb::MCParticle> > &, 
							   std::vector< art::Ptr<recob::PFParticle> > &, std::vector< anab::BackTrackerMatchingData > &);
	    
	    void FillShowerMCPVectors( std::unique_ptr<truth::IMCTruthMatching> &, art::FindManyP<recob::Hit> &, std::vector< art::Ptr<recob::Shower> > &,
			    			       art::Handle< std::vector<simb::MCParticle> > &, std::vector< art::Ptr<simb::MCParticle> > &, std::vector< art::Ptr<recob::Shower> > &, 
						       std::vector< anab::BackTrackerMatchingData > &);

	    void FillTrackMCPVectors( std::unique_ptr<truth::IMCTruthMatching> &, art::FindManyP<recob::Hit> &, std::vector< art::Ptr<recob::Track> > &,
			    			      art::Handle< std::vector<simb::MCParticle> > &, std::vector< art::Ptr<simb::MCParticle> > &, std::vector< art::Ptr<recob::Track> > &, 
						      std::vector< anab::BackTrackerMatchingData > &);
	
  };

}

#endif
