////////////////////////////////////////////////////////////////////////
// Class:       MCTruthProducer
// Plugin Type: producer (art v2_05_00)
// File:        MCTruthProducer_module.cc
//
// Generated at Mon Nov 05 14:43:37 2018 by Simone Marcocci using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "art/Utilities/make_tool.h"

#include "larcore/Geometry/Geometry.h"

// Data product include
#include "larcore/Geometry/Geometry.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Seed.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Wire.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

#include <memory>

// Algorithms include
#include "uboone/ProtonBenchmarker/Algos/MCTruthUtility.h"

class MCTruthProducer;


class MCTruthProducer : public art::EDProducer {
public:
  explicit MCTruthProducer(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MCTruthProducer(MCTruthProducer const &) = delete;
  MCTruthProducer(MCTruthProducer &&) = delete;
  MCTruthProducer & operator = (MCTruthProducer const &) = delete;
  MCTruthProducer & operator = (MCTruthProducer &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:

  std::string fTrackLabel;
  std::string fShowerLabel;
  std::string fPfpLabel;
  std::string fClusterLabel;
  std::string fG4TruthLabel;

  std::unique_ptr<truth::IMCTruthMatching> fMCTruthMatching;
  mctruthutil::MCTruthUtility mctruth_maker;

};


MCTruthProducer::MCTruthProducer(fhicl::ParameterSet const & p) {


  fTrackLabel = p.get<std::string> ("TrackLabel");
  fShowerLabel = p.get<std::string> ("ShowerLabel");
  fPfpLabel = p.get<std::string> ("PFParticleLabel");
  fClusterLabel = p.get<std::string> ("ClusterLabel");
  fG4TruthLabel = p.get<std::string> ("G4TruthLabel");
  
  // Get the tool for MC Truth matching
  fMCTruthMatching = art::make_tool<truth::IMCTruthMatching>(p.get<fhicl::ParameterSet>("MCTruthMatching"));
  
  produces< art::Assns<recob::Track , simb::MCParticle, anab::BackTrackerMatchingData > > ();
  produces< art::Assns<recob::Shower, simb::MCParticle, anab::BackTrackerMatchingData > > ();
  produces< art::Assns<recob::PFParticle , simb::MCParticle, anab::BackTrackerMatchingData > > ();

}

void MCTruthProducer::produce(art::Event & e)
{

  // Instantiate the output
  std::unique_ptr< art::Assns<recob::Track, simb::MCParticle, anab::BackTrackerMatchingData > > MCPartTrackassn( new art::Assns<recob::Track, simb::MCParticle, anab::BackTrackerMatchingData >);
  std::unique_ptr< art::Assns<recob::Shower, simb::MCParticle, anab::BackTrackerMatchingData> > MCPartShowerassn( new art::Assns<recob::Shower, simb::MCParticle, anab::BackTrackerMatchingData>);
  std::unique_ptr< art::Assns<recob::PFParticle, simb::MCParticle,anab::BackTrackerMatchingData> > MCPartPFParticleassn( new art::Assns<recob::PFParticle, simb::MCParticle, anab::BackTrackerMatchingData>);

  if ( e.isRealData() ) {
    std::cout << "[MCTruthProducer] Running on a real data file. No MC-PFP matching will be attempted." << std::endl;
    e.put(std::move(MCPartTrackassn));
    e.put(std::move(MCPartShowerassn));
    e.put(std::move(MCPartPFParticleassn));
    return;
  } 

  fMCTruthMatching->Rebuild(e);

  //get data products
  art::ValidHandle< std::vector<recob::Track> > trackHandle = e.getValidHandle< std::vector<recob::Track> > (fTrackLabel);
  std::vector< art::Ptr<recob::Track> > trackList;
  art::fill_ptr_vector(trackList, trackHandle); 
  art::FindManyP<recob::Hit> hitsFromTracks( trackHandle, e, fTrackLabel);

  art::Handle< std::vector<simb::MCParticle> > mcParticleHandle; 
  std::vector< art::Ptr<simb::MCParticle> > mcList;
  if (e.getByLabel( fG4TruthLabel, mcParticleHandle))
    art::fill_ptr_vector(mcList, mcParticleHandle); 

  art::ValidHandle< std::vector<recob::Shower> > showerHandle =
    e.getValidHandle< std::vector<recob::Shower> >(fShowerLabel);
  std::vector< art::Ptr<recob::Shower> >  showerList;  
  art::fill_ptr_vector(showerList, showerHandle); 
  art::FindManyP<recob::Hit> hitsFromShowers( showerHandle, e, fShowerLabel);

  art::ValidHandle< std::vector< recob::PFParticle > > pfpHandle = 
    e.getValidHandle< std::vector< recob::PFParticle > >(fPfpLabel);
  std::vector< art::Ptr<recob::PFParticle> > pfpList;
  art::fill_ptr_vector( pfpList, pfpHandle);
  art::FindManyP<recob::Cluster> clustersFromPfp( pfpHandle, e, fPfpLabel);
  
  art::ValidHandle< std::vector< recob::Cluster > > clusterHandle = 
    e.getValidHandle< std::vector< recob::Cluster > >(fClusterLabel);
  std::vector< art::Ptr<recob::Cluster> > clusterList;
  art::fill_ptr_vector( clusterList, clusterHandle);
  art::FindManyP<recob::Hit> hitsFromClusters(clusterHandle, e, fClusterLabel);

  //prepare vectors for making associations
  std::vector< art::Ptr<recob::Track> > track_v;
  std::vector< art::Ptr<simb::MCParticle> > track_mcp_v;
  std::vector< anab::BackTrackerMatchingData > track_btdata_v;

  std::vector< art::Ptr<recob::Shower> > shower_v;
  std::vector< art::Ptr<simb::MCParticle> > shower_mcp_v;
  std::vector< anab::BackTrackerMatchingData > shower_btdata_v;

  std::vector< art::Ptr<recob::PFParticle> > pfp_v;
  std::vector< art::Ptr<simb::MCParticle> > pfp_mcp_v;
  std::vector< anab::BackTrackerMatchingData > pfp_btdata_v;

  //make tracks
  mctruth_maker.FillTrackMCPVectors( fMCTruthMatching, hitsFromTracks, trackList, mcParticleHandle, track_mcp_v, track_v, track_btdata_v);
  //make showers
  mctruth_maker.FillShowerMCPVectors( fMCTruthMatching, hitsFromShowers, showerList, mcParticleHandle, shower_mcp_v, shower_v, shower_btdata_v);
  //make pfs
  mctruth_maker.FillPFParticleMCPVectors( fMCTruthMatching, clustersFromPfp, hitsFromClusters, pfpList, mcParticleHandle, pfp_mcp_v, pfp_v, pfp_btdata_v);

  std::cout << "VECTOR SIZE " << track_mcp_v.size() << std::endl;
  for ( unsigned ii=0; ii<track_mcp_v.size(); ii++ )
    	MCPartTrackassn->addSingle( track_v[ii], track_mcp_v[ii], track_btdata_v[ii] );

  for ( unsigned ii=0; ii<shower_mcp_v.size(); ii++ )
    	MCPartShowerassn->addSingle( shower_v[ii], shower_mcp_v[ii], shower_btdata_v[ii] );

  for ( unsigned ii=0; ii<pfp_mcp_v.size(); ii++ )
    	MCPartPFParticleassn->addSingle( pfp_v[ii], pfp_mcp_v[ii], pfp_btdata_v[ii] );


  e.put(std::move( MCPartTrackassn ));
  e.put(std::move( MCPartShowerassn ));
  e.put(std::move( MCPartPFParticleassn ));
}

DEFINE_ART_MODULE(MCTruthProducer)
