////////////////////////////////////////////////////////////////////////
// Class:       ProtonBenchmarker
// Plugin Type: analyzer (art v2_05_00)
// File:        ProtonBenchmarker_module.cc
//
// Generated at Tue Sep 19 17:11:42 2017 by Adam Lister using cetskelgen
// from cetlib version v1_21_00.
//
//
// /** 
// Benchmarker for BNB only events (no cosmics) 
// */
////////////////////////////////////////////////////////////////////////

// base includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
//#include "art/Framework/IO/Root/RootInputFile.h"
#include "art/Framework/IO/Catalog/InputFileCatalog.h"

#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "lardataobj/MCBase/MCHitCollection.h"
#include "larpandora/LArPandoraInterface/LArPandora.h"
//#include "larpandora/LArPandoraObjects/PFParticleMetadata.h"
#include "Pandora/PdgTable.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"
#include <algorithm>

// local includes
#include "uboone/ProtonBenchmarker/Algos/protonBenchmarkerUtility.h"
#include "uboone/ProtonBenchmarker/Algos/HistoMaker.h"
#include "uboone/ProtonBenchmarker/Datatypes/StoredEvent.h"

//UBXSec includes
#include "uboone/UBXSec/DataTypes/SelectionResult.h"

#define isDebug 1

#define isSingleParticle 0

namespace recohelper {
  class ProtonBenchmarker;
}


class recohelper::ProtonBenchmarker : public art::EDAnalyzer {
  public:
    typedef std::map< size_t, art::Ptr<recob::PFParticle> > PFParticleIdMap; 

    explicit ProtonBenchmarker(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    ProtonBenchmarker(ProtonBenchmarker const &) = delete;
    ProtonBenchmarker(ProtonBenchmarker &&) = delete;
    ProtonBenchmarker & operator = (ProtonBenchmarker const &) = delete;
    ProtonBenchmarker & operator = (ProtonBenchmarker &&) = delete;

    // Required functions.
    void analyze(art::Event const & e) override;

    // Selected optional functions.
    void beginJob() override;
    void endJob() override;
    void GetPFParticleIdMap(const art::ValidHandle< std::vector<recob::PFParticle> > &, PFParticleIdMap &);
    void GetFinalStatePFParticleVectors(const PFParticleIdMap &, std::vector< art::Ptr<recob::PFParticle> > &, std::vector< art::Ptr<recob::PFParticle> > &);
    void CollectTracksAndShowers(const std::vector< art::Ptr<recob::PFParticle> >&, const art::FindManyP<recob::Track>, const art::FindManyP<recob::Shower>, std::vector< art::Ptr<recob::Track> >&, std::vector< art::Ptr<recob::Shower> > &);

  private:
    // fcl input parameters
    std::string fTrackLabel;
    std::string fPfpLabel;
    std::string fPfpAssnLabel;
    std::string fCalorimetryLabel;
    std::string fTrackTruthLabel;
    std::string fHitAssnTruthLabel;
    std::string fShowerLabel;
    std::string fShowerTruthLabel;
    std::string fClusterLabel;
    std::string fVertexFitterLabel;
    bool 	fIsVertexFitter;
    std::string fHitLabel;
    std::string fMCTruthLabel;
    std::string fG4TruthLabel;
    std::string fMCHitLabel;
    bool 	fWriteHistograms;
    bool	fIsUBXSec;
    bool	fIsConsolidated;


    pbutil::protonBenchmarkerUtility _pbutilInstance;

    StoredEvent* event_store;
    reco_histo::HistoMaker* histo_maker;

    art::ServiceHandle< art::TFileService > tfs;
    TTree* recoTree;

   int n_events;

};


recohelper::ProtonBenchmarker::ProtonBenchmarker(fhicl::ParameterSet const & p)
  :
    EDAnalyzer(p)  // ,
    // More initializers here.
{

  fTrackLabel = p.get<std::string> ("TrackLabel");
  fPfpLabel = p.get<std::string> ("PfpLabel");
  fPfpAssnLabel = p.get<std::string> ("PfpAssnLabel");
  fCalorimetryLabel = p.get<std::string> ("CalorimetryLabel");
  fTrackTruthLabel = p.get<std::string> ("TrackTruthLabel");
  fHitAssnTruthLabel = p.get<std::string> ("HitAssnTruthLabel");
  fShowerLabel = p.get<std::string> ("ShowerLabel");
  fVertexFitterLabel = p.get<std::string> ("VertexFitterLabel");
  fIsVertexFitter = p.get<bool> ("IsVertexFitter");
  fWriteHistograms = p.get<bool> ("WriteHistograms");
  fMCTruthLabel = p.get<std::string> ("MCTruthLabel");
  fG4TruthLabel = p.get<std::string> ("G4TruthLabel");
  fShowerTruthLabel = p.get<std::string> ("ShowerTruthLabel");
  fClusterLabel = p.get<std::string> ("ClusterLabel");
  fHitLabel = p.get<std::string> ("HitLabel");
  fMCHitLabel = p.get<std::string> ("MCHitLabel");
  fIsUBXSec = p.get<bool> ("UBXSecInput");
  fIsConsolidated = p.get<bool> ("IsConsolidated");

}

void recohelper::ProtonBenchmarker::beginJob()
{


  //std::cout << ">>>>>>>" << art::InputFileCatalog::currentFile()->fileName() << std::endl;
  recoTree = tfs->make<TTree>("recotree", "recotree");
  event_store = new StoredEvent();
  histo_maker = new reco_histo::HistoMaker();
  if ( fWriteHistograms ) {
  histo_maker->Init( tfs );
  art::TFileDirectory hits_dir = tfs->mkdir("hits_dir");
  histo_maker->Init_Hit(  hits_dir );
  }

  int bufsize    = 16000;
  int splitlevel = 99;
  // define branches
  recoTree->Branch("stored", &event_store, bufsize, splitlevel);

  n_events = 0;
   
}

void recohelper::ProtonBenchmarker::analyze(art::Event const & e)
{
  if (e.isRealData() == true) return;
  
  //auto const* theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
  auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();  
  bool space_charge = true;

  event_store->fEvent  = e.id().event();
  event_store->fRun    = e.run();
  event_store->fSubRun = e.subRun();

  //pandora MAP of PFPs and IDs
 PFParticleIdMap pfParticleMap;

  // get handles to objects of interest
 //pfps 
 art::ValidHandle< std::vector< recob::PFParticle > > pfpHandle = 
    e.getValidHandle< std::vector< recob::PFParticle > >(fPfpLabel);
  std::vector< art::Ptr<recob::PFParticle> > pfpList;
  
  //get it only if is consolidated
  if (fIsConsolidated)
  this->GetPFParticleIdMap(pfpHandle, pfParticleMap);
  
  art::ValidHandle< std::vector<recob::Track> > trackHandle = e.getValidHandle< std::vector<recob::Track> > (fTrackLabel);
  std::vector< art::Ptr<recob::Track> > trackList;
  art::ValidHandle< std::vector<recob::Shower> > showerHandle = e.getValidHandle< std::vector<recob::Shower> >(fShowerLabel);
  std::vector< art::Ptr<recob::Shower> >  showerList;  
  
  
  std::vector< art::Ptr<recob::PFParticle> > crParticles;//cosmic ray particles - ignore for now
  if (fIsConsolidated) {
  	//get the nu PFPs
	this->GetFinalStatePFParticleVectors(pfParticleMap, crParticles, pfpList);
  	
	art::FindManyP<recob::Track> trackFromPfp(pfpList, e, fPfpAssnLabel);
	art::FindManyP<recob::Shower> showerFromPfp(pfpList, e, fPfpAssnLabel);
	this->CollectTracksAndShowers( pfpList, trackFromPfp, showerFromPfp, trackList, showerList);
	std::cout << "Getting CONSOLIDATED pandora data products!" << std::endl;
  } else {
  	art::fill_ptr_vector(trackList, trackHandle); 
  	art::fill_ptr_vector(showerList, showerHandle); 
  	art::fill_ptr_vector( pfpList, pfpHandle);
  }

  art::FindManyP<simb::MCParticle, anab::BackTrackerMatchingData> mcpsFromTracks(trackList, e, fTrackTruthLabel);
  art::FindManyP<anab::Calorimetry> caloFromTracks(trackList, e, fCalorimetryLabel);
  art::FindManyP<simb::MCParticle, anab::BackTrackerMatchingData> mcpsFromShowers(showerList, e, fShowerTruthLabel);
  art::FindManyP<recob::PFParticle> pfpFromTrack(trackList, e, fPfpAssnLabel);

  art::Handle< std::vector<simb::MCParticle> > mcParticleHandle; 
  std::vector< art::Ptr<simb::MCParticle> > mcList;
  if (e.getByLabel( fG4TruthLabel, mcParticleHandle))
    art::fill_ptr_vector(mcList, mcParticleHandle); 

  art::Handle< std::vector<simb::MCTruth> > mcTruthHandle; 
  std::vector< art::Ptr<simb::MCTruth> > mcTruth;
  if (e.getByLabel( fMCTruthLabel, mcTruthHandle))
    art::fill_ptr_vector(mcTruth, mcTruthHandle); 

  // cluster/hit handles
  art::ValidHandle< std::vector< recob::Cluster > > clusterHandle = 
    e.getValidHandle< std::vector< recob::Cluster > >(fClusterLabel);
  std::vector< art::Ptr<recob::Cluster> > clusterList;
  art::fill_ptr_vector( clusterList, clusterHandle);

  art::FindManyP<recob::Hit> hitsFromClusters(clusterHandle, e, fClusterLabel);
  
  art::ValidHandle< std::vector< recob::Hit > > hitHandle = 
    e.getValidHandle< std::vector< recob::Hit > >(fHitLabel);
  std::vector< art::Ptr<recob::Hit> > hitList;
  art::fill_ptr_vector( hitList, hitHandle);

  art::FindManyP<recob::Hit> hitsFromTracks( trackList, e, fTrackLabel);
  art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> MCPfromhits( hitHandle, e, fHitAssnTruthLabel);
  art::FindManyP<recob::Hit, anab::BackTrackerHitMatchingData> hitsFromMCP( mcParticleHandle, e, fHitAssnTruthLabel);
  
  art::Handle< std::vector<sim::MCHitCollection> > mcHitHandle; 
  std::vector< art::Ptr<sim::MCHitCollection> > mcHitList;
  if (e.getByLabel( fMCHitLabel, mcHitHandle) )
    art::fill_ptr_vector(mcHitList, mcHitHandle); 
 
  art::Handle< std::vector<sim::SimChannel> > simChannelHandle; 
  std::vector< art::Ptr<sim::SimChannel> > simChannelList;
  if (e.getByLabel( fG4TruthLabel, simChannelHandle) )
    art::fill_ptr_vector(simChannelList, simChannelHandle); 
 
  //---------------------------------
  // MCParticles
  //---------------------------------

  //get MCparticles associated to MCtruth
  art::FindManyP<simb::MCParticle> MCpFromMCtruth( mcTruth , e, fG4TruthLabel );

  //vertex stuff
  art::FindManyP<recob::Vertex> vertexFromPfp(pfpList, e, fPfpLabel);
  art::FindManyP<recob::Track> trackFromPfp(pfpList, e, fPfpAssnLabel);
  
//hits associations
  art::FindManyP<recob::Cluster> clustersFromHits(hitHandle, e, fClusterLabel);
  art::FindManyP<recob::PFParticle> pfpFromCluster(clusterHandle, e, fClusterLabel);
  art::FindManyP<recob::SpacePoint> spacepointFromHits( hitHandle, e, fTrackLabel );

  //loop on all MC truth frames (mostly 1 per event)
  for ( unsigned n_truth = 0; n_truth < mcTruth.size(); n_truth++ ) {
	  //std::cout << ">>>>>>>>>>>>>>>>>EVENT" << std::endl; 
  n_events++;
  event_store->ClearVectors();

#if isSingleParticle == 0 
  //check that the neutrino is in the active volume
  const simb::MCNeutrino thisNeutrino = mcTruth[n_truth]->GetNeutrino();
  const simb::MCParticle thisLepton = thisNeutrino.Lepton();
  if ( !_pbutilInstance.isInTPC( thisLepton ) ) continue;
  
  //store the interaction info
  event_store->fccnc = thisNeutrino.CCNC();
  event_store->finteraction = thisNeutrino.InteractionType();
 
  double xOffset = 0.7 - SCE->GetPosOffsets( geo::Point_t(thisNeutrino.Nu().Position().X(),  thisNeutrino.Nu().Position().Y(), thisNeutrino.Nu().Position().Z())  ).x();
  double yOffset = SCE->GetPosOffsets( geo::Point_t(thisNeutrino.Nu().Position().X(),  thisNeutrino.Nu().Position().Y(), thisNeutrino.Nu().Position().Z()) ).y();
  double zOffset = SCE->GetPosOffsets( geo::Point_t(thisNeutrino.Nu().Position().X(),  thisNeutrino.Nu().Position().Y(), thisNeutrino.Nu().Position().Z()) ).z();
  if (!space_charge) {
	  xOffset=0;
	  yOffset=0;  zOffset=0;
  }
  event_store->fneutrino_x = thisNeutrino.Nu().Position().X() + xOffset;
  event_store->fneutrino_y = thisNeutrino.Nu().Position().Y() + yOffset;
  event_store->fneutrino_z = thisNeutrino.Nu().Position().Z() + zOffset;

#endif
  // first loop muons to find true neutrino induced muon (NIM)
  double muon_p = 0;
  double muon_px = 0;
  double muon_py = 0;
  double muon_pz = 0;
  double muon_endx = 0;
  double muon_endy = 0;
  double muon_endz = 0;
  int muon_id = -1;
  for (unsigned i = 0; i < MCpFromMCtruth.at(n_truth).size();i++){
  //for ( auto const& thisMcp : mcList ){
    const art::Ptr< simb::MCParticle >& thisMcp = MCpFromMCtruth.at(n_truth).at(i);
    if (std::abs(thisMcp->PdgCode()) == 13 && thisMcp->StatusCode()==1 && thisMcp->Mother()==0 && thisMcp->Process() == "primary" && _pbutilInstance.isInTPC(thisMcp) == true) {
      muon_p = thisMcp->Momentum().Rho();
      muon_px = thisMcp->Momentum().X();
      muon_py = thisMcp->Momentum().Y();
      muon_pz = thisMcp->Momentum().Z();
      muon_endx = thisMcp->EndPosition().X();
      muon_endy = thisMcp->EndPosition().Y();
      muon_endz = thisMcp->EndPosition().Z();
      muon_id = thisMcp->TrackId();
    }
  }

  if (muon_p==0) {
	  mf::LogError(__FUNCTION__) << "Error! There must be at least 1 muon in this analysis! " << std::endl;
	  return;
  }
 
  //now other MC particles
  for (unsigned i = 0; i < MCpFromMCtruth.at(n_truth).size(); i++){
    const art::Ptr< simb::MCParticle >& thisMcp = MCpFromMCtruth.at(n_truth).at(i);
    
    //in case of NC, the neutrino is kept and thus we should drop it
    //if it's an electron, check if it's a Michel. If yes, keep it.
    if ( !(std::abs(thisMcp->PdgCode()) == 11 && thisMcp->Mother() == muon_id && std::abs(thisMcp->Position().X()-muon_endx) < DBL_EPSILON &&
		    std::abs(thisMcp->Position().Y()-muon_endy) < DBL_EPSILON && std::abs(thisMcp->Position().Z()-muon_endz) < DBL_EPSILON ) )
    if ( std::abs(thisMcp->PdgCode()) == 14 || thisMcp->Process()!="primary" || thisMcp->StatusCode()!=1 || thisMcp->Mother()>0 ) continue; //we only want primaries

    if ( (std::abs(thisMcp->PdgCode()) == 11 && thisMcp->Mother() == muon_id && std::abs(thisMcp->Position().X()-muon_endx) < DBL_EPSILON &&
		    std::abs(thisMcp->Position().Y()-muon_endy) < DBL_EPSILON && std::abs(thisMcp->Position().Z()-muon_endz) < DBL_EPSILON ) )
	    std::cout << ">>>>>>>>>>>>>>>>>FOUND A MICHEL!!!!!" << std::endl;

#if isDebug == 1
      std::cout << "---- MCParticle Information ----"
        << "\nTrack ID   : " << thisMcp->TrackId() 
        << "\nPdgCode    : " << thisMcp->PdgCode()
        << "\nProcess    : " << thisMcp->Process()
        << "\nStatusCode : " << thisMcp->StatusCode()
        << "\nMother Pdg : " << thisMcp->Mother()
        << "\nPx, Py, Pz : " << thisMcp->Px() << ", " << thisMcp->Py() << ", " << thisMcp->Pz()
        << std::endl;
#endif

	//count particle types
	if (event_store->fparticle_count.find( thisMcp->PdgCode() ) == event_store->fparticle_count.end()) //not found
	     event_store->fparticle_count[ thisMcp->PdgCode() ] = 1;
	else
	     event_store->fparticle_count[thisMcp->PdgCode()] = event_store->fparticle_count[thisMcp->PdgCode()] +1;
	
	event_store->flength.push_back( thisMcp->Trajectory().TotalLength() );
	event_store->fn_steps.push_back ( thisMcp->Trajectory().size() );
	event_store->fstartT.push_back( thisMcp->T() );

	//anatree recipe:
  	xOffset = 0.7 - SCE->GetPosOffsets( geo::Point_t(thisNeutrino.Nu().Position().X(),  thisNeutrino.Nu().Position().Y(), thisNeutrino.Nu().Position().Z())  ).x();
 	yOffset = SCE->GetPosOffsets( geo::Point_t(thisNeutrino.Nu().Position().X(),  thisNeutrino.Nu().Position().Y(), thisNeutrino.Nu().Position().Z()) ).y();
  	zOffset = SCE->GetPosOffsets( geo::Point_t(thisNeutrino.Nu().Position().X(),  thisNeutrino.Nu().Position().Y(), thisNeutrino.Nu().Position().Z()) ).z();
	if (!space_charge) {
		xOffset = 0;
		yOffset = 0;
		zOffset = 0;
	}
	event_store->fstart_x.push_back ( thisMcp->Position().X() + xOffset );
	event_store->fstart_y.push_back ( thisMcp->Position().Y() + yOffset );
	event_store->fstart_z.push_back ( thisMcp->Position().Z() + zOffset );
	
	if (space_charge) {
  	xOffset = 0.7 - SCE->GetPosOffsets( geo::Point_t(thisNeutrino.Nu().Position().X(),  thisNeutrino.Nu().Position().Y(), thisNeutrino.Nu().Position().Z())  ).x();
 	yOffset = SCE->GetPosOffsets( geo::Point_t(thisNeutrino.Nu().Position().X(),  thisNeutrino.Nu().Position().Y(), thisNeutrino.Nu().Position().Z()) ).y();
  	zOffset = SCE->GetPosOffsets( geo::Point_t(thisNeutrino.Nu().Position().X(),  thisNeutrino.Nu().Position().Y(), thisNeutrino.Nu().Position().Z()) ).z();
	}
	event_store->fend_x.push_back ( thisMcp->EndPosition().X() + xOffset );
	event_store->fend_y.push_back ( thisMcp->EndPosition().Y() + yOffset );
	event_store->fend_z.push_back ( thisMcp->EndPosition().Z() + zOffset );
	event_store->fmother_id.push_back ( thisMcp->Mother() );
	event_store->fpdg.push_back( thisMcp->PdgCode() );	
	event_store->fg4_id.push_back(thisMcp->TrackId());
	event_store->fp0.push_back( thisMcp->Momentum().Rho());
	event_store->fp0x.push_back( thisMcp->Momentum().X() );
	event_store->fp0y.push_back( thisMcp->Momentum().Y() );
	event_store->fp0z.push_back( thisMcp->Momentum().Z() );
	
	//loop on other tracked_particles and decide if this is the leading based on initial kinetic energy. This must be done before filling the kinE vector for the current particle!
	bool is_leading = true;
	float current_kinE = thisMcp->E() - thisMcp->Mass() ;
	for (unsigned jj=0; jj< event_store->fkinE.size(); jj++) { //loop on previous particles
	if ( event_store->fpdg[jj] == thisMcp->PdgCode() ) {
	    if ( event_store->fkinE[jj] >= current_kinE )
	       is_leading = false;
	    else 
	       event_store->fis_leading[jj] = false;
	    } 
	}
	
	event_store->fis_leading.push_back( is_leading );

	//now write current kin E
	event_store->fkinE.push_back( thisMcp->E() - thisMcp->Mass() );

	if ( thisMcp->Momentum().Rho() != 0) {
	     event_store->fcostheta_muon.push_back( muon_px*thisMcp->Momentum().X()/thisMcp->Momentum().Rho()/muon_p
			     	+ muon_py*thisMcp->Momentum().Y()/thisMcp->Momentum().Rho()/muon_p
				+ muon_pz*thisMcp->Momentum().Z()/thisMcp->Momentum().Rho()/muon_p );
	} else {
	     event_store->fcostheta_muon.push_back(-2);
	}
	
	
	//add dummy entries for the reco variables. They will be possibly updated later, if a matching reco track is found
	event_store->AllocateRecoVectors();
      
  } // MCParticles
  
  //this stays out of the previous loop but it is somewhat connected
   for ( auto const& mcp : mcList ) {
	   auto itt = std::find ( event_store->fg4_id.begin(), event_store->fg4_id.end(),  mcp->TrackId());
	   if ( itt == event_store->fg4_id.end() ) continue;
  //nhits from MCP
	std::vector< art::Ptr < recob::Hit> > hits_mcp = hitsFromMCP.at( &mcp - &mcList[0] );
	//std::cout <<  itt - fg4_id.begin() << " " << freco_mcp_hits.size() << std::endl;
	event_store->freco_mcp_hits[ itt - event_store->fg4_id.begin() ] = hits_mcp.size() ;
	
	int n_collection_hits = 0;
	int total_hit_charge = 0;
	//std::cout << "NUMBER OF HITS " << hits_mcp.size() << std::endl;
	event_store->fnhits[ itt - event_store->fg4_id.begin() ] = 0;
	if ( hits_mcp.size() ) {
	//mcp hits
	for ( auto const& iter_hit : hits_mcp ){
		if ( iter_hit->View() == geo::kW ){	
			n_collection_hits++;
			total_hit_charge += iter_hit->Integral();
		}
		event_store->fnhits[ itt - event_store->fg4_id.begin() ] = event_store->fnhits[ itt - event_store->fg4_id.begin() ] + 1;
	}
	event_store->freco_mcp_collection_hits[itt - event_store->fg4_id.begin() ] =  n_collection_hits ;
	event_store->freco_mcp_collection_charge[itt - event_store->fg4_id.begin()] = total_hit_charge ;
	}
	hits_mcp.clear();
   }

   if ( fWriteHistograms ) 
   histo_maker->Fill_Truth_Histos( event_store );

  //----------------------------
  // Tracks
  //----------------------------

  //I can't use thisTrack->ID() as index, because in case of Kalman Fitter failure, there could be mismatches between the number
  //of elements in the associations and the track index. So I need to index manually. 
  //There is an additional crosscheck for which the trackID must be increasing in the loop, just to be sure that nothing nasty is happening.
  // loop tracks and do truth matching 
  for (auto const& thisTrack : trackList) {
    
    std::vector< art::Ptr<simb::MCParticle> > mcps = mcpsFromTracks.at( &thisTrack - &trackList[0] );
    std::vector< art::Ptr<anab::Calorimetry> > calos = caloFromTracks.at( &thisTrack - &trackList[0] );

    if (mcps.size() >1 ) mf::LogWarning(__FUNCTION__) << "Warning !!! More than 1 MCparticle associated to the same track!" << std::endl;
    if (calos.size() != 3 ) mf::LogWarning(__FUNCTION__) << "Warning !!! Calorimetry info size " << calos.size() << " != 3 associated to tracks!" << std::endl;
	
    for (auto const& thisMcp : mcps){
    
    //this if is necessary, since sometimes a particle is matched to a secondary (electron, etc) -> to be checked
    if ( !(std::abs(thisMcp->PdgCode()) == 11 && thisMcp->Mother() == muon_id && std::abs(thisMcp->Position().X()-muon_endx) < DBL_EPSILON &&
		     std::abs(thisMcp->Position().Y()-muon_endy) < DBL_EPSILON && std::abs(thisMcp->Position().Z()-muon_endz) < DBL_EPSILON ) )
    if ( std::abs(thisMcp->PdgCode()) == 14 || thisMcp->Process()!="primary" || thisMcp->StatusCode()!=1 || thisMcp->Mother()>0 ) continue; //we only want primaries
      
#if isDebug == 1
      std::cout << "MATCHED PARTICLE: " << std::endl;
      std::cout << "---- MCParticle Information ----"
        << "\nTrack ID   : " << thisMcp->TrackId() 
        << "\nPdgCode    : " << thisMcp->PdgCode()
        << "\nProcess    : " << thisMcp->Process()
        << "\nStatusCode : " << thisMcp->StatusCode()
        << "\nMother Pdg : " << thisMcp->Mother()
        << "\nPx, Py, Pz : " << thisMcp->Px() << ", " << thisMcp->Py() << ", " << thisMcp->Pz()
        << std::endl;
#endif
    auto it_found = std::find( event_store->fg4_id.begin(), event_store->fg4_id.end(), thisMcp->TrackId() ) ;
    if (it_found==event_store->fg4_id.end()) {
	    mf::LogWarning(__FUNCTION__) << "ERROR!!! Matched particle not found!" << std::endl;
	    continue;
    }
    size_t pos = it_found - event_store->fg4_id.begin();

    //save information on reco track
    if (event_store->fis_tracked[pos] == true ) {
	    std::cout << "Probably broken track!" << std::endl;
	    event_store->fmatch_multiplicity[pos] = event_store->fmatch_multiplicity[pos] + 1;
    	    //decide if keeping the old info or the new one - based on the minimum distance between real start and recoed end/start
	    float start_distance = sqrt( pow( thisTrack->Vertex().X() - event_store->fstart_x[pos] , 2 ) 
		      +  pow( thisTrack->Vertex().Y() - event_store->fstart_y[pos] , 2 )
		      +   pow( thisTrack->Vertex().Z() - event_store->fstart_z[pos] , 2 ) );
	    float end_distance = sqrt( pow( thisTrack->End().X() - event_store->fstart_x[pos] , 2 ) 
		      +  pow( thisTrack->End().Y() - event_store->fstart_y[pos] , 2 )
		      +   pow( thisTrack->End().Z() - event_store->fstart_z[pos] , 2 ) );
	    float new_distance = min( start_distance, end_distance );
	    start_distance = sqrt( pow( event_store->freco_startx[pos]  - event_store->fstart_x[pos] , 2 ) 
		      +  pow( event_store->freco_starty[pos] - event_store->fstart_y[pos] , 2 )
		      +   pow( event_store->freco_startz[pos] - event_store->fstart_z[pos] , 2 ) );
	    end_distance = sqrt( pow( event_store->freco_endx[pos] - event_store->fstart_x[pos] , 2 ) 
		      +  pow( event_store->freco_endy[pos] - event_store->fstart_y[pos] , 2 )
		      +  pow(event_store->freco_endz[pos] - event_store->fstart_z[pos] , 2 ) );
	    float old_distance = min( start_distance, end_distance );
	    if ( old_distance < new_distance )
		    continue;
    }
        
        event_store->fis_tracked[pos] = true;
	event_store->fmatch_multiplicity[pos] = event_store->fmatch_multiplicity[pos] + 1;
	event_store->flength_reco[pos] = thisTrack->Length();
	trkf::TrackMomentumCalculator trkm; //track momentum calculator
	trkm.SetMinLength(0); //minimum track length for momentum calculation
	event_store->freco_momentum_mcs[pos] = trkm.GetMomentumMultiScatterChi2( thisTrack );
	event_store->freco_momentum_mcs_llhd[pos] = trkm.GetMomentumMultiScatterLLHD( thisTrack ) ;
	event_store->freco_momentum_range[pos] = trkm.GetTrackMomentum( thisTrack->Length(), event_store->fpdg[pos] ); //use info on pdg
	event_store->freco_startx[pos] = thisTrack->Vertex().X(); 
	event_store->freco_starty[pos] = thisTrack->Vertex().Y(); 
	event_store->freco_startz[pos] = thisTrack->Vertex().Z(); 
	event_store->freco_endx[pos] = thisTrack->End().X(); 
	event_store->freco_endy[pos] = thisTrack->End().Y(); 
	event_store->freco_endz[pos] = thisTrack->End().Z(); 
	event_store->freco_trackid[pos] = thisTrack->ID();
	
//	std::cout << "TRACK ID " << freco_trackid[pos] << " pos=" << pos << std::endl;
 	
	if (thisMcp->PdgCode() == 13 ) { //for the muon, only when there is the info
		if (event_store->fmuon_dqdx.size()!=0) {
			mf::LogError(__FUNCTION__) << "Calorimetry should be filled only once!!!! Clearing..." << std::endl;
			event_store->fmuon_dqdx.clear();
			event_store->fmuon_dedx.clear();
			event_store->fmuon_residual.clear();
		}

		for (size_t position=0; position<calos.at( geo::kCollection )->dQdx().size(); position++) {
		event_store->fmuon_dqdx.push_back(calos.at( geo::kCollection )->dQdx()[position]); //look only at collection
		event_store->fmuon_dedx.push_back(calos.at( geo::kCollection )->dEdx()[position]); //look only at collection
		event_store->fmuon_residual.push_back(calos.at( geo::kCollection )->ResidualRange()[position]); //look only at collection
		}
	event_store->fmuon_range = calos.at( geo::kCollection )->Range();
	}
    
	auto thisMcpCleanliness = mcpsFromTracks.data( &thisTrack - &trackList[0] ).at(0)->cleanliness;
	auto thisMcpCompleteness = mcpsFromTracks.data( &thisTrack - &trackList[0] ).at(0)->completeness;
	event_store->fpurity[pos] = thisMcpCleanliness;
	event_store->fcompleteness[pos] = thisMcpCompleteness;

        //check association of hits <-> tracks
	//nhits from track
	std::vector< art::Ptr < recob::Hit> > hits_tracks = hitsFromTracks.at(  &thisTrack - &trackList[0]  );
	event_store->freco_track_hits[pos] = hits_tracks.size() ;

	//loop over hits, and save some information for collection hits only (for now)
	int n_collection_hits = 0;
	float total_hit_charge = 0;
	//track hits
	for ( auto const& iter_hit : hits_tracks ){
		if ( iter_hit->View() == geo::kW ){	
			n_collection_hits++;
			total_hit_charge += iter_hit->Integral();
		}
	}
	event_store->freco_track_collection_hits[pos] = n_collection_hits;
	event_store->freco_track_collection_charge[pos] = total_hit_charge;
	hits_tracks.clear();	
    }//MCParticle
        
	mcps.clear();
	calos.clear();
  }//Tracks


  //checks on showers
  //I am indexing showers manually because I am lazy and I am copying what I did for tracks
  std::vector<int> mcp_showers_ids;
  std::cout << "LOOPING over " << showerList.size() << " showers " << std::endl;
  for (auto const& thisShower : showerList) {
    
    std::vector< art::Ptr<simb::MCParticle> > mcps = mcpsFromShowers.at( &thisShower - &showerList[0] );

    if (mcps.size() >1 ) mf::LogWarning(__FUNCTION__) << "Warning !!! More than 1 MCparticle associated to the same shower!" << std::endl;
    std::cout << "Found " << mcps.size() << " MCPs for this shower " << &thisShower - &showerList[0] << std::endl;
    for (auto const& thisMcp : mcps){
    
    //this if is necessary, since sometimes a particle is matched to a secondary (electron, etc) -> to be checked
    if ( !(std::abs(thisMcp->PdgCode()) == 11 && thisMcp->Mother() == muon_id && std::abs(thisMcp->Position().X()-muon_endx) < DBL_EPSILON &&
		    std::abs(thisMcp->Position().Y()-muon_endy) < DBL_EPSILON && std::abs(thisMcp->Position().Z()-muon_endz) < DBL_EPSILON ) )
    if ( std::abs(thisMcp->PdgCode()) == 14 || thisMcp->Process()!="primary" || thisMcp->StatusCode()!=1 || thisMcp->Mother()>0 ) continue; //we only want primaries
    
#if isDebug == 1
      std::cout << "SHOWER MATCHED PARTICLE: " << std::endl;
      std::cout << "---- MCParticle Information ----"
        << "\nTrack ID   : " << thisMcp->TrackId() 
        << "\nPdgCode    : " << thisMcp->PdgCode()
        << "\nProcess    : " << thisMcp->Process()
        << "\nStatusCode : " << thisMcp->StatusCode()
        << "\nMother Pdg : " << thisMcp->Mother()
        << "\nPx, Py, Pz : " << thisMcp->Px() << ", " << thisMcp->Py() << ", " << thisMcp->Pz()
        << std::endl;
#endif
    auto it_found = std::find( event_store->fg4_id.begin(), event_store->fg4_id.end(), thisMcp->TrackId() ) ;
    if (it_found==event_store->fg4_id.end()) {
	    mf::LogWarning(__FUNCTION__) << "ERROR!!! Matched particle not found!" << std::endl;
	    continue;
    }
    size_t pos = it_found - event_store->fg4_id.begin();
	

    //save information on reco track
    std::cout << "FOUND A SHOWER!!! >>>>>>>>> pdg=" << thisMcp->PdgCode() << std::endl;
    event_store->fshower_pdg[pos] = thisMcp->PdgCode();
    
    if (event_store->fis_shower_matched[pos] == true ) {
	    std::cout << "Already matched to a shower!!!! Skipping." <<  std::endl;
	    continue;
    }
    
    event_store->fis_shower_matched[pos] = true;
    if ( event_store->fpdg[pos] == 2212 ) { //protons
	    if ( std::find( mcp_showers_ids.begin(), mcp_showers_ids.end(), event_store->fg4_id[pos] ) == mcp_showers_ids.end() ) event_store->fcount_proton_showers = event_store->fcount_proton_showers + 1;
	    std::cout << event_store->fcount_proton_showers << std::endl;
	    std::cout << ">>>>PROTON SHOWER" << std::endl;
	    //std::cout << ">>>>>>>FILE " << art::RootInputFile::fileName << std::endl;
	    std::cout << "Event=" << event_store->fEvent << " Run=" << event_store->fRun << " SubRun=" << event_store->fSubRun << std::endl;
    }
    mcp_showers_ids.push_back(event_store->fg4_id[pos]);
 	
    }//MCParticle
        
	mcps.clear();
  }//Shower

 	mcp_showers_ids.clear();	


  //----------------------------
  // Vertex
  //----------------------------
  
  //look for neutrino pfp
  bool neutrino_set = false;
  for ( auto const& pfp : pfpList ) {
	  std::vector< art::Ptr <recob::Vertex> > vertex_pfp = vertexFromPfp.at( &pfp - &pfpList[0] );
	  
	  if ( vertex_pfp.size()==0 ) continue; //no vertex
	  if ( vertex_pfp.size()>1 ) mf::LogError(__FUNCTION__) << "ERROR! There should be only 1 vertex per PFP! Considering only the first one..." <<std::endl;
	  double xyzz[3];
	  vertex_pfp[0]->XYZ(xyzz);
	  
	  if ( lar_pandora::LArPandoraHelper::IsNeutrino(pfp) ) {
	  	  if (neutrino_set) mf::LogError(__FUNCTION__) << "There should be only 1 neutrino per PFP!!!!" << std::endl;
		  neutrino_set = true;
	  //save info on neutrino reco'ed vertex
	  event_store->fnu_reco_x = xyzz[0];
	  event_store->fnu_reco_y = xyzz[1];
	  event_store->fnu_reco_z = xyzz[2];
	  }
	  
	  //investigate matching particles and tracks
	  std::vector< art::Ptr <recob::Track> > track_pfp = trackFromPfp.at( &pfp - &pfpList[0] );
	  if (track_pfp.size()==0) continue; //no track! (shower?)
	  int trackid = track_pfp[0]->ID();
	  auto iter = std::find( event_store->freco_trackid.begin(), event_store->freco_trackid.end(), trackid );
	  if ( iter == event_store->freco_trackid.end() ) continue; //not found
	  unsigned mc_pos = iter - event_store->freco_trackid.begin();
          
	  event_store->freco_vertex_x[ mc_pos ] = xyzz[0];
	  event_store->freco_vertex_y[ mc_pos ] = xyzz[1];
	  event_store->freco_vertex_z[ mc_pos ] = xyzz[2];
	  vertex_pfp.clear();
	  track_pfp.clear();
  }

  
  //----------------------------
  // Vertex fitter (Giuseppe)
  //----------------------------
 
  if (fIsVertexFitter) {
  art::FindManyP<recob::Vertex> vertexfitterFromPfp(pfpList, e, fVertexFitterLabel);
  
  //look for neutrino pfp
  neutrino_set = false;
  for ( auto const& pfp : pfpList ) {
	  std::vector< art::Ptr <recob::Vertex> > vertexfitter_pfp = vertexfitterFromPfp.at( &pfp - &pfpList[0] );
	  
	  if ( vertexfitter_pfp.size()==0 ) continue; //no vertex
	  if ( vertexfitter_pfp.size()>1 ) mf::LogError(__FUNCTION__) << "Vertex Fitter: ERROR! There should be only 1 vertex per PFP! Considering only the first one..." <<std::endl;
	  double xyzz[3];
	  vertexfitter_pfp[0]->XYZ(xyzz);
	  
	  if ( lar_pandora::LArPandoraHelper::IsNeutrino(pfp) ) {
	  	  if (neutrino_set) mf::LogError(__FUNCTION__) << "Vertex Fitter: There should be only 1 neutrino per PFP!!!!" << std::endl;
		  neutrino_set = true;
	  //save info on neutrino reco'ed vertex
	  event_store->fnu_reco_fitter_x = xyzz[0];
	  event_store->fnu_reco_fitter_y = xyzz[1];
	  event_store->fnu_reco_fitter_z = xyzz[2];
	  //fnu_reco_fitter_chi2ndf = vertexfitter_pfp[0]->chi2PerNdof();
	  //fnu_reco_fitter_chi2 = vertexfitter_pfp[0]->chi2();
	  event_store->fnu_reco_fitter_chi2ndf = -1;
	  event_store->fnu_reco_fitter_chi2 = -1;
	  }
	  
	  //investigate matching particles and tracks
	  std::vector< art::Ptr <recob::Track> > track_pfp = trackFromPfp.at( &pfp - &pfpList[0] );
	  if (track_pfp.size()==0) continue; //no track! (shower?)
	  int trackid = track_pfp[0]->ID();
	  auto iter = std::find( event_store->freco_trackid.begin(), event_store->freco_trackid.end(), trackid );
	  if ( iter == event_store->freco_trackid.end() ) continue; //not found
	  unsigned mc_pos = iter - event_store->freco_trackid.begin();
          
	  event_store->freco_vertexfitter_x[ mc_pos ] = xyzz[0];
	  event_store->freco_vertexfitter_y[ mc_pos ] = xyzz[1];
	  event_store->freco_vertexfitter_z[ mc_pos ] = xyzz[2];
	  //freco_vertexfitter_chi2ndf[ mc_pos ] = vertexfitter_pfp[0]->chi2PerNdof(); 
	  //freco_vertexfitter_chi2[ mc_pos ] = vertexfitter_pfp[0]->chi2(); 
	  event_store->freco_vertexfitter_chi2ndf[ mc_pos ] = -1; 
	  event_store->freco_vertexfitter_chi2[ mc_pos ] = -1; 
	  track_pfp.clear();
	  vertexfitter_pfp.clear();
  }
  } //is vertex fitter


  //std::cout << "NEUTRINO " <<  fnu_reco_x << " " << fnu_reco_y << " " << fnu_reco_z << std::endl;

  int muon_pos;
  histo_maker->Fill_Analysis_Histos( event_store, muon_pos, fWriteHistograms); //plots tracking and vertexing information
  
  if (muon_pos != -1) {
  //now loop over hits and check:
  //-which MCP the hit belongs to
  //-is the MCP tracked?
  //-if is tracked, check if the hit is clustered and attributed properly to the track
  //-if it is not tracked, save which other particle/track got the hit or if it was not clustered
  
  for ( auto const& hit: hitList ) {
	  //only care at collection for now
	if ( hit->View() != geo::kW ) continue;
	std::vector< art::Ptr <simb::MCParticle> > mcparticle_hits = MCPfromhits.at( &hit - &hitList[0] );
	if (  mcparticle_hits.size() == 0 ) continue; 
	std::vector< art::Ptr <recob::SpacePoint> > spacepoint_hits = spacepointFromHits.at( &hit - &hitList[0] );
	art::ServiceHandle<geo::Geometry> geom;
	std::vector<double> xyz = { 0, 0, 0};
	double xyz2[3] = {0,0,0};
	bool is_spacepoint = false;
	if ( spacepoint_hits.size() != 1) {
		geom->WireIDToWireGeo( hit->WireID() ).GetCenter( xyz2, 0 );
		xyz[0] = xyz2[0];
		xyz[1] = xyz2[1];
		xyz[2] = xyz2[2];
//		std::cout << "NO SPACE POINT!" << std::endl;
	} else {
		is_spacepoint = true;
		xyz[0] = spacepoint_hits[0]->XYZ()[0];
		xyz[1] = spacepoint_hits[0]->XYZ()[1];
		xyz[2] = spacepoint_hits[0]->XYZ()[2];
	}
   	 
	//std::cout << "size2 " <<  mcparticle_hits.size() << std::endl;
	int index = -1;
	for (unsigned jj=0; jj<mcparticle_hits.size(); jj++) {
			bool is_max = MCPfromhits.data( &hit - &hitList[0] ).at( jj )->isMaxIDE;
			if ( is_max ) {
				index = jj;
				//std::cout << "WARNING! Can't be already !=0" << std::endl;
			}
	}

//	std::cout << "INDEX " << index << std::endl;

	if ( index == -1 ) {
		std::cout << ">>>>>ERROR!!!!! There should always be 1 MCP associated to the hit! " << mcparticle_hits.size() << std::endl;
		continue;
	}

	int trackid = mcparticle_hits[index]->TrackId();
	//lookup the originating MC particle
	bool found = false;
	auto iter = std::find (event_store->fg4_id.begin(), event_store->fg4_id.end(), trackid );
	if ( iter != event_store->fg4_id.end() ) found = true;
	else {
		iter = std::find (event_store->fg4_id.begin(), event_store->fg4_id.end(), mcparticle_hits[index]->Mother() );
		if ( iter != event_store->fg4_id.end() ) found = true;
	}
	if (!found) { 
#if DEBUG == 1
		std::cout<< "WARNING!!!!!!! MCParticle NOT found! " << std::endl;
      std::cout << "---- MCParticle Information ----"
        << "\nTrack ID   : " << mcparticle_hits[index]->TrackId() 
        << "\nPdgCode    : " << mcparticle_hits[index]->PdgCode()
        << "\nProcess    : " << mcparticle_hits[index]->Process()
        << "\nStatusCode : " << mcparticle_hits[index]->StatusCode()
        << "\nMother Pdg : " << mcparticle_hits[index]->Mother()
        << "\nPx, Py, Pz : " << mcparticle_hits[index]->Px() << ", " << mcparticle_hits[index]->Py() << ", " << mcparticle_hits[index]->Pz()
        << std::endl;
#endif
      continue;
	}

	

	int mcp_index = iter - event_store->fg4_id.begin();
	event_store->fis_spacepoint[mcp_index].push_back(is_spacepoint);
	event_store->fspacepoint_xyz[mcp_index].push_back(xyz);

	//freco_mcp_hits	
  //mcp index allows me to understand if this particle is tracked or not
      // std::cout << "clusters from hits " << clustersFromHits.size() << std::endl;
//	       std::cout << "pos " << &hit - &hitList[0] << std::endl;
	       std::vector< art::Ptr <recob::Cluster>> cluster_hits = clustersFromHits.at(  &hit - &hitList[0] );
	       if ( cluster_hits.size() > 3) {
		       //std::cout << ">>>>>>>>>>>>>>>>>>>>>>>" << clustersFromHits.size() << std::endl;
		       throw cet::exception("Error") << "Number of clusters must be maximum one!!!" << std::endl;
	       }


	       if ( cluster_hits.size() == 0 ) { //this hit is not clustered!
		      
		       event_store->fnot_clustered_hit_index[mcp_index].push_back( event_store->fis_spacepoint[mcp_index].size() - 1); //save index in terms of all hits
		       event_store->fnot_clustered[mcp_index] = event_store->fnot_clustered[mcp_index] + 1;
		       event_store->fnot_clustered_charge[mcp_index].push_back( hit->Integral() );
	   		if ( event_store->fis_tracked[mcp_index] ) //check and fill histos for a MCP which is tracked
				event_store->fnot_clustered_tracked_charge[mcp_index].push_back( hit->Integral() );
			else
				event_store->fnot_clustered_not_tracked_charge[mcp_index].push_back( hit->Integral() );

	       } else { //in this case the hit is clustered
			unsigned cluster_index = 0;
			for ( unsigned zz=0; zz<cluster_hits.size(); zz++) {
				if ( cluster_hits[zz]->View() == geo::kW ) cluster_index = zz;  //get collection plane ones
			}
		       
			event_store->fclustered[mcp_index] = event_store->fclustered[mcp_index] + 1;
		        event_store->fclustered_charge[mcp_index].push_back( hit->Integral() );
		        event_store->fclustered_hit_index[mcp_index].push_back( event_store->fis_spacepoint[mcp_index].size() - 1); //save index in terms of all hits
			
			auto cluster_ID = cluster_hits.at(cluster_index)->ID();
			//loop on clusters and find PFP
			size_t pfp_id = -1;
			for ( auto const& cluster : clusterList) {
			if ( cluster->ID() != cluster_ID ) continue; //select the cluster I am interested in
			std::vector< art::Ptr <recob::PFParticle> > pfp_cluster = pfpFromCluster.at( &cluster - &clusterList[0] );
			if (pfp_cluster.size()==0) continue;
			if (pfp_cluster.size()!=1) std::cout << "MORE THAN 1 PFP! size=" << pfp_cluster.size() << std::endl;
			else pfp_id = pfp_cluster[0]->Self();
			pfp_cluster.clear();
			}
			
			//associate the pfp to the tracks
			for ( auto const& track : trackList) {
			std::vector< art::Ptr <recob::PFParticle> > pfp_track = pfpFromTrack.at( &track - &trackList[0] );
			std::vector< art::Ptr <simb::MCParticle> > mcp_track = mcpsFromTracks.at( &track - &trackList[0] );
			unsigned mm = 0;
			float purityy = -1;
			for ( unsigned nn=0; nn<mcp_track.size(); nn++ ) { //select highest purity
				if ( mcpsFromTracks.data( &track - &trackList[0] ).at( nn )->cleanliness > purityy ) {
					purityy = mcpsFromTracks.data( &track - &trackList[0] ).at( nn )->cleanliness;
					mm = nn;
				}
			}
			if (pfp_track.size()==0 || mcp_track.size()==0) {
			 	continue;
			}

			if (pfp_track.size()!=1) std::cout << "MORE THAN 1 PFP per track!" << std::endl;
			if ( pfp_track[0]->Self() != pfp_id ) continue;
			//std::cout << "FOUND PFP" << std::endl;
		
			//fpdg[mcp_index] is the particle truth matched to the hit
			//mcp_track[mm] is the particle matched to the track associated to the cluster (to which the hit is clustered)
			if ( !event_store->fis_tracked[mcp_index] && track->ID() == event_store->freco_trackid[mcp_index] )
				throw cet::exception("Line916") << "ERROR! This should never happen..." << std::endl;
			
			if ( track->ID() == event_store->freco_trackid[mcp_index] ) {
					//the hit and track matching agree
				event_store->fclustered_matched[mcp_index] = event_store->fclustered_matched[mcp_index] + 1;
		                event_store->fclustered_matched_charge[mcp_index].push_back( hit->Integral() );
			} else { //mismatch
				if ( event_store->fpdg[mcp_index] == 13 && mcp_track[mm]->PdgCode()==13) {
					std::cout << "Broken track!" << std::endl;
					//std::cout << "trackid1 " << track->ID() << " trackid2 " << freco_trackid[mcp_index] << std::endl; 
					//std::cout << "mcp_index " << mcp_index << " muon_index " << muon_pos << std::endl; 
					//std::cout << "numero track " << &track - &trackList[0] << std::endl;
				} else if (event_store->fpdg[mcp_index] == 2212 && event_store->fp0[mcp_index] < 0.2) {
					std::cout << "Low proton!!!" << std:: endl;
				} else {
				      event_store->fhit_mismatch_pdg[mcp_index].push_back ( mcp_track[mm]->PdgCode() );
				      event_store->fclustered_mismatched[mcp_index] = event_store->fclustered_mismatched[mcp_index] + 1;
		                      event_store->fclustered_mismatched_charge[mcp_index].push_back( hit->Integral() );
				}
			}
			pfp_track.clear();
			mcp_track.clear();
		} //tracklist
	       }//is clustered
	//is it clustered?
	//is it attached to a track?
	//are those the right track and the right cluster?
	//check protons and merged protons
	cluster_hits.clear();
	mcparticle_hits.clear();
	spacepoint_hits.clear();

  } // loop on hits
  
  //make hit plots
  if (fWriteHistograms)
  histo_maker->Fill_Hit_Histos( event_store );

  } // good muon

  recoTree->Fill();
  } //MCtruth

  //clear memory
  pfParticleMap.clear();
  pfpList.clear();
  trackList.clear();
  showerList.clear();
  mcList.clear();
  mcTruth.clear();
  clusterList.clear();
  hitList.clear();
  mcHitList.clear();
  simChannelList.clear();
  

}

void recohelper::ProtonBenchmarker::GetPFParticleIdMap(const art::ValidHandle< std::vector<recob::PFParticle> > &pfParticleHandle, recohelper::ProtonBenchmarker::PFParticleIdMap &pfParticleMap) {
	
	for (unsigned int i = 0; i < pfParticleHandle->size(); ++i){
		const art::Ptr<recob::PFParticle> pParticle(pfParticleHandle, i);
			if (!pfParticleMap.insert(PFParticleIdMap::value_type(pParticle->Self(), pParticle)).second){
				throw cet::exception("ProtonBenchmarker::GetPFParticleIdMap") << "  Unable to get PFParticle ID map, the input PFParticle collection has repeat IDs!";
			}
	}
}


void recohelper::ProtonBenchmarker::GetFinalStatePFParticleVectors(const PFParticleIdMap &pfParticleMap, std::vector< art::Ptr<recob::PFParticle> > &crParticles, std::vector< art::Ptr<recob::PFParticle> > &nuParticles) {
    for (PFParticleIdMap::const_iterator it = pfParticleMap.begin(); it != pfParticleMap.end(); ++it)  {
        const art::Ptr<recob::PFParticle> pParticle(it->second);
        // Only look for primary particles
        if (!pParticle->IsPrimary()) continue;
        // Check if this particle is identified as the neutrino
        const int pdg(pParticle->PdgCode());
        const bool isNeutrino(std::abs(pdg) == pandora::NU_E || std::abs(pdg) == pandora::NU_MU || std::abs(pdg) == pandora::NU_TAU);
        // All non-neutrino primary particles are reconstructed under the cosmic hypothesis
        if (!isNeutrino) {
            crParticles.push_back(pParticle);
            continue;
        }

        // ATTN. We are filling nuParticles under the assumption that there is only one reconstructed neutrino identified per event.
        //       If this is not the case please handle accordingly
        if (!nuParticles.empty()){
            throw cet::exception("ProtonBenchmarker::GetFinalStatePFParticleVectors") << "  This event contains multiple reconstructed neutrinos!";
        }

        // Add the daughters of the neutrino PFParticle to the nuPFParticles vector
        for (const size_t daughterId : pParticle->Daughters()){
            if (pfParticleMap.find(daughterId) == pfParticleMap.end())
                throw cet::exception("ProtonBenchmarker::GetFinalStatePFParticleVectors") << "  Invalid PFParticle collection!";

            nuParticles.push_back(pfParticleMap.at(daughterId));
        }
    }
}

    
void recohelper::ProtonBenchmarker::CollectTracksAndShowers(const std::vector< art::Ptr<recob::PFParticle> > &particles, const art::FindManyP<recob::Track> pfPartToTrackAssoc, const art::FindManyP<recob::Shower> pfPartToShowerAssoc, std::vector< art::Ptr<recob::Track> > &tracks, std::vector< art::Ptr<recob::Shower> > &showers) {
   
    for (const art::Ptr<recob::PFParticle> &pParticle : particles) {
        const std::vector< art::Ptr<recob::Track> > associatedTracks( pfPartToTrackAssoc.at( &pParticle - &particles[0]) );
        const std::vector< art::Ptr<recob::Shower> > associatedShowers( pfPartToShowerAssoc.at( &pParticle - &particles[0] ) );
        const unsigned int nTracks(associatedTracks.size());
        const unsigned int nShowers(associatedShowers.size());

        // Check if the PFParticle has no associated tracks or showers
        if (nTracks == 0 && nShowers == 0)        {
            mf::LogDebug("ConsolidatedPFParticleAnalysisTemplate") << "  No tracks or showers were associated to PFParticle " << pParticle->Self() << std::endl;
            continue;
        }

        // Check if there is an associated track
        if (nTracks == 1 && nShowers == 0){
            tracks.push_back(associatedTracks.front());
            continue;
        }

        // Check if there is an associated shower
        if (nTracks == 0 && nShowers == 1){
            showers.push_back(associatedShowers.front());
            continue;
        }

        throw cet::exception("ConsolidatedPFParticleAnalysisTemplate") << "  There were " << nTracks << " tracks and " << nShowers << " showers associated with PFParticle " << pParticle->Self();
    }
}

    
void recohelper::ProtonBenchmarker::endJob()
{

  if (fWriteHistograms) {
  histo_maker->ScalePlots( n_events );
  histo_maker->FillCumulativeHistograms();
  }

}

DEFINE_ART_MODULE(recohelper::ProtonBenchmarker)
