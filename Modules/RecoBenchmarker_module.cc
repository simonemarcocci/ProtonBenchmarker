////////////////////////////////////////////////////////////////////////
// Class:       RecoBenchmarker
// Plugin Type: analyzer (art v2_05_00)
// File:        RecoBenchmarker_module.cc
//
// Generated at Tue Sep 19 17:11:42 2017 by Adam Lister using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

// base includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

// LArSoft includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

// ROOT includes
#include "TTree.h"

// local includes
#include "uboone/RecoBenchmarker/Algos/recoBenchmarkerUtility.h"

namespace recohelper {
  class RecoBenchmarker;
}


class recohelper::RecoBenchmarker : public art::EDAnalyzer {
  public:
    explicit RecoBenchmarker(fhicl::ParameterSet const & p);
    // The compiler-generated destructor is fine for non-base
    // classes without bare pointers or other resource use.

    // Plugins should not be copied or assigned.
    RecoBenchmarker(RecoBenchmarker const &) = delete;
    RecoBenchmarker(RecoBenchmarker &&) = delete;
    RecoBenchmarker & operator = (RecoBenchmarker const &) = delete;
    RecoBenchmarker & operator = (RecoBenchmarker &&) = delete;

    // Required functions.
    void analyze(art::Event const & e) override;

    // Selected optional functions.
    void beginJob() override;
    void endJob() override;

  private:
    // fcl input parameters
    std::string fTrackLabel;
    std::string fShowerLabel;

    bool isDebug = false;

    rbutil::recoBenchmarkerUtility _rbutilInstance;

    TLorentzVector thisMcpMomentum;
    TLorentzVector nextMcpMomentum;
   
    TTree* recoTree;

    int fEvent;
    int fRun;
    int fSubRun;

    std::vector<float> thisNextMcpAngles;
    float thisNextMcpAngle;
    std::vector<float> thisNextMcpAnglesXZ;
    float thisNextMcpAngleXZ;
    

};


recohelper::RecoBenchmarker::RecoBenchmarker(fhicl::ParameterSet const & p)
  :
    EDAnalyzer(p)  // ,
    // More initializers here.
{

  fTrackLabel = p.get<std::string> ("TrackLabel");
  fShowerLabel = p.get<std::string> ("ShowerLabel");

}

void recohelper::RecoBenchmarker::beginJob()
{

  art::ServiceHandle< art::TFileService > tfs;

  recoTree = tfs->make<TTree>("recotree", "recotree");

  // define branches
  recoTree->Branch("Event", &fEvent, "Event/I");
  recoTree->Branch("SubRun", &fSubRun, "SubRun/I");
  recoTree->Branch("Run", &fRun, "Run/I");
  recoTree->Branch("thisNextMcpAngles", &thisNextMcpAngles);
  recoTree->Branch("thisNextMcpAnglesXZ", &thisNextMcpAnglesXZ);

}

void recohelper::RecoBenchmarker::analyze(art::Event const & e)
{

  if (e.isRealData() == true) return;

  fEvent  = e.id().event();
  fRun    = e.run();
  fSubRun = e.subRun();

  // get handles to objects of interest

  art::ValidHandle< std::vector<recob::Track> > trackHandle = 
    e.getValidHandle< std::vector<recob::Track> >(fTrackLabel);

  art::FindManyP<simb::MCParticle, anab::BackTrackerMatchingData> mcpsFromTracks(trackHandle, e, "pandoraNuTruthMatch");

  art::ValidHandle< std::vector<recob::Shower> > showerHandle =
    e.getValidHandle< std::vector<recob::Shower> >(fShowerLabel);

  art::Handle< std::vector<simb::MCParticle> > mcParticleHandle; 
  std::vector< art::Ptr<simb::MCParticle> > mcList;
  if (e.getByLabel("largeant", mcParticleHandle))
      art::fill_ptr_vector(mcList, mcParticleHandle); 

  for (auto const& track : (*trackHandle)) {

    std::cout << "I found a track with length " << track.Length() << std::endl;

    std::vector< art::Ptr<simb::MCParticle> > mcps = mcpsFromTracks.at(track.ID());
    for (auto const& mcp: mcps){

      std::cout << "matched mcp with ID " << mcp.get()->TrackId() << std::endl;
    
    }
  }

  for (auto const& shower : (*showerHandle)){

    std::cout << "I found a shower with length " << shower.Length() << std::endl;

  }

  thisNextMcpAngles.clear();
  thisNextMcpAnglesXZ.clear();

  for (size_t i = 0; i < mcList.size(); i++){
    
    const art::Ptr<simb::MCParticle>& thisMcp = mcList.at(i);
    if (thisMcp->Process() != "primary") continue;

    if (isDebug == true){
      std::cout << "---- MCParticle Information ----"
        << "\nTrack ID   : " << thisMcp->TrackId() 
        << "\nPdgCode    : " << thisMcp->PdgCode()
        << "\nProcess    : " << thisMcp->Process()
        << "\nStatusCode : " << thisMcp->StatusCode()
        << "\nMother Pdg : " << thisMcp->Mother()
        << "\nPx, Py, Pz : " << thisMcp->Px() << ", " << thisMcp->Py() << ", " << thisMcp->Pz()
        << std::endl;
    }

    thisMcpMomentum = thisMcp->Momentum();

    for (size_t j = i+1; j < mcList.size(); j++){

      const art::Ptr<simb::MCParticle>& nextMcp = mcList.at(j); 
      if (nextMcp->Process() != "primary") continue;

      nextMcpMomentum = nextMcp->Momentum();

      thisNextMcpAngle = _rbutilInstance.getAngle(thisMcp, nextMcp, _rbutilInstance, "no");
      thisNextMcpAngles.push_back(thisNextMcpAngle);
      
      thisNextMcpAngleXZ = _rbutilInstance.getAngle(thisMcp, nextMcp, _rbutilInstance, "xz");
      thisNextMcpAnglesXZ.push_back(thisNextMcpAngleXZ);

    }


  }

  recoTree->Fill();
}

void recohelper::RecoBenchmarker::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(recohelper::RecoBenchmarker)
