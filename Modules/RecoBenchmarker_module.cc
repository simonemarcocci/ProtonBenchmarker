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
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// LArSoft includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Shower.h"
#include "nusimdata/SimulationBase/MCParticle.h"

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
   
    float thisNextMcpAngle;
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
  // Implementation of optional member function here.
}

void recohelper::RecoBenchmarker::analyze(art::Event const & e)
{

  if (e.isRealData() == true) return;

  // get handles to objects of interest

  art::ValidHandle< std::vector<recob::Track> > trackHandle = 
    e.getValidHandle< std::vector<recob::Track> >(fTrackLabel);

  art::ValidHandle< std::vector<recob::Shower> > showerHandle =
    e.getValidHandle< std::vector<recob::Shower> >(fShowerLabel);

  art::Handle< std::vector<simb::MCParticle> > mcParticleHandle; 
  std::vector< art::Ptr<simb::MCParticle> > mcList;
  if (e.getByLabel("largeant", mcParticleHandle))
      art::fill_ptr_vector(mcList, mcParticleHandle); 

  for (auto const& track : (*trackHandle)) {

    std::cout << "I found a track with length " << track.Length() << std::endl;

  }

  for (auto const& shower : (*showerHandle)){

    std::cout << "I found a shower with length " << shower.Length() << std::endl;

  }

  for (size_t i = 0; i < mcList.size(); i++){
    
    auto const & thisMcp = mcList.at(i);
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

      auto const & nextMcp = mcList.at(j); 
      
      nextMcpMomentum = nextMcp->Momentum();

      thisNextMcpAngle = _rbutilInstance.getAngle(thisMcp, nextMcp);

      std::cout << thisNextMcpAngle << std::endl;
    }

  }

}

void recohelper::RecoBenchmarker::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(recohelper::RecoBenchmarker)
