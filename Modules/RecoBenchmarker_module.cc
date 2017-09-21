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
    void reconfigure(fhicl::ParameterSet const & p) override;

  private:
    // fcl input parameters
    std::string fTrackLabel;
 //   std::string fShowerLabel;

};


recohelper::RecoBenchmarker::RecoBenchmarker(fhicl::ParameterSet const & p)
  :
    EDAnalyzer(p)  // ,
    // More initializers here.
{}

void recohelper::RecoBenchmarker::reconfigure(fhicl::ParameterSet const & p)
{

  fTrackLabel = p.get<std::string> ("TrackLabel");
 // fShowerLabel = p.get<std::string> ("ShowerLabel");

}

void recohelper::RecoBenchmarker::beginJob()
{
  // Implementation of optional member function here.
}

void recohelper::RecoBenchmarker::analyze(art::Event const & e)
{

  if (e.isRealData() == false) return;

  art::ValidHandle< std::vector<recob::Track> > trackHandle = 
    e.getValidHandle< std::vector< recob::Track> >(fTrackLabel);
/*
  art::Handle< std::vector<recob::Shower> > showerHandle =
    e.getValidHandle< std::vector< recob::Shower> >(fShowerLabel);

  art::Handle< std::vector<simb::MCParticle> > mcParticleHandle =
    e.getValidHandle< std::vector< simb::MCParticle> >("largeant"); 
*/
  for (auto const& track : (*trackHandle)) {

    std::cout << "I found a track with length << " << track.Length() << std::endl;

  }

}

void recohelper::RecoBenchmarker::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(recohelper::RecoBenchmarker)
