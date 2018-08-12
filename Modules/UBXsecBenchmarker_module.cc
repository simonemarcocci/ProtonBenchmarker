////////////////////////////////////////////////////////////////////////
// Class:       UBXsecBenchmarker
// Plugin Type: analyzer (art v2_05_01)
// File:        UBXsecBenchmarker_module.cc
//
// Generated at Wed Aug  8 14:20:05 2018 by Simone Marcocci using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"
#include <stdexcept>


#include "uboone/UBXSec/DataTypes/SelectionResult.h"


namespace recobenchmarker {
  class UBXsecBenchmarker;
}


class recobenchmarker::UBXsecBenchmarker : public art::EDAnalyzer {
public:
  explicit UBXsecBenchmarker(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  UBXsecBenchmarker(UBXsecBenchmarker const &) = delete;
  UBXsecBenchmarker(UBXsecBenchmarker &&) = delete;
  UBXsecBenchmarker & operator = (UBXsecBenchmarker const &) = delete;
  UBXsecBenchmarker & operator = (UBXsecBenchmarker &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  // Declare member data here.

};


recobenchmarker::UBXsecBenchmarker::UBXsecBenchmarker(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
/*  fTrackLabel = p.get<std::string> ("TrackLabel");
  fPfpLabel = p.get<std::string> ("PfpLabel");
  fPfpAssnLabel = p.get<std::string> ("PfpAssnLabel");
  fCalorimetryLabel = p.get<std::string> ("CalorimetryLabel");
  fTrackTruthLabel = p.get<std::string> ("TrackTruthLabel");
  fHitAssnTruthLabel = p.get<std::string> ("HitAssnTruthLabel");
  fShowerLabel = p.get<std::string> ("ShowerLabel");
  fVertexFitterLabel = p.get<std::string> ("VertexFitterLabel");
  fIsVertexFitter = p.get<bool> ("IsVertexFitter");
  fMCTruthLabel = p.get<std::string> ("MCTruthLabel");
  fG4TruthLabel = p.get<std::string> ("G4TruthLabel");
  fShowerTruthLabel = p.get<std::string> ("ShowerTruthLabel");
  fClusterLabel = p.get<std::string> ("ClusterLabel");
  fHitLabel = p.get<std::string> ("HitLabel");
  fMCHitLabel = p.get<std::string> ("MCHitLabel");
  */

}

void recobenchmarker::UBXsecBenchmarker::analyze(art::Event const & e)
{
  // Implementation of required member function here.
 	art::Handle<std::vector<ubana::SelectionResult> > selection_h;
	e.getByLabel("UBXSec", selection_h);
	if(!selection_h.isValid()){
		  mf::LogError(__PRETTY_FUNCTION__) << "SelectionResult product not found." << std::endl;
		    throw cet::exception("SkipEvent");
	}

	std::vector<art::Ptr<ubana::SelectionResult>> selection_v;
	art::fill_ptr_vector(selection_v, selection_h);
	// The selection result vector will always only contain one entry (at least 1 neutrino per event)
	
	


}

DEFINE_ART_MODULE(recobenchmarker::UBXsecBenchmarker)
