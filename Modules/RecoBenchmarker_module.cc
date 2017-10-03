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
#include "TH2D.h"
#include "TH1D.h"

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
    std::string fTrackTruthLabel;
    std::string fShowerLabel;
    std::string fShowerTruthLabel;

    bool isDebug = false;

    rbutil::recoBenchmarkerUtility _rbutilInstance;

    TTree* recoTree;

    // auxiliary 

    int fEvent;
    int fRun;
    int fSubRun;

    // particle counters

    int nParticlesAboveThreshold;
    int nProtonsAboveThreshold;
    int nMuons;
    int nPions;
    int nElectrons;

    // matching information
    std::vector<int> isRecoTrackTruthMatched;
    std::vector<int> isRecoShowerTruthMatched;

    // length information
    std::vector<double> thisRecoLength;
    std::vector<double> thisMcpLength;

    // angular information

    std::vector<double> thisMcpMomentum;
    std::vector<double> nextMcpMomentum;
    std::vector<double> thisRecoMomentum;
    std::vector<double> nextRecoMomentum;

    std::vector<float> thisNextMcpAngles;
    float thisNextMcpAngle;
    std::vector<float> thisNextMcpAnglesXZ;
    float thisNextMcpAngleXZ;
    std::vector<float> thisNextRecoAngles;
    float thisNextRecoAngle;
    std::vector<float> thisNextRecoAnglesXZ;
    float thisNextRecoAngleXZ;

    bool isMuPEvent;
    TH2D* mupAngleVersusPLength;
    TH2D* matchedShowerPdgVersusCleanliness;

    TH1D* showerCleanlinessPrimaryProton;
    TH1D* showerCleanlinessPrimaryMuonOrPion;
    TH1D* showerCleanlinessPrimaryElectron;
    TH1D* showerCleanlinessConv;
    TH1D* showerCleanlinessInelastic;
    TH1D* showerCleanlinessMuIoni;
    TH1D* showerCleanlinessOther;
};


recohelper::RecoBenchmarker::RecoBenchmarker(fhicl::ParameterSet const & p)
  :
    EDAnalyzer(p)  // ,
    // More initializers here.
{

  fTrackLabel = p.get<std::string> ("TrackLabel");
  fTrackTruthLabel = p.get<std::string> ("TrackTruthLabel");
  fShowerLabel = p.get<std::string> ("ShowerLabel");
  fShowerTruthLabel = p.get<std::string> ("ShowerTruthLabel");

}

void recohelper::RecoBenchmarker::beginJob()
{

  art::ServiceHandle< art::TFileService > tfs;

  recoTree = tfs->make<TTree>("recotree", "recotree");

  // define branches
  recoTree->Branch("Event", &fEvent, "Event/I");
  recoTree->Branch("SubRun", &fSubRun, "SubRun/I");
  recoTree->Branch("Run", &fRun, "Run/I");

  recoTree->Branch("isRecoTrackTruthMatched", &isRecoTrackTruthMatched);
  recoTree->Branch("isRecoShowerTruthMatched", &isRecoShowerTruthMatched);

  recoTree->Branch("thisNextRecoAngles", &thisNextRecoAngles);
  recoTree->Branch("thisNextRecoAnglesXZ", &thisNextRecoAnglesXZ);
  recoTree->Branch("thisRecoLength", &thisRecoLength);

  recoTree->Branch("thisNextMcpAngles", &thisNextMcpAngles);
  recoTree->Branch("thisNextMcpAnglesXZ", &thisNextMcpAnglesXZ);
  recoTree->Branch("thisMcpLength", &thisMcpLength);

  mupAngleVersusPLength = tfs->make<TH2D>("mupAngleVersusPLength", ";#theta_{#mup};length_{p}", 50, 0, 180, 50, 0, 10);

  matchedShowerPdgVersusCleanliness = tfs->make<TH2D>("matchedShowerPdgVersusCleanliness", ";Pdg code; Cleanliness", 3000, 0, 3000, 1, 0, 1);

  showerCleanlinessPrimaryProton = tfs->make<TH1D>("showerCleanlinessPrimaryProton", ";cleanliness;", 20, 0, 1);
  showerCleanlinessPrimaryMuonOrPion = tfs->make<TH1D>("showerCleanlinessPrimaryMuonOrPion", ";cleanliness;", 20, 0, 1);
  showerCleanlinessPrimaryElectron = tfs->make<TH1D>("showerCleanlinessPrimaryElectron", ";cleanliness;", 20, 0, 1);
  showerCleanlinessConv = tfs->make<TH1D>("showerCleanlinessConv", ";cleanliness;", 20, 0, 1);
  showerCleanlinessInelastic = tfs->make<TH1D>("showerCleanlinessInelastic", ";cleanliness;", 20, 0, 1);
  showerCleanlinessMuIoni = tfs->make<TH1D>("showerCleanlinessMuIoni", ";cleanliness;", 20, 0, 1);
  showerCleanlinessOther = tfs->make<TH1D>("showerCleanlinessOther", ";cleanliness;", 20, 0, 1);

}

void recohelper::RecoBenchmarker::analyze(art::Event const & e)
{

  if (e.isRealData() == true) return;

  fEvent  = e.id().event();
  fRun    = e.run();
  fSubRun = e.subRun();

  // setup variables

  isRecoTrackTruthMatched.clear(); 
  isRecoShowerTruthMatched.clear();

  nParticlesAboveThreshold = 0;
  nProtonsAboveThreshold = 0;
  nMuons = 0;
  nPions = 0;
  nElectrons = 0;
  thisNextMcpAngles.clear();
  thisNextMcpAnglesXZ.clear();
  thisMcpMomentum.clear();
  nextMcpMomentum.clear();
  thisMcpLength.clear();

  thisNextRecoAngles.clear();
  thisNextRecoAnglesXZ.clear();
  thisRecoLength.clear();
  thisRecoMomentum.clear();
  nextRecoMomentum.clear();

  isMuPEvent = false;

  // get handles to objects of interest

  art::ValidHandle< std::vector<recob::Track> > trackHandle = 
    e.getValidHandle< std::vector<recob::Track> >(fTrackLabel);

  art::FindManyP<simb::MCParticle, anab::BackTrackerMatchingData> mcpsFromTracks(trackHandle, e, fTrackTruthLabel);

  art::ValidHandle< std::vector<recob::Shower> > showerHandle =
    e.getValidHandle< std::vector<recob::Shower> >(fShowerLabel);

  art::FindManyP<simb::MCParticle, anab::BackTrackerMatchingData> mcpsFromShowers(showerHandle, e, fShowerTruthLabel);

  art::Handle< std::vector<simb::MCParticle> > mcParticleHandle; 
  std::vector< art::Ptr<simb::MCParticle> > mcList;
  if (e.getByLabel("largeant", mcParticleHandle))
    art::fill_ptr_vector(mcList, mcParticleHandle); 

  for (size_t i = 0; i < mcList.size(); i++){

    const art::Ptr<simb::MCParticle>& thisMcp = mcList.at(i);
    if ((thisMcp->Process() != "primary")
        || (std::abs(thisMcp->PdgCode()) == 2112)
        || (std::abs(thisMcp->PdgCode()) == 14)
        || (std::abs(thisMcp->PdgCode()) == 12)
        || (std::abs(thisMcp->PdgCode()) == 22)
        || (std::abs(thisMcp->PdgCode()) == 111)
        || (std::abs(thisMcp->PdgCode()) == 2212 && thisMcp->P() < 0.2)) continue;

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

    // update counter information
    if (std::abs(thisMcp->PdgCode()) == 2212)
      nProtonsAboveThreshold++;
    if (std::abs(thisMcp->PdgCode()) == 13)
      nMuons++;
    if (std::abs(thisMcp->PdgCode()) == 211)
      nPions++;
    if (std::abs(thisMcp->PdgCode()) == 11)
      nElectrons++;

    // angular stuff
    thisMcpMomentum.clear();
    thisMcpMomentum = _rbutilInstance.getMomentumVector(thisMcp);
    thisMcpLength.push_back(thisMcp->Trajectory().TotalLength());

    for (size_t j = i+1; j < mcList.size(); j++){

      const art::Ptr<simb::MCParticle>& nextMcp = mcList.at(j); 
      if ((nextMcp->Process() != "primary") 
          ||(std::abs(nextMcp->PdgCode()) == 2212 && nextMcp->P() < 0.2)) continue;

      nextMcpMomentum.clear();
      nextMcpMomentum = _rbutilInstance.getMomentumVector(nextMcp);

      thisNextMcpAngle = _rbutilInstance.getAngle(thisMcpMomentum, nextMcpMomentum, _rbutilInstance, "no");
      thisNextMcpAngles.push_back(thisNextMcpAngle);

      thisNextMcpAngleXZ = _rbutilInstance.getAngle(thisMcpMomentum, nextMcpMomentum, _rbutilInstance, "xz");
      thisNextMcpAnglesXZ.push_back(thisNextMcpAngleXZ);

    }

  } // MCParticles

  nParticlesAboveThreshold = nProtonsAboveThreshold + nMuons + nPions + nElectrons;

  if (nProtonsAboveThreshold == 1 && nMuons == 1 && nParticlesAboveThreshold == 2){
    std::cout << ">> Found classic CCQE signature"
      << "\n>> Saving for angle analysis" << std::endl;
    isMuPEvent = true;
  }

  int it = 0;
  for (auto const& thisTrack : (*trackHandle)) {

    thisRecoLength.push_back(thisTrack.Length());

    thisRecoMomentum.clear();
    thisRecoMomentum = _rbutilInstance.getMomentumVector(thisTrack);

    // pick up all other reconstructed particles in the event for angular studies

    for (size_t i = 1; i < trackHandle->size() - it; i++){ 
      auto const& nextTrack = *(&thisTrack + i);

      nextRecoMomentum.clear();
      nextRecoMomentum = _rbutilInstance.getMomentumVector(nextTrack);

      thisNextRecoAngle = _rbutilInstance.getAngle(thisRecoMomentum, nextRecoMomentum, _rbutilInstance, "no");
      thisNextRecoAngles.push_back(thisNextRecoAngle);

      thisNextRecoAngleXZ = _rbutilInstance.getAngle(thisRecoMomentum, nextRecoMomentum, _rbutilInstance, "xz");
      thisNextRecoAnglesXZ.push_back(thisNextRecoAngleXZ);

    }

    std::cout << "Found Track!" << std::endl;

    std::vector< art::Ptr<simb::MCParticle> > mcps = mcpsFromTracks.at(thisTrack.ID());

    isRecoTrackTruthMatched.push_back(mcps.size());

    for (auto const& thisMcp: mcps){
    
      auto thisMcpCleanliness = mcpsFromTracks.data(0).at(0)->cleanliness;
      auto thisMcpCompleteness = mcpsFromTracks.data(0).at(0)->completeness;

      std::cout << "MATCH CLEANLINESS: " << thisMcpCleanliness << std::endl;
      std::cout << "MATCH COMPLETENESS: " << thisMcpCompleteness << std::endl;

      std::cout << "---- MCParticle Information ----"
        << "\nTrack ID   : " << thisMcp->TrackId() 
        << "\nPdgCode    : " << thisMcp->PdgCode()
        << "\nProcess    : " << thisMcp->Process()
        << "\nStatusCode : " << thisMcp->StatusCode()
        << "\nMother Pdg : " << thisMcp->Mother()
        << "\nPx, Py, Pz : " << thisMcp->Px() << ", " << thisMcp->Py() << ", " << thisMcp->Pz()
        << std::endl;

    }
    it++;
  } // trackHandle

  std::sort(thisRecoLength.begin(), thisRecoLength.end());

  if (thisNextRecoAnglesXZ.size() == 1 && isMuPEvent == true){

    std::cout << ">> Filling histogram with "
      << "\n>>>> Angle:  " << thisNextRecoAnglesXZ.at(0)
      << "\n>>>> Length: " << thisRecoLength.at(0) 
      << std::endl;
    mupAngleVersusPLength->Fill(thisNextRecoAnglesXZ.at(0), thisRecoLength.at(0));

  }

  for (auto const& thisShower : (*showerHandle)){

    std::cout << "Found Shower!" << std::endl;

    std::vector< art::Ptr<simb::MCParticle> > mcps = mcpsFromShowers.at(thisShower.ID());

    isRecoShowerTruthMatched.push_back(mcps.size());

    for (auto const& thisMcp: mcps){
    
      auto thisMcpCleanliness = mcpsFromShowers.data(0).at(0)->cleanliness;
      auto thisMcpCompleteness = mcpsFromShowers.data(0).at(0)->completeness;

      std::cout << "MATCH CLEANLINESS: " << thisMcpCleanliness << std::endl;
      std::cout << "MATCH COMPLETENESS: " << thisMcpCompleteness << std::endl;
      std::cout << "---- MCParticle Information ----"
        << "\nTrack ID   : " << thisMcp->TrackId() 
        << "\nPdgCode    : " << thisMcp->PdgCode()
        << "\nProcess    : " << thisMcp->Process()
        << "\nStatusCode : " << thisMcp->StatusCode()
        << "\nMother Pdg : " << thisMcp->Mother()
        << "\nPx, Py, Pz : " << thisMcp->Px() << ", " << thisMcp->Py() << ", " << thisMcp->Pz()
        << std::endl;
   
      matchedShowerPdgVersusCleanliness->Fill(std::abs(thisMcp->PdgCode()), thisMcpCleanliness);

      if (thisMcp->Process() == "primary" && std::abs(thisMcp->PdgCode()) == 2212)
        showerCleanlinessPrimaryProton->Fill(thisMcpCleanliness);
      else if (thisMcp->Process() == "primary" 
          && ((std::abs(thisMcp->PdgCode()) == 13) 
              || (std::abs(thisMcp->PdgCode()) == 211)))
        showerCleanlinessPrimaryMuonOrPion->Fill(thisMcpCleanliness);
      else if (thisMcp->Process() == "primary" && std::abs(thisMcp->PdgCode()) == 11)
        showerCleanlinessPrimaryElectron->Fill(thisMcpCleanliness);
      else if (thisMcp->Process() == "conv")
        showerCleanlinessConv->Fill(thisMcpCleanliness);
      else if ((thisMcp->Process() == "neutronInelastic")
          || (thisMcp->Process() == "protonInelastic")
          || (thisMcp->Process() == "pi+Inelastic"))
        showerCleanlinessInelastic->Fill(thisMcpCleanliness);
      else if (thisMcp->Process() == "muIoni")
        showerCleanlinessMuIoni->Fill(thisMcpCleanliness);
      else showerCleanlinessOther->Fill(thisMcpCleanliness);

    
    }

  } // showerHandle

  recoTree->Fill();
}

void recohelper::RecoBenchmarker::endJob()
{
  // Implementation of optional member function here.
}

DEFINE_ART_MODULE(recohelper::RecoBenchmarker)
