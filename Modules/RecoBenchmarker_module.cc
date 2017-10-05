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
#include "TFile.h"
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

    art::ServiceHandle< art::TFileService > tfs;
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

    int nimID;

    std::vector<int> matchIDChecker;

    // matching information
    std::vector<int> isRecoTrackTruthMatched;
    std::vector<int> isRecoShowerTruthMatched;

    // length information
    std::vector<double> thisRecoLength;
    std::vector<double> thisMcpLength;
    float thisMatchedLength;

    // angular information

    std::vector<double> thisMcpMomentum;
    std::vector<double> nextMcpMomentum;
    std::vector<double> nimMcpMomentum;
    std::vector<double> thisMatchedMomentum;
    std::vector<double> nextRecoMomentum;
    std::vector<double> nimMatchedMomentum;

    std::vector<float> thisNimMcpAngles;
    float thisNimMcpAngle;
    std::vector<float> thisNimMcpAnglesXZ;
    float thisNimMcpAngleXZ;
    std::vector<float> thisNimMatchedAngles;
    float thisNimMatchedAngle;
    std::vector<float> thisNimMatchedAnglesXZ;
    float thisNimMatchedAngleXZ;
    float thisZmatchedAngleYZ;
    float thisZMcpAngleYZ;

    TH2D* matchedLengthVersusAngle;
    TH2D* mcpLengthVersusAngle;
    TH2D* lengthVersusAngleEfficiency;

    TH1D* matchedProjectedLength;
    TH1D* mcpProjectedLength;
    TH1D* trackLengthEfficiency;

    TH1D* matchedProjectedAngle;
    TH1D* mcpProjectedAngle;
    TH1D* trackAngleEfficiency;

    TH2D* mupMatchedAngleVersusPMom;
    TH2D* mupMcpAngleVersusPMom;
    TH2D* mupAngleMomentumEfficiency;

    TH1D* pMatchedProjectedMomentum;
    TH1D* pMcpProjectedMomentum;
    TH1D* pMomentumEfficiency;

    TH1D* pMatchedProjectedAngle;
    TH1D* pMcpProjectedAngle;
    TH1D* pAngleEfficiency;

    TH1D* showerCleanlinessPrimaryProton;
    TH1D* showerCleanlinessPrimaryMuonOrPion;
    TH1D* showerCleanlinessPrimaryElectron;
    TH1D* showerCleanlinessConv;
    TH1D* showerCleanlinessInelastic;
    TH1D* showerCleanlinessMuIoni;
    TH1D* showerCleanlinessOther;

    TH1D* trackCleanlinessPrimaryProton;
    TH1D* trackCleanlinessPrimaryMuonOrPion;
    TH1D* trackCleanlinessPrimaryElectron;
    TH1D* trackCleanlinessConv;
    TH1D* trackCleanlinessInelastic;
    TH1D* trackCleanlinessMuIoni;
    TH1D* trackCleanlinessOther;

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

  recoTree = tfs->make<TTree>("recotree", "recotree");

  // define branches
  recoTree->Branch("Event", &fEvent, "Event/I");
  recoTree->Branch("SubRun", &fSubRun, "SubRun/I");
  recoTree->Branch("Run", &fRun, "Run/I");

  recoTree->Branch("isRecoTrackTruthMatched", &isRecoTrackTruthMatched);
  recoTree->Branch("isRecoShowerTruthMatched", &isRecoShowerTruthMatched);

  recoTree->Branch("thisNimMatchedAngles", &thisNimMatchedAngles);
  recoTree->Branch("thisNimMatchedAnglesXZ", &thisNimMatchedAnglesXZ);
  recoTree->Branch("thisRecoLength", &thisRecoLength);

  recoTree->Branch("thisNimMcpAngles", &thisNimMcpAngles);
  recoTree->Branch("thisNimMcpAnglesXZ", &thisNimMcpAnglesXZ);
  recoTree->Branch("thisMcpLength", &thisMcpLength);

  matchedLengthVersusAngle = tfs->make<TH2D>("matchedLengthVersusAngle", ";#theta_{YZ} (degrees); Length (cm)", 20, 0, 180, 10, 0, 10);
  mcpLengthVersusAngle = tfs->make<TH2D>("mcpLengthVersusAngle", ";#theta_{YZ} (degrees); Length (cm)", 20, 0, 180, 10, 0, 10);

  mupMatchedAngleVersusPMom = tfs->make<TH2D>("mupMatchedAngleVersusPMom", ";#theta^{#mup}_{XZ} (degrees) ;P_{p} (Gev)", 20, 0, 180, 20, 0, 1);
  mupMcpAngleVersusPMom = tfs->make<TH2D>("mupMcpAngleVersusPMom", ";#theta^{#mup}_{XZ};P_{p}", 20, 0, 180, 20, 0, 1);

  // cleanliness plots
  showerCleanlinessPrimaryProton = tfs->make<TH1D>("showerCleanlinessPrimaryProton", ";cleanliness;", 20, 0, 1);
  showerCleanlinessPrimaryMuonOrPion = tfs->make<TH1D>("showerCleanlinessPrimaryMuonOrPion", ";cleanliness;", 20, 0, 1);
  showerCleanlinessPrimaryElectron = tfs->make<TH1D>("showerCleanlinessPrimaryElectron", ";cleanliness;", 20, 0, 1);
  showerCleanlinessConv = tfs->make<TH1D>("showerCleanlinessConv", ";cleanliness;", 20, 0, 1);
  showerCleanlinessInelastic = tfs->make<TH1D>("showerCleanlinessInelastic", ";cleanliness;", 20, 0, 1);
  showerCleanlinessMuIoni = tfs->make<TH1D>("showerCleanlinessMuIoni", ";cleanliness;", 20, 0, 1);
  showerCleanlinessOther = tfs->make<TH1D>("showerCleanlinessOther", ";cleanliness;", 20, 0, 1);

  trackCleanlinessPrimaryProton = tfs->make<TH1D>("trackCleanlinessPrimaryProton", ";cleanliness;", 20, 0, 1);
  trackCleanlinessPrimaryMuonOrPion = tfs->make<TH1D>("trackCleanlinessPrimaryMuonOrPion", ";cleanliness;", 20, 0, 1);
  trackCleanlinessPrimaryElectron = tfs->make<TH1D>("trackCleanlinessPrimaryElectron", ";cleanliness;", 20, 0, 1);
  trackCleanlinessConv = tfs->make<TH1D>("trackCleanlinessConv", ";cleanliness;", 20, 0, 1);
  trackCleanlinessInelastic = tfs->make<TH1D>("trackCleanlinessInelastic", ";cleanliness;", 20, 0, 1);
  trackCleanlinessMuIoni = tfs->make<TH1D>("trackCleanlinessMuIoni", ";cleanliness;", 20, 0, 1);
  trackCleanlinessOther = tfs->make<TH1D>("trackCleanlinessOther", ";cleanliness;", 20, 0, 1);


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

  matchIDChecker.clear();

  nParticlesAboveThreshold = 0;
  nProtonsAboveThreshold = 0;
  nMuons = 0;
  nPions = 0;
  nElectrons = 0;
  thisNimMcpAngles.clear();
  thisNimMcpAnglesXZ.clear();
  thisMcpMomentum.clear();
  nimMcpMomentum.clear();
  thisMcpLength.clear();

  thisNimMatchedAngles.clear();
  thisNimMatchedAnglesXZ.clear();
  nimMatchedMomentum.clear();
  thisRecoLength.clear();

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

  nimID = -1;
  // loop mcps to find neutrino induced muon id
  for (size_t i = 0; i < mcList.size();i++){

    const art::Ptr< simb::MCParticle >& thisMcp = mcList.at(i);
    if (thisMcp->PdgCode() == 13 && thisMcp->Process() == "primary") {
      nimID = thisMcp->TrackId();
      nimMcpMomentum = _rbutilInstance.getMomentumVector(thisMcp);
    }
  }

  for (size_t i = 0; i < mcList.size(); i++){

    const art::Ptr<simb::MCParticle>& thisMcp = mcList.at(i);
    if ((thisMcp->Process() != "primary")
        || (std::abs(thisMcp->PdgCode()) == 2112)
        || (std::abs(thisMcp->PdgCode()) == 14)
        || (std::abs(thisMcp->PdgCode()) == 12)
        || (std::abs(thisMcp->PdgCode()) == 22)
        || (std::abs(thisMcp->PdgCode()) == 111)
        || (std::abs(thisMcp->PdgCode()) == 11) // remove anything which is shower like for now
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

    thisMcpLength.push_back(thisMcp->Trajectory().TotalLength());
    thisMcpMomentum = _rbutilInstance.getMomentumVector(thisMcp);

    // all tracks
    if (thisMcp->Process() == "primary"){

      std::vector<double> zDir = {0,0,1};

      thisZMcpAngleYZ = _rbutilInstance.getAngle(zDir, thisMcpMomentum, _rbutilInstance, "yz");
      mcpLengthVersusAngle->Fill(thisZMcpAngleYZ, thisMcp->Trajectory().TotalLength());

    }

    // specifically protons
    if ((thisMcp->PdgCode() == 2212) && (thisMcp->Process() == "primary") && (nimID !=-1)){

      thisNimMcpAngle = 
        _rbutilInstance.getAngle(thisMcpMomentum, nimMcpMomentum, _rbutilInstance, "no");
      thisNimMcpAngles.push_back(thisNimMcpAngle);

      thisNimMcpAngleXZ =
        _rbutilInstance.getAngle(thisMcpMomentum, nimMcpMomentum, _rbutilInstance, "xz");
      thisNimMcpAnglesXZ.push_back(thisNimMcpAngleXZ);

      mupMcpAngleVersusPMom->Fill(thisNimMcpAngleXZ, thisMcp->P());

    }

  } // MCParticles

  nParticlesAboveThreshold = nProtonsAboveThreshold + nMuons + nPions + nElectrons;

  // loop tracks and do truth matching to find neutrino induced muon
  nimID = -1;
  for (auto const& thisTrack : (*trackHandle)) {

    std::vector< art::Ptr<simb::MCParticle> > mcps = mcpsFromTracks.at(thisTrack.ID());

    for (auto const& thisMcp : mcps){

      if (thisMcp->PdgCode() == 13 && thisMcp->Process() == "primary"){
        nimID = thisMcp->TrackId();
        nimMatchedMomentum = _rbutilInstance.getMomentumVector(thisMcp);

      }

    }

  }

  int it = 0;
  for (auto const& thisTrack : (*trackHandle)) {

    thisRecoLength.push_back(thisTrack.Length());

    std::vector< art::Ptr<simb::MCParticle> > mcps = mcpsFromTracks.at(it);

    isRecoTrackTruthMatched.push_back(mcps.size());

    bool isMatched = false;
    for (auto const& thisMcp : mcps){

      for (size_t i = 0; i < matchIDChecker.size(); i++){

        if (thisMcp->TrackId() == matchIDChecker.at(i))
          isMatched = true;

      }
      matchIDChecker.push_back(thisMcp->TrackId());

      if (isMatched == true) continue;

      if ((thisMcp->Process() != "primary")
          || (std::abs(thisMcp->PdgCode()) == 2112)
          || (std::abs(thisMcp->PdgCode()) == 14)
          || (std::abs(thisMcp->PdgCode()) == 12)
          || (std::abs(thisMcp->PdgCode()) == 22)
          || (std::abs(thisMcp->PdgCode()) == 111)
          || (std::abs(thisMcp->PdgCode()) == 11) // remove anything which is shower like for now
          || (std::abs(thisMcp->PdgCode()) == 2212 && thisMcp->P() < 0.2)) continue;

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

      if (thisMcp->Process() == "primary" && std::abs(thisMcp->PdgCode()) == 2212)
        trackCleanlinessPrimaryProton->Fill(thisMcpCleanliness);
      else if (thisMcp->Process() == "primary" 
          && ((std::abs(thisMcp->PdgCode()) == 13) 
            || (std::abs(thisMcp->PdgCode()) == 211)))
        trackCleanlinessPrimaryMuonOrPion->Fill(thisMcpCleanliness);
      else if (thisMcp->Process() == "primary" && std::abs(thisMcp->PdgCode()) == 11)
        trackCleanlinessPrimaryElectron->Fill(thisMcpCleanliness);
      else if (thisMcp->Process() == "conv")
        trackCleanlinessConv->Fill(thisMcpCleanliness);
      else if ((thisMcp->Process() == "neutronInelastic")
          || (thisMcp->Process() == "protonInelastic")
          || (thisMcp->Process() == "pi+Inelastic"))
        trackCleanlinessInelastic->Fill(thisMcpCleanliness);
      else if (thisMcp->Process() == "muIoni")
        trackCleanlinessMuIoni->Fill(thisMcpCleanliness);
      else trackCleanlinessOther->Fill(thisMcpCleanliness);

      // angular information

      thisMatchedMomentum = _rbutilInstance.getMomentumVector(thisMcp); 
      thisMatchedLength   = thisMcp->Trajectory().TotalLength();

      // all tracks
      if (thisMcp->Process() == "primary"){

        std::vector<double> zDir = {0,0,1};

        thisZmatchedAngleYZ = _rbutilInstance.getAngle(zDir, thisMatchedMomentum, _rbutilInstance, "yz");

        matchedLengthVersusAngle->Fill(thisZmatchedAngleYZ, thisMatchedLength);


      }

      // specifically protons
      if (thisMcp->PdgCode() == 2212 && thisMcp->Process() == "primary" && nimID != -1){


        thisNimMatchedAngle = _rbutilInstance.getAngle(nimMatchedMomentum, thisMatchedMomentum, _rbutilInstance, "no");
        thisNimMatchedAngles.push_back(thisNimMatchedAngle);

        thisNimMatchedAngleXZ = _rbutilInstance.getAngle(nimMatchedMomentum, thisMatchedMomentum, _rbutilInstance, "xz");
        thisNimMatchedAnglesXZ.push_back(thisNimMatchedAngleXZ);

        mupMatchedAngleVersusPMom->Fill(thisNimMatchedAngleXZ, thisMcp->P());
      }

    }
    it++;
  } // trackHandle

  it = 0;
  for (auto const& thisShower : (*showerHandle)){

    std::cout << "Found shower with ID " << thisShower.ID() << std::endl; 
    std::vector< art::Ptr<simb::MCParticle> > mcps = mcpsFromShowers.at(it);

    isRecoShowerTruthMatched.push_back(mcps.size());

    for (auto const& thisMcp : mcps){

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
    it++;
  } // showerHandle

  recoTree->Fill();
}

void recohelper::RecoBenchmarker::endJob()
{

  TFile& file = tfs->file();
  file.cd();

  // mupangleMomentumEfficiencies
  mupAngleMomentumEfficiency = 
    (TH2D*)mupMatchedAngleVersusPMom->Clone("mupAngleMomentumEfficiency");
  mupAngleMomentumEfficiency->Divide(mupMcpAngleVersusPMom);

  // pMomentumEfficiencies
  TH2D* h = (TH2D*)mupMatchedAngleVersusPMom->Clone("h");
  pMatchedProjectedMomentum = (TH1D*)h->ProjectionY("pMatchedProjectedMomentum");

  TH2D* h2 = (TH2D*)mupMcpAngleVersusPMom->Clone("h2");
  pMcpProjectedMomentum = (TH1D*)h2->ProjectionY("pMcpProjectedMomentum");

  pMomentumEfficiency = ((TH1D*)pMatchedProjectedMomentum->Clone("pMomentumEfficiency"));
  pMomentumEfficiency->Divide(pMcpProjectedMomentum);

  // pAngleEfficiencies
  pMatchedProjectedAngle = (TH1D*)h->ProjectionX("pMatchedProjectedAngle");
  pMcpProjectedAngle = (TH1D*)h2->ProjectionX("pMcpProjectedAngle");

  pAngleEfficiency = 
    ((TH1D*)pMatchedProjectedAngle->Clone("pAngleEfficiency"));
  pAngleEfficiency->Divide(pMcpProjectedAngle);

  h->Delete();
  h2->Delete();

  //trackAngleLengthEfficiencies
  lengthVersusAngleEfficiency 
    = (TH2D*)matchedLengthVersusAngle->Clone("lengthVersusAngleEfficiency");
  lengthVersusAngleEfficiency->Divide(mcpLengthVersusAngle);

  matchedProjectedLength = (TH1D*)matchedLengthVersusAngle->ProjectionY();
  mcpProjectedLength = (TH1D*)mcpLengthVersusAngle->ProjectionY();

  trackLengthEfficiency = 
    (TH1D*)matchedProjectedLength->Clone("trackLengthEfficiency");
  trackLengthEfficiency->Divide(mcpProjectedLength);

  matchedProjectedAngle = (TH1D*)matchedLengthVersusAngle->ProjectionX();
  mcpProjectedAngle = (TH1D*)mcpLengthVersusAngle->ProjectionX();

  trackAngleEfficiency =
    (TH1D*)matchedProjectedAngle->Clone("trackAngleEfficiency");
  trackAngleEfficiency->Divide(mcpProjectedAngle);

  mupAngleMomentumEfficiency->Write();
  pMatchedProjectedMomentum->Write();
  pMcpProjectedMomentum->Write();
  pMomentumEfficiency->Write();
  pMatchedProjectedAngle->Write();
  pMcpProjectedAngle->Write();
  pAngleEfficiency->Write();
  lengthVersusAngleEfficiency->Write();
  matchedProjectedLength->Write();
  mcpProjectedLength->Write();
  trackLengthEfficiency->Write();
  matchedProjectedAngle->Write();
  mcpProjectedAngle->Write();
  trackAngleEfficiency->Write();
}

DEFINE_ART_MODULE(recohelper::RecoBenchmarker)
