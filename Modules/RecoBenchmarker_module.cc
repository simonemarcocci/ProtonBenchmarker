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
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
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
    std::string fClusterLabel;
    std::string fHitLabel;

    bool isDebug = false;

    rbutil::recoBenchmarkerUtility _rbutilInstance;

    art::ServiceHandle< art::TFileService > tfs;
    TTree* recoTree;

    // auxiliary 

    int fEvent;
    int fRun;
    int fSubRun;

    int nimID;
    int recoNimID;

    // hit list stuff
    std::vector< std::pair<int, float> >allHitPositions;
    TH1D* hitMatchScore;

    std::vector<int> matchIDChecker;

    // trueVertexXZPosition
    std::vector<float> trueVertexXZPosition;

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

    // Nim == Neutrino Induced Muon
    std::vector<float> thisNimMcpAngles;
    float thisNimMcpAngle;
    std::vector<float> thisNimMcpAnglesXZ;
    float thisNimMcpAngleXZ;
    std::vector<float> thisNimMatchedMcpAngles;
    float thisNimMatchedMcpAngle;
    std::vector<float> thisNimMatchedMcpAnglesXZ;
    float thisNimMatchedMcpAngleXZ;
    float thisZmatchedAngleYZ;
    float thisZDirMcpAngleYZ;

    // efficiency plots and prerequisites
    TH1D* muMatchedMcpMomentum;
    TH1D* muMcpMomentum;
    TH1D* muMomentumEfficiency;

    TH2D* allMatchedMcpLengthAngle;
    TH2D* allMcpLengthAngleYZ;
    TH2D* allLengthAngleEfficiency;

    TH1D* allMatchedMcpProjectedLength;
    TH1D* allMcpProjectedLength;
    TH1D* allProjectedLengthEfficiency;

    TH1D* allMatchedMcpProjectedAngle;
    TH1D* allMcpProjectedAngle;
    TH1D* allProjectedAngleEfficiency;

    TH2D* mupMatchedMcpAnglePMom;
    TH2D* mupMcpAnglePMom;
    TH2D* mupAngleMomentumEfficiency;

    TH1D* pMatchedMcpProjectedMomentum;
    TH1D* pMcpProjectedMomentum;
    TH1D* pProjectedMomentumEfficiency;

    TH1D* pMatchedMcpProjectedAngle;
    TH1D* pMcpProjectedAngle;
    TH1D* pProjectedAngleEfficiency;

    // cleanliness plots
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
  fClusterLabel = p.get<std::string> ("ClusterLabel");
  fHitLabel = p.get<std::string> ("HitLabel");

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

  // mcp information
  recoTree->Branch("thisNimMcpAngles", &thisNimMcpAngles);
  recoTree->Branch("thisNimMcpAnglesXZ", &thisNimMcpAnglesXZ);
  recoTree->Branch("thisMcpLength", &thisMcpLength);

  // matched mcp information
  recoTree->Branch("thisNimMatchedMcpAngles", &thisNimMatchedMcpAngles);
  recoTree->Branch("thisNimMatchedMcpAnglesXZ", &thisNimMatchedMcpAnglesXZ);
  
  // reconstructed information
  recoTree->Branch("thisRecoLength", &thisRecoLength);

  hitMatchScore = tfs->make<TH1D>("hitMatchScore", ";hitMatchScore Score;", 2, 0, 2);

  muMatchedMcpMomentum = tfs->make<TH1D>("muMatchedMcpMomentum", ";#mu Momentum;", 20, 0, 3.5);
  muMcpMomentum = tfs->make<TH1D>("muMcpMomentum", ";#mu Momentum;", 20, 0, 3.5);

  allMatchedMcpLengthAngle = tfs->make<TH2D>("allMatchedMcpLengthAngle", ";#theta_{YZ} (degrees); Length (cm)", 20, 0, 180, 10, 0, 10);
  allMcpLengthAngleYZ = tfs->make<TH2D>("allMcpLengthAngleYZ", ";#theta_{YZ} (degrees); Length (cm)", 20, 0, 180, 10, 0, 10);

  mupMatchedMcpAnglePMom = tfs->make<TH2D>("mupMatchedMcpAnglePMom", ";#theta^{#mup}_{XZ} (degrees) ;P_{p} (Gev)", 20, 0, 180, 20, 0, 1);
  mupMcpAnglePMom = tfs->make<TH2D>("mupMcpAnglePMom", ";#theta^{#mup}_{XZ};P_{p}", 20, 0, 180, 20, 0, 1);

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
  thisNimMcpAngles.clear();
  thisNimMcpAnglesXZ.clear();
  thisMcpMomentum.clear();
  nimMcpMomentum.clear();
  thisMcpLength.clear();
  thisNimMatchedMcpAngles.clear();
  thisNimMatchedMcpAnglesXZ.clear();
  nimMatchedMomentum.clear();
  thisRecoLength.clear();
  trueVertexXZPosition.clear();
  allHitPositions.clear();

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

  art::ValidHandle< std::vector< recob::Cluster > > clusterHandle = 
    e.getValidHandle< std::vector< recob::Cluster > >(fClusterLabel);

  art::FindManyP<recob::Hit> hitsFromClusters(clusterHandle, e, fClusterLabel);

  art::ValidHandle< std::vector< recob::Hit > > hitHandle = 
    e.getValidHandle< std::vector< recob::Hit > >(fHitLabel);

  //---------------------------------
  // MCParticles
  //---------------------------------


  // first loop muons to find true neutrino induced muon (NIM)
  nimID = -1;
  for (size_t i = 0; i < mcList.size();i++){

    const art::Ptr< simb::MCParticle >& thisMcp = mcList.at(i);
    if (std::abs(thisMcp->PdgCode()) == 13 && thisMcp->Process() == "primary" && _rbutilInstance.isInTPC(thisMcp) == true) {
      nimID = thisMcp->TrackId();
      nimMcpMomentum = _rbutilInstance.getMomentumVector(thisMcp);
      muMcpMomentum->Fill(thisMcp->P());
      trueVertexXZPosition = {(float)thisMcp->Vx(), (float)thisMcp->Vz()};
    }
  }

  // then the other particles...
  for (size_t i = 0; i < mcList.size(); i++){

    const art::Ptr<simb::MCParticle>& thisMcp = mcList.at(i);
    
    // only interested in tracks for now...
    if ((thisMcp->Process() != "primary")
        || (std::abs(thisMcp->PdgCode()) == 2112)
        || (std::abs(thisMcp->PdgCode()) == 14)
        || (std::abs(thisMcp->PdgCode()) == 12)
        || (std::abs(thisMcp->PdgCode()) == 22)
        || (std::abs(thisMcp->PdgCode()) == 111)
        || (std::abs(thisMcp->PdgCode()) == 11) 
        || (std::abs(thisMcp->PdgCode()) == 2212 && thisMcp->P() < 0.2)
        || (_rbutilInstance.isInTPC(thisMcp) == false)) continue;

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

    thisMcpLength.push_back(thisMcp->Trajectory().TotalLength());
    thisMcpMomentum = _rbutilInstance.getMomentumVector(thisMcp);

    // all tracks
    if (thisMcp->Process() == "primary"){

      std::vector<double> zDir = {0,0,1};

      thisZDirMcpAngleYZ = _rbutilInstance.getAngle(zDir, thisMcpMomentum, _rbutilInstance, "yz");
      allMcpLengthAngleYZ->Fill(thisZDirMcpAngleYZ, thisMcp->Trajectory().TotalLength());

    }

    // specifically protons
    if ((thisMcp->PdgCode() == 2212) && (thisMcp->Process() == "primary") && (nimID !=-1)){

      thisNimMcpAngle = 
        _rbutilInstance.getAngle(thisMcpMomentum, nimMcpMomentum, _rbutilInstance, "no");
      thisNimMcpAngles.push_back(thisNimMcpAngle);

      thisNimMcpAngleXZ =
        _rbutilInstance.getAngle(thisMcpMomentum, nimMcpMomentum, _rbutilInstance, "xz");
      thisNimMcpAnglesXZ.push_back(thisNimMcpAngleXZ);

      mupMcpAnglePMom->Fill(thisNimMcpAngleXZ, thisMcp->P());

    }

  } // MCParticles


  //----------------------------
  // Tracks
  //----------------------------

  // loop tracks and do truth matching to find neutrino induced muon
  recoNimID = -1;
  bool isMatched = false; 
  for (auto const& thisTrack : (*trackHandle)) {

    std::vector< art::Ptr<simb::MCParticle> > mcps = mcpsFromTracks.at(thisTrack.ID());

    for (auto const& thisMcp : mcps){

      if (thisMcp->PdgCode() == 13 && thisMcp->Process() == "primary" && _rbutilInstance.isInTPC(thisMcp) == true && isMatched == false){
        recoNimID = thisMcp->TrackId();
        nimMatchedMomentum = _rbutilInstance.getMomentumVector(thisMcp);
        muMatchedMcpMomentum->Fill(thisMcp->P());
        isMatched = true;
      }

    }

  }

  int it = 0;
  for (auto const& thisTrack : (*trackHandle)) {

    thisRecoLength.push_back(thisTrack.Length());

    std::vector< art::Ptr<simb::MCParticle> > mcps = mcpsFromTracks.at(it);
    isRecoTrackTruthMatched.push_back(mcps.size());

    // check to make sure two reconstructed tracks aren't matching to the same mcp
    // i.e. should remove broken tracks
    bool isMatched = false;
    for (auto const& thisMcp : mcps){

      for (size_t i = 0; i < matchIDChecker.size(); i++){

        if (thisMcp->TrackId() == matchIDChecker.at(i))
          isMatched = true;

      }
      matchIDChecker.push_back(thisMcp->TrackId());

      if (isMatched == true) continue;

      auto thisMcpCleanliness = mcpsFromTracks.data(0).at(0)->cleanliness;
      auto thisMcpCompleteness = mcpsFromTracks.data(0).at(0)->completeness; // gives nonsense

      // filling cleanliness plots
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

      if ((thisMcp->Process() != "primary")
          || (std::abs(thisMcp->PdgCode()) == 2112)
          || (std::abs(thisMcp->PdgCode()) == 14)
          || (std::abs(thisMcp->PdgCode()) == 12)
          || (std::abs(thisMcp->PdgCode()) == 22)
          || (std::abs(thisMcp->PdgCode()) == 111)
          || (std::abs(thisMcp->PdgCode()) == 11) // remove anything which is shower like for now
          || (std::abs(thisMcp->PdgCode()) == 2212 && thisMcp->P() < 0.2)
          || (_rbutilInstance.isInTPC(thisMcp)) == false) continue;

      if (isDebug == true){
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


      // angular information

      thisMatchedMomentum = _rbutilInstance.getMomentumVector(thisMcp); 
      thisMatchedLength   = thisMcp->Trajectory().TotalLength();

      // all tracks
      if (thisMcp->Process() == "primary"){

        std::vector<double> zDir = {0,0,1};

        thisZmatchedAngleYZ = _rbutilInstance.getAngle(zDir, thisMatchedMomentum, _rbutilInstance, "yz");

        allMatchedMcpLengthAngle->Fill(thisZmatchedAngleYZ, thisMatchedLength);


      }

      // specifically protons
      if (thisMcp->PdgCode() == 2212 && thisMcp->Process() == "primary" && recoNimID != -1){

        thisNimMatchedMcpAngle = _rbutilInstance.getAngle(nimMatchedMomentum, thisMatchedMomentum, _rbutilInstance, "no");
        thisNimMatchedMcpAngles.push_back(thisNimMatchedMcpAngle);

        thisNimMatchedMcpAngleXZ = _rbutilInstance.getAngle(nimMatchedMomentum, thisMatchedMomentum, _rbutilInstance, "xz");
        thisNimMatchedMcpAnglesXZ.push_back(thisNimMatchedMcpAngleXZ);

        mupMatchedMcpAnglePMom->Fill(thisNimMatchedMcpAngleXZ, thisMcp->P());
      }


    }
    it++;
  } // trackHandle

  //------------------------------
  // Showers
  //------------------------------

  it = 0;
  for (auto const& thisShower : (*showerHandle)){

    // this is to stop LArSoft complaining...
    std::cout << "Found shower with ID " << thisShower.ID() << std::endl; 
    std::vector< art::Ptr<simb::MCParticle> > mcps = mcpsFromShowers.at(it);

    isRecoShowerTruthMatched.push_back(mcps.size());

    for (auto const& thisMcp : mcps){

      auto thisMcpCleanliness = mcpsFromShowers.data(0).at(0)->cleanliness;
      auto thisMcpCompleteness = mcpsFromShowers.data(0).at(0)->completeness;

      if (isDebug == true){
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

  //---------------------------
  // Hits
  //---------------------------
  
  for (auto const& thisHit : (*hitHandle)){

    if (trueVertexXZPosition.size() == 0) break;

    // put hit into hitlist
    std::pair<int, float> hitPair;
    hitPair.first = (int)thisHit.Channel();
    hitPair.second = (float)thisHit.PeakTime();

    std::vector<float> hitXZpos = _rbutilInstance.getHitXZPosition(thisHit, _rbutilInstance);
  
    bool isHitInRange = _rbutilInstance.isHitNearVertex(trueVertexXZPosition, hitXZpos);

    if (isHitInRange == true)
      allHitPositions.push_back(hitPair);

  }

  //---------------------------
  // Clusters
  //---------------------------
 
  it = 0;
  for (auto const& thisCluster : (*clusterHandle)){
    if (trueVertexXZPosition.size() == 0) break;

    std::cout << thisCluster.ID() << std::endl;

    // get associated hits
    std::vector< art::Ptr< recob::Hit > > hits = hitsFromClusters.at(it);

    for (auto const& thisHit : hits){

      int isMatched = 0;
      int hitChannel = thisHit->Channel();
      int hitPeakTime = thisHit->PeakTime();

      std::vector<float> hitXZpos = _rbutilInstance.getHitXZPosition(*thisHit, _rbutilInstance);

      bool isHitInRange = _rbutilInstance.isHitNearVertex(trueVertexXZPosition, hitXZpos);

      if (isHitInRange == false) continue;

      for (size_t i = 0; i < allHitPositions.size(); i++){

        if (hitChannel == allHitPositions.at(i).first && hitPeakTime == allHitPositions.at(i).second){

          isMatched = 1;
          hitMatchScore->Fill(isMatched);
        
        }

      }


    }

    hitMatchScore->SetBinContent(0, (int)allHitPositions.size() - hitMatchScore->GetBinContent(1));

    it++;
  }
  

  recoTree->Fill();
}

void recohelper::RecoBenchmarker::endJob()
{

  TFile& file = tfs->file();
  file.cd();

  //--------------------------
  // Produce eff. plots
  //--------------------------

  // mu efficiencies
  muMomentumEfficiency =
    (TH1D*)muMatchedMcpMomentum->Clone("muMomentumEfficiency");
  muMomentumEfficiency->Divide(muMcpMomentum);

  // mupangleMomentumEfficiencies
  mupAngleMomentumEfficiency = 
    (TH2D*)mupMatchedMcpAnglePMom->Clone("mupAngleMomentumEfficiency");
  mupAngleMomentumEfficiency->Divide(mupMcpAnglePMom);

  // pMomentumEfficiencies
  mupAngleMomentumEfficiency 
    =(TH2D*)mupMatchedMcpAnglePMom->Clone("mupAngleMomentumEfficiency");
  mupAngleMomentumEfficiency->Divide(mupMcpAnglePMom);

  pMatchedMcpProjectedAngle = (TH1D*)mupMatchedMcpAnglePMom->ProjectionY();
  pMcpProjectedAngle = (TH1D*)mupMcpAnglePMom->ProjectionY();

  pProjectedAngleEfficiency =
    (TH1D*)pMatchedMcpProjectedAngle->Clone("pProjectedAngleEfficiency");
  pProjectedAngleEfficiency->Divide(pMcpProjectedAngle);

  pMatchedMcpProjectedMomentum = (TH1D*)mupMatchedMcpAnglePMom->ProjectionX();
  pMcpProjectedMomentum = (TH1D*)mupMcpAnglePMom->ProjectionX();

  pProjectedMomentumEfficiency =
    (TH1D*)pMatchedMcpProjectedMomentum->Clone("pProjectedMomentumEfficiency");
  pProjectedMomentumEfficiency->Divide(pMcpProjectedMomentum);

  //trackAngleLengthEfficiencies
  allLengthAngleEfficiency 
    = (TH2D*)allMatchedMcpLengthAngle->Clone("allLengthAngleEfficiency");
  allLengthAngleEfficiency->Divide(allMcpLengthAngleYZ);

  allMatchedMcpProjectedLength = (TH1D*)allMatchedMcpLengthAngle->ProjectionY();
  allMcpProjectedLength = (TH1D*)allMcpLengthAngleYZ->ProjectionY();

  allProjectedLengthEfficiency = 
    (TH1D*)allMatchedMcpProjectedLength->Clone("allProjectedLengthEfficiency");
  allProjectedLengthEfficiency->Divide(allMcpProjectedLength);

  allMatchedMcpProjectedAngle = (TH1D*)allMatchedMcpLengthAngle->ProjectionX();
  allMcpProjectedAngle = (TH1D*)allMcpLengthAngleYZ->ProjectionX();

  allProjectedAngleEfficiency =
    (TH1D*)allMatchedMcpProjectedAngle->Clone("allProjectedAngleEfficiency");
  allProjectedAngleEfficiency->Divide(allMcpProjectedAngle);

  // aaaaand write...
  muMomentumEfficiency->Write();
  mupAngleMomentumEfficiency->Write();
  pMatchedMcpProjectedMomentum->Write();
  pMcpProjectedMomentum->Write();
  pProjectedMomentumEfficiency->Write();
  pMatchedMcpProjectedAngle->Write();
  pMcpProjectedAngle->Write();
  pProjectedAngleEfficiency->Write();
  allLengthAngleEfficiency->Write();
  allMatchedMcpProjectedLength->Write();
  allMcpProjectedLength->Write();
  allProjectedLengthEfficiency->Write();
  allMatchedMcpProjectedAngle->Write();
  allMcpProjectedAngle->Write();
  allProjectedAngleEfficiency->Write();

}

DEFINE_ART_MODULE(recohelper::RecoBenchmarker)
