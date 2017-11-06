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
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"

// local includes
#include "uboone/RecoBenchmarker/Algos/recoBenchmarkerUtility.h"

#define isDebug 0

#define isSingleParticle 0

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
    std::string fCalorimetryLabel;
    std::string fTrackTruthLabel;
    std::string fShowerLabel;
    std::string fShowerTruthLabel;
    std::string fClusterLabel;
    std::string fHitLabel;
    std::string fMCTruthLabel;
    std::string fG4TruthLabel;


    rbutil::recoBenchmarkerUtility _rbutilInstance;

    art::ServiceHandle< art::TFileService > tfs;
    TTree* recoTree;

    //a few parameters for histo plotting
    long low_edge = 350;
    long high_edge = 1500;
    long length_cut = 20; //8cm

    // auxiliary 

    int fEvent;
    int fRun;
    int fSubRun;
    int fccnc;
    int finteraction;

    int nimID;
    int recoNimID;

    //info on the MC particle
    std::map<Int_t,Int_t> fparticle_count;
    std::vector<float> flength;
    std::vector<float> fstartT;
    std::vector<float> fstart_x;
    std::vector<float> fstart_y;
    std::vector<float> fstart_z;
    std::vector<float> fend_x;
    std::vector<float> fend_y;
    std::vector<float> fend_z;
    std::vector<int> fn_steps;
    std::vector<int> fpdg;
    std::vector<int> fmother_id;
   std::vector<int > fg4_id; 
   std::vector<float> fp0;//initial momentum
   std::vector<float> fp0x;
   std::vector<float> fp0y;
   std::vector<float> fp0z;
   std::vector<float> fkinE;
   std::vector<float> fcostheta_muon;
   std::vector<bool> fis_leading;
   
   std::vector<double> fmuon_dqdx;
   std::vector<double> fmuon_residual;
   std::vector<double> fmuon_dedx;
   double fmuon_range;

   //info coming from the tracking algorithm - when there is mc truth
   std::vector<bool> fis_tracked;
   std::vector<bool> fmatch_multiplicity;
   //std::vector<bool> fis_mismatched; //it says if the MC truth assignment in different planes is different (possible hint for wrong or problematic backtracking)
   std::vector<float> fcostheta_muon_reco; //MCtruth info
   std::vector<int> fpdg_reco;
   std::vector<float> flength_reco;
   std::vector<float> freco_momentum_mcs; //(GeV) MCS
   std::vector<float> freco_momentum_mcs_llhd; //(GeV) MCS LLHD
   std::vector<float> freco_momentum_range; //MeV
   std::vector<float> fpurity;
   std::vector<float> fcompleteness;
   std::vector<double> fnhits;
   std::vector<float> freco_kinE;
   std::vector<float> freco_startx;
   std::vector<float> freco_starty;
   std::vector<float> freco_startz;
   std::vector<float> freco_endx;
   std::vector<float> freco_endy;
   std::vector<float> freco_endz;
   
   //info coming from the tracking algorithm - when there is NO mc truth
   std::vector<bool> ffake_is_tracked;
   std::vector<bool> ffake_is_mismatched;
   std::vector<float> ffake_costheta_muon_reco; //MCtruth info
   std::vector<int> ffake_pdg_reco;
   std::vector<float> ffake_length_reco;
   std::vector<float> ffake_reco_momentum_mcs; //(GeV) MCS
   std::vector<float> ffake_reco_momentum_mcs_llhd; //(GeV) MCS LLHD
   std::vector<float> ffake_reco_momentum_range; //MeV
   std::vector<float> ffake_purity;
   std::vector<float> ffake_completeness;
   std::vector<double> ffake_nhits;
   std::vector<float> ffake_kinE;


   void AllocateAnalysisHistograms();
   void FillAnalysisHistograms();
   void FillCumulativeHistograms();

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

    // Simone's plots
   TH1D* hmuon_pos_res;
   TH1D* hmuon_pos_res_goodprotons;
   TH1D* hmuon_pos_res_badprotons;
   TH1D* hmuon_pos_res_lowprotons;
   TH1D* hproton_pos_res;
   TH1D* hproton_pos_res_goodprotons;
   TH1D* hproton_pos_res_badprotons;
   TH1D* hproton_pos_res_lowprotons;
   TH1D* hmuon_proton_tracked;
   TH1D* hmuon_spectrum;
   TH1D* hmuon_spectrum_all;
   TH1D* hmuon_length;
   TH1D* hmuon_length_all;
   TH1D* hproton_kinE;
   TH1D* hproton_kinE_all;
   TH1D* hproton_p;
   TH1D* hproton_p_all;
   TH1D* hproton_l;
   TH1D* hproton_l_all;
   TH1D* h_pmu_end_not_tracked;
   TH1D* h_pmu_end_tracked;
   TH1D* h_theta_mu_tracked;
   TH1D* h_theta_mu_not_tracked;
   TH1D* h_theta_mu;
   TH2D* h_theta_mu_length;
   TH2D* h_theta_mu_length_all;
   TH2D* h_dqdx_merged;
   TH2D* h_dqdx_not_merged;
   TH2D*  h_dqdx_low_protons;
   TH1D* h_dqdx_1d_merged;
   TH1D* h_dqdx_1d_not_merged;
   TH2D* h_dqdx_tailtotot_length_merged;
   TH2D* h_dqdx_tailtotot_length_not_merged;
   TH1D* htail_to_tot_low_protons;
   TH1D* htail_to_tot_merged;
   TH1D* htail_to_tot_not_merged;
   TH2D* h_dqdx_merged_service;
   TH2D* h_dqdx_not_merged_service;
   TH2D* h_dqdx_low_protons_service;
    

    void AllocateRecoVectors();
    void clear_vectors();
};


recohelper::RecoBenchmarker::RecoBenchmarker(fhicl::ParameterSet const & p)
  :
    EDAnalyzer(p)  // ,
    // More initializers here.
{

  fTrackLabel = p.get<std::string> ("TrackLabel");
  fCalorimetryLabel = p.get<std::string> ("CalorimetryLabel");
  fTrackTruthLabel = p.get<std::string> ("TrackTruthLabel");
  fShowerLabel = p.get<std::string> ("ShowerLabel");
  fMCTruthLabel = p.get<std::string> ("MCTruthLabel");
  fG4TruthLabel = p.get<std::string> ("G4TruthLabel");
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

  //mctruth info
  recoTree->Branch("ccnc", &fccnc);
  recoTree->Branch("interaction", &finteraction);

    //info on the MC particle
    recoTree->Branch("length",&flength);
    recoTree->Branch("startT",&fstartT);
    recoTree->Branch("start_x",&fstart_x);
    recoTree->Branch("start_y",&fstart_y);
    recoTree->Branch("start_z",&fstart_z);
    recoTree->Branch("end_x",&fend_x);
    recoTree->Branch("end_y",&fend_y);
    recoTree->Branch("end_z",&fend_z);
    recoTree->Branch("steps",&fn_steps);
    recoTree->Branch("pdg",&fpdg);
    recoTree->Branch("mother_id",&fmother_id);
   recoTree->Branch("g4_id",&fg4_id); //first index is t
   recoTree->Branch("p0",&fp0);//initial momentum
   recoTree->Branch("p0x",&fp0x);
   recoTree->Branch("p0y",&fp0y);
   recoTree->Branch("p0z",&fp0z);
   recoTree->Branch("kinE",&fkinE);
   recoTree->Branch("costheta_muon",&fcostheta_muon);
   recoTree->Branch("is_leading",&fis_leading);
   recoTree->Branch("particle_count",&fparticle_count);
   recoTree->Branch("muon_dedx",&fmuon_dedx);
   recoTree->Branch("muon_dqdx",&fmuon_dqdx);
   recoTree->Branch("muon_residual",&fmuon_residual);
   recoTree->Branch("muon_range",&fmuon_range);


  recoTree->Branch("isRecoTrackTruthMatched", &isRecoTrackTruthMatched);
  recoTree->Branch("isRecoShowerTruthMatched", &isRecoShowerTruthMatched);

  // mcp information
  recoTree->Branch("thisNimMcpAngles", &thisNimMcpAngles);
  recoTree->Branch("thisNimMcpAnglesXZ", &thisNimMcpAnglesXZ);
  recoTree->Branch("thisMcpLength", &thisMcpLength);

  // matched mcp information
  recoTree->Branch("thisNimMatchedMcpAngles", &thisNimMatchedMcpAngles);
  recoTree->Branch("thisNimMatchedMcpAnglesXZ", &thisNimMatchedMcpAnglesXZ);
   
  //info coming from the tracking algorithm - when there is mc truth
  recoTree->Branch("is_tracked",&fis_tracked);
  recoTree->Branch("match_multiplicity",&fmatch_multiplicity);
  recoTree->Branch("length_reco",&flength_reco);
  recoTree->Branch("reco_momentum_mcs",&freco_momentum_mcs); //(GeV) MCS
  recoTree->Branch("reco_momentum_mcs_llhd",&freco_momentum_mcs_llhd); //(GeV) MCS LLHD
  recoTree->Branch("reco_momentum_range",&freco_momentum_range); //MeV
  recoTree->Branch("purity",&fpurity);
  recoTree->Branch("completeness",&fcompleteness);
  recoTree->Branch("nhits",&fnhits);
  recoTree->Branch("reco_kinE",&freco_kinE);
  recoTree->Branch("reco_startx",&freco_startx);
  recoTree->Branch("reco_starty",&freco_starty);
  recoTree->Branch("reco_startz",&freco_startz);
  recoTree->Branch("reco_endx",&freco_endx);
  recoTree->Branch("reco_endy",&freco_endy);
  recoTree->Branch("reco_endz",&freco_endz);

  AllocateAnalysisHistograms();
   
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
 
  //auto const* theDetector = lar::providerFrom<detinfo::DetectorPropertiesService>();
  auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();  
  bool space_charge = true;

  fEvent  = e.id().event();
  fRun    = e.run();
  fSubRun = e.subRun();


  // get handles to objects of interest

  art::ValidHandle< std::vector<recob::Track> > trackHandle = 
    e.getValidHandle< std::vector<recob::Track> > (fTrackLabel);
  std::vector< art::Ptr<recob::Track> > trackList;
  //if (e.getByLabel( fTrackLabel, trackHandle))
    art::fill_ptr_vector(trackList, trackHandle); 

  art::FindManyP<simb::MCParticle, anab::BackTrackerMatchingData> mcpsFromTracks(trackHandle, e, fTrackTruthLabel);
  
  art::FindManyP<anab::Calorimetry> caloFromTracks(trackHandle, e, fCalorimetryLabel);

  art::ValidHandle< std::vector<recob::Shower> > showerHandle =
    e.getValidHandle< std::vector<recob::Shower> >(fShowerLabel);

  art::FindManyP<simb::MCParticle, anab::BackTrackerMatchingData> mcpsFromShowers(showerHandle, e, fShowerTruthLabel);

  art::Handle< std::vector<simb::MCParticle> > mcParticleHandle; 
  std::vector< art::Ptr<simb::MCParticle> > mcList;
  if (e.getByLabel( fG4TruthLabel, mcParticleHandle))
    art::fill_ptr_vector(mcList, mcParticleHandle); 

  art::Handle< std::vector<simb::MCTruth> > mcTruthHandle; 
  std::vector< art::Ptr<simb::MCTruth> > mcTruth;
  if (e.getByLabel( fMCTruthLabel, mcTruthHandle))
    art::fill_ptr_vector(mcTruth, mcTruthHandle); 
  
  art::ValidHandle< std::vector< recob::Cluster > > clusterHandle = 
    e.getValidHandle< std::vector< recob::Cluster > >(fClusterLabel);

  art::FindManyP<recob::Hit> hitsFromClusters(clusterHandle, e, fClusterLabel);

  //art::ValidHandle< std::vector< recob::Hit > > hitHandle = 
  //  e.getValidHandle< std::vector< recob::Hit > >(fHitLabel);

  //---------------------------------
  // MCParticles
  //---------------------------------

  //get MCparticles associated to MCtruth
  art::FindManyP<simb::MCParticle> MCpFromMCtruth( mcTruth , e, fG4TruthLabel );

  //loop on all MC truth frames (mostly 1 per event)
  for ( unsigned n_truth = 0; n_truth < mcTruth.size(); n_truth++ ) {
  
  clear_vectors();

#if isSingleParticle == 0 
  //check that the neutrino is in the active volume
  const simb::MCNeutrino thisNeutrino = mcTruth[n_truth]->GetNeutrino();
  const simb::MCParticle thisLepton = thisNeutrino.Lepton();
  if ( !_rbutilInstance.isInTPC( thisLepton ) ) continue;
  
  //store the interaction info
  fccnc = thisNeutrino.CCNC();
  finteraction = thisNeutrino.InteractionType();

#endif
  // first loop muons to find true neutrino induced muon (NIM)
  nimID = -1;
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
    if (std::abs(thisMcp->PdgCode()) == 13 && thisMcp->StatusCode()==1 && thisMcp->Mother()==0 && thisMcp->Process() == "primary" && _rbutilInstance.isInTPC(thisMcp) == true) {
      nimID = thisMcp->TrackId();
      nimMcpMomentum = _rbutilInstance.getMomentumVector(thisMcp);
      muMcpMomentum->Fill(thisMcp->P());
      trueVertexXZPosition = {(float)thisMcp->Vx(), (float)thisMcp->Vz()};
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

  if (muon_p==0) mf::LogError(__FUNCTION__) << "Error! There must be at least 1 muon in this analysis! " << std::endl;
 
  //now other MC particles
  for (unsigned i = 0; i < MCpFromMCtruth.at(n_truth).size(); i++){
    const art::Ptr< simb::MCParticle >& thisMcp = MCpFromMCtruth.at(n_truth).at(i);
    
    //in case of NC, the neutrino is kept and thus we should drop it
    //if it's an electron, check if it's a Michel. If yes, keep it.
    if ( !(std::abs(thisMcp->PdgCode()) == 11 && thisMcp->Mother() == muon_id && abs(thisMcp->Position().X()-muon_endx) < DBL_EPSILON &&
		    abs(thisMcp->Position().Y()-muon_endy) < DBL_EPSILON && abs(thisMcp->Position().Z()-muon_endz) < DBL_EPSILON ) )
    if ( std::abs(thisMcp->PdgCode()) == 14 || thisMcp->Process()!="primary" || thisMcp->StatusCode()!=1 || thisMcp->Mother()>0 ) continue; //we only want primaries

    if ( (std::abs(thisMcp->PdgCode()) == 11 && thisMcp->Mother() == muon_id && abs(thisMcp->Position().X()-muon_endx) < DBL_EPSILON &&
		    abs(thisMcp->Position().Y()-muon_endy) < DBL_EPSILON && abs(thisMcp->Position().Z()-muon_endz) < DBL_EPSILON ) )
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

	
	//count particle types
	if (fparticle_count.find( thisMcp->PdgCode() ) == fparticle_count.end()) //not found
	     fparticle_count[ thisMcp->PdgCode() ] = 1;
	else
	     fparticle_count[thisMcp->PdgCode()] = fparticle_count[thisMcp->PdgCode()] +1;
	
	flength.push_back( thisMcp->Trajectory().TotalLength() );
	fn_steps.push_back ( thisMcp->Trajectory().size() );
	fstartT.push_back( thisMcp->T() );

	//if space charge is ON, we should correct the MC truth position to perform true-reco comparisons (Giuseppe)
	/*auto scecorr = SCE->GetPosOffsets ( geo::Point_t( thisMcp->Position().X(), thisMcp->Position().Y(), thisMcp->Position().Z()) );
	double g4Ticks = detClocks->TPCG4Time2Tick( thisMcp->T() )+theDetector->GetXTicksOffset(0,0,0)-theDetector->TriggerOffset();
	double xOffset = theDetector->ConvertTicksToX(g4Ticks, 0, 0, 0)-scecorr.X();
	double yOffset = scecorr.Y();
	double zOffset = scecorr.Z();*/
	//anatree recipe:
	double xOffset = 0.7 - SCE->GetPosOffsets( thisMcp->Position().X(),  thisMcp->Position().Y(), thisMcp->Position().Z())[0];
	double yOffset = SCE->GetPosOffsets( thisMcp->Position().X(),  thisMcp->Position().Y(), thisMcp->Position().Z())[1];
	double zOffset = SCE->GetPosOffsets( thisMcp->Position().X(),  thisMcp->Position().Y(), thisMcp->Position().Z())[2];
	if (!space_charge) {
		xOffset = 0;
		yOffset = 0;
		zOffset = 0;
	}
	fstart_x.push_back ( thisMcp->Position().X() + xOffset );
	fstart_y.push_back ( thisMcp->Position().Y() + yOffset );
	fstart_z.push_back ( thisMcp->Position().Z() + zOffset );
	
	if (space_charge) {
	xOffset = 0.7 - SCE->GetPosOffsets( thisMcp->EndPosition().X(),  thisMcp->EndPosition().Y(), thisMcp->EndPosition().Z())[0];
	yOffset = SCE->GetPosOffsets( thisMcp->EndPosition().X(),  thisMcp->EndPosition().Y(), thisMcp->EndPosition().Z())[1];
	zOffset = SCE->GetPosOffsets( thisMcp->EndPosition().X(),  thisMcp->EndPosition().Y(), thisMcp->EndPosition().Z())[2];
	}
	fend_x.push_back ( thisMcp->EndPosition().X() + xOffset );
	fend_y.push_back ( thisMcp->EndPosition().Y() + yOffset );
	fend_z.push_back ( thisMcp->EndPosition().Z() + zOffset );
	fmother_id.push_back ( thisMcp->Mother() );
	fpdg.push_back( thisMcp->PdgCode() );	
	fg4_id.push_back(thisMcp->TrackId());
	fp0.push_back( thisMcp->Momentum().Rho());
	fp0x.push_back( thisMcp->Momentum().X() );
	fp0y.push_back( thisMcp->Momentum().Y() );
	fp0z.push_back( thisMcp->Momentum().Z() );
	
	//loop on other tracked_particles and decide if this is the leading based on initial kinetic energy. This must be done before filling the kinE vector for the current particle!
	bool is_leading = true;
	float current_kinE = thisMcp->E() - thisMcp->Mass() ;
	for (unsigned jj=0; jj< fkinE.size(); jj++) { //loop on previous particles
	if ( fpdg[jj] == thisMcp->PdgCode() ) {
	    if ( fkinE[jj] > current_kinE )
	       is_leading = false;
		}
	}
	fis_leading.push_back( is_leading );

	//now write current kin E
	fkinE.push_back( thisMcp->E() - thisMcp->Mass() );

	if ( thisMcp->Momentum().Rho() != 0) {
	     fcostheta_muon.push_back( muon_px*thisMcp->Momentum().X()/thisMcp->Momentum().Rho()/muon_p
			     	+ muon_py*thisMcp->Momentum().Y()/thisMcp->Momentum().Rho()/muon_p
				+ muon_pz*thisMcp->Momentum().Z()/thisMcp->Momentum().Rho()/muon_p );
	} else {
	     fcostheta_muon.push_back(-2);
	}
	
	//add dummy entries for the reco variables. They will be possibly updated later, if a matching reco track is found
	AllocateRecoVectors();
      
  } // MCParticles


  //----------------------------
  // Tracks
  //----------------------------

  // loop tracks and do truth matching 
  for (auto const& thisTrack : trackList) {

    std::vector< art::Ptr<simb::MCParticle> > mcps = mcpsFromTracks.at(thisTrack->ID());
    std::vector< art::Ptr<anab::Calorimetry> > calos = caloFromTracks.at(thisTrack->ID());

    if (mcps.size() >1 ) mf::LogWarning(__FUNCTION__) << "Warning !!! More than 1 MCparticle associated to the same track!" << std::endl;
    if (calos.size() != 3 ) mf::LogWarning(__FUNCTION__) << "Warning !!! More than 1 Calorimetry info associated to the same track!" << std::endl;

    for (auto const& thisMcp : mcps){
    
    //this if is necessary, since sometimes a particle is matched to a secondary (electron, etc) -> to be checked
    if ( !(std::abs(thisMcp->PdgCode()) == 11 && thisMcp->Mother() == muon_id && abs(thisMcp->Position().X()-muon_endx) < DBL_EPSILON &&
		    abs(thisMcp->Position().Y()-muon_endy) < DBL_EPSILON && abs(thisMcp->Position().Z()-muon_endz) < DBL_EPSILON ) )
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
    auto it_found = std::find( fg4_id.begin(), fg4_id.end(), thisMcp->TrackId() ) ;
    if (it_found==fg4_id.end()) mf::LogError(__FUNCTION__) << "ERROR!!! Matched particle not found!" << std::endl;
    size_t pos = it_found - fg4_id.begin();

    //save information on reco track
    if (fis_tracked[pos]) {
	    mf::LogDebug() << "Probably broken track!" << std::endl;
	    fmatch_multiplicity[pos] = fmatch_multiplicity[pos] + 1;
    	    //decide if keeping the old info or the new one - based on the minimum track length difference
	    if ( abs(thisTrack->Length() - flength[pos]) > abs(flength_reco[pos] - flength[pos]) )
		    continue;
    } else {
        fis_tracked[pos] = true;
	fmatch_multiplicity[pos] = fmatch_multiplicity[pos] + 1;
    }
	flength_reco[pos] = thisTrack->Length();
	trkf::TrackMomentumCalculator trkm; //track momentum calculator
	trkm.SetMinLength(0); //minimum track length for momentum calculation
	freco_momentum_mcs[pos] = trkm.GetMomentumMultiScatterChi2( thisTrack );
	freco_momentum_mcs_llhd[pos] = trkm.GetMomentumMultiScatterLLHD( thisTrack ) ;
	freco_momentum_range[pos] = trkm.GetTrackMomentum( thisTrack->Length(), fpdg[pos] ); //use info on pdg
	freco_startx[pos] = thisTrack->Vertex().X(); 
	freco_starty[pos] = thisTrack->Vertex().Y(); 
	freco_startz[pos] = thisTrack->Vertex().Z(); 
	freco_endx[pos] = thisTrack->End().X(); 
	freco_endy[pos] = thisTrack->End().Y(); 
	freco_endz[pos] = thisTrack->End().Z(); 
 	
	if (thisMcp->PdgCode() == 13 ) { //for the muon, only when there is the info
		if (fmuon_dqdx.size()!=0) mf::LogError(__FUNCTION__) << "Calorimetry should be filled only once!!!!" << std::endl;
		for (size_t position=0; position<calos.at(2)->dQdx().size(); position++) {
		fmuon_dqdx.push_back(calos.at(2)->dQdx()[position]); //look only at collection
		fmuon_dedx.push_back(calos.at(2)->dEdx()[position]); //look only at collection
		fmuon_residual.push_back(calos.at(2)->ResidualRange()[position]); //look only at collection
		}
	fmuon_range = calos.at(2)->Range();
	}

	auto thisMcpCleanliness = mcpsFromTracks.data(thisTrack->ID()).at(0)->cleanliness;
	auto thisMcpCompleteness = mcpsFromTracks.data(thisTrack->ID()).at(0)->completeness;
	fpurity[pos] = thisMcpCleanliness;
	fcompleteness[pos] = thisMcpCompleteness;
    }//MCParticle
  }//Tracks

  /*
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
      //auto thisMcpCompleteness = mcpsFromTracks.data(0).at(0)->completeness; // gives nonsense

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

#if isDebug == 1
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
#endif

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
      //auto thisMcpCompleteness = mcpsFromShowers.data(0).at(0)->completeness;
#if isDebug == 1
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
#endif 

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
  */

  recoTree->Fill();

  FillAnalysisHistograms();

  } //MCtruth
  

}

void recohelper::RecoBenchmarker::AllocateAnalysisHistograms() {
   
   hmuon_pos_res = tfs->make<TH1D>("muon_pos_res","All reco muons;Reco - True Muon start position (cm);",1000,0,50); //displacement between real muon start position and reco
   hmuon_pos_res_goodprotons = tfs->make<TH1D>("muon_pos_res_goodprotons","All reco muons with good proton reco (>20MeV);Reco - True Muon start position (cm);",1000,0,50); //displacement between real muon start position and reco for events with all reco protons
   hmuon_pos_res_badprotons = tfs->make<TH1D>("muon_pos_res_badprotons","All reco muons with at least a bad proton reco (>20MeV);Reco - True Muon start position (cm);",1000,0,50); //displacement between real muon start position and reco for events with some non-reco protons
   hmuon_pos_res_lowprotons = tfs->make<TH1D>("muon_pos_res_lowprotons","All reco muons when there are protons <20MeV;Reco - True Muon start position (cm);",1000,0,50); //displacement between real muon start position and reco for events with some non-reco low energy protons
   hproton_pos_res = tfs->make<TH1D>("proton_pos_res","All reco protons;Reco - True Proton start position (cm);",1000,0,50); //displacement between real proton start position and reco
   hproton_pos_res_goodprotons = tfs->make<TH1D>("proton_pos_res_goodprotons","All reco protons with E>20MeV;Reco - True Proton start position (cm);",1000,0,50); //displacement between real proton start position and reco
   hproton_pos_res_badprotons = tfs->make<TH1D>("proton_pos_res_badprotons","All reco protons with at least a bad proton reco (>20MeV);Reco - True Proton start position (cm);",1000,0,50); //displacement between real proton start position and reco
   hproton_pos_res_lowprotons = tfs->make<TH1D>("proton_pos_res_lowprotons","All reco protons when there are protons <20MeV;Reco - True Proton start position (cm);",1000,0,50); //displacement between real proton start position and reco
   hmuon_proton_tracked = tfs->make<TH1D>("muon_proton_tracked","Reco Proton Start - Reco Muon Start;Proton-Muon (cm);",1000,0,50); //displacement between reco muon start position and proton reco start position
   hmuon_spectrum = tfs->make<TH1D>("muon_spectrum","Muon kinetic energy; Kinetic Energy (GeV);",1000,0,5); //reco muons
   hmuon_spectrum_all = tfs->make<TH1D>("muon_spectrum_all","Muon kinetic energy; Kinetic Energy (GeV);",1000,0,5); //all of them, not just reco
   hmuon_length = tfs->make<TH1D>("muon_length","Muon length; True Length (cm);",1000,0,1000); //reco muons
   hmuon_length_all = tfs->make<TH1D>("muon_length_all","Muon length; True Length (cm);",1000,0,1000); //all of them, not just reco
   hproton_kinE = tfs->make<TH1D>("proton_kinE","Proton reco efficiency; Kinetic Energy (GeV)",1000,0,2); //reco efficiency protons vs kin E
   hproton_kinE_all = tfs->make<TH1D>("proton_kinE_all","Proton reco efficiency; Kinetic Energy (GeV)",1000,0,2); //reco efficiency protons vs kin E
   hproton_p = tfs->make<TH1D>("proton_p","Proton reco efficiency; Momentum (MeV/c)",1000,0,10); //reco efficiency protons vs p
   hproton_p_all = tfs->make<TH1D>("proton_p_all","Proton reco efficiency; Momentum (MeV/c)",1000,0,10); //reco efficiency protons vs p
   hproton_l = tfs->make<TH1D>("proton_l","Proton reco efficiency; True length (cm)",1000,0,200); //reco efficiency protons vs p
   hproton_l_all = tfs->make<TH1D>("proton_l_all","Proton reco efficiency; True length (cm)",1000,0,200); //reco efficiency protons vs p
   h_pmu_end_not_tracked = tfs->make<TH1D>("pmu_end_not_tracked","Not Tracked protons;Distance (cm);",1000,0,100); //lateral distance between proton end and muon
   h_pmu_end_tracked = tfs->make<TH1D>("pmu_end_tracked","Tracked protons;Distance (cm);",1000,0,100); //lateral distance between proton end and muon
   h_theta_mu_tracked = tfs->make<TH1D>("theta_mu_tracked","Cos #theta between muon and tracked protons;cos #theta",1000,-1,1); //costheta between muon and tracked protons
   h_theta_mu_not_tracked = tfs->make<TH1D>("theta_mu_not_tracked","Cos #theta between muon and non tracked protons;cos #theta",1000,-1,1);
   h_theta_mu = tfs->make<TH1D>("theta_mu","Cos #theta between muon and protons;cos #theta",1000,-1,1); //true costheta between muon and protons
   h_theta_mu_length = tfs->make<TH2D>("theta_mu_length","Tracking efficiency vs (length, cos #theta_{p#mu});Cos #theta; l (cm)",1000,-1,1,1000,0,100);
   h_theta_mu_length_all = tfs->make<TH2D>("theta_mu_length_all","Tracking efficiency vs (length, cos #theta_{p#mu});Cos #theta; l (cm)",1000,-1,1,1000,0,100);
   h_dqdx_merged = tfs->make<TH2D>("dqdx_merged","dq/dx for events with at least a merged proton; Distance from vertex (cm); dq/dx (ADC)",2500,0,250,1500,0,1500);
   h_dqdx_not_merged = tfs->make<TH2D>("dqdx_not_merged","dq/dx for events with no merged proton; Distance from vertex (cm); dq/dx (ADC)",2500,0,250,1500,0,1500);
   h_dqdx_low_protons = tfs->make<TH2D>("dqdx_low_protons","dq/dx for events with low E proton; Distance from vertex (cm); dq/dx (ADC)",2500,0,250,1500,0,1500);
   h_dqdx_1d_merged = tfs->make<TH1D>("dqdx_1d_merged","dq/dx integrated in (0,8cm) when protons are merged; dq/dx (ADC);",1500,0,1500);
   h_dqdx_1d_not_merged = tfs->make<TH1D>("dqdx_1d_not_merged","dq/dx integrated in (0,8cm) when protons are not merged; dq/dx (ADC);",1500,0,1500);
   h_dqdx_tailtotot_length_merged = tfs->make<TH2D>("dqdx_tailtotot_length_merged","Tail to tot vs length for merged tracks; Integration Length (mm); Tail to tot",1000,0,1000,1000,0,1);
   h_dqdx_tailtotot_length_not_merged = tfs->make<TH2D>("dqdx_tailtotot_length_not_merged","Tail to tot vs length for merged tracks; Integration Length (mm); Tail to tot",1000,0,1000,1000,0,1);
   htail_to_tot_low_protons = tfs->make<TH1D>("tailtotot_low_protons","Tail to tot w/ low energy protons; Tail to tot;",1000,0,1); // tail to tot merged
   htail_to_tot_merged = tfs->make<TH1D>("tailtotot_merged","Tail to tot w/ merged protons; Tail to tot;",100,0,100); // tail to tot merged
   htail_to_tot_not_merged = tfs->make<TH1D>("tailtotot_not_merged","Tail to tot w/ good protons; Tail to tot;",100,0,100); // tail to tot not merged
   h_dqdx_merged_service = tfs->make<TH2D>("dqdx_merged_service","dqdx_merged_service",2500,0,250,1500,0,1500);
   h_dqdx_not_merged_service = tfs->make<TH2D>("dqdx_not_merged_service","dqdx_not_merged_service",2500,0,250,1500,0,1500);
   h_dqdx_low_protons_service = tfs->make<TH2D>("dqdx_low_protons_service","dqdx_low_protons_service",2500,0,250,1500,0,1500);
}

void recohelper::RecoBenchmarker::FillAnalysisHistograms() {
  
    if (fccnc!=0 ) return; //comment if you're running over a std sample with NC interactions and you're interested in them
    unsigned muon_pos = -1;
    //check that the muon is reco
    bool reco_muon = false;
    bool is_pion=false;
    bool is_lowmomentum_p = false;
    for (unsigned i = 0; i < fpdg.size(); i++) {
	    if ( fpdg[i] == 13 ) { //ismuon
	    hmuon_length_all->Fill( flength[i] );
	    hmuon_spectrum_all->Fill( fkinE[i] );
	    if ( fis_tracked[i] &&  fmuon_dqdx.size() != 0) { //want that the muon is well matched and w/ dqdx info
		muon_pos = i; 
	    	hmuon_length->Fill( flength[i] );
	    	hmuon_spectrum->Fill( fkinE[i] );
		reco_muon = true;
		hmuon_pos_res->Fill( sqrt( pow( freco_startx[i]- fstart_x[i] ,2) + pow( freco_starty[i]- fstart_y[i] ,2) + pow( freco_startz[i]- fstart_z[i] ,2) ));
		if ( sqrt( pow( freco_startx[i]- fstart_x[i] ,2) + pow( freco_starty[i]- fstart_y[i] ,2) + pow( freco_startz[i]- fstart_z[i] ,2) ) > 50 ) {
			reco_muon = false; //skip these events, might have direction flipped
	    	}
    		}
	    } //record info on the muon
	    if ( fpdg[i] == 211 || fpdg[i] == -211 || fpdg[i] == 111 ) is_pion=true; //record if there are any pions
	    if (  fpdg[i] == 2212 && fp0[i] <= 0.2 ) is_lowmomentum_p = true;
    }
    
    long count_not_tracked = 0;
    long count_tracked = 0;
    for (unsigned j=0; j<fpdg.size(); j++) {
	    if ( !reco_muon ) break; //select events with a reco muon
	    if ( fpdg[j]!=2212 ) continue; //watch only protons
	    
	    hproton_p_all->Fill( fp0[j] );
	    hproton_l_all->Fill( flength[j] );
	    hproton_kinE_all->Fill( fkinE[j] );
	    h_theta_mu_length_all->Fill( fcostheta_muon[j], flength[j]);
	    if ( fis_tracked[j] ) { //efficiency plots for all protons (but don't divide yet)
	    	hproton_p->Fill( fp0[j] );
	    	hproton_l->Fill( flength[j] );
	    	hproton_kinE->Fill( fkinE[j] );
	        h_theta_mu_length->Fill( fcostheta_muon[j], flength[j]);
	    }


    //now study all the protons with momentum > 200MeV
    if (fp0[j] > 0.2) {
  	 h_theta_mu->Fill( fcostheta_muon[j] );	
	 if ( !fis_tracked[j] ) {
		count_not_tracked++;
		h_pmu_end_not_tracked->Fill( flength[j] * sqrt( 1 - pow( fcostheta_muon[j] ,2 ) ) ) ;
		h_theta_mu_not_tracked->Fill ( fcostheta_muon[j] ) ;
		}
    	else if ( fis_tracked[j] ) {
		count_tracked++;
		h_pmu_end_tracked->Fill( flength[j] * sqrt( 1 - pow( fcostheta_muon[j] ,2 ) ) ) ;
		h_theta_mu_tracked->Fill ( fcostheta_muon[j] ) ;
		hmuon_proton_tracked->Fill ( sqrt( pow( freco_startx[j]- freco_startx[muon_pos] ,2) + pow( freco_starty[j]- freco_starty[muon_pos] ,2) + pow( freco_startz[j]- freco_startz[muon_pos] ,2) ) );
		hproton_pos_res->Fill ( sqrt( pow( freco_startx[j]- fstart_x[j] ,2) + pow( freco_starty[j]- fstart_y[j] ,2) + pow( freco_startz[j]- fstart_z[j] ,2) ) );
		}
   	}
}

	if ( fmuon_residual.size() != fmuon_dqdx.size()) cout << "ERROR on calorimetry vector sizes!!!" << endl;
	if ( count_not_tracked == 0 && count_tracked > 0 ) { //all protons are tracked 
		for (unsigned jj=1; jj<fmuon_dqdx.size()-1; jj++) {
				if ( fmuon_range - fmuon_residual[jj] < 8 ) { //look at 8cm only
				h_dqdx_1d_not_merged->Fill(fmuon_dqdx[jj]);
				}

				h_dqdx_not_merged->Fill( fmuon_range - fmuon_residual[jj], fmuon_dqdx[jj] );
				h_dqdx_not_merged_service->Fill( fmuon_range -  fmuon_residual[jj], fmuon_dqdx[jj] );
		}
		TH1D* h1 = NULL;
		h1 = h_dqdx_not_merged_service->ProjectionY("h",1, length_cut);
		if (h1) htail_to_tot_not_merged->Fill( h1->Integral(low_edge,high_edge) ); 
		h_dqdx_not_merged_service->Reset();

		hmuon_pos_res_goodprotons->Fill( sqrt( pow( freco_startx[muon_pos]- fstart_x[muon_pos] ,2) + pow( freco_starty[muon_pos]- fstart_y[muon_pos] ,2) + pow( freco_startz[muon_pos]- fstart_z[muon_pos] ,2) ) );
    		for (unsigned j=0; j < fpdg.size(); j++) {
	    	if ( !reco_muon ) break; //select events with a reco muon
	    	if ( fpdg[j]!=2212 && fp0[j] <= 0.2 && !fis_tracked[j]) continue; //watch only reco protons with p>0.2Gev/c
   		hproton_pos_res_goodprotons->Fill(  sqrt( pow( freco_startx[j]- fstart_x[j] ,2) + pow( freco_starty[j]- fstart_y[j] ,2) + pow( freco_startz[j]- fstart_z[j] ,2) ) ); 
   		if (is_lowmomentum_p) hproton_pos_res_lowprotons->Fill(  sqrt( pow( freco_startx[j]- fstart_x[j] ,2) + pow( freco_starty[j]- fstart_y[j] ,2) + pow( freco_startz[j]- fstart_z[j] ,2) ) ); 
		}
	} 

	if ( count_not_tracked > 0 ) { //at least 1 proton is not tracked
		for (unsigned jj=1; jj<fmuon_dqdx.size()-1; jj++) {
				if ( fmuon_range - fmuon_residual[jj] < 8 ) //look at 8cm only
				h_dqdx_1d_merged->Fill(fmuon_dqdx[jj]);

				h_dqdx_merged->Fill(  fmuon_range - fmuon_residual[jj], fmuon_dqdx[jj] );
				h_dqdx_merged_service->Fill( fmuon_range -  fmuon_residual[jj], fmuon_dqdx[jj] );
		}
		TH1D* h1 = NULL;
		h1 = h_dqdx_merged_service->ProjectionY("h",1,length_cut);
		if (h1)
		htail_to_tot_merged->Fill( h1->Integral(low_edge,high_edge) );
		h_dqdx_merged_service->Reset();
		    
		hmuon_pos_res_badprotons->Fill( sqrt( pow( freco_startx[muon_pos]- fstart_x[muon_pos] ,2) + pow( freco_starty[muon_pos]- fstart_y[muon_pos] ,2) + pow( freco_startz[muon_pos]- fstart_z[muon_pos] ,2) ) );
    		for (unsigned j=0; j<fpdg.size(); j++) {
	    	if ( !reco_muon ) break; //select events with a reco muon
	    	if ( fpdg[j]!=2212 && fp0[j] <= 0.2 && !fis_tracked[j]) continue; //watch only reco protons with p>0.2Gev/c
   		hproton_pos_res_badprotons->Fill(  sqrt( pow( freco_startx[j]- fstart_x[j] ,2) + pow( freco_starty[j]- fstart_y[j] ,2) + pow( freco_startz[j]- fstart_z[j] ,2) ) ); 
   		if (is_lowmomentum_p) hproton_pos_res_lowprotons->Fill(  sqrt( pow( freco_startx[j]- fstart_x[j] ,2) + pow( freco_starty[j]- fstart_y[j] ,2) + pow( freco_startz[j]- fstart_z[j] ,2) ) ); 
		}
	}

	if ( reco_muon && count_tracked == 0 && count_not_tracked== 0 && !is_pion && is_lowmomentum_p ) { //no protons tracked but low energy ones
		for (unsigned jj=1; jj<fmuon_dqdx.size()-1; jj++) {
				h_dqdx_low_protons->Fill( fmuon_range -  fmuon_residual[jj], fmuon_dqdx[jj] );
				h_dqdx_low_protons_service->Fill( fmuon_range -  fmuon_residual[jj], fmuon_dqdx[jj] );
		}
		TH1D* h1 = NULL;
		h1 = h_dqdx_low_protons_service->ProjectionY("h",1,length_cut);
		if (h1)
		htail_to_tot_low_protons->Fill(h1->Integral(low_edge,high_edge)/h1->Integral(1,high_edge));
		h_dqdx_low_protons_service->Reset();
		
		hmuon_pos_res_lowprotons->Fill( sqrt( pow( freco_startx[muon_pos]- fstart_x[muon_pos] ,2) + pow( freco_starty[muon_pos]- fstart_y[muon_pos] ,2) + pow( freco_startz[muon_pos]- fstart_z[muon_pos] ,2) ) );
	}


}



void recohelper::RecoBenchmarker::FillCumulativeHistograms() {

    for (long ii=1; ii<1000;ii++) {
    TH1D* h = h_dqdx_not_merged->ProjectionY("h",1,ii);
    h_dqdx_tailtotot_length_not_merged->Fill(  ii, h->Integral(low_edge,high_edge)/h->Integral(1,high_edge) );
    h = h_dqdx_merged->ProjectionY("h",1,ii);
    h_dqdx_tailtotot_length_merged->Fill(  ii,  h->Integral(low_edge,high_edge)/h->Integral(1,high_edge));
    }

}


void recohelper::RecoBenchmarker::AllocateRecoVectors() {
	fis_tracked.push_back(0);
	fmatch_multiplicity.push_back(0);
	flength_reco.push_back(-1);
	freco_momentum_mcs.push_back(-1);
	freco_momentum_mcs_llhd.push_back(-1);
	freco_momentum_range.push_back(-1);
	fpurity.push_back(-1);
	fcompleteness.push_back(-1);
	fnhits.push_back(-1);
	freco_kinE.push_back(-1);
	freco_startx.push_back(-1);
	freco_starty.push_back(-1);
	freco_startz.push_back(-1);
	freco_endx.push_back(-1);
	freco_endy.push_back(-1);
	freco_endz.push_back(-1);
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

  FillCumulativeHistograms();

}

void recohelper::RecoBenchmarker::clear_vectors(){
    //info on the MC particle
    fmuon_dqdx.clear();
    fmuon_dedx.clear();
    fmuon_residual.clear();
    flength.clear();
    fstartT.clear();
    fstart_x.clear();
    fstart_y.clear();
    fstart_z.clear();
    fend_x.clear();
    fend_y.clear();
    fend_z.clear();
  fn_steps.clear();
  fpdg.clear();
  fmother_id.clear();
  fg4_id.clear(); 
   fp0.clear();//initial momentum
   fp0x.clear();
   fp0y.clear();
   fp0z.clear();
   fkinE.clear();
   fcostheta_muon.clear();
  fis_leading.clear();
  fparticle_count.clear();


   //info coming from the tracking algorithm - when there is mc truth
  fis_tracked.clear();
  fmatch_multiplicity.clear();
   fcostheta_muon_reco.clear(); //MCtruth info
   freco_momentum_mcs.clear(); //(GeV) MCS
   freco_momentum_mcs_llhd.clear(); //(GeV) MCS LLHD
   freco_momentum_range.clear(); //MeV
   freco_startx.clear();
   freco_starty.clear();
   freco_startz.clear();
   freco_endx.clear();
   freco_endy.clear();
   freco_endz.clear();
    flength_reco.clear();
   fpurity.clear();
   fcompleteness.clear();
   fnhits.clear();
   fkinE.clear();
   
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
  
  /* //info coming from the tracking algorithm - when there is NO mc truth
  ffake_is_tracked.clear();
  ffake_is_mismatched.clear();
   ffake_costheta_muon_reco.clear(); //MCtruth info
 ffake_pdg_reco.clear();
   ffake_length_reco.clear();
   ffake_reco_momentum_mcs.clear(); //(GeV) MCS
   ffake_reco_momentum_mcs_llhd.clear(); //(GeV) MCS LLHD
   ffake_reco_momentum_range.clear(); //MeV
   ffake_purity.clear();
   ffake_completeness.clear();
   ffake_nhits.clear();
   ffake_kinE.clear();
*/
}

DEFINE_ART_MODULE(recohelper::RecoBenchmarker)
