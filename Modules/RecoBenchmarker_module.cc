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
//#include "art/Framework/IO/Root/RootInputFile.h"
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
#include "larcore/Geometry/WireGeo.h"
#include "lardataobj/MCBase/MCHitCollection.h"
#include "larpandora/LArPandoraInterface/LArPandora.h"

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"
#include <algorithm>

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
    float fneutrino_x;
    float fneutrino_y;
    float fneutrino_z;

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
   float fnu_reco_x;
   float fnu_reco_y;
   float fnu_reco_z;
   float fnu_reco_fitter_x;
   float fnu_reco_fitter_y;
   float fnu_reco_fitter_z;
   float fnu_reco_fitter_chi2ndf;
   float fnu_reco_fitter_chi2;
   std::vector<bool> fis_tracked;
   std::vector<bool> fis_shower_matched;
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
   std::vector<int> freco_trackid;
   std::vector<float> freco_vertex_x;
   std::vector<float> freco_vertex_y;
   std::vector<float> freco_vertex_z;
   std::vector<float> freco_vertexfitter_x;
   std::vector<float> freco_vertexfitter_y;
   std::vector<float> freco_vertexfitter_z;
   std::vector<float> freco_vertexfitter_chi2ndf;
   std::vector<float> freco_vertexfitter_chi2;
   
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

   //hit analysis
   std::vector<int> freco_track_hits;
   std::vector<float> freco_track_collection_hits;
   std::vector<float> freco_track_collection_charge;
   std::vector<int> freco_mcp_hits;
   std::vector<float> freco_mcp_collection_hits;
   std::vector<float> freco_mcp_collection_charge;
   std::vector<float> fnot_clustered;
   std::vector< std::vector<float>> fnot_clustered_charge;
   std::vector<float> fclustered;
   std::vector< std::vector<float>> fclustered_charge;
   std::vector<float> fclustered_matched;
   std::vector< std::vector<float>> fclustered_matched_charge;
   std::vector<float> fclustered_mismatched;
   std::vector< std::vector<float>> fclustered_mismatched_charge;
   std::vector< std::vector<int>> fhit_mismatch_pdg;

   //service variables
   long tot_n_protons; //protons > 20MeV
   long n_all_protons; //all protons
   long low_protons; // protons <20MeV
   long n_muons; // total reco muons

   void AllocateAnalysisHistograms();
   int FillAnalysisHistograms(int&, int&, bool&);
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
    TH1D* hproton_multi_all;
    TH1D* hproton_leading_kinE;
    TH1D* hproton_multi_above20MeV;
    TH1D* hproton_multi_below20MeV;
    TH1D* hproton_merged_not_merged;
   
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
   TH1D* hmuon_spectrum_eff;
   TH1D* hmuon_length;
   TH1D* hmuon_length_all;
   TH1D* hmuon_length_eff;
   TH1D* hproton_kinE;
   TH1D* hproton_kinE_all;
   TH1D* hproton_kinE_eff;
   TH1D* hproton_p;
   TH1D* hproton_p_all;
   TH1D* hproton_p_eff;
   TH1D* hproton_l;
   TH1D* hproton_l_all;
   TH1D* hproton_l_eff;
   TH1D* hproton_kinE_tracked_angle1;
   TH1D* hproton_kinE_all_angle1;
   TH1D* hproton_kinE_tracked_angle2;
   TH1D* hproton_kinE_all_angle2;
   TH1D* hproton_kinE_tracked_angle3;
   TH1D* hproton_kinE_all_angle3;
   TH1D* hproton_l_tracked_angle1;
   TH1D* hproton_l_all_angle1;
   TH1D* hproton_l_tracked_angle2;
   TH1D* hproton_l_all_angle2;
   TH1D* hproton_l_tracked_angle3;
   TH1D* hproton_l_all_angle3;
   TH1D* hproton_nhits_tracked_angle1;
   TH1D* hproton_nhits_all_angle1;
   TH1D* hproton_nhits_tracked_angle2;
   TH1D* hproton_nhits_all_angle2;
   TH1D* hproton_nhits_tracked_angle3;
   TH1D* hproton_nhits_all_angle3;
   
   TH1D* hproton_nhits;
   TH1D* hproton_nhits_CP;
   TH2D* hproton_nhits_theta_mu;
   TH2D* hproton_nhits_CP_theta_mu;
   
   TH1D* h_pmu_end_not_tracked;
   TH1D* h_pmu_end_tracked;
   TH1D* h_theta_mu_tracked;
   TH1D* h_theta_mu_not_tracked;
   TH1D* h_theta_mu;
   TH1D* hproton_theta_mu;
   TH1D* hproton_theta_mu_eff;
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
   //vertexes
   TH1D* h_vertex_resolution_neutrino;
   TH1D* h_vertex_resolution_proton;
   TH1D* h_vertex_resolution_muon;
   TH1D* h_vertex_resolution_neutrino_not_merged;
   TH1D* h_vertex_resolution_proton_not_merged;
   TH1D* h_vertex_resolution_muon_not_merged;
   TH1D* h_vertex_resolution_neutrino_merged;
   TH1D* h_vertex_resolution_proton_merged;
   TH1D* h_vertex_resolution_muon_merged;
   TH2D* h_vertex_resolution_vs_not_tracked_above20MeV;
   TH2D* h_vertex_resolution_vs_not_tracked_below20MeV;
   TH2D* h_vertex_resolution_vs_not_tracked;
   
   //vertex fitter
   TH1D* h_vertexfitter_resolution_neutrino;
   TH1D* h_vertexfitter_chi2ndf_neutrino;
   TH1D* h_vertexfitter_resolution_proton;
   TH1D* h_vertexfitter_chi2ndf_proton;
   TH1D* h_vertexfitter_resolution_muon;
   TH1D* h_vertexfitter_chi2ndf_muon;
   TH1D* h_vertexfitter_resolution_neutrino_not_merged;
   TH1D* h_vertexfitter_chi2ndf_neutrino_not_merged;
   TH1D* h_vertexfitter_chi2ndf_proton_not_merged;
   TH1D* h_vertexfitter_resolution_muon_not_merged;
   TH1D* h_vertexfitter_chi2ndf_muon_not_merged;
   TH1D* h_vertexfitter_resolution_neutrino_merged;
   TH1D* h_vertexfitter_chi2ndf_neutrino_merged;
   TH1D* h_vertexfitter_chi2ndf_proton_merged;
   TH1D* h_vertexfitter_resolution_muon_merged;
   TH1D* h_vertexfitter_chi2ndf_muon_merged;
   TH1D* h_vertexfitter_resolution_proton_not_merged;
   TH1D* h_vertexfitter_resolution_proton_merged;
   TH2D* h_vertexfitter_resolution_vs_not_tracked_above20MeV;
   TH2D* h_vertexfitter_resolution_vs_not_tracked_below20MeV;
   TH2D* h_vertexfitter_resolution_vs_not_tracked;
  
   //shower monitoring
   TH1D* h_shower_pdg;
   TH1D* h_n_proton_showers;
   TH1D* h_shower_proton_kinE;
   TH1D* h_shower_proton_l;
   TH1D* h_shower_proton_nhits;
   TH1D* h_shower_proton_costheta_muon;
  
   //hits studies
   TH1D* h_tracked_not_clustered_distance_nuvtx;
   TH1D* h_tracked_not_clustered_muon_start;
   TH1D* h_tracked_not_clustered_muon_end;
   TH2D* h_tracked_not_clustered_proton_start;
   TH2D* h_tracked_not_clustered_proton_end;
   TH1D* h_hits_not_clustered_tracked_charge;
   TH1D* h_hits_not_clustered_tracked_charge_proton;
   TH1D* h_hits_not_clustered_tracked_charge_muon;
   TH1D* h_fraction_pdgs_not_tracked_not_clustered;
   TH1D* h_not_tracked_not_clustered_muon_start;
   TH1D* h_not_tracked_not_clustered_muon_end;
   TH2D* h_not_tracked_not_clustered_proton_start;
   TH2D* h_not_tracked_not_clustered_proton_end;
   TH1D* h_hits_not_clustered_not_tracked_charge_muon;
   TH1D* h_hits_not_clustered_not_tracked_charge_proton;
   
   TH1D* h_muon_clustering_prob_good_protons;
   TH1D* h_muon_clustering_mismatch_pdg_good_protons;
   TH2D* h_muon_not_clustered_reco_hits_good_protons;
   TH2D* h_muon_NC_lateral_hits_good_protons;
   TH2D* h_muon_NC_costheta_hits_good_protons;
   TH2D* h_muon_clustered_matched_reco_hits_good_protons;
   TH2D* h_muon_CMA_lateral_hits_good_protons;
   TH2D* h_muon_CMA_costheta_hits_good_protons;
   TH2D* h_muon_clustered_mismatched_reco_hits_good_protons;
   TH2D* h_muon_CMI_lateral_hits_good_protons;
   TH2D* h_muon_CMI_costheta_hits_good_protons;
   TH2D* h_muon_not_clustered_reco_charge_good_protons;
   TH2D* h_muon_NC_lateral_charge_good_protons;
   TH2D* h_muon_NC_costheta_charge_good_protons;
   TH2D* h_muon_clustered_matched_reco_charge_good_protons;
   TH2D* h_muon_CMA_lateral_charge_good_protons;
   TH2D* h_muon_CMA_costheta_charge_good_protons;
   TH2D* h_muon_clustered_mismatched_reco_charge_good_protons;
   TH2D* h_muon_CMI_lateral_charge_good_protons;
   TH2D* h_muon_CMI_costheta_charge_good_protons;

   TH1D* h_proton_clustering_prob_good_protons;
   TH1D* h_proton_clustering_mismatch_pdg_good_protons;
   TH2D* h_proton_not_clustered_reco_hits_good_protons;
   TH2D* h_proton_NC_lateral_hits_good_protons;
   TH2D* h_proton_NC_costheta_hits_good_protons;
   TH2D* h_proton_clustered_matched_reco_hits_good_protons;
   TH2D* h_proton_CMA_lateral_hits_good_protons;
   TH2D* h_proton_CMA_costheta_hits_good_protons;
   TH2D* h_proton_clustered_mismatched_reco_hits_good_protons;
   TH2D* h_proton_CMI_lateral_hits_good_protons;
   TH2D* h_proton_CMI_costheta_hits_good_protons;
   TH2D* h_proton_not_clustered_reco_charge_good_protons;
   TH2D* h_proton_NC_lateral_charge_good_protons;
   TH2D* h_proton_NC_costheta_charge_good_protons;
   TH2D* h_proton_clustered_matched_reco_charge_good_protons;
   TH2D* h_proton_CMA_lateral_charge_good_protons;
   TH2D* h_proton_CMA_costheta_charge_good_protons;
   TH2D* h_proton_clustered_mismatched_reco_charge_good_protons;
   TH2D* h_proton_CMI_lateral_charge_good_protons;
   TH2D* h_proton_CMI_costheta_charge_good_protons;

   TH1D* h_muon_clustering_prob_bad_protons;
   TH1D* h_muon_clustering_mismatch_pdg_bad_protons;
   TH2D* h_muon_not_clustered_reco_hits_bad_protons;
   TH2D* h_muon_NC_lateral_hits_bad_protons;
   TH2D* h_muon_NC_costheta_hits_bad_protons;
   TH2D* h_muon_clustered_matched_reco_hits_bad_protons;
   TH2D* h_muon_CMA_lateral_hits_bad_protons;
   TH2D* h_muon_CMA_costheta_hits_bad_protons;
   TH2D* h_muon_clustered_mismatched_reco_hits_bad_protons;
   TH2D* h_muon_CMI_lateral_hits_bad_protons;
   TH2D* h_muon_CMI_costheta_hits_bad_protons;
   TH2D* h_muon_not_clustered_reco_charge_bad_protons;
   TH2D* h_muon_NC_lateral_charge_bad_protons;
   TH2D* h_muon_NC_costheta_charge_bad_protons;
   TH2D* h_muon_clustered_matched_reco_charge_bad_protons;
   TH2D* h_muon_CMA_lateral_charge_bad_protons;
   TH2D* h_muon_CMA_costheta_charge_bad_protons;
   TH2D* h_muon_clustered_mismatched_reco_charge_bad_protons;
   TH2D* h_muon_CMI_lateral_charge_bad_protons;
   TH2D* h_muon_CMI_costheta_charge_bad_protons;

   TH1D* h_proton_clustering_prob_bad_protons;
   TH1D* h_proton_clustering_mismatch_pdg_bad_protons;
   TH2D* h_proton_not_clustered_reco_hits_bad_protons;
   TH2D* h_proton_NC_lateral_hits_bad_protons;
   TH2D* h_proton_NC_costheta_hits_bad_protons;
   TH2D* h_proton_clustered_matched_reco_hits_bad_protons;
   TH2D* h_proton_CMA_lateral_hits_bad_protons;
   TH2D* h_proton_CMA_costheta_hits_bad_protons;
   TH2D* h_proton_clustered_mismatched_reco_hits_bad_protons;
   TH2D* h_proton_CMI_lateral_hits_bad_protons;
   TH2D* h_proton_CMI_costheta_hits_bad_protons;
   TH2D* h_proton_not_clustered_reco_charge_bad_protons;
   TH2D* h_proton_NC_lateral_charge_bad_protons;
   TH2D* h_proton_NC_costheta_charge_bad_protons;
   TH2D* h_proton_clustered_matched_reco_charge_bad_protons;
   TH2D* h_proton_CMA_lateral_charge_bad_protons;
   TH2D* h_proton_CMA_costheta_charge_bad_protons;
   TH2D* h_proton_clustered_mismatched_reco_charge_bad_protons;
   TH2D* h_proton_CMI_lateral_charge_bad_protons;
   TH2D* h_proton_CMI_costheta_charge_bad_protons;

   TH1D* h_proton_clustering_prob_bad_protons_not_tracked;
   TH1D* h_proton_clustering_mismatch_pdg_bad_protons_not_tracked;
   TH2D* h_proton_not_clustered_reco_hits_bad_protons_not_tracked;
   TH2D* h_proton_NC_lateral_hits_bad_protons_not_tracked;
   TH2D* h_proton_NC_costheta_hits_bad_protons_not_tracked;
   TH2D* h_proton_clustered_matched_reco_hits_bad_protons_not_tracked;
   TH2D* h_proton_CMA_lateral_hits_bad_protons_not_tracked;
   TH2D* h_proton_CMA_costheta_hits_bad_protons_not_tracked;
   TH2D* h_proton_clustered_mismatched_reco_hits_bad_protons_not_tracked;
   TH2D* h_proton_CMI_lateral_hits_bad_protons_not_tracked;
   TH2D* h_proton_CMI_costheta_hits_bad_protons_not_tracked;
   TH2D* h_proton_not_clustered_reco_charge_bad_protons_not_tracked;
   TH2D* h_proton_NC_lateral_charge_bad_protons_not_tracked;
   TH2D* h_proton_NC_costheta_charge_bad_protons_not_tracked;
   TH2D* h_proton_clustered_matched_reco_charge_bad_protons_not_tracked;
   TH2D* h_proton_CMA_lateral_charge_bad_protons_not_tracked;
   TH2D* h_proton_CMA_costheta_charge_bad_protons_not_tracked;
   TH2D* h_proton_clustered_mismatched_reco_charge_bad_protons_not_tracked;
   TH2D* h_proton_CMI_lateral_charge_bad_protons_not_tracked;
   TH2D* h_proton_CMI_costheta_charge_bad_protons_not_tracked;


   TH1D* h_muon_clustering_prob_low_protons;
   TH1D* h_muon_clustering_mismatch_pdg_low_protons;
   TH2D* h_muon_not_clustered_reco_hits_low_protons;
   TH2D* h_muon_NC_lateral_hits_low_protons;
   TH2D* h_muon_NC_costheta_hits_low_protons;
   TH2D* h_muon_clustered_matched_reco_hits_low_protons;
   TH2D* h_muon_CMA_lateral_hits_low_protons;
   TH2D* h_muon_CMA_costheta_hits_low_protons;
   TH2D* h_muon_clustered_mismatched_reco_hits_low_protons;
   TH2D* h_muon_CMI_lateral_hits_low_protons;
   TH2D* h_muon_CMI_costheta_hits_low_protons;
   TH2D* h_muon_not_clustered_reco_charge_low_protons;
   TH2D* h_muon_NC_lateral_charge_low_protons;
   TH2D* h_muon_NC_costheta_charge_low_protons;
   TH2D* h_muon_clustered_matched_reco_charge_low_protons;
   TH2D* h_muon_CMA_lateral_charge_low_protons;
   TH2D* h_muon_CMA_costheta_charge_low_protons;
   TH2D* h_muon_clustered_mismatched_reco_charge_low_protons;
   TH2D* h_muon_CMI_lateral_charge_low_protons;
   TH2D* h_muon_CMI_costheta_charge_low_protons;


   int n_events;

    void AllocateRecoVectors();
    void clear_vectors();
};


recohelper::RecoBenchmarker::RecoBenchmarker(fhicl::ParameterSet const & p)
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
  fMCTruthLabel = p.get<std::string> ("MCTruthLabel");
  fG4TruthLabel = p.get<std::string> ("G4TruthLabel");
  fShowerTruthLabel = p.get<std::string> ("ShowerTruthLabel");
  fClusterLabel = p.get<std::string> ("ClusterLabel");
  fHitLabel = p.get<std::string> ("HitLabel");
  fMCHitLabel = p.get<std::string> ("MCHitLabel");

}

void recohelper::RecoBenchmarker::beginJob()
{


	//std::cout << ">>>>>>>" << art::RootInputFile::fileName() << std::endl;
  recoTree = tfs->make<TTree>("recotree", "recotree");

  // define branches
  recoTree->Branch("Event", &fEvent, "Event/I");
  recoTree->Branch("SubRun", &fSubRun, "SubRun/I");
  recoTree->Branch("Run", &fRun, "Run/I");

  //mctruth info
  recoTree->Branch("ccnc", &fccnc);
  recoTree->Branch("interaction", &finteraction);
  recoTree->Branch("neutrino_x", &fneutrino_x);
  recoTree->Branch("neutrino_y", &fneutrino_y);
  recoTree->Branch("neutrino_z", &fneutrino_z);

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
  recoTree->Branch("nu_reco_x",&fnu_reco_x);
  recoTree->Branch("nu_reco_y",&fnu_reco_y);
  recoTree->Branch("nu_reco_z",&fnu_reco_z);
  recoTree->Branch("nu_reco_fitter_x",&fnu_reco_fitter_x);
  recoTree->Branch("nu_reco_fitter_y",&fnu_reco_fitter_y);
  recoTree->Branch("nu_reco_fitter_z",&fnu_reco_fitter_z);
  recoTree->Branch("nu_reco_fitter_chi2ndf",&fnu_reco_fitter_chi2ndf);
  recoTree->Branch("nu_reco_fitter_chi2",&fnu_reco_fitter_chi2);
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
  recoTree->Branch("reco_trackid",&freco_trackid);
  recoTree->Branch("reco_vertex_x",&freco_vertex_x);
  recoTree->Branch("reco_vertex_y",&freco_vertex_y);
  recoTree->Branch("reco_vertex_z",&freco_vertex_z);
  recoTree->Branch("reco_vertexfitter_x",&freco_vertexfitter_x);
  recoTree->Branch("reco_vertexfitter_y",&freco_vertexfitter_y);
  recoTree->Branch("reco_vertexfitter_z",&freco_vertexfitter_z);
  recoTree->Branch("reco_vertexfitter_chi2ndf",&freco_vertexfitter_chi2ndf);
  recoTree->Branch("reco_vertexfitter_chi2",&freco_vertexfitter_chi2);
  recoTree->Branch("is_shower_matched",&fis_shower_matched);

  //hits analysis
  recoTree->Branch("reco_track_hits",&freco_track_hits);
  recoTree->Branch("reco_track_collection_hits",&freco_track_collection_hits);
  recoTree->Branch("reco_track_collection_charge",&freco_track_collection_charge);
  recoTree->Branch("reco_mcp_hits",&freco_mcp_hits);
  recoTree->Branch("reco_mcp_collection_hits",&freco_mcp_collection_hits);
  recoTree->Branch("reco_mcp_collection_charge",&freco_mcp_collection_charge);
  recoTree->Branch("not_clustered",&fnot_clustered);
  recoTree->Branch("not_clustered_charge",&fnot_clustered_charge);
  recoTree->Branch("clustered",&fclustered);
  recoTree->Branch("clustered_charge",&fclustered_charge);
  recoTree->Branch("clustered_matched",&fclustered_matched);
  recoTree->Branch("clustered_matched_charge",&fclustered_matched_charge);
  recoTree->Branch("clustered_mismatched",&fclustered_mismatched);
  recoTree->Branch("clustered_mismatched_charge",&fclustered_mismatched_charge);
  recoTree->Branch("hit_mismatch_pdg",&fhit_mismatch_pdg);


  AllocateAnalysisHistograms();

  n_events = 0;
  tot_n_protons = 0;
  n_muons = 0;
  n_all_protons = 0;
  low_protons = 0;
   
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
  std::vector< art::Ptr<recob::Shower> >  showerList;  
  art::fill_ptr_vector(showerList, showerHandle); 

  art::FindManyP<simb::MCParticle, anab::BackTrackerMatchingData> mcpsFromShowers(showerHandle, e, fShowerTruthLabel);

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

  art::ValidHandle< std::vector< recob::PFParticle > > pfpHandle = 
    e.getValidHandle< std::vector< recob::PFParticle > >(fPfpLabel);
  std::vector< art::Ptr<recob::PFParticle> > pfpList;
  art::fill_ptr_vector( pfpList, pfpHandle);

  art::FindManyP<recob::Hit> hitsFromTracks( trackHandle, e, fTrackLabel);
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

  //loop on all MC truth frames (mostly 1 per event)
  for ( unsigned n_truth = 0; n_truth < mcTruth.size(); n_truth++ ) {
	  //std::cout << ">>>>>>>>>>>>>>>>>EVENT" << std::endl; 
  n_events++;
  clear_vectors();

#if isSingleParticle == 0 
  //check that the neutrino is in the active volume
  const simb::MCNeutrino thisNeutrino = mcTruth[n_truth]->GetNeutrino();
  const simb::MCParticle thisLepton = thisNeutrino.Lepton();
  if ( !_rbutilInstance.isInTPC( thisLepton ) ) continue;
  
  //store the interaction info
  fccnc = thisNeutrino.CCNC();
  finteraction = thisNeutrino.InteractionType();
 
  double xOffset = 0.7 - SCE->GetPosOffsets( thisNeutrino.Nu().Position().X(),  thisNeutrino.Nu().Position().Y(), thisNeutrino.Nu().Position().Z()  )[0];
  double yOffset = SCE->GetPosOffsets( thisNeutrino.Nu().Position().X(),  thisNeutrino.Nu().Position().Y(), thisNeutrino.Nu().Position().Z() )[1];
  double zOffset = SCE->GetPosOffsets( thisNeutrino.Nu().Position().X(),  thisNeutrino.Nu().Position().Y(), thisNeutrino.Nu().Position().Z() )[2];
  if (!space_charge) {
	  xOffset=0;
	  yOffset=0;
	  zOffset=0;
  }
  fneutrino_x = thisNeutrino.Nu().Position().X() + xOffset;
  fneutrino_y = thisNeutrino.Nu().Position().Y() + yOffset;
  fneutrino_z = thisNeutrino.Nu().Position().Z() + zOffset;

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
        xOffset = 0.7 - SCE->GetPosOffsets( thisMcp->Position().X(),  thisMcp->Position().Y(), thisMcp->Position().Z() )[0];
        yOffset = SCE->GetPosOffsets( thisMcp->Position().X(),  thisMcp->Position().Y(), thisMcp->Position().Z() )[1];
        zOffset = SCE->GetPosOffsets( thisMcp->Position().X(),  thisMcp->Position().Y(), thisMcp->Position().Z() )[2];
	if (!space_charge) {
		xOffset = 0;
		yOffset = 0;
		zOffset = 0;
	}
	fstart_x.push_back ( thisMcp->Position().X() + xOffset );
	fstart_y.push_back ( thisMcp->Position().Y() + yOffset );
	fstart_z.push_back ( thisMcp->Position().Z() + zOffset );
	
	if (space_charge) {
        xOffset = 0.7 - SCE->GetPosOffsets( thisMcp->Position().X(),  thisMcp->Position().Y(), thisMcp->Position().Z() )[0];
        yOffset = SCE->GetPosOffsets( thisMcp->Position().X(),  thisMcp->Position().Y(), thisMcp->Position().Z() )[1];
        zOffset = SCE->GetPosOffsets( thisMcp->Position().X(),  thisMcp->Position().Y(), thisMcp->Position().Z() )[2];
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
	    if ( fkinE[jj] >= current_kinE )
	       is_leading = false;
	    else 
	       fis_leading[jj] = false;
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
  
  //this stays out of the previous loop but it is somewhat connected
   for ( auto const& mcp : mcList ) {
	   auto itt = std::find ( fg4_id.begin(), fg4_id.end(),  mcp->TrackId());
	   if ( itt == fg4_id.end() ) continue;
  //nhits from MCP
	std::vector< art::Ptr < recob::Hit> > hits_mcp = hitsFromMCP.at( &mcp - &mcList[0] );
	//std::cout <<  itt - fg4_id.begin() << " " << freco_mcp_hits.size() << std::endl;
	freco_mcp_hits[ itt - fg4_id.begin() ] = hits_mcp.size() ;
	
	int n_collection_hits = 0;
	int total_hit_charge = 0;
	//std::cout << "NUMBER OF HITS " << hits_mcp.size() << std::endl;
	fnhits[ itt - fg4_id.begin() ] = 0;
	if ( hits_mcp.size() ) {
	//mcp hits
	for ( auto const& iter_hit : hits_mcp ){
		if ( iter_hit->View() == geo::kW ){	
			n_collection_hits++;
			total_hit_charge += iter_hit->Integral();
		}
		fnhits[ itt - fg4_id.begin() ] = fnhits[ itt - fg4_id.begin() ] + 1;
	}
	freco_mcp_collection_hits[itt - fg4_id.begin() ] =  n_collection_hits ;
	freco_mcp_collection_charge[itt - fg4_id.begin()] = total_hit_charge ;
	}
   }



  bool lleading = false;
  //loop on the just saved truth info and fill a bunch of summary histograms
  int count_protons = 0;
  int count_protons_above = 0;
  int count_protons_below = 0;
  for ( unsigned ii = 0; ii < fpdg.size(); ii++ ) {
	if ( fpdg[ii] == 2212 && fis_leading[ii] ) {
		hproton_leading_kinE->Fill( fkinE[ii] * 1000. ); //bring it to MeV
		if (lleading==true) std::cout << "ERROR!!!! CAN'T BE ALREADY TRUE" << std::endl;
		lleading = true;
	}
	if ( fpdg[ii] == 2212 ) count_protons++;
	if ( fpdg[ii] == 2212 && fkinE[ii] >= 0.02 ) count_protons_above++;
	if ( fpdg[ii] == 2212 && fkinE[ii] < 0.02 ) count_protons_below++;
  }
	hproton_multi_all->Fill( count_protons );
	hproton_multi_above20MeV->Fill( count_protons_above );
	hproton_multi_below20MeV->Fill( count_protons_below );

  //----------------------------
  // Tracks
  //----------------------------

  //I can't use thisTrack->ID() as index, because in case of Kalman Fitter failure, there could be mismatches between the number
  //of elements in the associations and the track index. So I need to index manually. 
  //There is an additional crosscheck for which the trackID must be increasing in the loop, just to be sure that nothing nasty is happening.
  int track_id_counter = 0;
  int previous_track_id = -1;

  //int all_hits_tracks = 0;
  // loop tracks and do truth matching 
  for (auto const& thisTrack : trackList) {
    
    std::vector< art::Ptr<simb::MCParticle> > mcps = mcpsFromTracks.at(track_id_counter);
    std::vector< art::Ptr<anab::Calorimetry> > calos = caloFromTracks.at(track_id_counter);

    if (mcps.size() >1 ) mf::LogWarning(__FUNCTION__) << "Warning !!! More than 1 MCparticle associated to the same track!" << std::endl;
    if (calos.size() != 3 ) mf::LogWarning(__FUNCTION__) << "Warning !!! Calorimetry info size " << calos.size() << " != 3 associated to tracks!" << std::endl;
    if ( previous_track_id >= thisTrack->ID() ) mf::LogError(__FUNCTION__) << "ERROR! The Track ID's are not in ascending order! " << std::endl;

	
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
    if (fis_tracked[pos] == true ) {
	    std::cout << "Probably broken track!" << std::endl;
	    fmatch_multiplicity[pos] = fmatch_multiplicity[pos] + 1;
    	    //decide if keeping the old info or the new one - based on the minimum distance between real start and recoed end/start
	    float start_distance = sqrt( pow( thisTrack->Vertex().X() - fstart_x[pos] , 2 ) 
		      +  pow( thisTrack->Vertex().Y() - fstart_y[pos] , 2 )
		      +   pow( thisTrack->Vertex().Z() - fstart_z[pos] , 2 ) );
	    float end_distance = sqrt( pow( thisTrack->End().X() - fstart_x[pos] , 2 ) 
		      +  pow( thisTrack->End().Y() - fstart_y[pos] , 2 )
		      +   pow( thisTrack->End().Z() - fstart_z[pos] , 2 ) );
	    float new_distance = min( start_distance, end_distance );
	    start_distance = sqrt( pow( freco_startx[pos]  - fstart_x[pos] , 2 ) 
		      +  pow( freco_starty[pos] - fstart_y[pos] , 2 )
		      +   pow( freco_startz[pos] - fstart_z[pos] , 2 ) );
	    end_distance = sqrt( pow( freco_endx[pos] - fstart_x[pos] , 2 ) 
		      +  pow( freco_endy[pos] - fstart_y[pos] , 2 )
		      +  pow(freco_endz[pos] - fstart_z[pos] , 2 ) );
	    float old_distance = min( start_distance, end_distance );
	    if ( old_distance < new_distance )
		    continue;
    }
        
        fis_tracked[pos] = true;
	fmatch_multiplicity[pos] = fmatch_multiplicity[pos] + 1;
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
	freco_trackid[pos] = thisTrack->ID();
	
//	std::cout << "TRACK ID " << freco_trackid[pos] << " pos=" << pos << std::endl;
 	
	if (thisMcp->PdgCode() == 13 ) { //for the muon, only when there is the info
		if (fmuon_dqdx.size()!=0) {
			mf::LogError(__FUNCTION__) << "Calorimetry should be filled only once!!!! Clearing..." << std::endl;
			fmuon_dqdx.clear();
			fmuon_dedx.clear();
			fmuon_residual.clear();
		}

		for (size_t position=0; position<calos.at( geo::kCollection )->dQdx().size(); position++) {
		fmuon_dqdx.push_back(calos.at( geo::kCollection )->dQdx()[position]); //look only at collection
		fmuon_dedx.push_back(calos.at( geo::kCollection )->dEdx()[position]); //look only at collection
		fmuon_residual.push_back(calos.at( geo::kCollection )->ResidualRange()[position]); //look only at collection
		}
	fmuon_range = calos.at( geo::kCollection )->Range();
	}
    
	auto thisMcpCleanliness = mcpsFromTracks.data(track_id_counter).at(0)->cleanliness;
	auto thisMcpCompleteness = mcpsFromTracks.data(track_id_counter).at(0)->completeness;
	fpurity[pos] = thisMcpCleanliness;
	fcompleteness[pos] = thisMcpCompleteness;

        //check association of hits <-> tracks
	//nhits from track
	std::vector< art::Ptr < recob::Hit> > hits_tracks = hitsFromTracks.at( track_id_counter );
	freco_track_hits[pos] = hits_tracks.size() ;

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
	freco_track_collection_hits[pos] = n_collection_hits;
	freco_track_collection_charge[pos] = total_hit_charge;
	
    }//MCParticle
        
 
        //get clusters from the reco track, and track down which MCP's contribute to the clustered hits
    
	/*int nhits = 0;
	//std::cout << ">>>>>>>Hits from track " << hits_tracks.size() << " hits from MCP " << nhits << std::endl;
	std::cout << ">>>>>>>Hits from track " << hits_tracks.size() << " hits from MCP " << hits_mcp.size() << std::endl;

	std::vector< art::Ptr < recob::Hit> > hits_tracks = hitsFromTracks.at( track_id_counter );
	all_hits_tracks+=hits_tracks.size();
  
    //study the hits and the clustering
  art::ValidHandle< std::vector< recob::Cluster > > clusterHandle = 
    e.getValidHandle< std::vector< recob::Cluster > >(fClusterLabel);

  art::FindManyP<recob::Hit> hitsFromClusters(clusterHandle, e, fClusterLabel);
  
  art::ValidHandle< std::vector< recob::Hit > > hitHandle = 
    e.getValidHandle< std::vector< recob::Hit > >(fHitLabel);
  std::vector< recob::Hit > hitList;
  art::fill_ptr_vector( hitList, hitHandle);

  art::FindManyP<recob::Hit> hitsFromTracks( trackHandle, e, fTrackLabel);
  art::FindManyP<recob::Hit, anab::BackTrackerMatchingData> hitsFromMCP( mcParticleHandle, e, fTrackTruthLabel);
*/
    
    
        track_id_counter++;
  }//Tracks


  //checks on showers
  //I am indexing showers manually because I am lazy and I am copying what I did for tracks
  int shower_id_counter = 0;
  int previous_shower_id = -1;
  int count_proton_showers = 0;
  std::vector<int> mcp_showers_ids;

  //int all_hits_tracks = 0;
  // loop tracks and do truth matching 
  for (auto const& thisShower : showerList) {
    
    std::vector< art::Ptr<simb::MCParticle> > mcps = mcpsFromShowers.at(shower_id_counter);

    if (mcps.size() >1 ) mf::LogWarning(__FUNCTION__) << "Warning !!! More than 1 MCparticle associated to the same shower!" << std::endl;
    if ( previous_shower_id >= thisShower->ID() ) mf::LogError(__FUNCTION__) << "ERROR! The Shower ID's are not in ascending order! " << std::endl;

	
    for (auto const& thisMcp : mcps){
    
    //this if is necessary, since sometimes a particle is matched to a secondary (electron, etc) -> to be checked
    if ( !(std::abs(thisMcp->PdgCode()) == 11 && thisMcp->Mother() == muon_id && abs(thisMcp->Position().X()-muon_endx) < DBL_EPSILON &&
		    abs(thisMcp->Position().Y()-muon_endy) < DBL_EPSILON && abs(thisMcp->Position().Z()-muon_endz) < DBL_EPSILON ) )
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
    auto it_found = std::find( fg4_id.begin(), fg4_id.end(), thisMcp->TrackId() ) ;
    if (it_found==fg4_id.end()) mf::LogError(__FUNCTION__) << "ERROR!!! Matched particle not found!" << std::endl;
    size_t pos = it_found - fg4_id.begin();
	

    //save information on reco track
    if (fis_tracked[pos] == true ) std::cout << "Tracked particle matched to a shower?!?! >>>>>> SEEMS WRONG!!!" << std::endl;
    
    std::cout << "FOUND A SHOWER!!! >>>>>>>>> pdg=" << thisMcp->PdgCode() << std::endl;
    fis_shower_matched[pos] = true;
    h_shower_pdg->Fill( thisMcp->PdgCode() );
	
    if ( fpdg[pos] == 2212 ) { //protons
	    h_shower_proton_l->Fill( flength[pos] );
	    h_shower_proton_kinE->Fill( fkinE[pos] );
	    h_shower_proton_nhits->Fill( fnhits[pos] );
	    h_shower_proton_costheta_muon->Fill( fcostheta_muon[pos] );
	    if ( std::find( mcp_showers_ids.begin(), mcp_showers_ids.end(), fg4_id[pos] ) == mcp_showers_ids.end() )count_proton_showers++;
	    std::cout << ">>>>PROTON SHOWER" << std::endl;
	    //std::cout << ">>>>>>>FILE " << art::RootInputFile::fileName << std::endl;
	    std::cout << "Event=" << fEvent << " Run=" << fRun << " SubRun=" << fSubRun << std::endl;

    }
    mcp_showers_ids.push_back(fg4_id[pos]);
 	
    }//MCParticle
        
        shower_id_counter++;
  }//Shower

	
  h_n_proton_showers->Fill(count_proton_showers);


  //----------------------------
  // Vertex
  //----------------------------
  
  art::FindManyP<recob::Vertex> vertexFromPfp(pfpHandle, e, fPfpLabel);
  art::FindManyP<recob::Track> trackFromPfp(pfpHandle, e, fPfpAssnLabel);
  
  //look for neutrino pfp
  bool neutrino_set = false;
  for ( auto const& pfp : pfpList ) {
	  std::vector< art::Ptr <recob::Vertex> > vertex_pfp = vertexFromPfp.at( &pfp - &pfpList[0] );
	  
	  if ( vertex_pfp.size()==0 ) continue; //no vertex
	  if ( vertex_pfp.size()>1 ) mf::LogError(__FUNCTION__) << "ERROR! There should be only 1 vertex per PFP!" <<std::endl;
	  double xyzz[3];
	  vertex_pfp[0]->XYZ(xyzz);
	  
	  if ( lar_pandora::LArPandoraHelper::IsNeutrino(pfp) ) {
	  	  if (neutrino_set) mf::LogError(__FUNCTION__) << "There should be only 1 neutrino per PFP!!!!" << std::endl;
		  neutrino_set = true;
	  //save info on neutrino reco'ed vertex
	  fnu_reco_x = xyzz[0];
	  fnu_reco_y = xyzz[1];
	  fnu_reco_z = xyzz[2];
	  }
	  
	  //investigate matching particles and tracks
	  std::vector< art::Ptr <recob::Track> > track_pfp = trackFromPfp.at( &pfp - &pfpList[0] );
	  if (track_pfp.size()==0) continue; //no track! (shower?)
	  int trackid = track_pfp[0]->ID();
	  auto iter = std::find( freco_trackid.begin(), freco_trackid.end(), trackid );
	  if ( iter == freco_trackid.end() ) continue; //not found
	  unsigned mc_pos = iter - freco_trackid.begin();
          
	  freco_vertex_x[ mc_pos ] = xyzz[0];
	  freco_vertex_y[ mc_pos ] = xyzz[1];
	  freco_vertex_z[ mc_pos ] = xyzz[2];
  }

  
  //----------------------------
  // Vertex fitter (Giuseppe)
  //----------------------------
 
  if (fIsVertexFitter) {
  art::FindManyP<recob::Vertex> vertexfitterFromPfp(pfpHandle, e, fVertexFitterLabel);
  
  //look for neutrino pfp
  neutrino_set = false;
  for ( auto const& pfp : pfpList ) {
	  std::vector< art::Ptr <recob::Vertex> > vertexfitter_pfp = vertexfitterFromPfp.at( &pfp - &pfpList[0] );
	  
	  if ( vertexfitter_pfp.size()==0 ) continue; //no vertex
	  if ( vertexfitter_pfp.size()>1 ) mf::LogError(__FUNCTION__) << "Vertex Fitter: ERROR! There should be only 1 vertex per PFP!" <<std::endl;
	  double xyzz[3];
	  vertexfitter_pfp[0]->XYZ(xyzz);
	  
	  if ( lar_pandora::LArPandoraHelper::IsNeutrino(pfp) ) {
	  	  if (neutrino_set) mf::LogError(__FUNCTION__) << "Vertex Fitter: There should be only 1 neutrino per PFP!!!!" << std::endl;
		  neutrino_set = true;
	  //save info on neutrino reco'ed vertex
	  fnu_reco_fitter_x = xyzz[0];
	  fnu_reco_fitter_y = xyzz[1];
	  fnu_reco_fitter_z = xyzz[2];
	  //fnu_reco_fitter_chi2ndf = vertexfitter_pfp[0]->chi2PerNdof();
	  //fnu_reco_fitter_chi2 = vertexfitter_pfp[0]->chi2();
	  fnu_reco_fitter_chi2ndf = -1;
	  fnu_reco_fitter_chi2 = -1;
	  }
	  
	  //investigate matching particles and tracks
	  std::vector< art::Ptr <recob::Track> > track_pfp = trackFromPfp.at( &pfp - &pfpList[0] );
	  if (track_pfp.size()==0) continue; //no track! (shower?)
	  int trackid = track_pfp[0]->ID();
	  auto iter = std::find( freco_trackid.begin(), freco_trackid.end(), trackid );
	  if ( iter == freco_trackid.end() ) continue; //not found
	  unsigned mc_pos = iter - freco_trackid.begin();
          
	  freco_vertexfitter_x[ mc_pos ] = xyzz[0];
	  freco_vertexfitter_y[ mc_pos ] = xyzz[1];
	  freco_vertexfitter_z[ mc_pos ] = xyzz[2];
	  //freco_vertexfitter_chi2ndf[ mc_pos ] = vertexfitter_pfp[0]->chi2PerNdof(); 
	  //freco_vertexfitter_chi2[ mc_pos ] = vertexfitter_pfp[0]->chi2(); 
	  freco_vertexfitter_chi2ndf[ mc_pos ] = -1; 
	  freco_vertexfitter_chi2[ mc_pos ] = -1; 
  }
  } //is vertex fitter


  //std::cout << "NEUTRINO " <<  fnu_reco_x << " " << fnu_reco_y << " " << fnu_reco_z << std::endl;

  int count_tracked = 0;
  int count_not_tracked = 0;
  bool is_lowmomentum_p = false;
  int muon_pos = FillAnalysisHistograms( count_tracked, count_not_tracked, is_lowmomentum_p ); //plots tracking and vertexing information
  
  if (count_not_tracked) {
	  std::cout << ">>>>NOT TRACKED!" << std::endl;
	  std::cout << "Number of not tracked " << count_not_tracked << ", number of tracked " << count_tracked << std::endl;
	  std::cout<< "Run: " << fRun << " SubRun: " << fSubRun << " Event: " << fEvent << std::endl;
  }

  art::FindManyP<recob::Cluster> clustersFromHits(hitHandle, e, fClusterLabel);
  art::FindManyP<recob::PFParticle> pfpFromCluster(clusterHandle, e, fClusterLabel);
  art::FindManyP<recob::PFParticle> pfpFromTrack(trackHandle, e, fPfpAssnLabel);
  art::FindManyP<recob::SpacePoint> spacepointFromHits( hitHandle, e, fTrackLabel );
  
  if (muon_pos != -1) {
  //now loop over hits and check:
  //-which MCP the hit belongs to
  //-is the MCP tracked?
  //-if is tracked, check if the hit is clustered and attributed properly to the track
  //-if it is not tracked, save which other particle/track got the hit or if it was not clustered
  
  //std::cout << "SIZE: " << MCPfromhits.size() << std::endl;
  //std::cout << "SIZE hits: " << hitList.size() << std::endl;
  for ( auto const& hit: hitList ) {
	  //only care at collection for now
	if ( hit->View() != geo::kW ) continue;
	std::vector< art::Ptr <simb::MCParticle> > mcparticle_hits = MCPfromhits.at( &hit - &hitList[0] );
	std::vector< art::Ptr <recob::SpacePoint> > spacepoint_hits = spacepointFromHits.at( &hit - &hitList[0] );
	art::ServiceHandle<geo::Geometry> geom;
	double xyz[3] = { 0, 0, 0};
	bool is_spacepoint = false;
	if ( spacepoint_hits.size() != 1) {
		geom->WireIDToWireGeo( hit->WireID() ).GetCenter( xyz, 0 );
//		std::cout << "NO SPACE POINT!" << std::endl;
	} else {
		is_spacepoint = true;
		xyz[0] = spacepoint_hits[0]->XYZ()[0];
		xyz[1] = spacepoint_hits[0]->XYZ()[1];
		xyz[2] = spacepoint_hits[0]->XYZ()[2];
	}
   	 
	//std::cout << "size2 " <<  mcparticle_hits.size() << std::endl;
	if (  mcparticle_hits.size() == 0 ) continue; 
	unsigned index = 0;
	for (unsigned jj=0; jj<mcparticle_hits.size(); jj++) {
			bool is_max = MCPfromhits.data( &hit - &hitList[0] ).at( jj )->isMaxIDE;
			if ( is_max ) {
				index = jj;
				//std::cout << "WARNING! Can't be already !=0" << std::endl;
			}
	}

//	std::cout << "INDEX " << index << std::endl;

	if ( index >= mcparticle_hits.size() ) {
		mf::LogError(__FUNCTION__) << ">>>>>ERROR!!!!! There should always be 1 MCP associated to the hit! " << mcparticle_hits.size() << std::endl;
		exit (-1);
	}

	int trackid = mcparticle_hits[index]->TrackId();
	//lookup the originating MC particle
	bool found = false;
	auto iter = std::find (fg4_id.begin(), fg4_id.end(), trackid );
	if ( iter != fg4_id.end() ) found = true;
	else {
		iter = std::find (fg4_id.begin(), fg4_id.end(), mcparticle_hits[index]->Mother() );
		if ( iter != fg4_id.end() ) found = true;
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

	

	int mcp_index = iter - fg4_id.begin();

	//freco_mcp_hits	
  //mcp index allows me to understand if this particle is tracked or not
      // std::cout << "clusters from hits " << clustersFromHits.size() << std::endl;
//	       std::cout << "pos " << &hit - &hitList[0] << std::endl;
	       std::vector< art::Ptr <recob::Cluster>> cluster_hits = clustersFromHits.at(  &hit - &hitList[0] );
	       if ( cluster_hits.size() > 3) {
		       //std::cout << ">>>>>>>>>>>>>>>>>>>>>>>" << clustersFromHits.size() << std::endl;
		       mf::LogError(__FUNCTION__) << "Number of clusters must be maximum one!!!" << std::endl;
		       exit(-1);
	       }


	       if ( cluster_hits.size() == 0 ) { //this hit is not clustered!
		       
		       fnot_clustered[mcp_index] = fnot_clustered[mcp_index] + 1;
		       fnot_clustered_charge[mcp_index].push_back( hit->Integral() );
		       //tracked particle, with a non-clustered hit
			h_hits_not_clustered_tracked_charge->Fill ( hit->Integral() );
		       if ( fis_tracked[mcp_index] ) { //check and fill histos for a MCP which is tracked
		        if (mcp_index==muon_pos) //is muon
			h_hits_not_clustered_tracked_charge_muon->Fill ( hit->Integral() );
			else if (fpdg[mcp_index] == 2212) //is proton
			h_hits_not_clustered_tracked_charge_proton->Fill ( hit->Integral() );

		        
		        if (is_spacepoint) {
			//	std::cout << "IS SPACEPOINT" << std::endl;
			h_tracked_not_clustered_distance_nuvtx->Fill( sqrt( pow( xyz[0] - fneutrino_x,2) + pow( xyz[1] - fneutrino_y,2) + pow( xyz[2] - fneutrino_z,2)) ); //distance between neutrino vertex z and hit z
			if ( mcp_index==muon_pos ) { //if is muon
			h_tracked_not_clustered_muon_start->Fill( sqrt( pow( xyz[0] - fstart_x[mcp_index],2) 
									+ pow( fstart_y[mcp_index] - xyz[1],2) 
									+ pow( xyz[2] - fstart_z[mcp_index],2)) ); //z distance between muon starting point and not clustered hit
			h_tracked_not_clustered_muon_end->Fill( sqrt( pow( xyz[0] - fend_x[mcp_index],2) 
								      + pow( fend_y[mcp_index] - xyz[1],2) 
								      + pow( xyz[2] - fend_z[mcp_index],2)) ); //z distance between muon ending point and not clustered hit
				} //is muon
			if ( fpdg[mcp_index] == 2212 ) { //if is proton
			h_tracked_not_clustered_proton_start->Fill( sqrt( pow( xyz[0] - fstart_x[mcp_index],2) 
									  + pow( fstart_y[mcp_index] - xyz[1],2) 
									  + pow( xyz[2] - fstart_z[mcp_index],2)), flength[mcp_index] ); //z distance between muon starting point and not clustered hit
			h_tracked_not_clustered_proton_end->Fill( sqrt( pow( xyz[0] - fend_x[mcp_index],2) 
									+ pow( fend_y[mcp_index] - xyz[1],2) 
									+ pow( xyz[2] - fend_z[mcp_index],2)), flength[mcp_index] ); //z distance between muon ending point and not clustered hit
			 }//is proton
			}//is spacepoint
		       } else { //the particle is not tracked
			
				if ( fpdg[mcp_index] ==13 ) {
				h_fraction_pdgs_not_tracked_not_clustered->Fill(0.,1./freco_mcp_collection_hits[mcp_index]); //muon
			  	h_hits_not_clustered_not_tracked_charge_muon->Fill ( hit->Integral() );
				if (is_spacepoint){
			h_not_tracked_not_clustered_muon_start->Fill( sqrt( pow( xyz[0] - fstart_x[mcp_index],2) 
									+ pow( fstart_y[mcp_index] - xyz[1],2) 
									+ pow( xyz[2] - fstart_z[mcp_index],2)) ); //z distance between muon starting point and not clustered hit
			h_not_tracked_not_clustered_muon_end->Fill( sqrt( pow( xyz[0] - fend_x[mcp_index],2) 
								      + pow( fend_y[mcp_index] - xyz[1],2) 
								      + pow( xyz[2] - fend_z[mcp_index],2)) ); //z distance between muon ending point and not clustered hit
				} //is spacepoint
			} else if ( fpdg[mcp_index] == 2212 ) {
				h_fraction_pdgs_not_tracked_not_clustered->Fill(1.,1./freco_mcp_collection_hits[mcp_index]); //proton all
				if ( fp0[mcp_index] > 0.2)
				h_fraction_pdgs_not_tracked_not_clustered->Fill(2.,1./freco_mcp_collection_hits[mcp_index]); //protons > 20MeV
				else
				h_fraction_pdgs_not_tracked_not_clustered->Fill(3.,1./freco_mcp_collection_hits[mcp_index]); //protons <20 MeV
			  	
				h_hits_not_clustered_not_tracked_charge_proton->Fill ( hit->Integral() );
				if (is_spacepoint){
			h_not_tracked_not_clustered_proton_start->Fill( sqrt( pow( xyz[0] - fstart_x[mcp_index],2) 
									+ pow( fstart_y[mcp_index] - xyz[1],2) 
									+ pow( xyz[2] - fstart_z[mcp_index],2)), flength[mcp_index] ); //z distance between muon starting point and not clustered hit
			h_not_tracked_not_clustered_proton_end->Fill( sqrt( pow( xyz[0] - fend_x[mcp_index],2) 
								      + pow( fend_y[mcp_index] - xyz[1],2) 
								      + pow( xyz[2] - fend_z[mcp_index],2)), flength[mcp_index] ); //z distance between muon ending point and not clustered hit
				} //is spacepoint
		       		} //is proton
		       } //not tracked
		
	       } else { //in this case the hit is clustered
			unsigned cluster_index = 0;
			for ( unsigned zz=0; zz<cluster_hits.size(); zz++) {
				if ( cluster_hits[zz]->View() == geo::kW ) cluster_index = zz;
			}
		       
			fclustered[mcp_index] = fclustered[mcp_index] + 1;
			//std::cout << "SIZE1 " << fclustered.size() << " index " << mcp_index << std::endl;
			//std::cout << "SIZE2 " << fclustered_charge.size() << " index " << mcp_index << std::endl;
		        fclustered_charge[mcp_index].push_back( hit->Integral() );
			
			auto cluster_ID = cluster_hits.at(cluster_index)->ID();
			//loop on clusters and find PFP
			size_t pfp_id = -1;
			for ( auto const& cluster : clusterList) {
			if ( cluster->ID() != cluster_ID ) continue; //select the cluster I am interested in
				//std::cout << "SONO QUA" << std::endl;
			std::vector< art::Ptr <recob::PFParticle> > pfp_cluster = pfpFromCluster.at( &cluster - &clusterList[0] );
			if (pfp_cluster.size()!=1) std::cout << "MORE THAN 1 PFP! size=" << pfp_cluster.size() << std::endl;
			if ( pfp_cluster.size()!=0 ) pfp_id = pfp_cluster[0]->Self();
			}
			
			//associate the pfp to the tracks
			for ( auto const& track : trackList) {
				//std::cout << "SONO QUA 1" << std::endl;
			std::vector< art::Ptr <recob::PFParticle> > pfp_track = pfpFromTrack.at( &track - &trackList[0] );
			std::vector< art::Ptr <simb::MCParticle> > mcp_track = mcpsFromTracks.at( &track - &trackList[0] );
			unsigned mm = 0;
			float purityy = -1;
			for ( unsigned nn=0; nn<mcp_track.size(); nn++ ) { //select highest purity
				if ( mcpsFromTracks.data( &track - &trackList[0] ).at( nn )->cleanliness > purityy ) {
					purityy = mcpsFromTracks.data( &track - &trackList[0] ).at( nn )->cleanliness;
					mm = nn;
				}
				//std::cout << "SONO QUA 2" << std::endl;
			}
			if (pfp_track.size()==0) std::cout << ">>>>>>>>SHOWER!!!!!!!!!!!!!!" << std::endl;
			if (pfp_track.size()!=1) std::cout << "MORE THAN 1 PFP per track!" << std::endl;
			if ( pfp_track[0]->Self() != pfp_id ) continue;
			//std::cout << "FOUND PFP" << std::endl;
		
			//fpdg[mcp_index] is the particle truth matched to the hit
			//mcp_track[mm] is the particle matched to the track associated to the cluster (to which the hit is clustered)
			if ( !fis_tracked[mcp_index] && track->ID() == freco_trackid[mcp_index] )
				mf::LogError(__FUNCTION__) << "ERROR! This should never happen..." << std::endl;
			
			if ( track->ID() == freco_trackid[mcp_index] ) {
					//the hit and track matching agree
				fclustered_matched[mcp_index] = fclustered_matched[mcp_index] + 1;
		                fclustered_matched_charge[mcp_index].push_back( hit->Integral() );
			} else { //mismatch
				if ( fpdg[mcp_index] == 13 && mcp_track[mm]->PdgCode()==13) {
					std::cout << "Broken track!" << std::endl;
					//std::cout << "trackid1 " << track->ID() << " trackid2 " << freco_trackid[mcp_index] << std::endl; 
					//std::cout << "mcp_index " << mcp_index << " muon_index " << muon_pos << std::endl; 
					//std::cout << "numero track " << &track - &trackList[0] << std::endl;
				} else if (fpdg[mcp_index] == 2212 && fp0[mcp_index] < 0.2) {
					std::cout << "Low proton!!!" << std:: endl;
				} else {
				      fhit_mismatch_pdg[mcp_index].push_back ( mcp_track[mm]->PdgCode() );
				      fclustered_mismatched[mcp_index] = fclustered_mismatched[mcp_index] + 1;
		                      fclustered_mismatched_charge[mcp_index].push_back( hit->Integral() );
				}
			}
			
		} //tracklist
	       }//is clustered
	//is it clustered?
	//is it attached to a track?
	//are those the right track and the right cluster?
	//check protons and merged protons

  } // loop on hits

  //now fill histograms on hit merging
  for (unsigned jj=0; jj<fpdg.size(); jj++) {

	  if ( count_not_tracked == 0 ) { //all protons above 20MeV are tracked
	  if ( fpdg[jj]==13 ) {
		  if ( freco_mcp_collection_hits[jj] ) {
		  h_muon_clustering_prob_good_protons->Fill(0., float(fnot_clustered[jj])/float(freco_mcp_collection_hits[jj]) );	
		  h_muon_clustering_prob_good_protons->Fill(1., float(fclustered_matched[jj])/float(freco_mcp_collection_hits[jj]) );	
		  h_muon_clustering_prob_good_protons->Fill(2., float(fclustered_mismatched[jj])/float(freco_mcp_collection_hits[jj]) );	
		  }
		
		  for (unsigned ii=0; ii<fhit_mismatch_pdg[jj].size(); ii++)
		  h_muon_clustering_mismatch_pdg_good_protons->Fill(fhit_mismatch_pdg[jj][ii]);	

		  if ( fnot_clustered[jj] > 0) 
		  h_muon_not_clustered_reco_hits_good_protons->Fill( fnot_clustered[jj] , freco_mcp_collection_hits[jj] );
		  if ( fclustered_matched[jj] > 0) 
		  h_muon_clustered_matched_reco_hits_good_protons->Fill( fclustered_matched[jj] , freco_mcp_collection_hits[jj] );
		  if ( fclustered_mismatched[jj] > 0) 
		  h_muon_clustered_mismatched_reco_hits_good_protons->Fill( fclustered_mismatched[jj] , freco_mcp_collection_hits[jj] );
                  
		  for ( unsigned ii=0; ii< fnot_clustered_charge[jj].size(); ii++ ) 
		  h_muon_not_clustered_reco_charge_good_protons->Fill( fnot_clustered_charge[jj][ii], float(fnot_clustered[jj])/float(freco_mcp_collection_hits[jj]) );
		  for ( unsigned ii=0; ii< fclustered_matched_charge[jj].size(); ii++ ) 
		  h_muon_clustered_matched_reco_charge_good_protons->Fill( fclustered_matched_charge[jj][ii], float(fclustered_matched[jj])/float(freco_mcp_collection_hits[jj]) );
		  for ( unsigned ii=0; ii< fclustered_mismatched_charge[jj].size(); ii++ ) 
		  h_muon_clustered_mismatched_reco_charge_good_protons->Fill( fclustered_mismatched_charge[jj][ii], float(fclustered_mismatched[jj])/float(freco_mcp_collection_hits[jj]) );
		  
		  for ( unsigned zz=0; zz<fpdg.size(); zz++) {
			  if (fpdg[zz]!=2212) continue;
		  if (fnot_clustered[jj]) {
		  h_muon_NC_lateral_hits_good_protons->Fill( float(fnot_clustered[jj])/float(freco_mcp_collection_hits[jj]),  flength[zz] * sqrt( 1 - pow( fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_NC_costheta_hits_good_protons->Fill( float(fnot_clustered[jj])/float(freco_mcp_collection_hits[jj]),  fcostheta_muon[zz]  ); //lateral distance between muon and proton
		  }
		  if (fclustered_matched[jj]) {
		  h_muon_CMA_lateral_hits_good_protons->Fill( float(fclustered_matched[jj])/float(freco_mcp_collection_hits[jj]),  flength[zz] * sqrt( 1 - pow( fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_CMA_costheta_hits_good_protons->Fill( float(fclustered_matched[jj])/float(freco_mcp_collection_hits[jj]),  fcostheta_muon[zz]  ); //lateral distance between muon and proton
		  }
		  if (fclustered_mismatched[jj]){
		  h_muon_CMI_lateral_hits_good_protons->Fill( float(fclustered_mismatched[jj])/float(freco_mcp_collection_hits[jj]),  flength[zz] * sqrt( 1 - pow( fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_CMI_costheta_hits_good_protons->Fill( float(fclustered_mismatched[jj])/float(freco_mcp_collection_hits[jj]),  fcostheta_muon[zz]  ); //lateral distance between muon and proton
		  }

		  for ( unsigned ii=0; ii< fnot_clustered_charge[jj].size(); ii++ ) {
                  h_muon_NC_lateral_charge_good_protons->Fill( fnot_clustered_charge[jj][ii],  flength[zz] * sqrt( 1 - pow( fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_NC_costheta_charge_good_protons->Fill( fnot_clustered_charge[jj][ii], fcostheta_muon[zz] ); //lateral distance between muon and proton
		  }
		  for ( unsigned ii=0; ii< fclustered_matched_charge[jj].size(); ii++ ) {
                  h_muon_CMA_lateral_charge_good_protons->Fill( fclustered_matched_charge[jj][ii],  flength[zz] * sqrt( 1 - pow( fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_CMA_costheta_charge_good_protons->Fill( fclustered_matched_charge[jj][ii], fcostheta_muon[zz] ); //lateral distance between muon and proton
		  }
		  for ( unsigned ii=0; ii< fclustered_mismatched_charge[jj].size(); ii++ ) {
                  h_muon_CMI_lateral_charge_good_protons->Fill( fclustered_mismatched_charge[jj][ii],  flength[zz] * sqrt( 1 - pow( fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_CMI_costheta_charge_good_protons->Fill( fclustered_mismatched_charge[jj][ii], fcostheta_muon[zz] ); //lateral distance between muon and proton
		  }
		  } //fpdg
		  } //muon
	  
	  if (fpdg[jj] ==2212 && fis_tracked[jj]) {

		  if ( freco_mcp_collection_hits[jj] ) {
		  h_proton_clustering_prob_good_protons->Fill(0., float(fnot_clustered[jj])/float(freco_mcp_collection_hits[jj]) );	
		  h_proton_clustering_prob_good_protons->Fill(1., float(fclustered_matched[jj])/float(freco_mcp_collection_hits[jj]) );	
		  h_proton_clustering_prob_good_protons->Fill(2., float(fclustered_mismatched[jj])/float(freco_mcp_collection_hits[jj]) );	
		  }

		  for (unsigned ii=0; ii<fhit_mismatch_pdg[jj].size(); ii++)
		  h_proton_clustering_mismatch_pdg_good_protons->Fill(fhit_mismatch_pdg[jj][ii]);	

		  if (fnot_clustered[jj]) {
		  h_proton_not_clustered_reco_hits_good_protons->Fill( fnot_clustered[jj] , freco_mcp_collection_hits[jj] );
                  h_proton_NC_lateral_hits_good_protons->Fill( float(fnot_clustered[jj])/float(freco_mcp_collection_hits[jj]),  flength[jj] * sqrt( 1 - pow( fcostheta_muon[jj] ,2 ) ) ); //lateral distance between proton and proton
                  h_proton_NC_costheta_hits_good_protons->Fill( float(fnot_clustered[jj])/float(freco_mcp_collection_hits[jj]),  fcostheta_muon[jj]  ); //lateral distance between muon and proton
		  }
		  if (fclustered_matched[jj]) {
		  h_proton_clustered_matched_reco_hits_good_protons->Fill( fclustered_matched[jj] , freco_mcp_collection_hits[jj] );
                  h_proton_CMA_lateral_hits_good_protons->Fill( float(fclustered_matched[jj])/float(freco_mcp_collection_hits[jj]),  flength[jj] * sqrt( 1 - pow( fcostheta_muon[jj] ,2 ) ) ); //lateral distance between proton and proton
                  h_proton_CMA_costheta_hits_good_protons->Fill( float(fclustered_matched[jj])/float(freco_mcp_collection_hits[jj]),  fcostheta_muon[jj]  ); //lateral distance between muon and proton
		  }
		  if (fclustered_mismatched[jj]) {
		  h_proton_clustered_mismatched_reco_hits_good_protons->Fill( fclustered_mismatched[jj] , freco_mcp_collection_hits[jj] );
                  h_proton_CMI_lateral_hits_good_protons->Fill( float(fclustered_mismatched[jj])/float(freco_mcp_collection_hits[jj]),  flength[jj] * sqrt( 1 - pow( fcostheta_muon[jj] ,2 ) ) ); //lateral distance between proton and proton
                  h_proton_CMI_costheta_hits_good_protons->Fill( float(fclustered_mismatched[jj])/float(freco_mcp_collection_hits[jj]),  fcostheta_muon[jj]  ); //lateral distance between muon and proton
		  }
		 
		  for ( unsigned ii=0; ii< fnot_clustered_charge[jj].size(); ii++ ) {
		  h_proton_not_clustered_reco_charge_good_protons->Fill( fnot_clustered_charge[jj][ii], float(fnot_clustered[jj])/float(freco_mcp_collection_hits[jj]) );
                  h_proton_NC_lateral_charge_good_protons->Fill( fnot_clustered_charge[jj][ii],  flength[jj] * sqrt( 1 - pow( fcostheta_muon[jj] ,2 ) ) ); //lateral distance between muon and proton
                  h_proton_NC_costheta_charge_good_protons->Fill( fnot_clustered_charge[jj][ii], fcostheta_muon[jj] ); //lateral distance between muon and proton
		  }
		  for ( unsigned ii=0; ii< fclustered_matched_charge[jj].size(); ii++ ) {
		  h_proton_clustered_matched_reco_charge_good_protons->Fill( fclustered_matched_charge[jj][ii], float(fclustered_matched[jj])/float(freco_mcp_collection_hits[jj]) );
                  h_proton_CMA_lateral_charge_good_protons->Fill( fclustered_matched_charge[jj][ii],  flength[jj] * sqrt( 1 - pow( fcostheta_muon[jj] ,2 ) ) ); //lateral distance between muon and proton
                  h_proton_CMA_costheta_charge_good_protons->Fill( fclustered_matched_charge[jj][ii], fcostheta_muon[jj] ); //lateral distance between muon and proton
		  }
		  for ( unsigned ii=0; ii< fclustered_mismatched_charge[jj].size(); ii++ ) {
		  h_proton_clustered_mismatched_reco_charge_good_protons->Fill( fclustered_mismatched_charge[jj][ii], float(fclustered_mismatched[jj])/float(freco_mcp_collection_hits[jj]) );
                  h_proton_CMI_lateral_charge_good_protons->Fill( fclustered_mismatched_charge[jj][ii],  flength[jj] * sqrt( 1 - pow( fcostheta_muon[jj] ,2 ) ) ); //lateral distance between muon and proton
                  h_proton_CMI_costheta_charge_good_protons->Fill( fclustered_mismatched_charge[jj][ii], fcostheta_muon[jj] ); //lateral distance between muon and proton
		  }
	     }//proton

	  } //if all are tracked
	  
if ( count_not_tracked > 0 ) { //some are not tracked
	  
	if ( fpdg[jj]==13 ) { //the muon is tracked by definition
		  if ( freco_mcp_collection_hits[jj] ) {
		  h_muon_clustering_prob_bad_protons->Fill(0., float(fnot_clustered[jj])/float(freco_mcp_collection_hits[jj]) );	
		  h_muon_clustering_prob_bad_protons->Fill(1., float(fclustered_matched[jj])/float(freco_mcp_collection_hits[jj]) );	
		  h_muon_clustering_prob_bad_protons->Fill(2., float(fclustered_mismatched[jj])/float(freco_mcp_collection_hits[jj]) );	
		  }

		  for (unsigned ii=0; ii<fhit_mismatch_pdg[jj].size(); ii++)
		  h_muon_clustering_mismatch_pdg_bad_protons->Fill(fhit_mismatch_pdg[jj][ii]);	
		 
		  if (fnot_clustered[jj])
		  h_muon_not_clustered_reco_hits_bad_protons->Fill( fnot_clustered[jj] , freco_mcp_collection_hits[jj] );
		  if (fclustered_matched[jj])
		  h_muon_clustered_matched_reco_hits_bad_protons->Fill( fclustered_matched[jj] , freco_mcp_collection_hits[jj] );
		  if (fclustered_mismatched[jj])
		  h_muon_clustered_mismatched_reco_hits_bad_protons->Fill( fclustered_mismatched[jj] , freco_mcp_collection_hits[jj] );
		  
		  for ( unsigned ii=0; ii< fnot_clustered_charge[jj].size(); ii++ ) 
		  h_muon_not_clustered_reco_charge_bad_protons->Fill( fnot_clustered_charge[jj][ii], float(fnot_clustered[jj])/float(freco_mcp_collection_hits[jj]) );
		  for ( unsigned ii=0; ii< fclustered_matched_charge[jj].size(); ii++ ) 
		  h_muon_clustered_matched_reco_charge_bad_protons->Fill( fclustered_matched_charge[jj][ii], float(fclustered_matched[jj])/float(freco_mcp_collection_hits[jj]) );
		  for ( unsigned ii=0; ii< fclustered_mismatched_charge[jj].size(); ii++ ) 
		  h_muon_clustered_mismatched_reco_charge_bad_protons->Fill( fclustered_mismatched_charge[jj][ii], float(fclustered_mismatched[jj])/float(freco_mcp_collection_hits[jj]) );

		  for ( unsigned zz=0; zz<fpdg.size(); zz++) {
			  if (fpdg[zz]!=2212) continue;
                  if (fnot_clustered[jj]) {
		  h_muon_NC_lateral_hits_bad_protons->Fill( float(fnot_clustered[jj])/float(freco_mcp_collection_hits[jj]),  flength[zz] * sqrt( 1 - pow( fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_NC_costheta_hits_bad_protons->Fill( float(fnot_clustered[jj])/float(freco_mcp_collection_hits[jj]),  fcostheta_muon[zz]  ); //lateral distance between muon and proton
		  }
		  if (fclustered_matched[jj]) {
                  h_muon_CMA_lateral_hits_bad_protons->Fill( float(fclustered_matched[jj])/float(freco_mcp_collection_hits[jj]),  flength[zz] * sqrt( 1 - pow( fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_CMA_costheta_hits_bad_protons->Fill( float(fclustered_matched[jj])/float(freco_mcp_collection_hits[jj]),  fcostheta_muon[zz]  ); //lateral distance between muon and proton
		  }
		  if (fclustered_mismatched[jj]){
                  h_muon_CMI_lateral_hits_bad_protons->Fill( float(fclustered_mismatched[jj])/float(freco_mcp_collection_hits[jj]),  flength[zz] * sqrt( 1 - pow( fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_CMI_costheta_hits_bad_protons->Fill( float(fclustered_mismatched[jj])/float(freco_mcp_collection_hits[jj]),  fcostheta_muon[zz]  ); //lateral distance between muon and proton
		  }
		 
		  for ( unsigned ii=0; ii< fnot_clustered_charge[jj].size(); ii++ ) {
                  h_muon_NC_lateral_charge_bad_protons->Fill( fnot_clustered_charge[jj][ii],  flength[zz] * sqrt( 1 - pow( fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_NC_costheta_charge_bad_protons->Fill( fnot_clustered_charge[jj][ii], fcostheta_muon[zz] ); //lateral distance between muon and proton
		  }
		  for ( unsigned ii=0; ii< fclustered_matched_charge[jj].size(); ii++ ) {
                  h_muon_CMA_lateral_charge_bad_protons->Fill( fclustered_matched_charge[jj][ii],  flength[zz] * sqrt( 1 - pow( fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_CMA_costheta_charge_bad_protons->Fill( fclustered_matched_charge[jj][ii], fcostheta_muon[zz] ); //lateral distance between muon and proton
		  }
		  for ( unsigned ii=0; ii< fclustered_mismatched_charge[jj].size(); ii++ ) {
                  h_muon_CMI_lateral_charge_bad_protons->Fill( fclustered_mismatched_charge[jj][ii],  flength[zz] * sqrt( 1 - pow( fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_CMI_costheta_charge_bad_protons->Fill( fclustered_mismatched_charge[jj][ii], fcostheta_muon[zz] ); //lateral distance between muon and proton
		  }
		  } //fpdg
		} //muon
	  
	if (fpdg[jj] ==2212 && fis_tracked[jj]) {
		  if ( freco_mcp_collection_hits[jj] ) {
		  h_proton_clustering_prob_bad_protons->Fill(0., float(fnot_clustered[jj])/float(freco_mcp_collection_hits[jj]) );	
		  h_proton_clustering_prob_bad_protons->Fill(1., float(fclustered_matched[jj])/float(freco_mcp_collection_hits[jj]) );	
		  h_proton_clustering_prob_bad_protons->Fill(2., float(fclustered_mismatched[jj])/float(freco_mcp_collection_hits[jj]) );	
		  } 
		  for (unsigned ii=0; ii<fhit_mismatch_pdg[jj].size(); ii++) 
		  h_proton_clustering_mismatch_pdg_bad_protons->Fill(fhit_mismatch_pdg[jj][ii]);	

		  if (fnot_clustered[jj]) {
		  h_proton_not_clustered_reco_hits_bad_protons->Fill( fnot_clustered[jj] , freco_mcp_collection_hits[jj] );
                  h_proton_NC_lateral_hits_bad_protons->Fill( float(fnot_clustered[jj])/float(freco_mcp_collection_hits[jj]),  flength[jj] * sqrt( 1 - pow( fcostheta_muon[jj] ,2 ) ) ); //lateral distance between proton and proton
                  h_proton_NC_costheta_hits_bad_protons->Fill( float(fnot_clustered[jj])/float(freco_mcp_collection_hits[jj]),  fcostheta_muon[jj]  ); //lateral distance between muon and proton
		  }
		  if (fclustered_matched[jj]){
		  h_proton_clustered_matched_reco_hits_bad_protons->Fill( fclustered_matched[jj] , freco_mcp_collection_hits[jj] );
                  h_proton_CMA_lateral_hits_bad_protons->Fill( float(fclustered_matched[jj])/float(freco_mcp_collection_hits[jj]),  flength[jj] * sqrt( 1 - pow( fcostheta_muon[jj] ,2 ) ) ); //lateral distance between proton and proton
                  h_proton_CMA_costheta_hits_bad_protons->Fill( float(fclustered_matched[jj])/float(freco_mcp_collection_hits[jj]),  fcostheta_muon[jj]  ); //lateral distance between muon and proton
		  }
		  if (fclustered_mismatched[jj]){
		  h_proton_clustered_mismatched_reco_hits_bad_protons->Fill( fclustered_mismatched[jj] , freco_mcp_collection_hits[jj] );
                  h_proton_CMI_lateral_hits_bad_protons->Fill( float(fclustered_mismatched[jj])/float(freco_mcp_collection_hits[jj]),  flength[jj] * sqrt( 1 - pow( fcostheta_muon[jj] ,2 ) ) ); //lateral distance between proton and proton
                  h_proton_CMI_costheta_hits_bad_protons->Fill( float(fclustered_mismatched[jj])/float(freco_mcp_collection_hits[jj]),  fcostheta_muon[jj]  ); //lateral distance between muon and proton
		  }

		  for ( unsigned ii=0; ii< fnot_clustered_charge[jj].size(); ii++ ) {
		  h_proton_not_clustered_reco_charge_bad_protons->Fill( fnot_clustered_charge[jj][ii], float(fnot_clustered[jj])/float(freco_mcp_collection_hits[jj]) );
                  h_proton_NC_lateral_charge_bad_protons->Fill( fnot_clustered_charge[jj][ii],  flength[jj] * sqrt( 1 - pow( fcostheta_muon[jj] ,2 ) ) ); //lateral distance between muon and proton
                  h_proton_NC_costheta_charge_bad_protons->Fill( fnot_clustered_charge[jj][ii], fcostheta_muon[jj] ); //lateral distance between muon and proton
		  }
		  for ( unsigned ii=0; ii< fclustered_matched_charge[jj].size(); ii++ ) {
		  h_proton_clustered_matched_reco_charge_bad_protons->Fill( fclustered_matched_charge[jj][ii], float(fclustered_matched[jj])/float(freco_mcp_collection_hits[jj]) );
                  h_proton_CMA_lateral_charge_bad_protons->Fill( fclustered_matched_charge[jj][ii],  flength[jj] * sqrt( 1 - pow( fcostheta_muon[jj] ,2 ) ) ); //lateral distance between muon and proton
                  h_proton_CMA_costheta_charge_bad_protons->Fill( fclustered_matched_charge[jj][ii], fcostheta_muon[jj] ); //lateral distance between muon and proton
		  }
		  for ( unsigned ii=0; ii< fclustered_mismatched_charge[jj].size(); ii++ ) {
		  h_proton_clustered_mismatched_reco_charge_bad_protons->Fill( fclustered_mismatched_charge[jj][ii], float(fclustered_mismatched[jj])/float(freco_mcp_collection_hits[jj]) );
                  h_proton_CMI_lateral_charge_bad_protons->Fill( fclustered_mismatched_charge[jj][ii],  flength[jj] * sqrt( 1 - pow( fcostheta_muon[jj] ,2 ) ) ); //lateral distance between muon and proton
                  h_proton_CMI_costheta_charge_bad_protons->Fill( fclustered_mismatched_charge[jj][ii], fcostheta_muon[jj] ); //lateral distance between muon and proton
		  }
	  } //protons && is tracked 
	
	if ( fpdg[jj] == 2212 && !fis_tracked[jj] ) {
		  if ( freco_mcp_collection_hits[jj] ) {
		  h_proton_clustering_prob_bad_protons_not_tracked->Fill(0., float(fnot_clustered[jj])/float(freco_mcp_collection_hits[jj]) );	
		  h_proton_clustering_prob_bad_protons_not_tracked->Fill(1., float(fclustered_matched[jj])/float(freco_mcp_collection_hits[jj]) );	
		  h_proton_clustering_prob_bad_protons_not_tracked->Fill(2., float(fclustered_mismatched[jj])/float(freco_mcp_collection_hits[jj]) );	
		  }
		  for (unsigned ii=0; ii<fhit_mismatch_pdg[jj].size(); ii++)
		  h_proton_clustering_mismatch_pdg_bad_protons_not_tracked->Fill(fhit_mismatch_pdg[jj][ii]);	
		 
		  if (fnot_clustered[jj]){
		  h_proton_not_clustered_reco_hits_bad_protons_not_tracked->Fill( fnot_clustered[jj] , freco_mcp_collection_hits[jj] );
                  h_proton_NC_lateral_hits_bad_protons_not_tracked->Fill( float(fnot_clustered[jj])/float(freco_mcp_collection_hits[jj]),  flength[jj] * sqrt( 1 - pow( fcostheta_muon[jj] ,2 ) ) ); //lateral distance between proton and proton
                  h_proton_NC_costheta_hits_bad_protons_not_tracked->Fill( float(fnot_clustered[jj])/float(freco_mcp_collection_hits[jj]),  fcostheta_muon[jj]  ); //lateral distance between muon and proton
		  }
		  if (fclustered_matched[jj]){
		  h_proton_clustered_matched_reco_hits_bad_protons_not_tracked->Fill( fclustered_matched[jj] , freco_mcp_collection_hits[jj] );
                  h_proton_CMA_lateral_hits_bad_protons_not_tracked->Fill( float(fclustered_matched[jj])/float(freco_mcp_collection_hits[jj]),  flength[jj] * sqrt( 1 - pow( fcostheta_muon[jj] ,2 ) ) ); //lateral distance between proton and proton
                  h_proton_CMA_costheta_hits_bad_protons_not_tracked->Fill( float(fclustered_matched[jj])/float(freco_mcp_collection_hits[jj]),  fcostheta_muon[jj]  ); //lateral distance between muon and proton
		  }
		  if (fclustered_mismatched[jj]){
		  h_proton_clustered_mismatched_reco_hits_bad_protons_not_tracked->Fill( fclustered_mismatched[jj] , freco_mcp_collection_hits[jj] );
                  h_proton_CMI_lateral_hits_bad_protons_not_tracked->Fill( float(fclustered_mismatched[jj])/float(freco_mcp_collection_hits[jj]),  flength[jj] * sqrt( 1 - pow( fcostheta_muon[jj] ,2 ) ) ); //lateral distance between proton and proton
                  h_proton_CMI_costheta_hits_bad_protons_not_tracked->Fill( float(fclustered_mismatched[jj])/float(freco_mcp_collection_hits[jj]),  fcostheta_muon[jj]  ); //lateral distance between muon and proton
		  }

		  for ( unsigned ii=0; ii< fnot_clustered_charge[jj].size(); ii++ ) {
		  h_proton_not_clustered_reco_charge_bad_protons_not_tracked->Fill( fnot_clustered_charge[jj][ii], float(fnot_clustered[jj])/float(freco_mcp_collection_hits[jj]) );
                  h_proton_NC_lateral_charge_bad_protons_not_tracked->Fill( fnot_clustered_charge[jj][ii],  flength[jj] * sqrt( 1 - pow( fcostheta_muon[jj] ,2 ) ) ); //lateral distance between muon and proton
                  h_proton_NC_costheta_charge_bad_protons_not_tracked->Fill( fnot_clustered_charge[jj][ii], fcostheta_muon[jj] ); //lateral distance between muon and proton
		  }
		  for ( unsigned ii=0; ii< fclustered_matched_charge[jj].size(); ii++ ) {
		  h_proton_clustered_matched_reco_charge_bad_protons_not_tracked->Fill( fclustered_matched_charge[jj][ii], float(fclustered_matched[jj])/float(freco_mcp_collection_hits[jj]) );
                  h_proton_CMA_lateral_charge_bad_protons_not_tracked->Fill( fclustered_matched_charge[jj][ii],  flength[jj] * sqrt( 1 - pow( fcostheta_muon[jj] ,2 ) ) ); //lateral distance between muon and proton
                  h_proton_CMA_costheta_charge_bad_protons_not_tracked->Fill( fclustered_matched_charge[jj][ii], fcostheta_muon[jj] ); //lateral distance between muon and proton
		  }
		  for ( unsigned ii=0; ii< fclustered_mismatched_charge[jj].size(); ii++ ) {
		  h_proton_clustered_mismatched_reco_charge_bad_protons_not_tracked->Fill( fclustered_mismatched_charge[jj][ii], float(fclustered_mismatched[jj])/float(freco_mcp_collection_hits[jj]) );
                  h_proton_CMI_lateral_charge_bad_protons_not_tracked->Fill( fclustered_mismatched_charge[jj][ii],  flength[jj] * sqrt( 1 - pow( fcostheta_muon[jj] ,2 ) ) ); //lateral distance between muon and proton
                  h_proton_CMI_costheta_charge_bad_protons_not_tracked->Fill( fclustered_mismatched_charge[jj][ii], fcostheta_muon[jj] ); //lateral distance between muon and proton
		  }
		} //proton && ! is tracked

}  //some are not tracked
	  
  if ( is_lowmomentum_p ) { //only protons below 20MeV, none is tracked
	if ( fpdg[jj]==13 ) { //the muon is tracked by definition
		  if ( freco_mcp_collection_hits[jj] ) {
		  h_muon_clustering_prob_low_protons->Fill(0., float(fnot_clustered[jj])/float(freco_mcp_collection_hits[jj]) );	
		  h_muon_clustering_prob_low_protons->Fill(1., float(fclustered_matched[jj])/float(freco_mcp_collection_hits[jj]) );	
		  h_muon_clustering_prob_low_protons->Fill(2., float(fclustered_mismatched[jj])/float(freco_mcp_collection_hits[jj]) );	
		  }
		  for (unsigned ii=0; ii<fhit_mismatch_pdg[jj].size(); ii++)
		  h_muon_clustering_mismatch_pdg_low_protons->Fill(fhit_mismatch_pdg[jj][ii]);	
		  
		  for ( unsigned ii=0; ii< fnot_clustered_charge[jj].size(); ii++ )
		  h_muon_not_clustered_reco_charge_low_protons->Fill( fnot_clustered_charge[jj][ii], float(fnot_clustered[jj])/float(freco_mcp_collection_hits[jj]) );
		  for ( unsigned ii=0; ii< fclustered_matched_charge[jj].size(); ii++ )
		  h_muon_clustered_matched_reco_charge_low_protons->Fill( fclustered_matched_charge[jj][ii], float(fclustered_matched[jj])/float(freco_mcp_collection_hits[jj]) );
		  for ( unsigned ii=0; ii< fclustered_mismatched_charge[jj].size(); ii++ )
		  h_muon_clustered_mismatched_reco_charge_low_protons->Fill( fclustered_mismatched_charge[jj][ii], float(fclustered_mismatched[jj])/float(freco_mcp_collection_hits[jj]) );
		  
		  for ( unsigned zz=0; zz<fpdg.size(); zz++) {
			  if (fpdg[zz]!=2212) continue;
		  if (fnot_clustered[jj]) {
		  h_muon_not_clustered_reco_hits_low_protons->Fill( fnot_clustered[jj] , freco_mcp_collection_hits[jj] );
                  h_muon_NC_lateral_hits_low_protons->Fill( float(fnot_clustered[jj])/float(freco_mcp_collection_hits[jj]),  flength[zz] * sqrt( 1 - pow( fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_NC_costheta_hits_low_protons->Fill( float(fnot_clustered[jj])/float(freco_mcp_collection_hits[jj]),  fcostheta_muon[zz]  ); //lateral distance between muon and proton
		  } 
		  if (fclustered_matched[jj]) {
		  h_muon_clustered_matched_reco_hits_low_protons->Fill( fclustered_matched[jj] , freco_mcp_collection_hits[jj] );
                  h_muon_CMA_lateral_hits_low_protons->Fill( float(fclustered_matched[jj])/float(freco_mcp_collection_hits[jj]),  flength[zz] * sqrt( 1 - pow( fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_CMA_costheta_hits_low_protons->Fill( float(fclustered_matched[jj])/float(freco_mcp_collection_hits[jj]),  fcostheta_muon[zz]  ); //lateral distance between muon and proton
		  }
		  if (fclustered_mismatched[jj]) {
		  h_muon_clustered_mismatched_reco_hits_low_protons->Fill( fclustered_mismatched[jj] , freco_mcp_collection_hits[jj] );
                  h_muon_CMI_lateral_hits_low_protons->Fill( float(fclustered_mismatched[jj])/float(freco_mcp_collection_hits[jj]),  flength[zz] * sqrt( 1 - pow( fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_CMI_costheta_hits_low_protons->Fill( float(fclustered_mismatched[jj])/float(freco_mcp_collection_hits[jj]),  fcostheta_muon[zz]  ); //lateral distance between muon and proton
		  }

		  for ( unsigned ii=0; ii< fnot_clustered_charge[jj].size(); ii++ ) {
                  h_muon_NC_lateral_charge_low_protons->Fill( fnot_clustered_charge[jj][ii],  flength[zz] * sqrt( 1 - pow( fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_NC_costheta_charge_low_protons->Fill( fnot_clustered_charge[jj][ii], fcostheta_muon[zz] ); //lateral distance between muon and proton
		  }
		  for ( unsigned ii=0; ii< fclustered_matched_charge[jj].size(); ii++ ) {
                  h_muon_CMA_lateral_charge_low_protons->Fill( fclustered_matched_charge[jj][ii],  flength[zz] * sqrt( 1 - pow( fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_CMA_costheta_charge_low_protons->Fill( fclustered_matched_charge[jj][ii], fcostheta_muon[zz] ); //lateral distance between muon and proton
		  }
		  for ( unsigned ii=0; ii< fclustered_mismatched_charge[jj].size(); ii++ ) {
                  h_muon_CMI_lateral_charge_low_protons->Fill( fclustered_mismatched_charge[jj][ii],  flength[zz] * sqrt( 1 - pow( fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_CMI_costheta_charge_low_protons->Fill( fclustered_mismatched_charge[jj][ii], fcostheta_muon[zz] ); //lateral distance between muon and proton
		  }
		    } //fpdg
		  } //muon
	  	} //low p

  } //loop on MCP for hit merging histograms
	
/*
  //MC total amount of hits (or something like that)
	  float tot_electrons = 0;
  for ( auto const& simch: simChannelList ) {
	  auto map = simch->TDCIDEMap();
	  for (auto const& it: map) {
		for (auto const& it2: it.second) 
		  tot_electrons+=it2.numElectrons;	
	  }
  }


    int all_hits_mcp = 0;
    for (auto const& mccp: mcList) {
    if ( !(std::abs(mccp->PdgCode()) == 11 && mccp->Mother() == muon_id && abs(mccp->Position().X()-muon_endx) < DBL_EPSILON &&
		    abs(mccp->Position().Y()-muon_endy) < DBL_EPSILON && abs(mccp->Position().Z()-muon_endz) < DBL_EPSILON ) )
    if ( std::abs(mccp->PdgCode()) == 14 || mccp->Process()!="primary" || mccp->StatusCode()!=1 || mccp->Mother()>0 ) continue; //we only want primaries
    all_hits_mcp += hitsFromMCP.at( &mccp - &mcList[0] ).size();
    }

    for (auto const& mccp: mcList) {
    if ( !(std::abs(mccp->PdgCode()) == 11 && mccp->Mother() == muon_id && abs(mccp->Position().X()-muon_endx) < DBL_EPSILON &&
		    abs(mccp->Position().Y()-muon_endy) < DBL_EPSILON && abs(mccp->Position().Z()-muon_endz) < DBL_EPSILON ) )
    if ( std::abs(mccp->PdgCode()) == 14 || mccp->Process()!="primary" || mccp->StatusCode()!=1 || mccp->Mother()>0 ) continue; //we only want primaries
    	float hits_this_mcp = 0;
	std::vector< art::Ptr < recob::Hit> > hits_mcp = hitsFromMCP.at( &mccp - &mcList[0] );
	for (auto const& hit: hits_mcp) {
	//bool max = hitsFromMCP.data( &mccp - &mcList[0] ).at( &hit - &hits_mcp[0] )->isMaxIDE;
	//bool max2 = hitsFromMCP.data( &mccp - &mcList[0] ).at( &hit - &hits_mcp[0] )->isMaxIDEN;
	//hits_this_mcp += hitsFromMCP.data( &mccp - &mcList[0] ).at( &hit - &hits_mcp[0] )->ideNFraction * hit->Integral() ;
	hits_this_mcp += hitsFromMCP.data( &mccp - &mcList[0] ).at( &hit - &hits_mcp[0] )->numElectrons ;
	cout << "hits " << hitsFromMCP.data( &mccp - &mcList[0] ).at( &hit - &hits_mcp[0] )->energy << endl; 
	cout << "size " << hitsFromMCP.data( &mccp - &mcList[0] ).size() << endl; 
	cout << "ide " << hitsFromMCP.data( &mccp - &mcList[0] ).at( &hit - &hits_mcp[0] )->ideFraction << endl; 
	cout << "electrons " << hitsFromMCP.data( &mccp - &mcList[0] ).at( &hit - &hits_mcp[0] )->numElectrons << endl; 
	//if ( hit->View() == geo::kCollection ) hits_this_mcp += hit->Integral() ;
	}
	float tmp = 0;
    for (auto const& aa: mcHitList) {
	    	for ( unsigned i=0; i<aa->size(); i++) {
	        if ( aa->at(i).PartTrackId() == mccp->TrackId())
			tmp+=aa->at(i).Charge();
		}
    }
    	//std::cout << "This MCP pdg=" << mccp->PdgCode() << " reco-matched hits: " << hits_this_mcp << " mctruth hits: " << tmp << std:: endl;

    tmp=0;
    for (auto const& simchannel: simChannelList) {
	    auto ides = simchannel->TrackIDEs(0,100000);
	    for (auto const& ide: ides) {
		if (ide.trackID == mccp->TrackId()) tmp+=ide.numElectrons;
	    }
    	}
    	std::cout << "This MCP pdg=" << mccp->PdgCode() << " reco-matched hits: " << hits_this_mcp/all_hits_mcp << " mctruth hits: " << tmp/tot_electrons << std:: endl;
    }

    for (auto const& aa: mcHitList) {
	    std::cout << "ch: " << aa->Channel() << " MCHitCollection size " << aa->size() << std::endl;
	    	for ( unsigned i=0; i<aa->size(); i++) {
	        std::cout << "charge " << aa->at(i).Charge() << " partID " << aa->at(i).PartTrackId() << std::endl;
		}
*/	    
    //for (auto const& bb: aa) {
	//    std::cout << "charge " << bb->Charge() << " partID " << bb->PartTrackId() << std::endl;
   // }
   // }
  //std::cout << ">>>><<<< TOTAL HITS: " << hitList.size() << " HITS FROM MCPS: " << all_hits_mcp << " TOTAL HITS FROM TRACKS: " << all_hits_tracks << " TOTAL MC HITS: " << mcHitList.size() << std::endl;

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
	
  } // good muon

  //make hit plots
		        
  //for (unsigned ii=0; ii<fnot_clustered_tracked.size(); ii++) {
//	  if ( fpdg[ii] == 13 ) h_hit_fraction_notclustered_tracked_muon->Fill( fnot_clustered_tracked[ii] );
//	  if ( fpdg[ii] == 2212 ) h_hit_fraction_notclustered_tracked_proton->Fill( fnot_clustered_tracked[ii] );
  //}

  recoTree->Fill();
  } //MCtruth
  

}

void recohelper::RecoBenchmarker::AllocateAnalysisHistograms() {
   hproton_multi_all = tfs->make<TH1D>("proton_multi_all","Proton multiplicity for all CC events;# protons;# of events normalized",20,-0.5,19.5); //proton multiplicity
   hproton_leading_kinE = tfs->make<TH1D>("proton_leading_kinE","Leading proton kinE;Kinetic Energy (MeV);",4000,0,4000); //leading proton kinE
   hproton_multi_above20MeV = tfs->make<TH1D>("proton_multi_above20MeV","Proton multiplicity above 20MeV;# protons;# of events normalized",20,-0.5,19.5); //proton multiplicity
   hproton_multi_below20MeV = tfs->make<TH1D>("proton_multi_below20MeV","Proton multiplicity below 20MeV;# protons;# of events normalized",20,-0.5,19.5); //proton multiplicity
   hproton_merged_not_merged = tfs->make<TH1D>("proton_merged_not_merged","Proton merged vs not merged (above 20MeV); true/false ;# of protons normalized",2,-0.5,1.5); //proton merged/not merged
   
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
   hproton_p = tfs->make<TH1D>("proton_p","Proton reco efficiency; Momentum (GeV/c)",1000,0,10); //reco efficiency protons vs p
   hproton_p_all = tfs->make<TH1D>("proton_p_all","Proton reco efficiency; Momentum (GeV/c)",1000,0,10); //reco efficiency protons vs p
   hproton_l = tfs->make<TH1D>("proton_l","Proton reco efficiency; True length (cm)",1000,0,200); //reco efficiency protons vs p
   hproton_l_all = tfs->make<TH1D>("proton_l_all","Proton reco efficiency; True length (cm)",1000,0,200); //reco efficiency protons vs p
   
   hproton_kinE_tracked_angle1 = tfs->make<TH1D>("proton_kinE_tracked_angle1","Proton reco efficiency when #mu-p angle is <30degree; Kinetic Energy (GeV)",1000,0,2); //reco efficiency protons vs kin E
   hproton_kinE_all_angle1 = tfs->make<TH1D>("proton_kinE_all_angle1","Proton reco efficiency when #mu-p angle is <30degree; Kinetic Energy (GeV)",1000,0,2); //reco efficiency protons vs kin E
   hproton_kinE_tracked_angle2 = tfs->make<TH1D>("proton_kinE_tracked_angle2","Proton reco efficiency when #mu-p angle is <60 && > 30degree; Kinetic Energy (GeV)",1000,0,2); //reco efficiency protons vs kin E
   hproton_kinE_all_angle2 = tfs->make<TH1D>("proton_kinE_all_angle2","Proton reco efficiency when #mu-p angle is <60 && > 30degree; Kinetic Energy (GeV)",1000,0,2); //reco efficiency protons vs kin E
   hproton_kinE_tracked_angle3 = tfs->make<TH1D>("proton_kinE_tracked_angle3","Proton reco efficiency when #mu-p angle is >60degree; Kinetic Energy (GeV)",1000,0,2); //reco efficiency protons vs kin E
   hproton_kinE_all_angle3 = tfs->make<TH1D>("proton_kinE_all_angle3","Proton reco efficiency when #mu-p angle is >60degree; Kinetic Energy (GeV)",1000,0,2); //reco efficiency protons vs kin E
   hproton_l_tracked_angle1 = tfs->make<TH1D>("proton_l_tracked_angle1","Proton reco efficiency when #mu-p angle is <30degree; Kinetic Energy (GeV)",1000,0,200); //reco efficiency protons vs kin E
   hproton_l_all_angle1 = tfs->make<TH1D>("proton_l_all_angle1","Proton reco efficiency when #mu-p angle is <30degree; True Length (cm)",1000,0,200); //reco efficiency protons vs kin E
   hproton_l_tracked_angle2 = tfs->make<TH1D>("proton_l_tracked_angle2","Proton reco efficiency when #mu-p angle is <60 && > 30degree; True Length (cm)",1000,0,200); //reco efficiency protons vs kin E
   hproton_l_all_angle2 = tfs->make<TH1D>("proton_l_all_angle2","Proton reco efficiency when #mu-p angle is <60 && > 30degree; True Length (cm)",1000,0,200); //reco efficiency protons vs kin E
   hproton_l_tracked_angle3 = tfs->make<TH1D>("proton_l_tracked_angle3","Proton reco efficiency when #mu-p angle is >60degree; True Length (cm)",1000,0,200); //reco efficiency protons vs kin E
   hproton_l_all_angle3 = tfs->make<TH1D>("proton_l_all_angle3","Proton reco efficiency when #mu-p angle is >60degree; True Length (cm)",1000,0,200); //reco efficiency protons vs kin E
   hproton_nhits_tracked_angle1 = tfs->make<TH1D>("proton_nhits_tracked_angle1","Proton reco efficiency when #mu-p angle is <30degree; nhits",1000,0,1000); //reco efficiency protons vs kin E
   hproton_nhits_all_angle1 = tfs->make<TH1D>("proton_nhits_all_angle1","Proton reco efficiency when #mu-p angle is <30degree; nhits",1000,0,1000); //reco efficiency protons vs kin E
   hproton_nhits_tracked_angle2 = tfs->make<TH1D>("proton_nhits_tracked_angle2","Proton reco efficiency when #mu-p angle is <60 && > 30degree; nhits",1000,0,1000); //reco efficiency protons vs kin E
   hproton_nhits_all_angle2 = tfs->make<TH1D>("proton_nhits_all_angle2","Proton reco efficiency when #mu-p angle is <60 && > 30degree; nhits",1000,0,1000); //reco efficiency protons vs kin E
   hproton_nhits_tracked_angle3 = tfs->make<TH1D>("proton_nhits_tracked_angle3","Proton reco efficiency when #mu-p angle is >60degree; nhits",1000,0,1000); //reco efficiency protons vs kin E
   hproton_nhits_all_angle3 = tfs->make<TH1D>("proton_nhits_all_angle3","Proton reco efficiency when #mu-p angle is >60degree; nhits",1000,0,1000); //reco efficiency protons vs kin E
  
   hproton_nhits = tfs->make<TH1D>("proton_nhits","nhits for reco'ed protons; nhits",1000,0,1000); //reco efficiency protons vs kin E
   hproton_nhits_CP = tfs->make<TH1D>("proton_nhits_CP","Collection Plane nhits for reco'ed protons; nhits",1000,0,1000); //reco efficiency protons vs kin E
   hproton_nhits_theta_mu = tfs->make<TH2D>("proton_nhits_theta_mu","Proton nhits vs angle with muon; nhits; #theta",1000,0,1000,100,0, 3.1415);
   hproton_nhits_CP_theta_mu = tfs->make<TH2D>("proton_nhits_CP_theta_mu","Proton collection plane nhits vs angle with muon; nhits; #theta",1000,0,1000,100,0, 3.1415);
   
   h_pmu_end_not_tracked = tfs->make<TH1D>("pmu_end_not_tracked","Not Tracked protons;Distance (cm);",1000,0,100); //lateral distance between proton end and muon
   h_pmu_end_tracked = tfs->make<TH1D>("pmu_end_tracked","Tracked protons;Distance (cm);",1000,0,100); //lateral distance between proton end and muon
   h_theta_mu_tracked = tfs->make<TH1D>("theta_mu_tracked","Cos #theta between muon and tracked protons;cos #theta",1000,-1,1); //costheta between muon and tracked protons
   hproton_theta_mu = tfs->make<TH1D>("proton_theta_mu","Cos #theta between muon and all protons >20MeV;cos #theta",1000,-1,1); //costheta between muon and tracked protons
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

   //vertex resolution
   h_vertex_resolution_neutrino = tfs->make<TH1D>("vertex_resolution_neutrino","Vertex resolution for neutrino vertexes",1000,0,100);
   h_vertex_resolution_proton = tfs->make<TH1D>("vertex_resolution_proton","Vertex resolution for proton vertexes",1000,0,100);
   h_vertex_resolution_muon = tfs->make<TH1D>("vertex_resolution_muon","Vertex resolution for muon vertexes",1000,0,100);
   h_vertexfitter_resolution_neutrino = tfs->make<TH1D>("vertexfitter_resolution_neutrino","Vertex (fitter) resolution for neutrino vertexes",1000,0,100);
   h_vertexfitter_chi2ndf_neutrino = tfs->make<TH1D>("vertexfitter_chi2ndf_neutrino","Vertex (fitter) chi2ndf for neutrino vertexes",1000,0,100);
   h_vertexfitter_resolution_proton = tfs->make<TH1D>("vertexfitter_resolution_proton","Vertex (fitter) resolution for proton vertexes",1000,0,100);
   h_vertexfitter_chi2ndf_proton = tfs->make<TH1D>("vertexfitter_chi2ndf_proton","Vertex (fitter) chi2ndf for proton vertexes",1000,0,100);
   h_vertexfitter_resolution_muon = tfs->make<TH1D>("vertexfitter_resolution_muon","Vertex (fitter) resolution for muon vertexes",1000,0,100);
   h_vertexfitter_chi2ndf_muon = tfs->make<TH1D>("vertexfitter_chi2ndf_muon","Vertex (fitter) chi2ndf for muon vertexes",1000,0,100);
   h_vertex_resolution_neutrino_not_merged = tfs->make<TH1D>("vertex_resolution_neutrino_not_merged","Vertex resolution for neutrino vertexes for NON merged protons",1000,0,100);
   h_vertexfitter_resolution_neutrino_not_merged = tfs->make<TH1D>("vertexfitter_resolution_neutrino_not_merged","Vertex (fitter) resolution for neutrino vertexes for NON merged protons",1000,0,100);
   h_vertexfitter_chi2ndf_neutrino_not_merged = tfs->make<TH1D>("vertexfitter_chi2_neutrino_not_merged","Vertex (fitter) chi2ndf for neutrino vertexes for NON merged protons",1000,0,100);
   h_vertex_resolution_proton_not_merged = tfs->make<TH1D>("vertex_resolution_proton_not_merged","Vertex resolution for proton vertexes for NON merged protons",1000,0,100);
   h_vertexfitter_resolution_proton_not_merged = tfs->make<TH1D>("vertexfitter_resolution_proton_not_merged","Vertex (fitter) resolution for proton vertexes for NON merged protons",1000,0,100);
   h_vertexfitter_chi2ndf_proton_not_merged = tfs->make<TH1D>("vertexfitter_chi2ndf_proton_not_merged","Vertex (fitter) chi2ndf for proton vertexes for NON merged protons",1000,0,100);
   h_vertex_resolution_muon_not_merged = tfs->make<TH1D>("vertex_resolution_muon_not_merged","Vertex resolution for muon vertexes for NON merged protons",1000,0,100);
   h_vertexfitter_resolution_muon_not_merged = tfs->make<TH1D>("vertexfitter_resolution_muon_not_merged","Vertex (fitter) resolution for muon vertexes for NON merged protons",1000,0,100);
   h_vertexfitter_chi2ndf_muon_not_merged = tfs->make<TH1D>("vertexfitter_chi2ndf_muon_not_merged","Vertex (fitter) chi2ndf for muon vertexes for NON merged protons",1000,0,100);
   h_vertex_resolution_neutrino_merged = tfs->make<TH1D>("vertex_resolution_neutrino_merged","Vertex resolution for neutrino vertexes for merged protons",1000,0,100);
   h_vertexfitter_resolution_neutrino_merged = tfs->make<TH1D>("vertexfitter_resolution_neutrino_merged","Vertex (fitter) resolution for neutrino vertexes for merged protons",1000,0,100);
   h_vertexfitter_chi2ndf_neutrino_merged = tfs->make<TH1D>("vertexfitter_chi2_neutrino_merged","Vertex (fitter) chi2ndf for neutrino vertexes for merged protons",1000,0,100);
   h_vertex_resolution_proton_merged = tfs->make<TH1D>("vertex_resolution_proton_merged","Vertex resolution for proton vertexes for merged protons",1000,0,100);
   h_vertexfitter_resolution_proton_merged = tfs->make<TH1D>("vertexfitter_resolution_proton_merged","Vertex (fitter) resolution for proton vertexes for merged protons",1000,0,100);
   h_vertexfitter_chi2ndf_proton_merged = tfs->make<TH1D>("vertexfitter_chi2ndf_proton_merged","Vertex (fitter) chi2ndf for proton vertexes for merged protons",1000,0,100);
   h_vertex_resolution_muon_merged = tfs->make<TH1D>("vertex_resolution_muon_merged","Vertex resolution for muon vertexes for merged protons",1000,0,100);
   h_vertexfitter_resolution_muon_merged = tfs->make<TH1D>("vertexfitter_resolution_muon_merged","Vertex (fitter) resolution for muon vertexes for merged protons",1000,0,100);
   h_vertexfitter_chi2ndf_muon_merged = tfs->make<TH1D>("vertexfitter_chi2ndf_muon_merged","Vertex (fitter) chi2ndf for muon vertexes for merged protons",1000,0,100);
   
   h_vertex_resolution_vs_not_tracked_above20MeV = tfs->make<TH2D>("h_vertex_resolution_vs_not_tracked_above20MeV","Vertex resolution for neutrino vertexes vs fraction of not tracked above 20MeV w.r.t to total # of protons",200,0,20,100,0,1);
   h_vertex_resolution_vs_not_tracked_below20MeV = tfs->make<TH2D>("h_vertex_resolution_vs_not_tracked_below","Vertex resolution for neutrino vertexes vs fraction of not tracked below 20MeV w.r.t to total # of protons",200,0,20,100,0,1);
   h_vertex_resolution_vs_not_tracked = tfs->make<TH2D>("h_vertex_resolution_vs_not_tracked","Vertex resolution for neutrino vertexes vs fraction of not tracked w.r.t to total # of protons",200,0,20,100,0,1);
   h_vertexfitter_resolution_vs_not_tracked_above20MeV = tfs->make<TH2D>("h_vertexfitter_resolution_vs_not_tracked_above20MeV","Vertex (fitter) resolution for neutrino vertexes vs fraction of not tracked above 20MeV w.r.t to total # of protons",200,0,20,100,0,1);
   h_vertexfitter_resolution_vs_not_tracked_below20MeV = tfs->make<TH2D>("h_vertexfitter_resolution_vs_not_tracked_below","Vertex (fitter) resolution for neutrino vertexes vs fraction of not tracked below 20MeV w.r.t to total # of protons",200,0,20,100,0,1);
   h_vertexfitter_resolution_vs_not_tracked = tfs->make<TH2D>("h_vertexfitter_resolution_vs_not_tracked","Vertex (fitter) resolution for neutrino vertexes vs fraction of not tracked w.r.t to total # of protons",200,0,20,100,0,1);
   
   h_shower_pdg = tfs->make<TH1D>("shower_pdg","PDG of a PFP recon as a shower",10000,0,10000);
   h_n_proton_showers = tfs->make<TH1D>("n_proton_showers","Number of protons recon as shower per event",20,0,20);
   h_shower_proton_kinE = tfs->make<TH1D>("shower_proton_kinE","Shower - Proton kinE; Kinetic Energy (GeV)",1000,0,2); 
   h_shower_proton_l = tfs->make<TH1D>("shower_proton_l","Shower - Proton Length; True length (cm)",1000,0,200); 
   h_shower_proton_nhits = tfs->make<TH1D>("shower_proton_nhits","Shower - Proton nhits; nhits",1000,0,1000);
   h_shower_proton_costheta_muon = tfs->make<TH1D>("shower_proton_costheta_muon","Cos #theta between muon and protons for SHOWER protons;cos #theta",1000,-1,1); 	    
   
   art::TFileDirectory hits_dir = tfs->mkdir("hits_dir");
   //hits analysis
   h_tracked_not_clustered_distance_nuvtx = hits_dir.make<TH1D>("tracked_not_clustered_distance_nuvtx", "Distance between hit and nu vertex for not clustered hits for tracked particles", 1000,0,100 );
   h_tracked_not_clustered_muon_start = hits_dir.make<TH1D>("tracked_not_clustered_muon_start", "Distance between hit and muon start position for not clustered hits for tracked muons", 1000,0,100 );
   h_tracked_not_clustered_muon_end = hits_dir.make<TH1D>("tracked_not_clustered_muon_end", "Distance between hit and muon end position for not clustered hits for tracked muons", 1000,0,100 );
   h_tracked_not_clustered_proton_start = hits_dir.make<TH2D>("tracked_not_clustered_proton_start", "Distance between hit and proton  start position for not clustered hits for tracked protons vs length", 1000,0,100, 1000, 0, 100 );
   h_tracked_not_clustered_proton_end = hits_dir.make<TH2D>("tracked_not_clustered_proton_end", "Distance between hit and proton end position for not clustered hits for tracked protons vs length", 1000,0,100, 1000,0,100 );
   h_hits_not_clustered_tracked_charge = hits_dir.make<TH1D>("hits_not_clustered_tracked_charge", "Charge of not clustered hits but tracked particles", 1000,0,1000 );
   h_hits_not_clustered_tracked_charge_proton = hits_dir.make<TH1D>("hits_not_clustered_tracked_charge_proton", "Charge of not clustered hits but tracked particles - protons", 1000,0,1000 );
   h_hits_not_clustered_tracked_charge_muon = hits_dir.make<TH1D>("hits_not_clustered_tracked_charge_muon", "Charge of not clustered hits but tracked particles - muons", 1000,0,1000 );
   h_fraction_pdgs_not_tracked_not_clustered = hits_dir.make<TH1D>("fraction_pdgs_not_tracked_not_clustered", "Fraction of not-tracked not-clustered hits for muon (0) and protons (1) (2for protons >20MeV, 3 for <20MeV)", 4,-0.5,3.5);
   h_not_tracked_not_clustered_muon_start = hits_dir.make<TH1D>("not_tracked_not_clustered_muon_start", "Distance between hit and muon start position for not clustered hits for NON tracked muons", 1000,0,100 );
   h_not_tracked_not_clustered_muon_end = hits_dir.make<TH1D>("not_tracked_not_clustered_muon_end", "Distance between hit and muon end position for not clustered hits for NON tracked muons", 1000,0,100 );
   h_not_tracked_not_clustered_proton_start = hits_dir.make<TH2D>("not_tracked_not_clustered_proton_start", "Distance between hit and proton  start position for not clustered hits for NON tracked protons vs length", 1000,0,100, 1000, 0, 100 );
   h_not_tracked_not_clustered_proton_end = hits_dir.make<TH2D>("not_tracked_not_clustered_proton_end", "Distance between hit and proton end position for not clustered hits for NON tracked protons vs length", 1000,0,100, 1000,0,100 );
   h_hits_not_clustered_not_tracked_charge_muon = hits_dir.make<TH1D>("hits_not_clustered_not_tracked_charge_muon", "Charge of not clustered hits and NON tracked particles - muons", 1000,0,1000 );
   h_hits_not_clustered_not_tracked_charge_proton = hits_dir.make<TH1D>("hits_not_clustered_not_tracked_charge_proton", "Charge of not clustered hits and NON tracked particles - protons", 1000,0,1000 );
   
   h_muon_clustering_prob_good_protons = hits_dir.make<TH1D>("muon_clustering_prob_good_protons","Probability of not-clustering (0), clustering correctly (1) and clustering wrongly (2) for muon in events w/ all reco'ed protons",3,-0.5,2.5);
   h_muon_clustering_mismatch_pdg_good_protons = hits_dir.make<TH1D>("muon_clustering_mismatch_pdg_good_protons","PDG of the particle the muon hits are wrongly assigned to",10000,0,10000);
   h_muon_not_clustered_reco_hits_good_protons = hits_dir.make<TH2D>("muon_not_clustered_reco_hits_good_protons","Number of not clustered hits vs number of hits for muon in \"good proton\" events", 500,0,500,1000,0,10000);
   h_muon_clustered_matched_reco_hits_good_protons = hits_dir.make<TH2D>("muon_clustered_matched_reco_hits_good_protons","Number of clustered and matched hits vs number of hits for muon in \"good proton\" events", 500,0,500,1000,0,10000);
   h_muon_clustered_mismatched_reco_hits_good_protons = hits_dir.make<TH2D>("muon_clustered_mismatched_reco_hits_good_protons","Number of clustered and mismatched hits vs number of hits for muon in \"good proton\" events", 500,0,500,1000,0,10000);
   h_muon_NC_lateral_hits_good_protons = hits_dir.make<TH2D>("muon_NC_lateral_hits_good_protons","Fraction of not clustered hits/total hits for muons in \"good proton\" events VS lateral distance between muon and proton",1000,0,1,800,0,200);
   h_muon_NC_costheta_hits_good_protons = hits_dir.make<TH2D>("muon_NC_costheta_hits_good_protons","Fraction of not clustered hits/total hits for muons in \"good proton\" events VS costheta between muon and proton",500,0,1,100,-1,1);
   h_muon_CMA_lateral_hits_good_protons = hits_dir.make<TH2D>("muon_CMA_lateral_hits_good_protons","Fraction of clusterend and matched hits/total hits for muons in \"good proton\" events VS lateral distance between muon and proton",1000,0,1,800,0,200);
   h_muon_CMA_costheta_hits_good_protons = hits_dir.make<TH2D>("muon_CMA_costheta_hits_good_protons","Fraction of clustered and matched hits/total hits for muons in \"good proton\" events VS costheta between muon and proton",500,0,1,100,-1,1);
   h_muon_CMI_lateral_hits_good_protons = hits_dir.make<TH2D>("muon_CMI_lateral_hits_good_protons","Fraction of clustered and mismatched hits/total hits for muons in \"good proton\" events VS lateral distance between muon and proton",1000,0,1,800,0,200);
   h_muon_CMI_costheta_hits_good_protons = hits_dir.make<TH2D>("muon_CMI_costheta_hits_good_protons","Fraction of clustered and mismatched hits/total hits for muons in \"good proton\" events VS costheta between muon and proton",500,0,1,100,-1,1);
   h_muon_not_clustered_reco_charge_good_protons = hits_dir.make<TH2D>("muon_not_clustered_reco_charge_good_protons", "Not clustered hit charge vs fraction of non clustered/total hits for muons in \"good proton\" events",1000,0,1000,500,0,1);
   h_muon_NC_lateral_charge_good_protons = hits_dir.make<TH2D>("muon_NC_lateral_charge_good_protons","Not clustered hit charge vs lateral distance muon - proton for muons in \"good proton\" events", 1000,0,1000,800,0,200);
   h_muon_NC_costheta_charge_good_protons = hits_dir.make<TH2D>("muon_NC_costheta_charge_good_protons","Not clustered hit charge vs costheta muon - proton for muons in \"good proton\" events", 1000,0,1000,100,-1,1);
   h_muon_clustered_matched_reco_charge_good_protons = hits_dir.make<TH2D>("muon_clustered_matched_reco_charge_good_protons", "Clustered and matched hit charge vs fraction of non clustered/total hits for muons in \"good proton\" events",1000,0,1000,500,0,1);
   h_muon_CMA_lateral_charge_good_protons = hits_dir.make<TH2D>("muon_CMA_lateral_charge_good_protons","Clustered and matched hit charge vs lateral distance muon - proton for muons in \"good proton\" events", 1000,0,1000,800,0,200);
   h_muon_CMA_costheta_charge_good_protons = hits_dir.make<TH2D>("muon_CMA_costheta_charge_good_protons","Clustered and matched hit charge vs costheta muon - proton for muons in \"good proton\" events", 1000,0,1000,100,-1,1);
   h_muon_clustered_mismatched_reco_charge_good_protons = hits_dir.make<TH2D>("muon_clustered_mismatched_reco_charge_good_protons", "Clustered and mismatched hit charge vs fraction of non clustered/total hits for muons in \"good proton\" events",1000,0,1000,500,0,1);
   h_muon_CMI_lateral_charge_good_protons = hits_dir.make<TH2D>("muon_CMI_lateral_charge_good_protons","Clustered and mismatched hit charge vs lateral distance muon - proton for muons in \"good proton\" events", 1000,0,1000,800,0,200);
   h_muon_CMI_costheta_charge_good_protons = hits_dir.make<TH2D>("muon_CMI_costheta_charge_good_protons","Clustered and mismatchedhit charge vs costheta muon - proton for muons in \"good proton\" events", 1000,0,1000,100,-1,1);

   h_proton_clustering_prob_good_protons = hits_dir.make<TH1D>("proton_clustering_prob_good_protons","Probability of not-clustering (0), clustering correctly (1) and clustering wrongly (2) for proton in events w/ all reco'ed protons",3,-0.5,2.5);
   h_proton_clustering_mismatch_pdg_good_protons = hits_dir.make<TH1D>("proton_clustering_mismatch_pdg_good_protons","PDG of the particle the proton hits are wrongly assigned to",10000,0,10000);
   h_proton_not_clustered_reco_hits_good_protons = hits_dir.make<TH2D>("proton_not_clustered_reco_hits_good_protons","Number of not clustered hits vs number of hits for proton in \"good proton\" events", 500,0,500,1000,0,10000);
   h_proton_NC_lateral_hits_good_protons = hits_dir.make<TH2D>("proton_NC_lateral_hits_good_protons","Fraction of not clustered hits/total hits for protons in \"good proton\" events VS lateral distance between muon and proton",1000,0,1,800,0,200);
   h_proton_NC_costheta_hits_good_protons = hits_dir.make<TH2D>("proton_NC_costheta_hits_good_protons","Fraction of not clustered hits/total hits for protons in \"good proton\" events VS costheta between muon and proton",500,0,1,100,-1,1);
   h_proton_clustered_matched_reco_hits_good_protons = hits_dir.make<TH2D>("proton_clustered_matched_reco_hits_good_protons","Number of clustered and matched hits vs number of hits for proton in \"good proton\" events", 500,0,500,1000,0,10000);
   h_proton_CMA_lateral_hits_good_protons = hits_dir.make<TH2D>("proton_CMA_lateral_hits_good_protons","Fraction of clustered and matched hits/total hits for protons in \"good proton\" events VS lateral distance between muon and proton",1000,0,1,800,0,200);
   h_proton_CMA_costheta_hits_good_protons = hits_dir.make<TH2D>("proton_CMA_costheta_hits_good_protons","Fraction of clustered and matched hits/total hits for protons in \"good proton\" events VS costheta between muon and proton",500,0,1,100,-1,1);
   h_proton_clustered_mismatched_reco_hits_good_protons = hits_dir.make<TH2D>("proton_clustered_mismatched_reco_hits_good_protons","Number of clustered and mismatched hits vs number of hits for proton in \"good proton\" events", 500,0,500,1000,0,10000);
   h_proton_CMI_lateral_hits_good_protons = hits_dir.make<TH2D>("proton_CMI_lateral_hits_good_protons","Fraction of clustered and mismatched hits/total hits for protons in \"good proton\" events VS lateral distance between muon and proton",1000,0,1,800,0,200);
   h_proton_CMI_costheta_hits_good_protons = hits_dir.make<TH2D>("proton_CMI_costheta_hits_good_protons","Fraction of clustered and mismatched hits/total hits for protons in \"good proton\" events VS costheta between muon and proton",500,0,1,100,-1,1);
   h_proton_not_clustered_reco_charge_good_protons = hits_dir.make<TH2D>("proton_not_clustered_reco_charge_good_protons", "Not clustered hit charge vs fraction of non clustered/total hits for protons in \"good proton\" events",1000,0,1000,500,0,1);
   h_proton_NC_lateral_charge_good_protons = hits_dir.make<TH2D>("proton_NC_lateral_charge_good_protons","Not clustered hit charge vs lateral distance muon - proton for protons in \"good proton\" events", 1000,0,1000,800,0,200);
   h_proton_NC_costheta_charge_good_protons = hits_dir.make<TH2D>("proton_NC_costheta_charge_good_protons","Not clustered hit charge vs costheta muon - proton for protons in \"good proton\" events", 1000,0,1000,100,-1,1);
   h_proton_clustered_matched_reco_charge_good_protons = hits_dir.make<TH2D>("proton_clustered_matched_reco_charge_good_protons", "Clustered and matched hit charge vs fraction of non clustered/total hits for protons in \"good proton\" events",1000,0,1000,500,0,1);
   h_proton_CMA_lateral_charge_good_protons = hits_dir.make<TH2D>("proton_CMA_lateral_charge_good_protons","Clustered and matched hit charge vs lateral distance muon - proton for protons in \"good proton\" events", 1000,0,1000,800,0,200);
   h_proton_CMA_costheta_charge_good_protons = hits_dir.make<TH2D>("proton_CMA_costheta_charge_good_protons","Clustered and matched hit charge vs costheta muon - proton for protons in \"good proton\" events", 1000,0,1000,100,-1,1);
   h_proton_clustered_mismatched_reco_charge_good_protons = hits_dir.make<TH2D>("proton_clustered_mismatched_reco_charge_good_protons", "Clustered and mismatched hit charge vs fraction of non clustered/total hits for protons in \"good proton\" events",1000,0,1000,500,0,1);
   h_proton_CMI_lateral_charge_good_protons = hits_dir.make<TH2D>("proton_CMI_lateral_charge_good_protons","Clustered and mismatched hit charge vs lateral distance muon - proton for protons in \"good proton\" events", 1000,0,1000,800,0,200);
   h_proton_CMI_costheta_charge_good_protons = hits_dir.make<TH2D>("proton_CMI_costheta_charge_good_protons","Clustered and mismatched hit charge vs costheta muon - proton for protons in \"good proton\" events", 1000,0,1000,100,-1,1);

   h_muon_clustering_prob_bad_protons = hits_dir.make<TH1D>("muon_clustering_prob_bad_protons","Probability of not-clustering (0), clustering correctly (1) and clustering wrongly (2) for muon in events w/ all reco'ed protons",3,-0.5,2.5);
   h_muon_clustering_mismatch_pdg_bad_protons = hits_dir.make<TH1D>("muon_clustering_mismatch_pdg_bad_protons","PDG of the particle the muon hits are wrongly assigned to",10000,0,10000);
   h_muon_not_clustered_reco_hits_bad_protons = hits_dir.make<TH2D>("muon_not_clustered_reco_hits_bad_protons","Number of not clustered hits vs number of hits for muon in \"bad proton\" events", 500,0,500,1000,0,10000);
   h_muon_NC_lateral_hits_bad_protons = hits_dir.make<TH2D>("muon_NC_lateral_hits_bad_protons","Fraction of not clustered hits/total hits for muons in \"bad proton\" events VS lateral distance between muon and proton",1000,0,1,800,0,200);
   h_muon_NC_costheta_hits_bad_protons = hits_dir.make<TH2D>("muon_NC_costheta_hits_bad_protons","Fraction of not clustered hits/total hits for muons in \"bad proton\" events VS costheta between muon and proton",500,0,1,100,-1,1);
   h_muon_clustered_matched_reco_hits_bad_protons = hits_dir.make<TH2D>("muon_clustered_matched_reco_hits_bad_protons","Number of clustered and matched hits vs number of hits for muon in \"bad proton\" events", 500,0,500,1000,0,10000);
   h_muon_CMA_lateral_hits_bad_protons = hits_dir.make<TH2D>("muon_CMA_lateral_hits_bad_protons","Fraction of clustered and matched hits/total hits for muons in \"bad proton\" events VS lateral distance between muon and proton",1000,0,1,800,0,200);
   h_muon_CMA_costheta_hits_bad_protons = hits_dir.make<TH2D>("muon_CMA_costheta_hits_bad_protons","Fraction of clustered and matched hits/total hits for muons in \"bad proton\" events VS costheta between muon and proton",500,0,1,100,-1,1);
   h_muon_clustered_mismatched_reco_hits_bad_protons = hits_dir.make<TH2D>("muon_clustered_mismatched_reco_hits_bad_protons","Number of clustered and mismatched hits vs number of hits for muon in \"bad proton\" events", 500,0,500,1000,0,10000);
   h_muon_CMI_lateral_hits_bad_protons = hits_dir.make<TH2D>("muon_CMI_lateral_hits_bad_protons","Fraction of clustered and mismatched hits/total hits for muons in \"bad proton\" events VS lateral distance between muon and proton",1000,0,1,800,0,200);
   h_muon_CMI_costheta_hits_bad_protons = hits_dir.make<TH2D>("muon_CMI_costheta_hits_bad_protons","Fraction of clustered and mismatched hits/total hits for muons in \"bad proton\" events VS costheta between muon and proton",500,0,1,100,-1,1);
   h_muon_not_clustered_reco_charge_bad_protons = hits_dir.make<TH2D>("muon_not_clustered_reco_charge_bad_protons", "Not clustered hit charge vs fraction of non clustered/total hits for muons in \"bad proton\" events",1000,0,1000,500,0,1);
   h_muon_NC_lateral_charge_bad_protons = hits_dir.make<TH2D>("muon_NC_lateral_charge_bad_protons","Not clustered hit charge vs lateral distance muon - proton for muons in \"bad proton\" events", 1000,0,1000,800,0,200);
   h_muon_NC_costheta_charge_bad_protons = hits_dir.make<TH2D>("muon_NC_costheta_charge_bad_protons","Not clustered hit charge vs costheta muon - proton for muons in \"bad proton\" events", 1000,0,1000,100,-1,1);
   h_muon_clustered_matched_reco_charge_bad_protons = hits_dir.make<TH2D>("muon_clustered_matched_reco_charge_bad_protons", "Clustered and matched hit charge vs fraction of non clustered/total hits for muons in \"bad proton\" events",1000,0,1000,500,0,1);
   h_muon_CMA_lateral_charge_bad_protons = hits_dir.make<TH2D>("muon_CMA_lateral_charge_bad_protons","Clustered and matched hit charge vs lateral distance muon - proton for muons in \"bad proton\" events", 1000,0,1000,800,0,200);
   h_muon_CMA_costheta_charge_bad_protons = hits_dir.make<TH2D>("muon_CMA_costheta_charge_bad_protons","Clustered and matched hit charge vs costheta muon - proton for muons in \"bad proton\" events", 1000,0,1000,100,-1,1);
   h_muon_clustered_mismatched_reco_charge_bad_protons = hits_dir.make<TH2D>("muon_clustered_mismatched_reco_charge_bad_protons", "Clustered and mismatched hit charge vs fraction of non clustered/total hits for muons in \"bad proton\" events",1000,0,1000,500,0,1);
   h_muon_CMI_lateral_charge_bad_protons = hits_dir.make<TH2D>("muon_CMI_lateral_charge_bad_protons","Clustered and mismatched hit charge vs lateral distance muon - proton for muons in \"bad proton\" events", 1000,0,1000,800,0,200);
   h_muon_CMI_costheta_charge_bad_protons = hits_dir.make<TH2D>("muon_CMI_costheta_charge_bad_protons","Clustered and mismatched hit charge vs costheta muon - proton for muons in \"bad proton\" events", 1000,0,1000,100,-1,1);

   h_proton_clustering_prob_bad_protons = hits_dir.make<TH1D>("proton_clustering_prob_bad_protons","Probability of not-clustering (0), clustering correctly (1) and clustering wrongly (2) for proton in events w/ all reco'ed protons",3,-0.5,2.5);
   h_proton_clustering_mismatch_pdg_bad_protons = hits_dir.make<TH1D>("proton_clustering_mismatch_pdg_bad_protons","PDG of the particle the proton hits are wrongly assigned to",10000,0,10000);
   h_proton_not_clustered_reco_hits_bad_protons = hits_dir.make<TH2D>("proton_not_clustered_reco_hits_bad_protons","Number of not clustered hits vs number of hits for proton in \"bad proton\" events", 500,0,500,1000,0,10000);
   h_proton_NC_lateral_hits_bad_protons = hits_dir.make<TH2D>("proton_NC_lateral_hits_bad_protons","Fraction of not clustered hits/total hits for protons in \"bad proton\" events VS lateral distance between muon and proton",100,0,1,800,0,200);
   h_proton_NC_costheta_hits_bad_protons = hits_dir.make<TH2D>("proton_NC_costheta_hits_bad_protons","Fraction of not clustered hits/total hits for protons in \"bad proton\" events VS costheta between muon and proton",500,0,1,100,-1,1);
   h_proton_clustered_matched_reco_hits_bad_protons = hits_dir.make<TH2D>("proton_clustered_matched_reco_hits_bad_protons","Number of clustered and matched hits vs number of hits for proton in \"bad proton\" events", 500,0,500,1000,0,10000);
   h_proton_CMA_lateral_hits_bad_protons = hits_dir.make<TH2D>("proton_CMA_lateral_hits_bad_protons","Fraction of clustered and matched hits/total hits for protons in \"bad proton\" events VS lateral distance between muon and proton",100,0,1,800,0,200);
   h_proton_CMA_costheta_hits_bad_protons = hits_dir.make<TH2D>("proton_CMA_costheta_hits_bad_protons","Fraction of clustered and matched hits/total hits for protons in \"bad proton\" events VS costheta between muon and proton",500,0,1,100,-1,1);
   h_proton_clustered_mismatched_reco_hits_bad_protons = hits_dir.make<TH2D>("proton_clustered_mismatched_reco_hits_bad_protons","Number of clustered and mismatched hits vs number of hits for proton in \"bad proton\" events", 500,0,500,1000,0,10000);
   h_proton_CMI_lateral_hits_bad_protons = hits_dir.make<TH2D>("proton_CMI_lateral_hits_bad_protons","Fraction of clustered and mismatched hits/total hits for protons in \"bad proton\" events VS lateral distance between muon and proton",100,0,1,800,0,200);
   h_proton_CMI_costheta_hits_bad_protons = hits_dir.make<TH2D>("proton_CMI_costheta_hits_bad_protons","Fraction of clustered and mismatched hits/total hits for protons in \"bad proton\" events VS costheta between muon and proton",500,0,1,100,-1,1);
   h_proton_not_clustered_reco_charge_bad_protons = hits_dir.make<TH2D>("proton_not_clustered_reco_charge_bad_protons", "Not clustered hit charge vs fraction of non clustered/total hits for protons in \"bad proton\" events",1000,0,1000,500,0,1);
   h_proton_NC_lateral_charge_bad_protons = hits_dir.make<TH2D>("proton_NC_lateral_charge_bad_protons","Not clustered hit charge vs lateral distance muon - proton for protons in \"bad proton\" events", 1000,0,1000,800,0,200);
   h_proton_NC_costheta_charge_bad_protons = hits_dir.make<TH2D>("proton_NC_costheta_charge_bad_protons","Not clustered hit charge vs costheta muon - proton for protons in \"bad proton\" events", 1000,0,1000,100,-1,1);
   h_proton_clustered_matched_reco_charge_bad_protons = hits_dir.make<TH2D>("proton_clustered_matched_reco_charge_bad_protons", "Clustered and matched hit charge vs fraction of non clustered/total hits for protons in \"bad proton\" events",1000,0,1000,500,0,1);
   h_proton_CMA_lateral_charge_bad_protons = hits_dir.make<TH2D>("proton_CMA_lateral_charge_bad_protons","Clustered and matched hit charge vs lateral distance muon - proton for protons in \"bad proton\" events", 1000,0,1000,800,0,200);
   h_proton_CMA_costheta_charge_bad_protons = hits_dir.make<TH2D>("proton_CMA_costheta_charge_bad_protons","Clustered and matched hit charge vs costheta muon - proton for protons in \"bad proton\" events", 1000,0,1000,100,-1,1);
   h_proton_clustered_mismatched_reco_charge_bad_protons = hits_dir.make<TH2D>("proton_clustered_mismatched_reco_charge_bad_protons", "Clustered and mismatched hit charge vs fraction of non clustered/total hits for protons in \"bad proton\" events",1000,0,1000,500,0,1);
   h_proton_CMI_lateral_charge_bad_protons = hits_dir.make<TH2D>("proton_CMI_lateral_charge_bad_protons","Clustered and mismatched hit charge vs lateral distance muon - proton for protons in \"bad proton\" events", 1000,0,1000,800,0,200);
   h_proton_CMI_costheta_charge_bad_protons = hits_dir.make<TH2D>("proton_CMI_costheta_charge_bad_protons","Clustered and mismatched hit charge vs costheta muon - proton for protons in \"bad proton\" events", 1000,0,1000,100,-1,1);


   h_proton_clustering_prob_bad_protons_not_tracked = hits_dir.make<TH1D>("proton_clustering_prob_bad_protons_not_tracked","Probability of not-clustering (0), clustering correctly (1) and clustering wrongly (2) for proton in events w/ all reco'ed protons",3,-0.5,2.5);
   h_proton_clustering_mismatch_pdg_bad_protons_not_tracked = hits_dir.make<TH1D>("proton_clustering_mismatch_pdg_bad_protons_not_tracked","PDG of the particle the proton hits are wrongly assigned to",10000,0,10000);
   h_proton_not_clustered_reco_hits_bad_protons_not_tracked = hits_dir.make<TH2D>("proton_not_clustered_reco_hits_bad_protons_not_tracked","Number of not clustered hits vs number of hits for proton in \"bad proton\" events", 500,0,500,1000,0,10000);
   h_proton_NC_lateral_hits_bad_protons_not_tracked = hits_dir.make<TH2D>("proton_NC_lateral_hits_bad_protons_not_tracked","Fraction of not clustered hits/total hits for protons in \"bad proton\" events VS lateral distance between muon and proton",100,0,1,800,0,200);
   h_proton_NC_costheta_hits_bad_protons_not_tracked = hits_dir.make<TH2D>("proton_NC_costheta_hits_bad_protons_not_tracked","Fraction of not clustered hits/total hits for protons in \"bad proton\" events VS costheta between muon and proton",500,0,1,100,-1,1);
   h_proton_clustered_matched_reco_hits_bad_protons_not_tracked = hits_dir.make<TH2D>("proton_clustered_matched_reco_hits_bad_protons_not_tracked","Number of clustered and matched hits vs number of hits for proton in \"bad proton\" events", 500,0,500,1000,0,10000);
   h_proton_CMA_lateral_hits_bad_protons_not_tracked = hits_dir.make<TH2D>("proton_CMA_lateral_hits_bad_protons_not_tracked","Fraction of clustered and matched hits/total hits for protons in \"bad proton\" events VS lateral distance between muon and proton",100,0,1,800,0,200);
   h_proton_CMA_costheta_hits_bad_protons_not_tracked = hits_dir.make<TH2D>("proton_CMA_costheta_hits_bad_protons_not_tracked","Fraction of clustered and matched hits/total hits for protons in \"bad proton\" events VS costheta between muon and proton",500,0,1,100,-1,1);
   h_proton_clustered_mismatched_reco_hits_bad_protons_not_tracked = hits_dir.make<TH2D>("proton_clustered_mismatched_reco_hits_bad_protons_not_tracked","Number of clustered and mismatched hits vs number of hits for proton in \"bad proton\" events", 500,0,500,1000,0,10000);
   h_proton_CMI_lateral_hits_bad_protons_not_tracked = hits_dir.make<TH2D>("proton_CMI_lateral_hits_bad_protons_not_tracked","Fraction of clustered and mismatched hits/total hits for protons in \"bad proton\" events VS lateral distance between muon and proton",100,0,1,800,0,200);
   h_proton_CMI_costheta_hits_bad_protons_not_tracked = hits_dir.make<TH2D>("proton_CMI_costheta_hits_bad_protons_not_tracked","Fraction of clustered and mismatched hits/total hits for protons in \"bad proton\" events VS costheta between muon and proton",500,0,1,100,-1,1);
   h_proton_not_clustered_reco_charge_bad_protons_not_tracked = hits_dir.make<TH2D>("proton_not_clustered_reco_charge_bad_protons_not_tracked", "Not clustered hit charge vs fraction of non clustered/total hits for protons in \"bad proton\" events",1000,0,1000,500,0,1);
   h_proton_NC_lateral_charge_bad_protons_not_tracked = hits_dir.make<TH2D>("proton_NC_lateral_charge_bad_protons_not_tracked","Not clustered hit charge vs lateral distance muon - proton for protons in \"bad proton\" events", 1000,0,1000,800,0,200);
   h_proton_NC_costheta_charge_bad_protons_not_tracked = hits_dir.make<TH2D>("proton_NC_costheta_charge_bad_protons_not_tracked","Not clustered hit charge vs costheta muon - proton for protons in \"bad proton\" events", 1000,0,1000,100,-1,1);
   h_proton_clustered_matched_reco_charge_bad_protons_not_tracked = hits_dir.make<TH2D>("proton_clustered_matched_reco_charge_bad_protons_not_tracked", "Clustered and matched hit charge vs fraction of non clustered/total hits for protons in \"bad proton\" events",1000,0,1000,500,0,1);
   h_proton_CMA_lateral_charge_bad_protons_not_tracked = hits_dir.make<TH2D>("proton_CMA_lateral_charge_bad_protons_not_tracked","Clustered and matched hit charge vs lateral distance muon - proton for protons in \"bad proton\" events", 1000,0,1000,800,0,200);
   h_proton_CMA_costheta_charge_bad_protons_not_tracked = hits_dir.make<TH2D>("proton_CMA_costheta_charge_bad_protons_not_tracked","Clustered and matched hit charge vs costheta muon - proton for protons in \"bad proton\" events", 1000,0,1000,100,-1,1);
   h_proton_clustered_mismatched_reco_charge_bad_protons_not_tracked = hits_dir.make<TH2D>("proton_clustered_mismatched_reco_charge_bad_protons_not_tracked", "Clustered and mismatched hit charge vs fraction of non clustered/total hits for protons in \"bad proton\" events",1000,0,1000,500,0,1);
   h_proton_CMI_lateral_charge_bad_protons_not_tracked = hits_dir.make<TH2D>("proton_CMI_lateral_charge_bad_protons_not_tracked","Clustered and mismatched hit charge vs lateral distance muon - proton for protons in \"bad proton\" events", 1000,0,1000,800,0,200);
   h_proton_CMI_costheta_charge_bad_protons_not_tracked = hits_dir.make<TH2D>("proton_CMI_costheta_charge_bad_protons_not_tracked","Clustered and mismatched hit charge vs costheta muon - proton for protons in \"bad proton\" events", 1000,0,1000,100,-1,1);


   h_muon_clustering_prob_low_protons = hits_dir.make<TH1D>("muon_clustering_prob_low_protons","Probability of not-clustering (0), clustering correctly (1) and clustering wrongly (2) for muon in events w/ all reco'ed protons",3,-0.5,2.5);
   h_muon_clustering_mismatch_pdg_low_protons = hits_dir.make<TH1D>("muon_clustering_mismatch_pdg_low_protons","PDG of the particle the muon hits are wrongly assigned to",10000,0,10000);
   h_muon_not_clustered_reco_hits_low_protons = hits_dir.make<TH2D>("muon_not_clustered_reco_hits_low_protons","Number of not clustered hits vs number of hits for muon in \"bad proton\" events", 500,0,500,1000,0,10000);
   h_muon_NC_lateral_hits_low_protons = hits_dir.make<TH2D>("muon_NC_lateral_hits_low_protons","Fraction of not clustered hits/total hits for muons in \"bad proton\" events VS lateral distance between muon and proton",100,0,1,800,0,200);
   h_muon_NC_costheta_hits_low_protons = hits_dir.make<TH2D>("muon_NC_costheta_hits_low_protons","Fraction of not clustered hits/total hits for muons in \"bad proton\" events VS costheta between muon and proton",500,0,1,100,-1,1);
   h_muon_clustered_matched_reco_hits_low_protons = hits_dir.make<TH2D>("muon_clustered_matched_reco_hits_low_protons","Number of clustered and matched hits vs number of hits for muon in \"bad proton\" events", 500,0,500,1000,0,10000);
   h_muon_CMA_lateral_hits_low_protons = hits_dir.make<TH2D>("muon_CMA_lateral_hits_low_protons","Fraction of clustered and matched hits/total hits for muons in \"bad proton\" events VS lateral distance between muon and proton",100,0,1,800,0,200);
   h_muon_CMA_costheta_hits_low_protons = hits_dir.make<TH2D>("muon_CMA_costheta_hits_low_protons","Fraction of clustered and matched hits/total hits for muons in \"bad proton\" events VS costheta between muon and proton",500,0,1,100,-1,1);
   h_muon_clustered_mismatched_reco_hits_low_protons = hits_dir.make<TH2D>("muon_clustered_mismatched_reco_hits_low_protons","Number of clustered and mismatched hits vs number of hits for muon in \"bad proton\" events", 500,0,500,1000,0,10000);
   h_muon_CMI_lateral_hits_low_protons = hits_dir.make<TH2D>("muon_CMI_lateral_hits_low_protons","Fraction of clustered and mismatched hits/total hits for muons in \"bad proton\" events VS lateral distance between muon and proton",100,0,1,800,0,200);
   h_muon_CMI_costheta_hits_low_protons = hits_dir.make<TH2D>("muon_CMI_costheta_hits_low_protons","Fraction of clustered and mismatched hits/total hits for muons in \"bad proton\" events VS costheta between muon and proton",500,0,1,100,-1,1);
   h_muon_not_clustered_reco_charge_low_protons = hits_dir.make<TH2D>("muon_not_clustered_reco_charge_low_protons", "Not clustered hit charge vs fraction of non clustered/total hits for muons in \"bad proton\" events",1000,0,1000,500,0,1);
   h_muon_NC_lateral_charge_low_protons = hits_dir.make<TH2D>("muon_NC_lateral_charge_low_protons","Not clustered hit charge vs lateral distance muon - proton for muons in \"bad proton\" events", 1000,0,1000,800,0,200);
   h_muon_NC_costheta_charge_low_protons = hits_dir.make<TH2D>("muon_NC_costheta_charge_low_protons","Not clustered hit charge vs costheta muon - proton for muons in \"bad proton\" events", 1000,0,1000,100,-1,1);
   h_muon_clustered_matched_reco_charge_low_protons = hits_dir.make<TH2D>("muon_clustered_matched_reco_charge_low_protons", "Clustered and matched hit charge vs fraction of non clustered/total hits for muons in \"bad proton\" events",1000,0,1000,500,0,1);
   h_muon_CMA_lateral_charge_low_protons = hits_dir.make<TH2D>("muon_CMA_lateral_charge_low_protons","Clustered and matched hit charge vs lateral distance muon - proton for muons in \"bad proton\" events", 1000,0,1000,800,0,200);
   h_muon_CMA_costheta_charge_low_protons = hits_dir.make<TH2D>("muon_CMA_costheta_charge_low_protons","Clustered and matched hit charge vs costheta muon - proton for muons in \"bad proton\" events", 1000,0,1000,100,-1,1);
   h_muon_clustered_mismatched_reco_charge_low_protons = hits_dir.make<TH2D>("muon_clustered_mismatched_reco_charge_low_protons", "Clustered and mismatched hit charge vs fraction of non clustered/total hits for muons in \"bad proton\" events",1000,0,1000,500,0,1);
   h_muon_CMI_lateral_charge_low_protons = hits_dir.make<TH2D>("muon_CMI_lateral_charge_low_protons","Clustered and mismatched hit charge vs lateral distance muon - proton for muons in \"bad proton\" events", 1000,0,1000,800,0,200);
   h_muon_CMI_costheta_charge_low_protons = hits_dir.make<TH2D>("muon_CMI_costheta_charge_low_protons","Clustered and mismatched hit charge vs costheta muon - proton for muons in \"bad proton\" events", 1000,0,1000,100,-1,1);

}

int recohelper::RecoBenchmarker::FillAnalysisHistograms( int& count_tracked, int& count_not_tracked, bool & is_lowmomentum_p ) {
  
    if (fccnc!=0 ) return -1; //comment if you're running over a std sample with NC interactions and you're interested in them
    unsigned muon_pos = -1;
    //check that the muon is reco
    bool reco_muon = false;
    bool is_pion=false;
    bool lowmomentum_p = false;
    for (unsigned i = 0; i < fpdg.size(); i++) {
	    if ( fpdg[i] == 13 ) { //ismuon
	    hmuon_length_all->Fill( flength[i] );
	    hmuon_spectrum_all->Fill( fkinE[i] );
	    if ( fis_tracked[i] &&  fmuon_dqdx.size() != 0 ) { //want that the muon is well matched and w/ dqdx info
		muon_pos = i; 
	    	hmuon_length->Fill( flength[i] );
	    	hmuon_spectrum->Fill( fkinE[i] );
		reco_muon = true;
		hmuon_pos_res->Fill( sqrt( pow( freco_startx[i]- fstart_x[i] ,2) + pow( freco_starty[i]- fstart_y[i] ,2) + pow( freco_startz[i]- fstart_z[i] ,2) ));
		float res = min ( sqrt( pow(freco_vertex_x[i] - fstart_x[i],2) + pow(freco_vertex_y[i] - fstart_y[i],2) + pow(freco_vertex_z[i] - fstart_z[i],2)) ,
					sqrt( pow(freco_vertex_x[i] - fend_x[i],2) + pow(freco_vertex_y[i] - fend_y[i],2) + pow(freco_vertex_z[i] - fend_z[i],2)) );
		h_vertex_resolution_muon->Fill( res );
		float res_fitter = min ( sqrt( pow(freco_vertexfitter_x[i] - fstart_x[i],2) + pow(freco_vertexfitter_y[i] - fstart_y[i],2) + pow(freco_vertexfitter_z[i] - fstart_z[i],2)) ,
					sqrt( pow(freco_vertexfitter_x[i] - fend_x[i],2) + pow(freco_vertexfitter_y[i] - fend_y[i],2) + pow(freco_vertexfitter_z[i] - fend_z[i],2)) );
		h_vertexfitter_resolution_muon->Fill( res_fitter );
		h_vertexfitter_chi2ndf_muon->Fill( freco_vertexfitter_chi2ndf[i] );
		if ( sqrt( pow( freco_startx[i]- fstart_x[i] ,2) + pow( freco_starty[i]- fstart_y[i] ,2) + pow( freco_startz[i]- fstart_z[i] ,2) ) > 50 ) {
			reco_muon = false; //skip these events, might have direction flipped
	    	}
    		}
	    } //record info on the muon
	    if ( fpdg[i] == 211 || fpdg[i] == -211 || fpdg[i] == 111 ) is_pion=true; //record if there are any pions
	    if (  fpdg[i] == 2212 && fp0[i] <= 0.2 ) lowmomentum_p = true;
	    if (  fpdg[i] == 2212 && fp0[i] >= 0.2 ) hproton_theta_mu->Fill(fcostheta_muon[i]);
    }
    
    if ( reco_muon == false ) return muon_pos; //select events with a reco muon
    
    n_muons++;

    //check if the neutrino reco'ed vertex is there 
    if ( !(abs(fnu_reco_x+1)<DBL_EPSILON && abs(fnu_reco_y+1)<DBL_EPSILON && abs(fnu_reco_z+1)<DBL_EPSILON) ) 
	h_vertex_resolution_neutrino->Fill( sqrt( pow(fnu_reco_x - fneutrino_x,2) + pow(fnu_reco_y - fneutrino_y,2) + pow(fnu_reco_z - fneutrino_z,2) ) );	    
    if ( !(abs(fnu_reco_fitter_x+1)<DBL_EPSILON && abs(fnu_reco_fitter_y+1)<DBL_EPSILON && abs(fnu_reco_fitter_z+1)<DBL_EPSILON) ) {
	h_vertexfitter_resolution_neutrino->Fill( sqrt( pow(fnu_reco_fitter_x - fneutrino_x,2) + pow(fnu_reco_fitter_y - fneutrino_y,2) + pow(fnu_reco_fitter_z - fneutrino_z,2) ) );	    
	h_vertexfitter_chi2ndf_neutrino->Fill( fnu_reco_fitter_chi2ndf );
    }

    count_not_tracked = 0;
    count_tracked = 0;
    is_lowmomentum_p = false;
    for (unsigned j=0; j<fpdg.size(); j++) {
	    if ( !reco_muon ) break; //select events with a reco muon
	    if ( fpdg[j]!=2212 ) continue; //watch only protons
	 
	    n_all_protons++;
	    hproton_p_all->Fill( fp0[j] );
	    hproton_l_all->Fill( flength[j] );
	    hproton_kinE_all->Fill( fkinE[j] );
	    
	    if ( abs( fcostheta_muon[j] ) > sqrt(3)/2. ) {
	    hproton_l_all_angle1->Fill( flength[j] );
	    hproton_kinE_all_angle1->Fill( fkinE[j] );
	    hproton_nhits_all_angle1->Fill( freco_mcp_collection_hits[j] ); //FIXME
	    } else if ( abs( fcostheta_muon[j] ) >= 0.5 && abs( fcostheta_muon[j] ) <= sqrt(3)/2. ) { 
	    hproton_l_all_angle2->Fill( flength[j] );
	    hproton_kinE_all_angle2->Fill( fkinE[j] );
	    hproton_nhits_all_angle2->Fill( freco_mcp_collection_hits[j] ); //FIXME
	    } else {
	    hproton_l_all_angle3->Fill( flength[j] );
	    hproton_kinE_all_angle3->Fill( fkinE[j] );
	    hproton_nhits_all_angle3->Fill( freco_mcp_collection_hits[j] ); //FIXME
	    }

	    h_theta_mu_length_all->Fill( fcostheta_muon[j], flength[j]);
	    float res = min ( sqrt( pow(freco_vertex_x[j] - fstart_x[j],2) + pow(freco_vertex_y[j] - fstart_y[j],2) + pow(freco_vertex_z[j] - fstart_z[j],2)) ,
					sqrt( pow(freco_vertex_x[j] - fend_x[j],2) + pow(freco_vertex_y[j] - fend_y[j],2) + pow(freco_vertex_z[j] - fend_z[j],2)) );
	    h_vertex_resolution_proton->Fill( res );
	    float res_fitter = min ( sqrt( pow(freco_vertexfitter_x[j] - fstart_x[j],2) + pow(freco_vertexfitter_y[j] - fstart_y[j],2) + pow(freco_vertexfitter_z[j] - fstart_z[j],2)) ,
					sqrt( pow(freco_vertexfitter_x[j] - fend_x[j],2) + pow(freco_vertexfitter_y[j] - fend_y[j],2) + pow(freco_vertexfitter_z[j] - fend_z[j],2)) );
	    h_vertexfitter_resolution_proton->Fill( res_fitter );
	    h_vertexfitter_chi2ndf_proton->Fill( freco_vertexfitter_chi2ndf[j] );

	    if ( fis_tracked[j] ) { //efficiency plots for all protons (but don't divide yet)
	    	hproton_p->Fill( fp0[j] );
	    	hproton_l->Fill( flength[j] );
	    	hproton_kinE->Fill( fkinE[j] );
	        h_theta_mu_length->Fill( fcostheta_muon[j], flength[j]);
   		hproton_nhits->Fill( fnhits[j] );
   		hproton_nhits_CP->Fill( freco_mcp_collection_hits[j] );
   		hproton_nhits_theta_mu->Fill ( fnhits[j], acos(fcostheta_muon[j]) ) ;
   		hproton_nhits_CP_theta_mu->Fill( freco_mcp_collection_hits[j], acos(fcostheta_muon[j]) );
	    
		if ( abs( fcostheta_muon[j] ) > sqrt(3)/2. ) { //theta < 30degrees
	    		hproton_l_tracked_angle1->Fill( flength[j] );
	    		hproton_kinE_tracked_angle1->Fill( fkinE[j] );
	    		hproton_nhits_tracked_angle1->Fill( freco_mcp_collection_hits[j] ); //FIXME
	    	} else if ( abs( fcostheta_muon[j] ) >= 0.5 && abs( fcostheta_muon[j] ) <= sqrt(3)/2. ) { 
	    		hproton_l_tracked_angle2->Fill( flength[j] );
	    		hproton_kinE_tracked_angle2->Fill( fkinE[j] );
	    		hproton_nhits_tracked_angle2->Fill( freco_mcp_collection_hits[j] ); //FIXME
	    	} else {
	    		hproton_l_tracked_angle3->Fill( flength[j] );
	    		hproton_kinE_tracked_angle3->Fill( fkinE[j] );
	    		hproton_nhits_tracked_angle3->Fill( freco_mcp_collection_hits[j] ); //FIXME
	  	}
	    }


    //now study all the protons with momentum > 200MeV
    if (fp0[j] > 0.2) {
	 tot_n_protons++;
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
   	} else {
		low_protons++;
	}
}
	
	hproton_merged_not_merged->Fill(0., double(count_tracked) );
	hproton_merged_not_merged->Fill(1., double(count_not_tracked) );
        
 	//neutrino vertex resolution vs proton recon performance:
    	if ( !(abs(fnu_reco_x+1)<DBL_EPSILON && abs(fnu_reco_y+1)<DBL_EPSILON && abs(fnu_reco_z+1)<DBL_EPSILON) ) {
	int tot_protons = count_not_tracked + count_tracked ;
	h_vertex_resolution_vs_not_tracked_above20MeV->Fill( sqrt( pow(fnu_reco_x - fneutrino_x,2) + pow(fnu_reco_y - fneutrino_y,2) + pow(fnu_reco_z - fneutrino_z,2) ), double(count_not_tracked)/tot_protons );
	h_vertex_resolution_vs_not_tracked->Fill( sqrt( pow(fnu_reco_x - fneutrino_x,2) + pow(fnu_reco_y - fneutrino_y,2) + pow(fnu_reco_z - fneutrino_z,2) ), double(count_not_tracked + low_protons)/ (tot_protons+low_protons) );
	h_vertex_resolution_vs_not_tracked_below20MeV->Fill( sqrt( pow(fnu_reco_x - fneutrino_x,2) + pow(fnu_reco_y - fneutrino_y,2) + pow(fnu_reco_z - fneutrino_z,2) ), double(low_protons)/ (tot_protons+low_protons) );
	h_vertexfitter_resolution_vs_not_tracked_above20MeV->Fill( sqrt( pow(fnu_reco_fitter_x - fneutrino_x,2) + pow(fnu_reco_fitter_y - fneutrino_y,2) + pow(fnu_reco_fitter_z - fneutrino_z,2) ), double(count_not_tracked)/tot_protons );
	h_vertexfitter_resolution_vs_not_tracked->Fill( sqrt( pow(fnu_reco_fitter_x - fneutrino_x,2) + pow(fnu_reco_fitter_y - fneutrino_y,2) + pow(fnu_reco_fitter_z - fneutrino_z,2) ), double(count_not_tracked + low_protons)/ (tot_protons+low_protons) );
	h_vertexfitter_resolution_vs_not_tracked_below20MeV->Fill( sqrt( pow(fnu_reco_fitter_x - fneutrino_x,2) + pow(fnu_reco_fitter_y - fneutrino_y,2) + pow(fnu_reco_fitter_z - fneutrino_z,2) ), double(low_protons)/ (tot_protons+low_protons) );
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
    		if ( !(abs(fnu_reco_x+1)<DBL_EPSILON && abs(fnu_reco_y+1)<DBL_EPSILON && abs(fnu_reco_z+1)<DBL_EPSILON) ) 
	        h_vertex_resolution_neutrino_not_merged->Fill( sqrt( pow(fnu_reco_x - fneutrino_x,2) + pow(fnu_reco_y - fneutrino_y,2) + pow(fnu_reco_z - fneutrino_z,2) ) );	    
	        h_vertexfitter_resolution_neutrino_not_merged->Fill( sqrt( pow(fnu_reco_fitter_x - fneutrino_x,2) + pow(fnu_reco_fitter_y - fneutrino_y,2) + pow(fnu_reco_fitter_z - fneutrino_z,2) ) );	    
	        h_vertexfitter_chi2ndf_neutrino_not_merged->Fill( fnu_reco_fitter_chi2ndf  );	    
		
    		for (unsigned j=0; j < fpdg.size(); j++) {
	    	if ( !reco_muon ) break; //select events with a reco muon
	    	float res = min ( sqrt( pow(freco_vertex_x[j] - fstart_x[j],2) + pow(freco_vertex_y[j] - fstart_y[j],2) + pow(freco_vertex_z[j] - fstart_z[j],2)) ,
					sqrt( pow(freco_vertex_x[j] - fend_x[j],2) + pow(freco_vertex_y[j] - fend_y[j],2) + pow(freco_vertex_z[j] - fend_z[j],2)) );
	    	float res_fitter = min ( sqrt( pow(freco_vertexfitter_x[j] - fstart_x[j],2) + pow(freco_vertexfitter_y[j] - fstart_y[j],2) + pow(freco_vertexfitter_z[j] - fstart_z[j],2)) ,
					sqrt( pow(freco_vertexfitter_x[j] - fend_x[j],2) + pow(freco_vertexfitter_y[j] - fend_y[j],2) + pow(freco_vertexfitter_z[j] - fend_z[j],2)) );
		
		if ( fpdg[j]==13 ) {
			h_vertex_resolution_muon_not_merged->Fill( res );
			h_vertexfitter_resolution_muon_not_merged->Fill( res_fitter );
			h_vertexfitter_chi2ndf_muon_not_merged->Fill( freco_vertexfitter_chi2ndf[j] );
		}
	    	if ( fpdg[j]!=2212 && fp0[j] <= 0.2 && !fis_tracked[j]) continue; //watch only reco protons with p>0.2Gev/c
   		hproton_pos_res_goodprotons->Fill(  sqrt( pow( freco_startx[j]- fstart_x[j] ,2) + pow( freco_starty[j]- fstart_y[j] ,2) + pow( freco_startz[j]- fstart_z[j] ,2) ) ); 
	    	h_vertex_resolution_proton_not_merged->Fill( res );
	    	h_vertexfitter_resolution_proton_not_merged->Fill( res_fitter );
	    	h_vertexfitter_chi2ndf_proton_not_merged->Fill( freco_vertexfitter_chi2ndf[j] );

   		if (lowmomentum_p) hproton_pos_res_lowprotons->Fill(  sqrt( pow( freco_startx[j]- fstart_x[j] ,2) + pow( freco_starty[j]- fstart_y[j] ,2) + pow( freco_startz[j]- fstart_z[j] ,2) ) ); 
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
    		if ( !(abs(fnu_reco_x+1)<DBL_EPSILON && abs(fnu_reco_y+1)<DBL_EPSILON && abs(fnu_reco_z+1)<DBL_EPSILON) ) {
	        h_vertex_resolution_neutrino_merged->Fill( sqrt( pow(fnu_reco_x - fneutrino_x,2) + pow(fnu_reco_y - fneutrino_y,2) + pow(fnu_reco_z - fneutrino_z,2) ) );	    
	        h_vertexfitter_resolution_neutrino_merged->Fill( sqrt( pow(fnu_reco_fitter_x - fneutrino_x,2) + pow(fnu_reco_fitter_y - fneutrino_y,2) + pow(fnu_reco_fitter_z - fneutrino_z,2) ) );	    
	        h_vertexfitter_chi2ndf_neutrino_merged->Fill( fnu_reco_fitter_chi2ndf );	    
		}
    		for (unsigned j=0; j<fpdg.size(); j++) {
	    	if ( !reco_muon ) break; //select events with a reco muon
	    	float res = min ( sqrt( pow(freco_vertex_x[j] - fstart_x[j],2) + pow(freco_vertex_y[j] - fstart_y[j],2) + pow(freco_vertex_z[j] - fstart_z[j],2)) ,
					sqrt( pow(freco_vertex_x[j] - fend_x[j],2) + pow(freco_vertex_y[j] - fend_y[j],2) + pow(freco_vertex_z[j] - fend_z[j],2)) );
	    	float res_fitter = min ( sqrt( pow(freco_vertexfitter_x[j] - fstart_x[j],2) + pow(freco_vertexfitter_y[j] - fstart_y[j],2) + pow(freco_vertexfitter_z[j] - fstart_z[j],2)) ,
					sqrt( pow(freco_vertexfitter_x[j] - fend_x[j],2) + pow(freco_vertexfitter_y[j] - fend_y[j],2) + pow(freco_vertexfitter_z[j] - fend_z[j],2)) );
		
		if ( fpdg[j]==13 ) {
			h_vertex_resolution_muon_merged->Fill( res );
			h_vertexfitter_resolution_muon_merged->Fill( res_fitter );
			h_vertexfitter_chi2ndf_muon_merged->Fill( freco_vertexfitter_chi2ndf[j] );
		}
	    	if ( fpdg[j]!=2212 && fp0[j] <= 0.2 && !fis_tracked[j]) continue; //watch only reco protons with p>0.2Gev/c
   		hproton_pos_res_badprotons->Fill(  sqrt( pow( freco_startx[j]- fstart_x[j] ,2) + pow( freco_starty[j]- fstart_y[j] ,2) + pow( freco_startz[j]- fstart_z[j] ,2) ) ); 
	    	h_vertex_resolution_proton_merged->Fill( res );
	    	h_vertexfitter_resolution_proton_merged->Fill( res_fitter );
	    	h_vertexfitter_chi2ndf_proton_merged->Fill( freco_vertexfitter_chi2ndf[j] );
   		if (lowmomentum_p) hproton_pos_res_lowprotons->Fill(  sqrt( pow( freco_startx[j]- fstart_x[j] ,2) + pow( freco_starty[j]- fstart_y[j] ,2) + pow( freco_startz[j]- fstart_z[j] ,2) ) ); 
		}
	}

	if ( reco_muon && count_tracked == 0 && count_not_tracked== 0 && !is_pion && lowmomentum_p ) { //no protons tracked but low energy ones
		is_lowmomentum_p = true;
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


	return muon_pos;

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
	fis_tracked.push_back(false);
	fis_shower_matched.push_back(false);
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
	freco_trackid.push_back(-99);
	freco_track_hits.push_back(0);
	freco_track_collection_hits.push_back(0);
	freco_track_collection_charge.push_back(0);
	freco_mcp_hits.push_back(-1);
	freco_mcp_collection_hits.push_back(-1);
	freco_mcp_collection_charge.push_back(-1);
	freco_vertex_x.push_back(-1);
	freco_vertex_y.push_back(-1);
	freco_vertex_z.push_back(-1);
	freco_vertexfitter_x.push_back(-1);
	freco_vertexfitter_y.push_back(-1);
	freco_vertexfitter_z.push_back(-1);
	freco_vertexfitter_chi2ndf.push_back(-1);
	freco_vertexfitter_chi2.push_back(-1);
	fnot_clustered.push_back(0);
	fclustered.push_back(0);
	fclustered_matched.push_back(0);
	fclustered_mismatched.push_back(0);
	fnot_clustered_charge.push_back( std::vector<float>() );
	fclustered_charge.push_back( std::vector<float>() );
	fclustered_matched_charge.push_back( std::vector<float>() );
	fclustered_mismatched_charge.push_back( std::vector<float>() );
	fhit_mismatch_pdg.push_back( std::vector<int>() );

}

void recohelper::RecoBenchmarker::endJob()
{

  //scale some histos	
  hproton_multi_all->Scale(1./n_events);
  hproton_multi_above20MeV->Scale(1./n_events);
  hproton_multi_below20MeV->Scale(1./n_events);
 
  hproton_merged_not_merged->Scale(1./tot_n_protons);

  h_fraction_pdgs_not_tracked_not_clustered->SetBinContent( 2 , h_fraction_pdgs_not_tracked_not_clustered->GetBinContent(2)/float(n_all_protons) ); //all protons
  h_fraction_pdgs_not_tracked_not_clustered->SetBinError( 2 , h_fraction_pdgs_not_tracked_not_clustered->GetBinContent(2)/sqrt(h_fraction_pdgs_not_tracked_not_clustered->GetEntries()) ); //all protons
  h_fraction_pdgs_not_tracked_not_clustered->SetBinContent( 3 , h_fraction_pdgs_not_tracked_not_clustered->GetBinContent(3)/float(tot_n_protons) ); //>20MeV
  h_fraction_pdgs_not_tracked_not_clustered->SetBinError( 3 , h_fraction_pdgs_not_tracked_not_clustered->GetBinContent(3)/sqrt(h_fraction_pdgs_not_tracked_not_clustered->GetEntries()) ); //>20MeV
  h_fraction_pdgs_not_tracked_not_clustered->SetBinContent( 4 , h_fraction_pdgs_not_tracked_not_clustered->GetBinContent(4)/float(low_protons) ); //<20MeV
  h_fraction_pdgs_not_tracked_not_clustered->SetBinError( 4 , h_fraction_pdgs_not_tracked_not_clustered->GetBinContent(4)/sqrt(h_fraction_pdgs_not_tracked_not_clustered->GetEntries()) ); //<20MeV
   
  /*h_muon_clustering_prob_good_protons->Scale( 1./h_muon_clustering_prob_good_protons->Integral() );
  h_proton_clustering_prob_good_protons->Scale( 1./h_proton_clustering_prob_good_protons->Integral() );
  h_muon_clustering_prob_bad_protons->Scale( 1./h_muon_clustering_prob_bad_protons->Integral() );
  h_proton_clustering_prob_bad_protons->Scale( 1./h_proton_clustering_prob_bad_protons->Integral() );*/

  TFile& file = tfs->file();
  file.cd();

  //--------------------------
  // Produce eff. plots
  //--------------------------
   hmuon_spectrum_eff = (TH1D*) hmuon_spectrum->Clone("muon_spectrum_eff");
   hmuon_spectrum_eff->Divide(hmuon_spectrum_all);
   hmuon_spectrum_eff->Write();
   
   hmuon_length_eff = (TH1D*) hmuon_length->Clone("muon_length_eff");
   hmuon_length_eff->Divide(hmuon_length_all);
   hmuon_length_eff->Write();
   
   hproton_kinE_eff = (TH1D*) hproton_kinE->Clone("proton_kinE_eff");
   hproton_kinE_eff->Divide(hproton_kinE_all);
   hproton_kinE_eff->Write();
   
   hproton_p_eff = (TH1D*) hproton_p->Clone("proton_p_eff");
   hproton_p_eff->Divide(hproton_p_all);
   hproton_p_eff->Write();
   
   hproton_l_eff = (TH1D*) hproton_l->Clone("proton_l_eff");
   hproton_l_eff->Divide(hproton_l_all);
   hproton_l_eff->Write();
   
   hproton_l_eff = (TH1D*) hproton_l->Clone("proton_l_eff");
   hproton_l_eff->Divide(hproton_l_all);
   hproton_l_eff->Write();

   hproton_theta_mu_eff = (TH1D*) h_theta_mu_tracked->Clone("proton_theta_mu_eff");
   hproton_theta_mu_eff->Divide(hproton_theta_mu);
   hproton_theta_mu_eff->Write();
/*
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
*/
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
  fis_shower_matched.clear();
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
   freco_trackid.clear();
    flength_reco.clear();
   fpurity.clear();
   fcompleteness.clear();
   fnhits.clear();
   fkinE.clear();
  
   //hits analysis
   freco_track_hits.clear();
   freco_track_collection_hits.clear();
   freco_track_collection_charge.clear();
   freco_mcp_hits.clear();
   freco_mcp_collection_hits.clear();
   freco_mcp_collection_charge.clear();
   fnot_clustered.clear();
   fclustered.clear();
   fclustered_matched.clear();
   fclustered_mismatched.clear();

   //delicate vectors
   for (unsigned ii=0; ii<fnot_clustered_charge.size(); ii++) {
   fnot_clustered_charge[ii].clear();
   fclustered_charge[ii].clear();
   fclustered_matched_charge[ii].clear();
   fclustered_mismatched_charge[ii].clear();
   fhit_mismatch_pdg[ii].clear();
   }
   fnot_clustered_charge.clear();
   fclustered_charge.clear();
   fclustered_matched_charge.clear();
   fclustered_mismatched_charge.clear();
   fhit_mismatch_pdg.clear();



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
