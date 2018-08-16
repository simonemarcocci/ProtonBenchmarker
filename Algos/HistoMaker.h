#ifndef HistoMaker_h
#define HistoMaker_h

// larsoft includes
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// ROOT includes
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"

// cc includes
#include <numeric>
#include <vector>
#include <functional>
#include <stdexcept>
#include <cmath>
#include <algorithm>

//ProtonBenchmarker includes
#include "uboone/ProtonBenchmarker/Datatypes/StoredEvent.h"

namespace reco_histo{

  class HistoMaker{

    public:
	
	    HistoMaker( bool isdata = false, bool isverbose = false ){
		    is_init = false;
		    is_init_hit = false;
    	 	    low_edge = 350;
    		    high_edge = 1500;
    		    length_cut = 20; //8cm
		    tot_n_protons = 0;
		    n_all_protons = 0;
		    low_protons = 0;
		    count_tracked = 0;
		    count_not_tracked = 0;
		    is_lowmomentum_p = false;
		    fmuon_pos = -1;
		    isVerbose = isverbose;
		    isData = isdata;
	    }
	    ~HistoMaker(){}

	    void Init( art::ServiceHandle< art::TFileService > ); //init histos for tracks/showers/pfp
	    void Init_Hit(  art::TFileDirectory ); //init histos for hit studies

	    void Fill_Truth_Histos( StoredEvent* );
	    void Fill_Hit_Histos( StoredEvent* );
	    void Fill_Analysis_Histos( StoredEvent*, int&, bool );
	    void Fill_Analysis_Histos_Data( StoredEvent*, int&, bool );
	    void FillCumulativeHistograms();
	    void ScalePlots(int);

   private:

   bool isData;
   bool isVerbose;
   bool is_init;
   bool is_init_hit;
   //a few parameters for histo plotting
   long low_edge;
   long high_edge;
   long length_cut;
   
   //service variables
   long tot_n_protons; //protons > 20MeV
   long n_all_protons; //all protons
   long low_protons; // protons <20MeV
   int  count_tracked;
   int count_not_tracked;
   bool is_lowmomentum_p;	
   int fmuon_pos;

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
   TH1D* hproton_nhits_all;
   TH1D* hproton_nhits_CP_all;
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

  };

}

#endif
