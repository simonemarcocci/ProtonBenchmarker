////////////////////////////////////////////////////////////////////////
// Class:       StoredEvent
// File:        StoredEvent.h/cc
//
// Generated at Wed Aug  8 16:15:28 2018 by Simone Marcocci using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#ifndef StoredEvent_h
#define StoredEvent_h

#include <vector>
#include <map>
#include <string>
#include "TVector3.h"


class StoredEvent {

public:
   
    StoredEvent();
       ~StoredEvent();

    int fEvent;
    int fRun;
    int fSubRun;
    int fccnc;
    int finteraction;

    float fneutrino_x;
    float fneutrino_y;
    float fneutrino_z;

    int fcount_proton_showers;

    //info on the MC particle
    std::map<int,int> fparticle_count;
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
   std::vector<int> fshower_pdg;
   std::vector<bool> fmatch_multiplicity;
   //std::vector<bool> fis_mismatched; //it says if the MC truth assignment in different planes is different (possible hint for wrong or problematic backtracking)
   std::vector<float> fcostheta_muon_reco; //MCtruth info
   std::vector<int> fpdg_reco;
   std::vector<float> flength_reco;
   std::vector<float> freco_costheta_muon;
   std::vector<TVector3> freco_start_direction;
   std::vector<float> freco_momentum_mcs; //(GeV) MCS
   std::vector<float> freco_momentum_mcs_llhd; //(GeV) MCS LLHD
   std::vector<float> freco_momentum_range; //MeV
   std::vector<float> fpurity;
   std::vector<float> fcompleteness;
   std::vector<double> fnhits;
   std::vector<bool> freco_ismuon;
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
   std::vector< std::vector<bool> > fis_spacepoint;
   std::vector< std::vector< std::vector<double> > > fspacepoint_xyz;
   std::vector< std::vector<int> > fclustered_hit_index;
   std::vector< std::vector<int> > fnot_clustered_hit_index;
   std::vector< std::vector<float>> fnot_clustered_charge;
   std::vector< std::vector<float>> fnot_clustered_tracked_charge;
   std::vector< std::vector<float>> fnot_clustered_not_tracked_charge;
   std::vector<float> fclustered;
   std::vector< std::vector<float>> fclustered_charge;
   std::vector<float> fclustered_matched;
   std::vector< std::vector<float>> fclustered_matched_charge;
   std::vector<float> fclustered_mismatched;
   std::vector< std::vector<float>> fclustered_mismatched_charge;
   std::vector< std::vector<int>> fhit_mismatch_pdg;
  
   void    AllocateRecoVectors();
   void    ClearVectors();

};
#endif
