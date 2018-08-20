////////////////////////////////////////////////////////////////////////
// Class:       StoredEvent
// File:        StoredEvent_module.cc
//
// Generated at Wed Aug  8 16:15:28 2018 by Simone Marcocci using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#ifndef StoredEvent_cxx
#define StoredEvent_cxx

#include "StoredEvent.h"


StoredEvent::StoredEvent() {
	AllocateRecoVectors();
	ClearVectors();
	fcount_proton_showers = 0;
}

StoredEvent::~StoredEvent() {
	ClearVectors();
}

void StoredEvent::AllocateRecoVectors() {
	fis_tracked.push_back(false);
	fis_shower_matched.push_back(false);
	fshower_pdg.push_back(-1);
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
	fis_spacepoint.push_back( std::vector<bool>() );
	fnot_clustered_charge.push_back( std::vector<float>() );
	fnot_clustered_hit_index.push_back( std::vector<int>() );
	fclustered_hit_index.push_back( std::vector<int>() );
	fnot_clustered_tracked_charge.push_back( std::vector<float>() );
	fnot_clustered_not_tracked_charge.push_back( std::vector<float>() );
	fspacepoint_xyz.push_back( std::vector< std::vector<double> >() );
	fclustered_charge.push_back( std::vector<float>() );
	fclustered_matched_charge.push_back( std::vector<float>() );
	fclustered_mismatched_charge.push_back( std::vector<float>() );
	fhit_mismatch_pdg.push_back( std::vector<int>() );

}


void StoredEvent::ClearVectors(){
    fmuon_dqdx.clear();
    fmuon_dedx.clear();
    fmuon_residual.clear();
    flength.clear();
    fshower_pdg.clear();
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
   for (unsigned ii=0; ii<fis_spacepoint.size(); ii++){
   	fis_spacepoint[ii].clear();
   	fspacepoint_xyz[ii].clear();
   }
   for (unsigned ii=0; ii<fnot_clustered_charge.size(); ii++) {
   fnot_clustered_hit_index[ii].clear();
   fnot_clustered_charge[ii].clear();
   fnot_clustered_tracked_charge[ii].clear();
   fnot_clustered_not_tracked_charge[ii].clear();
   }
   for (unsigned ii=0; ii<fclustered_charge.size(); ii++) {
   fclustered_hit_index[ii].clear();
   fclustered_charge[ii].clear();
   fclustered_matched_charge[ii].clear();
   fclustered_mismatched_charge[ii].clear();
   fhit_mismatch_pdg[ii].clear();
   }
   fnot_clustered_hit_index.clear();
   fclustered_hit_index.clear();
   fnot_clustered_charge.clear();
   fnot_clustered_tracked_charge.clear();
   fnot_clustered_not_tracked_charge.clear();
   fspacepoint_xyz.clear();
   fis_spacepoint.clear();
   fclustered_charge.clear();
   fclustered_matched_charge.clear();
   fclustered_mismatched_charge.clear();
   fhit_mismatch_pdg.clear();
}


#endif
