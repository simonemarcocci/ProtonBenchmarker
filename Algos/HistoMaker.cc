#ifndef HistoMaker_cxx
#define HistoMaker_cxx

#include "HistoMaker.h"

namespace reco_histo{

void HistoMaker::Init( art::ServiceHandle< art::TFileService > tfs ) {
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
   hproton_nhits_all = tfs->make<TH1D>("proton_nhits_all","nhits for ALL protons; nhits",1000,0,1000); //reco efficiency protons vs kin E
   hproton_nhits_CP_all = tfs->make<TH1D>("proton_nhits_CP_all","Collection Plane nhits for ALL protons; nhits",1000,0,1000); //reco efficiency protons vs kin E
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
   
   is_init = true;
}//end general init

void HistoMaker::Init_Hit(  art::TFileDirectory hits_dir ) {
   
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

   is_init_hit = true;
} //end hits init


void HistoMaker::Fill_Truth_Histos( StoredEvent* event_store ) {
  if (!is_init) {
	  std::cout << "Histo Maker NOT initialized" << std::endl;
	  throw cet::exception("Configuration");
  }
  bool lleading = false;
  int count_protons = 0;
  int count_protons_above = 0;
  int count_protons_below = 0;
  for ( unsigned ii = 0; ii < event_store->fpdg.size(); ii++ ) {
	if ( event_store->fpdg[ii] == 2212 && event_store->fis_leading[ii] ) {
		hproton_leading_kinE->Fill( event_store->fkinE[ii] * 1000. ); //bring it to MeV
		if (lleading==true) std::cout << "ERROR!!!! CAN'T BE ALREADY TRUE" << std::endl;
		lleading = true;
	}
	if ( event_store->fpdg[ii] == 2212 ) count_protons++;
	if ( event_store->fpdg[ii] == 2212 && event_store->fkinE[ii] >= 0.02 ) count_protons_above++;
	if ( event_store->fpdg[ii] == 2212 && event_store->fkinE[ii] < 0.02 ) count_protons_below++;

    	if ( event_store->fis_shower_matched[ii] == true ) {
    		h_shower_pdg->Fill( event_store->fshower_pdg[ii] );
    			if ( event_store->fpdg[ii] == 2212 ) { //protons
	    			h_shower_proton_l->Fill( event_store->flength[ii] );
	    			h_shower_proton_kinE->Fill( event_store->fkinE[ii] );
	    			h_shower_proton_nhits->Fill( event_store->fnhits[ii] );
	    			h_shower_proton_costheta_muon->Fill( event_store->fcostheta_muon[ii] );

    			}
 	}//is shower matched
        
  }
	hproton_multi_all->Fill( count_protons );
	hproton_multi_above20MeV->Fill( count_protons_above );
	hproton_multi_below20MeV->Fill( count_protons_below );
  	h_n_proton_showers->Fill( event_store->fcount_proton_showers );
}


void HistoMaker::Fill_Analysis_Histos( StoredEvent* event_store, int & muon_pos ) {
  
  if (!is_init) {
	  std::cout << "Histo Maker NOT initialized" << std::endl;
	  throw cet::exception("Configuration");
  }
  
    muon_pos = -1;
    fmuon_pos = muon_pos;
    if (event_store->fccnc!=0 ) {
	    return; //comment if you're running over a std sample with NC interactions and you're interested in them
    }
    //check that the muon is reco
    bool reco_muon = false;
    bool is_pion=false;
    bool lowmomentum_p = false;
    for (unsigned i = 0; i < event_store->fpdg.size(); i++) {
	    if ( event_store->fpdg[i] == 13 ) { //ismuon
	    hmuon_length_all->Fill( event_store->flength[i] );
	    hmuon_spectrum_all->Fill( event_store->fkinE[i] );
	    if ( event_store->fis_tracked[i] &&  event_store->fmuon_dqdx.size() != 0 ) { //want that the muon is well matched and w/ dqdx info
		muon_pos = i; 
	    	hmuon_length->Fill( event_store->flength[i] );
	    	hmuon_spectrum->Fill( event_store->fkinE[i] );
		reco_muon = true;
		hmuon_pos_res->Fill( sqrt( pow( event_store->freco_startx[i]- event_store->fstart_x[i] ,2) + pow( event_store->freco_starty[i]- event_store->fstart_y[i] ,2) + pow( event_store->freco_startz[i]- event_store->fstart_z[i] ,2) ));
		float res = std::min ( sqrt( pow(event_store->freco_vertex_x[i] - event_store->fstart_x[i],2) + pow(event_store->freco_vertex_y[i] - event_store->fstart_y[i],2) + pow(event_store->freco_vertex_z[i] - event_store->fstart_z[i],2)) ,
					sqrt( pow(event_store->freco_vertex_x[i] - event_store->fend_x[i],2) + pow(event_store->freco_vertex_y[i] - event_store->fend_y[i],2) + pow(event_store->freco_vertex_z[i] - event_store->fend_z[i],2)) );
		h_vertex_resolution_muon->Fill( res );
		float res_fitter = std::min ( sqrt( pow(event_store->freco_vertexfitter_x[i] - event_store->fstart_x[i],2) + pow(event_store->freco_vertexfitter_y[i] - event_store->fstart_y[i],2) + pow(event_store->freco_vertexfitter_z[i] - event_store->fstart_z[i],2)) ,
					sqrt( pow(event_store->freco_vertexfitter_x[i] - event_store->fend_x[i],2) + pow(event_store->freco_vertexfitter_y[i] - event_store->fend_y[i],2) + pow(event_store->freco_vertexfitter_z[i] - event_store->fend_z[i],2)) );
		h_vertexfitter_resolution_muon->Fill( res_fitter );
		h_vertexfitter_chi2ndf_muon->Fill( event_store->freco_vertexfitter_chi2ndf[i] );
		if ( sqrt( pow( event_store->freco_startx[i]- event_store->fstart_x[i] ,2) + pow( event_store->freco_starty[i]- event_store->fstart_y[i] ,2) + pow( event_store->freco_startz[i]- event_store->fstart_z[i] ,2) ) > 50 ) {
			reco_muon = false; //skip these events, might have direction flipped
	    	}
    		}
	    } //record info on the muon
	    if ( event_store->fpdg[i] == 211 || event_store->fpdg[i] == -211 || event_store->fpdg[i] == 111 ) is_pion=true; //record if there are any pions
	    if (  event_store->fpdg[i] == 2212 && event_store->fp0[i] <= 0.2 ) lowmomentum_p = true;
	    if (  event_store->fpdg[i] == 2212 && event_store->fp0[i] >= 0.2 ) hproton_theta_mu->Fill(event_store->fcostheta_muon[i]);
    }
    
    if ( reco_muon == false ) return; //select events with a reco muon
    
    //check if the neutrino reco'ed vertex is there 
    if ( !(std::abs(event_store->fnu_reco_x+1)<DBL_EPSILON && std::abs(event_store->fnu_reco_y+1)<DBL_EPSILON && std::abs(event_store->fnu_reco_z+1)<DBL_EPSILON) ) 
	h_vertex_resolution_neutrino->Fill( sqrt( pow(event_store->fnu_reco_x - event_store->fneutrino_x,2) + pow(event_store->fnu_reco_y - event_store->fneutrino_y,2) + pow(event_store->fnu_reco_z - event_store->fneutrino_z,2) ) );	    
    if ( !(std::abs(event_store->fnu_reco_fitter_x+1)<DBL_EPSILON && std::abs(event_store->fnu_reco_fitter_y+1)<DBL_EPSILON && std::abs(event_store->fnu_reco_fitter_z+1)<DBL_EPSILON) ) {
	h_vertexfitter_resolution_neutrino->Fill( sqrt( pow(event_store->fnu_reco_fitter_x - event_store->fneutrino_x,2) + pow(event_store->fnu_reco_fitter_y - event_store->fneutrino_y,2) + pow(event_store->fnu_reco_fitter_z - event_store->fneutrino_z,2) ) );	    
	h_vertexfitter_chi2ndf_neutrino->Fill( event_store->fnu_reco_fitter_chi2ndf );
    }

    count_not_tracked = 0;
    count_tracked = 0;
    is_lowmomentum_p = false;
    for (unsigned j=0; j<event_store->fpdg.size(); j++) {
	    if ( !reco_muon ) break; //select events with a reco muon
	    if ( event_store->fpdg[j]!=2212 ) continue; //watch only protons
	 
	    n_all_protons++;
	    hproton_p_all->Fill( event_store->fp0[j] );
	    hproton_l_all->Fill( event_store->flength[j] );
	    hproton_kinE_all->Fill( event_store->fkinE[j] );
   	    hproton_nhits_all->Fill( event_store->fnhits[j] );
   	    hproton_nhits_CP_all->Fill( event_store->freco_mcp_collection_hits[j] );
	   	
	    if ( std::abs( event_store->fcostheta_muon[j] ) > sqrt(3)/2. ) {
	    hproton_l_all_angle1->Fill( event_store->flength[j] );
	    hproton_kinE_all_angle1->Fill( event_store->fkinE[j] );
	    hproton_nhits_all_angle1->Fill( event_store->freco_mcp_collection_hits[j] ); //FIXME
	    } else if ( std::abs( event_store->fcostheta_muon[j] ) >= 0.5 && std::abs( event_store->fcostheta_muon[j] ) <= sqrt(3)/2. ) { 
	    hproton_l_all_angle2->Fill( event_store->flength[j] );
	    hproton_kinE_all_angle2->Fill( event_store->fkinE[j] );
	    hproton_nhits_all_angle2->Fill( event_store->freco_mcp_collection_hits[j] ); //FIXME
	    } else if ( std::abs( event_store->fcostheta_muon[j] ) >= 0. && std::abs( event_store->fcostheta_muon[j] ) < 0.5)  {
	    hproton_l_all_angle3->Fill( event_store->flength[j] );
	    hproton_kinE_all_angle3->Fill( event_store->fkinE[j] );
	    hproton_nhits_all_angle3->Fill( event_store->freco_mcp_collection_hits[j] ); //FIXME
	    } else {
	    throw cet::exception("LogicError") << "Crazy angle!";
	    }

	    h_theta_mu_length_all->Fill( event_store->fcostheta_muon[j], event_store->flength[j]);
	    float res = std::min ( sqrt( pow(event_store->freco_vertex_x[j] - event_store->fstart_x[j],2) + pow(event_store->freco_vertex_y[j] - event_store->fstart_y[j],2) + pow(event_store->freco_vertex_z[j] - event_store->fstart_z[j],2)) ,
					sqrt( pow(event_store->freco_vertex_x[j] - event_store->fend_x[j],2) + pow(event_store->freco_vertex_y[j] - event_store->fend_y[j],2) + pow(event_store->freco_vertex_z[j] - event_store->fend_z[j],2)) );
	    h_vertex_resolution_proton->Fill( res );
	    float res_fitter = std::min ( sqrt( pow(event_store->freco_vertexfitter_x[j] - event_store->fstart_x[j],2) + pow(event_store->freco_vertexfitter_y[j] - event_store->fstart_y[j],2) + pow(event_store->freco_vertexfitter_z[j] - event_store->fstart_z[j],2)) ,
					sqrt( pow(event_store->freco_vertexfitter_x[j] - event_store->fend_x[j],2) + pow(event_store->freco_vertexfitter_y[j] - event_store->fend_y[j],2) + pow(event_store->freco_vertexfitter_z[j] - event_store->fend_z[j],2)) );
	    h_vertexfitter_resolution_proton->Fill( res_fitter );
	    h_vertexfitter_chi2ndf_proton->Fill( event_store->freco_vertexfitter_chi2ndf[j] );

	    if ( event_store->fis_tracked[j] ) { //efficiency plots for all protons (but don't divide yet)
	    	hproton_p->Fill( event_store->fp0[j] );
	    	hproton_l->Fill( event_store->flength[j] );
	    	hproton_kinE->Fill( event_store->fkinE[j] );
	        h_theta_mu_length->Fill( event_store->fcostheta_muon[j], event_store->flength[j]);
   		hproton_nhits->Fill( event_store->fnhits[j] );
   		hproton_nhits_CP->Fill( event_store->freco_mcp_collection_hits[j] );
   		hproton_nhits_theta_mu->Fill ( event_store->fnhits[j], acos(event_store->fcostheta_muon[j]) ) ;
   		hproton_nhits_CP_theta_mu->Fill( event_store->freco_mcp_collection_hits[j], acos(event_store->fcostheta_muon[j]) );
	    
		if ( std::abs( event_store->fcostheta_muon[j] ) > sqrt(3)/2. ) { //theta < 30degrees
	    		hproton_l_tracked_angle1->Fill( event_store->flength[j] );
	    		hproton_kinE_tracked_angle1->Fill( event_store->fkinE[j] );
	    		hproton_nhits_tracked_angle1->Fill( event_store->freco_mcp_collection_hits[j] ); //FIXME
	    	} else if ( std::abs( event_store->fcostheta_muon[j] ) >= 0.5 && std::abs( event_store->fcostheta_muon[j] ) <= sqrt(3)/2. ) { 
	    		hproton_l_tracked_angle2->Fill( event_store->flength[j] );
	    		hproton_kinE_tracked_angle2->Fill( event_store->fkinE[j] );
	    		hproton_nhits_tracked_angle2->Fill( event_store->freco_mcp_collection_hits[j] ); //FIXME
	    	} else {
	    		hproton_l_tracked_angle3->Fill( event_store->flength[j] );
	    		hproton_kinE_tracked_angle3->Fill( event_store->fkinE[j] );
	    		hproton_nhits_tracked_angle3->Fill( event_store->freco_mcp_collection_hits[j] ); //FIXME
	  	}
	    }


    //now study all the protons with momentum > 200MeV
    if (event_store->fp0[j] > 0.2) {
	 tot_n_protons++;
  	 h_theta_mu->Fill( event_store->fcostheta_muon[j] );	
	 if ( !event_store->fis_tracked[j] ) {
		count_not_tracked++;
		h_pmu_end_not_tracked->Fill( event_store->flength[j] * sqrt( 1 - pow( event_store->fcostheta_muon[j] ,2 ) ) ) ;
		h_theta_mu_not_tracked->Fill ( event_store->fcostheta_muon[j] ) ;
		}
    	else if ( event_store->fis_tracked[j] ) {
		count_tracked++;
		h_pmu_end_tracked->Fill( event_store->flength[j] * sqrt( 1 - pow( event_store->fcostheta_muon[j] ,2 ) ) ) ;
		h_theta_mu_tracked->Fill ( event_store->fcostheta_muon[j] ) ;
		hmuon_proton_tracked->Fill ( sqrt( pow( event_store->freco_startx[j]- event_store->freco_startx[muon_pos] ,2) + pow( event_store->freco_starty[j]- event_store->freco_starty[muon_pos] ,2) + pow( event_store->freco_startz[j]- event_store->freco_startz[muon_pos] ,2) ) );
		hproton_pos_res->Fill ( sqrt( pow( event_store->freco_startx[j]- event_store->fstart_x[j] ,2) + pow( event_store->freco_starty[j]- event_store->fstart_y[j] ,2) + pow( event_store->freco_startz[j]- event_store->fstart_z[j] ,2) ) );
		}
   	} else {
		low_protons++;
	}
}
	
	hproton_merged_not_merged->Fill(0., double(count_tracked) );
	hproton_merged_not_merged->Fill(1., double(count_not_tracked) );
        
 	//neutrino vertex resolution vs proton recon performance:
    	if ( !(std::abs(event_store->fnu_reco_x+1)<DBL_EPSILON && std::abs(event_store->fnu_reco_y+1)<DBL_EPSILON && std::abs(event_store->fnu_reco_z+1)<DBL_EPSILON) ) {
	int tot_protons = count_not_tracked + count_tracked ;
	h_vertex_resolution_vs_not_tracked_above20MeV->Fill( sqrt( pow(event_store->fnu_reco_x - event_store->fneutrino_x,2) + pow(event_store->fnu_reco_y - event_store->fneutrino_y,2) + pow(event_store->fnu_reco_z - event_store->fneutrino_z,2) ), double(count_not_tracked)/tot_protons );
	h_vertex_resolution_vs_not_tracked->Fill( sqrt( pow(event_store->fnu_reco_x - event_store->fneutrino_x,2) + pow(event_store->fnu_reco_y - event_store->fneutrino_y,2) + pow(event_store->fnu_reco_z - event_store->fneutrino_z,2) ), double(count_not_tracked + low_protons)/ (tot_protons+low_protons) );
	h_vertex_resolution_vs_not_tracked_below20MeV->Fill( sqrt( pow(event_store->fnu_reco_x - event_store->fneutrino_x,2) + pow(event_store->fnu_reco_y - event_store->fneutrino_y,2) + pow(event_store->fnu_reco_z - event_store->fneutrino_z,2) ), double(low_protons)/ (tot_protons+low_protons) );
	h_vertexfitter_resolution_vs_not_tracked_above20MeV->Fill( sqrt( pow(event_store->fnu_reco_fitter_x - event_store->fneutrino_x,2) + pow(event_store->fnu_reco_fitter_y - event_store->fneutrino_y,2) + pow(event_store->fnu_reco_fitter_z - event_store->fneutrino_z,2) ), double(count_not_tracked)/tot_protons );
	h_vertexfitter_resolution_vs_not_tracked->Fill( sqrt( pow(event_store->fnu_reco_fitter_x - event_store->fneutrino_x,2) + pow(event_store->fnu_reco_fitter_y - event_store->fneutrino_y,2) + pow(event_store->fnu_reco_fitter_z - event_store->fneutrino_z,2) ), double(count_not_tracked + low_protons)/ (tot_protons+low_protons) );
	h_vertexfitter_resolution_vs_not_tracked_below20MeV->Fill( sqrt( pow(event_store->fnu_reco_fitter_x - event_store->fneutrino_x,2) + pow(event_store->fnu_reco_fitter_y - event_store->fneutrino_y,2) + pow(event_store->fnu_reco_fitter_z - event_store->fneutrino_z,2) ), double(low_protons)/ (tot_protons+low_protons) );
	}

	if ( event_store->fmuon_residual.size() != event_store->fmuon_dqdx.size()) std::cout << "ERROR on calorimetry vector sizes!!!" << std::endl;
	if ( count_not_tracked == 0 && count_tracked > 0 ) { //all protons are tracked 
		for (unsigned jj=1; jj<event_store->fmuon_dqdx.size()-1; jj++) {
				if ( event_store->fmuon_range - event_store->fmuon_residual[jj] < 8 ) { //look at 8cm only
				h_dqdx_1d_not_merged->Fill(event_store->fmuon_dqdx[jj]);
				}

				h_dqdx_not_merged->Fill( event_store->fmuon_range - event_store->fmuon_residual[jj], event_store->fmuon_dqdx[jj] );
				h_dqdx_not_merged_service->Fill( event_store->fmuon_range -  event_store->fmuon_residual[jj], event_store->fmuon_dqdx[jj] );
		}
		TH1D* h1 = NULL;
		h1 = h_dqdx_not_merged_service->ProjectionY("h",1, length_cut);
		if (h1) htail_to_tot_not_merged->Fill( h1->Integral(low_edge,high_edge) ); 
		h_dqdx_not_merged_service->Reset();

		hmuon_pos_res_goodprotons->Fill( sqrt( pow( event_store->freco_startx[muon_pos]- event_store->fstart_x[muon_pos] ,2) + pow( event_store->freco_starty[muon_pos]- event_store->fstart_y[muon_pos] ,2) + pow( event_store->freco_startz[muon_pos]- event_store->fstart_z[muon_pos] ,2) ) );
    		if ( !(std::abs(event_store->fnu_reco_x+1)<DBL_EPSILON && std::abs(event_store->fnu_reco_y+1)<DBL_EPSILON && std::abs(event_store->fnu_reco_z+1)<DBL_EPSILON) ) 
	        h_vertex_resolution_neutrino_not_merged->Fill( sqrt( pow(event_store->fnu_reco_x - event_store->fneutrino_x,2) + pow(event_store->fnu_reco_y - event_store->fneutrino_y,2) + pow(event_store->fnu_reco_z - event_store->fneutrino_z,2) ) );	    
	        h_vertexfitter_resolution_neutrino_not_merged->Fill( sqrt( pow(event_store->fnu_reco_fitter_x - event_store->fneutrino_x,2) + pow(event_store->fnu_reco_fitter_y - event_store->fneutrino_y,2) + pow(event_store->fnu_reco_fitter_z - event_store->fneutrino_z,2) ) );	    
	        h_vertexfitter_chi2ndf_neutrino_not_merged->Fill( event_store->fnu_reco_fitter_chi2ndf  );	    
		
    		for (unsigned j=0; j < event_store->fpdg.size(); j++) {
	    	if ( !reco_muon ) break; //select events with a reco muon
	    	float res = std::min ( sqrt( pow(event_store->freco_vertex_x[j] - event_store->fstart_x[j],2) + pow(event_store->freco_vertex_y[j] - event_store->fstart_y[j],2) + pow(event_store->freco_vertex_z[j] - event_store->fstart_z[j],2)) ,
					sqrt( pow(event_store->freco_vertex_x[j] - event_store->fend_x[j],2) + pow(event_store->freco_vertex_y[j] - event_store->fend_y[j],2) + pow(event_store->freco_vertex_z[j] - event_store->fend_z[j],2)) );
	    	float res_fitter = std::min ( sqrt( pow(event_store->freco_vertexfitter_x[j] - event_store->fstart_x[j],2) + pow(event_store->freco_vertexfitter_y[j] - event_store->fstart_y[j],2) + pow(event_store->freco_vertexfitter_z[j] - event_store->fstart_z[j],2)) ,
					sqrt( pow(event_store->freco_vertexfitter_x[j] - event_store->fend_x[j],2) + pow(event_store->freco_vertexfitter_y[j] - event_store->fend_y[j],2) + pow(event_store->freco_vertexfitter_z[j] - event_store->fend_z[j],2)) );
		
		if ( event_store->fpdg[j]==13 ) {
			h_vertex_resolution_muon_not_merged->Fill( res );
			h_vertexfitter_resolution_muon_not_merged->Fill( res_fitter );
			h_vertexfitter_chi2ndf_muon_not_merged->Fill( event_store->freco_vertexfitter_chi2ndf[j] );
		}
	    	if ( event_store->fpdg[j]!=2212 && event_store->fp0[j] <= 0.2 && !event_store->fis_tracked[j]) continue; //watch only reco protons with p>0.2Gev/c
   		hproton_pos_res_goodprotons->Fill(  sqrt( pow( event_store->freco_startx[j]- event_store->fstart_x[j] ,2) + pow( event_store->freco_starty[j]- event_store->fstart_y[j] ,2) + pow( event_store->freco_startz[j]- event_store->fstart_z[j] ,2) ) ); 
	    	h_vertex_resolution_proton_not_merged->Fill( res );
	    	h_vertexfitter_resolution_proton_not_merged->Fill( res_fitter );
	    	h_vertexfitter_chi2ndf_proton_not_merged->Fill( event_store->freco_vertexfitter_chi2ndf[j] );

   		if (lowmomentum_p) hproton_pos_res_lowprotons->Fill(  sqrt( pow( event_store->freco_startx[j]- event_store->fstart_x[j] ,2) + pow( event_store->freco_starty[j]- event_store->fstart_y[j] ,2) + pow( event_store->freco_startz[j]- event_store->fstart_z[j] ,2) ) ); 
		}
	} 

	if ( count_not_tracked > 0 ) { //at least 1 proton is not tracked
		for (unsigned jj=1; jj<event_store->fmuon_dqdx.size()-1; jj++) {
				if ( event_store->fmuon_range - event_store->fmuon_residual[jj] < 8 ) //look at 8cm only
				h_dqdx_1d_merged->Fill(event_store->fmuon_dqdx[jj]);

				h_dqdx_merged->Fill(  event_store->fmuon_range - event_store->fmuon_residual[jj], event_store->fmuon_dqdx[jj] );
				h_dqdx_merged_service->Fill( event_store->fmuon_range -  event_store->fmuon_residual[jj], event_store->fmuon_dqdx[jj] );
		}
		TH1D* h1 = NULL;
		h1 = h_dqdx_merged_service->ProjectionY("h",1,length_cut);
		if (h1)
		htail_to_tot_merged->Fill( h1->Integral(low_edge,high_edge) );
		h_dqdx_merged_service->Reset();
		    
		hmuon_pos_res_badprotons->Fill( sqrt( pow( event_store->freco_startx[muon_pos]- event_store->fstart_x[muon_pos] ,2) + pow( event_store->freco_starty[muon_pos]- event_store->fstart_y[muon_pos] ,2) + pow( event_store->freco_startz[muon_pos]- event_store->fstart_z[muon_pos] ,2) ) );
    		if ( !(std::abs(event_store->fnu_reco_x+1)<DBL_EPSILON && std::abs(event_store->fnu_reco_y+1)<DBL_EPSILON && std::abs(event_store->fnu_reco_z+1)<DBL_EPSILON) ) {
	        h_vertex_resolution_neutrino_merged->Fill( sqrt( pow(event_store->fnu_reco_x - event_store->fneutrino_x,2) + pow(event_store->fnu_reco_y - event_store->fneutrino_y,2) + pow(event_store->fnu_reco_z - event_store->fneutrino_z,2) ) );	    
	        h_vertexfitter_resolution_neutrino_merged->Fill( sqrt( pow(event_store->fnu_reco_fitter_x - event_store->fneutrino_x,2) + pow(event_store->fnu_reco_fitter_y - event_store->fneutrino_y,2) + pow(event_store->fnu_reco_fitter_z - event_store->fneutrino_z,2) ) );	    
	        h_vertexfitter_chi2ndf_neutrino_merged->Fill( event_store->fnu_reco_fitter_chi2ndf );	    
		}
    		for (unsigned j=0; j<event_store->fpdg.size(); j++) {
	    	if ( !reco_muon ) break; //select events with a reco muon
	    	float res = std::min ( sqrt( pow(event_store->freco_vertex_x[j] - event_store->fstart_x[j],2) + pow(event_store->freco_vertex_y[j] - event_store->fstart_y[j],2) + pow(event_store->freco_vertex_z[j] - event_store->fstart_z[j],2)) ,
					sqrt( pow(event_store->freco_vertex_x[j] - event_store->fend_x[j],2) + pow(event_store->freco_vertex_y[j] - event_store->fend_y[j],2) + pow(event_store->freco_vertex_z[j] - event_store->fend_z[j],2)) );
	    	float res_fitter = std::min ( sqrt( pow(event_store->freco_vertexfitter_x[j] - event_store->fstart_x[j],2) + pow(event_store->freco_vertexfitter_y[j] - event_store->fstart_y[j],2) + pow(event_store->freco_vertexfitter_z[j] - event_store->fstart_z[j],2)) ,
					sqrt( pow(event_store->freco_vertexfitter_x[j] - event_store->fend_x[j],2) + pow(event_store->freco_vertexfitter_y[j] - event_store->fend_y[j],2) + pow(event_store->freco_vertexfitter_z[j] - event_store->fend_z[j],2)) );
		
		if ( event_store->fpdg[j]==13 ) {
			h_vertex_resolution_muon_merged->Fill( res );
			h_vertexfitter_resolution_muon_merged->Fill( res_fitter );
			h_vertexfitter_chi2ndf_muon_merged->Fill( event_store->freco_vertexfitter_chi2ndf[j] );
		}
	    	if ( event_store->fpdg[j]!=2212 && event_store->fp0[j] <= 0.2 && !event_store->fis_tracked[j]) continue; //watch only reco protons with p>0.2Gev/c
   		hproton_pos_res_badprotons->Fill(  sqrt( pow( event_store->freco_startx[j]- event_store->fstart_x[j] ,2) + pow( event_store->freco_starty[j]- event_store->fstart_y[j] ,2) + pow( event_store->freco_startz[j]- event_store->fstart_z[j] ,2) ) ); 
	    	h_vertex_resolution_proton_merged->Fill( res );
	    	h_vertexfitter_resolution_proton_merged->Fill( res_fitter );
	    	h_vertexfitter_chi2ndf_proton_merged->Fill( event_store->freco_vertexfitter_chi2ndf[j] );
   		if (lowmomentum_p) hproton_pos_res_lowprotons->Fill(  sqrt( pow( event_store->freco_startx[j]- event_store->fstart_x[j] ,2) + pow( event_store->freco_starty[j]- event_store->fstart_y[j] ,2) + pow( event_store->freco_startz[j]- event_store->fstart_z[j] ,2) ) ); 
		}
	}

	if ( reco_muon && count_tracked == 0 && count_not_tracked== 0 && !is_pion && lowmomentum_p ) { //no protons tracked but low energy ones
		is_lowmomentum_p = true;
		for (unsigned jj=1; jj<event_store->fmuon_dqdx.size()-1; jj++) {
				h_dqdx_low_protons->Fill( event_store->fmuon_range -  event_store->fmuon_residual[jj], event_store->fmuon_dqdx[jj] );
				h_dqdx_low_protons_service->Fill( event_store->fmuon_range -  event_store->fmuon_residual[jj], event_store->fmuon_dqdx[jj] );
		}
		TH1D* h1 = NULL;
		h1 = h_dqdx_low_protons_service->ProjectionY("h",1,length_cut);
		if (h1)
		htail_to_tot_low_protons->Fill(h1->Integral(low_edge,high_edge)/h1->Integral(1,high_edge));
		h_dqdx_low_protons_service->Reset();
		
		hmuon_pos_res_lowprotons->Fill( sqrt( pow( event_store->freco_startx[muon_pos]- event_store->fstart_x[muon_pos] ,2) + pow( event_store->freco_starty[muon_pos]- event_store->fstart_y[muon_pos] ,2) + pow( event_store->freco_startz[muon_pos]- event_store->fstart_z[muon_pos] ,2) ) );
	}
  
if (count_not_tracked) {
	  std::cout << ">>>>NOT TRACKED!" << std::endl;
	  std::cout << "Number of not tracked " << count_not_tracked << ", number of tracked " << count_tracked << std::endl;
	  std::cout<< "Run: " << event_store->fRun << " SubRun: " << event_store->fSubRun << " Event: " << event_store->fEvent << std::endl;
  }

	fmuon_pos = muon_pos;


}//end fill analysis histos

void HistoMaker::FillCumulativeHistograms() {

    for (long ii=1; ii<1000;ii++) {
    TH1D* h = h_dqdx_not_merged->ProjectionY("h",1,ii);
    h_dqdx_tailtotot_length_not_merged->Fill(  ii, h->Integral(low_edge,high_edge)/h->Integral(1,high_edge) );
    h = h_dqdx_merged->ProjectionY("h",1,ii);
    h_dqdx_tailtotot_length_merged->Fill(  ii,  h->Integral(low_edge,high_edge)/h->Integral(1,high_edge));
    }

}


void HistoMaker::Fill_Hit_Histos( StoredEvent* event_store ) {
  
  if (!is_init_hit) {
	  std::cout << "Histo Maker NOT initialized" << std::endl;
	  throw cet::exception("Configuration");
  }

  //now fill histograms on hit merging
  for (unsigned jj=0; jj<event_store->fpdg.size(); jj++) {

	
	   for (unsigned z=0; z<event_store->fnot_clustered_tracked_charge[jj].size(); z++ ) {
	   //tracked particle, with a non-clustered hit
	   if ( event_store->fis_tracked[jj] ) {  //check and fill histos for a MCP which is tracked
	   	h_hits_not_clustered_tracked_charge->Fill ( event_store->fnot_clustered_tracked_charge[jj].at(z) );
		if (jj==unsigned(fmuon_pos)) //is muon
			h_hits_not_clustered_tracked_charge_muon->Fill ( event_store->fnot_clustered_tracked_charge[jj].at(z) );
		else if (event_store->fpdg[jj] == 2212) //is proton
			h_hits_not_clustered_tracked_charge_proton->Fill ( event_store->fnot_clustered_tracked_charge[jj].at(z) );
		        
		        if (event_store->fis_spacepoint[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )) { //remember: jj indexes the mcp and z the hit
			h_tracked_not_clustered_distance_nuvtx->Fill( sqrt( pow( event_store->fspacepoint_xyz[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )[0] - event_store->fneutrino_x,2) 
									+ pow( event_store->fspacepoint_xyz[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )[1] - event_store->fneutrino_y,2) 
									+ pow( event_store->fspacepoint_xyz[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )[2] - event_store->fneutrino_z,2)) ); //distance between neutrino vertex z and hit z
			if ( jj==unsigned(fmuon_pos) ) { //if is muon
			h_tracked_not_clustered_muon_start->Fill( sqrt( pow( event_store->fspacepoint_xyz[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )[0] - event_store->fstart_x[jj],2) 
									+ pow( event_store->fstart_y[jj] - event_store->fspacepoint_xyz[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )[1],2) 
									+ pow( event_store->fspacepoint_xyz[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )[2] - event_store->fstart_z[jj],2)) ); //z distance between muon starting point and not clustered hit
			h_tracked_not_clustered_muon_end->Fill( sqrt( pow( event_store->fspacepoint_xyz[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )[0] - event_store->fend_x[jj],2) 
								      + pow( event_store->fend_y[jj] - event_store->fspacepoint_xyz[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )[1],2) 
								      + pow( event_store->fspacepoint_xyz[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )[2] - event_store->fend_z[jj],2)) ); //z distance between muon ending point and not clustered hit
			} //is muon
			
			if ( event_store->fpdg[jj] == 2212 ) { //if is proton
			h_tracked_not_clustered_proton_start->Fill( sqrt( pow( event_store->fspacepoint_xyz[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )[0] - event_store->fstart_x[jj],2) 
									  + pow( event_store->fstart_y[jj] - event_store->fspacepoint_xyz[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )[1],2) 
									  + pow( event_store->fspacepoint_xyz[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )[2] - event_store->fstart_z[jj],2)), event_store->flength[jj] ); //z distance between muon starting point and not clustered hit
			h_tracked_not_clustered_proton_end->Fill( sqrt( pow( event_store->fspacepoint_xyz[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )[0] - event_store->fend_x[jj],2) 
									+ pow( event_store->fend_y[jj] - event_store->fspacepoint_xyz[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )[1],2) 
									+ pow( event_store->fspacepoint_xyz[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )[2] - event_store->fend_z[jj],2)), event_store->flength[jj] ); //z distance between muon ending point and not clustered hit
			 }//is proton
			}//is spacepoint

		       } else { //the particle is not tracked
				if ( jj==unsigned(fmuon_pos) ) { //is muon
				h_fraction_pdgs_not_tracked_not_clustered->Fill(0.,1./event_store->freco_mcp_collection_hits[jj]); //muon
			  	h_hits_not_clustered_not_tracked_charge_muon->Fill ( event_store->fnot_clustered_not_tracked_charge[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) ) );
				if (event_store->fis_spacepoint[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )){
			h_not_tracked_not_clustered_muon_start->Fill( sqrt( pow( event_store->fspacepoint_xyz[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )[0] - event_store->fstart_x[jj],2) 
									+ pow( event_store->fstart_y[jj] - event_store->fspacepoint_xyz[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )[1],2) 
									+ pow( event_store->fspacepoint_xyz[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )[2] - event_store->fstart_z[jj],2)) ); //z distance between muon starting point and not clustered hit
			h_not_tracked_not_clustered_muon_end->Fill( sqrt( pow( event_store->fspacepoint_xyz[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )[0] - event_store->fend_x[jj],2) 
								      + pow( event_store->fend_y[jj] - event_store->fspacepoint_xyz[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )[1],2) 
								      + pow( event_store->fspacepoint_xyz[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )[2] - event_store->fend_z[jj],2)) ); //z distance between muon ending point and not clustered hit
				} //is spacepoint
				
				} else if ( event_store->fpdg[jj] == 2212 ) { //is proton
				h_fraction_pdgs_not_tracked_not_clustered->Fill(1.,1./event_store->freco_mcp_collection_hits[jj]); //proton all
				if ( event_store->fp0[jj] > 0.2)
				h_fraction_pdgs_not_tracked_not_clustered->Fill(2.,1./event_store->freco_mcp_collection_hits[jj]); //protons > 20MeV
				else
				h_fraction_pdgs_not_tracked_not_clustered->Fill(3.,1./event_store->freco_mcp_collection_hits[jj]); //protons <20 MeV
			  	
				h_hits_not_clustered_not_tracked_charge_proton->Fill ( event_store->fnot_clustered_not_tracked_charge[jj].at(z) );
				if (event_store->fis_spacepoint[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )){
			h_not_tracked_not_clustered_proton_start->Fill( sqrt( pow( event_store->fspacepoint_xyz[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )[0] - event_store->fstart_x[jj],2) 
									+ pow( event_store->fstart_y[jj] - event_store->fspacepoint_xyz[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )[1],2) 
									+ pow( event_store->fspacepoint_xyz[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )[2] - event_store->fstart_z[jj],2)), event_store->flength[jj] ); //z distance between muon starting point and not clustered hit
			h_not_tracked_not_clustered_proton_end->Fill( sqrt( pow( event_store->fspacepoint_xyz[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )[0] - event_store->fend_x[jj],2) 
								      + pow( event_store->fend_y[jj] - event_store->fspacepoint_xyz[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )[1],2) 
								      + pow( event_store->fspacepoint_xyz[jj].at( event_store->fnot_clustered_hit_index[jj].at(z) )[2] - event_store->fend_z[jj],2)), event_store->flength[jj] ); //z distance between muon ending point and not clustered hit
				} //is spacepoint
		       		} //is proton
		       } //not tracked
	   }//end loop on not clustered hits
	   



	  if ( count_not_tracked == 0 ) { //all protons above 20MeV are tracked
	  if ( jj == unsigned(fmuon_pos) ) { //is muon
		  if ( event_store->freco_mcp_collection_hits[jj] ) {
		  h_muon_clustering_prob_good_protons->Fill(0., float(event_store->fnot_clustered[jj])/float(event_store->freco_mcp_collection_hits[jj]) );	
		  h_muon_clustering_prob_good_protons->Fill(1., float(event_store->fclustered_matched[jj])/float(event_store->freco_mcp_collection_hits[jj]) );	
		  h_muon_clustering_prob_good_protons->Fill(2., float(event_store->fclustered_mismatched[jj])/float(event_store->freco_mcp_collection_hits[jj]) );	
		  }
		
		  for (unsigned ii=0; ii<event_store->fhit_mismatch_pdg[jj].size(); ii++)
		  h_muon_clustering_mismatch_pdg_good_protons->Fill(event_store->fhit_mismatch_pdg[jj][ii]);	

		  if ( event_store->fnot_clustered[jj] > 0) 
		  h_muon_not_clustered_reco_hits_good_protons->Fill( event_store->fnot_clustered[jj] , event_store->freco_mcp_collection_hits[jj] );
		  if ( event_store->fclustered_matched[jj] > 0) 
		  h_muon_clustered_matched_reco_hits_good_protons->Fill( event_store->fclustered_matched[jj] , event_store->freco_mcp_collection_hits[jj] );
		  if ( event_store->fclustered_mismatched[jj] > 0) 
		  h_muon_clustered_mismatched_reco_hits_good_protons->Fill( event_store->fclustered_mismatched[jj] , event_store->freco_mcp_collection_hits[jj] );
                  
		  for ( unsigned ii=0; ii< event_store->fnot_clustered_charge[jj].size(); ii++ ) 
		  h_muon_not_clustered_reco_charge_good_protons->Fill( event_store->fnot_clustered_charge[jj][ii], float(event_store->fnot_clustered[jj])/float(event_store->freco_mcp_collection_hits[jj]) );
		  for ( unsigned ii=0; ii< event_store->fclustered_matched_charge[jj].size(); ii++ ) 
		  h_muon_clustered_matched_reco_charge_good_protons->Fill( event_store->fclustered_matched_charge[jj][ii], float(event_store->fclustered_matched[jj])/float(event_store->freco_mcp_collection_hits[jj]) );
		  for ( unsigned ii=0; ii< event_store->fclustered_mismatched_charge[jj].size(); ii++ ) 
		  h_muon_clustered_mismatched_reco_charge_good_protons->Fill( event_store->fclustered_mismatched_charge[jj][ii], float(event_store->fclustered_mismatched[jj])/float(event_store->freco_mcp_collection_hits[jj]) );
		  
		  for ( unsigned zz=0; zz<event_store->fpdg.size(); zz++) { //muon and proton plots
			  if (event_store->fpdg[zz]!=2212) continue;
		  if (event_store->fnot_clustered[jj]) {
		  h_muon_NC_lateral_hits_good_protons->Fill( float(event_store->fnot_clustered[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->flength[zz] * sqrt( 1 - pow( event_store->fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_NC_costheta_hits_good_protons->Fill( float(event_store->fnot_clustered[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->fcostheta_muon[zz]  ); //lateral distance between muon and proton
		  }
		  if (event_store->fclustered_matched[jj]) {
		  h_muon_CMA_lateral_hits_good_protons->Fill( float(event_store->fclustered_matched[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->flength[zz] * sqrt( 1 - pow( event_store->fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_CMA_costheta_hits_good_protons->Fill( float(event_store->fclustered_matched[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->fcostheta_muon[zz]  ); //lateral distance between muon and proton
		  }
		  if (event_store->fclustered_mismatched[jj]){
		  h_muon_CMI_lateral_hits_good_protons->Fill( float(event_store->fclustered_mismatched[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->flength[zz] * sqrt( 1 - pow( event_store->fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_CMI_costheta_hits_good_protons->Fill( float(event_store->fclustered_mismatched[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->fcostheta_muon[zz]  ); //lateral distance between muon and proton
		  }

		  for ( unsigned ii=0; ii< event_store->fnot_clustered_charge[jj].size(); ii++ ) {
                  h_muon_NC_lateral_charge_good_protons->Fill( event_store->fnot_clustered_charge[jj][ii],  event_store->flength[zz] * sqrt( 1 - pow( event_store->fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_NC_costheta_charge_good_protons->Fill( event_store->fnot_clustered_charge[jj][ii], event_store->fcostheta_muon[zz] ); //lateral distance between muon and proton
		  }
		  for ( unsigned ii=0; ii< event_store->fclustered_matched_charge[jj].size(); ii++ ) {
                  h_muon_CMA_lateral_charge_good_protons->Fill( event_store->fclustered_matched_charge[jj][ii],  event_store->flength[zz] * sqrt( 1 - pow( event_store->fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_CMA_costheta_charge_good_protons->Fill( event_store->fclustered_matched_charge[jj][ii], event_store->fcostheta_muon[zz] ); //lateral distance between muon and proton
		  }
		  for ( unsigned ii=0; ii< event_store->fclustered_mismatched_charge[jj].size(); ii++ ) {
                  h_muon_CMI_lateral_charge_good_protons->Fill( event_store->fclustered_mismatched_charge[jj][ii],  event_store->flength[zz] * sqrt( 1 - pow( event_store->fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_CMI_costheta_charge_good_protons->Fill( event_store->fclustered_mismatched_charge[jj][ii], event_store->fcostheta_muon[zz] ); //lateral distance between muon and proton
		  }
		  } //fpdg
		  } //muon
	  
	  if (event_store->fpdg[jj] ==2212 && event_store->fis_tracked[jj]) {

		  if ( event_store->freco_mcp_collection_hits[jj] ) {
		  h_proton_clustering_prob_good_protons->Fill(0., float(event_store->fnot_clustered[jj])/float(event_store->freco_mcp_collection_hits[jj]) );	
		  h_proton_clustering_prob_good_protons->Fill(1., float(event_store->fclustered_matched[jj])/float(event_store->freco_mcp_collection_hits[jj]) );	
		  h_proton_clustering_prob_good_protons->Fill(2., float(event_store->fclustered_mismatched[jj])/float(event_store->freco_mcp_collection_hits[jj]) );	
		  }

		  for (unsigned ii=0; ii<event_store->fhit_mismatch_pdg[jj].size(); ii++)
		  h_proton_clustering_mismatch_pdg_good_protons->Fill(event_store->fhit_mismatch_pdg[jj][ii]);	

		  if (event_store->fnot_clustered[jj]) {
		  h_proton_not_clustered_reco_hits_good_protons->Fill( event_store->fnot_clustered[jj] , event_store->freco_mcp_collection_hits[jj] );
                  h_proton_NC_lateral_hits_good_protons->Fill( float(event_store->fnot_clustered[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->flength[jj] * sqrt( 1 - pow( event_store->fcostheta_muon[jj] ,2 ) ) ); //lateral distance between proton and proton
                  h_proton_NC_costheta_hits_good_protons->Fill( float(event_store->fnot_clustered[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->fcostheta_muon[jj]  ); //lateral distance between muon and proton
		  }
		  if (event_store->fclustered_matched[jj]) {
		  h_proton_clustered_matched_reco_hits_good_protons->Fill( event_store->fclustered_matched[jj] , event_store->freco_mcp_collection_hits[jj] );
                  h_proton_CMA_lateral_hits_good_protons->Fill( float(event_store->fclustered_matched[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->flength[jj] * sqrt( 1 - pow( event_store->fcostheta_muon[jj] ,2 ) ) ); //lateral distance between proton and proton
                  h_proton_CMA_costheta_hits_good_protons->Fill( float(event_store->fclustered_matched[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->fcostheta_muon[jj]  ); //lateral distance between muon and proton
		  }
		  if (event_store->fclustered_mismatched[jj]) {
		  h_proton_clustered_mismatched_reco_hits_good_protons->Fill( event_store->fclustered_mismatched[jj] , event_store->freco_mcp_collection_hits[jj] );
                  h_proton_CMI_lateral_hits_good_protons->Fill( float(event_store->fclustered_mismatched[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->flength[jj] * sqrt( 1 - pow( event_store->fcostheta_muon[jj] ,2 ) ) ); //lateral distance between proton and proton
                  h_proton_CMI_costheta_hits_good_protons->Fill( float(event_store->fclustered_mismatched[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->fcostheta_muon[jj]  ); //lateral distance between muon and proton
		  }
		 
		  for ( unsigned ii=0; ii< event_store->fnot_clustered_charge[jj].size(); ii++ ) {
		  h_proton_not_clustered_reco_charge_good_protons->Fill( event_store->fnot_clustered_charge[jj][ii], float(event_store->fnot_clustered[jj])/float(event_store->freco_mcp_collection_hits[jj]) );
                  h_proton_NC_lateral_charge_good_protons->Fill( event_store->fnot_clustered_charge[jj][ii],  event_store->flength[jj] * sqrt( 1 - pow( event_store->fcostheta_muon[jj] ,2 ) ) ); //lateral distance between muon and proton
                  h_proton_NC_costheta_charge_good_protons->Fill( event_store->fnot_clustered_charge[jj][ii], event_store->fcostheta_muon[jj] ); //lateral distance between muon and proton
		  }
		  for ( unsigned ii=0; ii< event_store->fclustered_matched_charge[jj].size(); ii++ ) {
		  h_proton_clustered_matched_reco_charge_good_protons->Fill( event_store->fclustered_matched_charge[jj][ii], float(event_store->fclustered_matched[jj])/float(event_store->freco_mcp_collection_hits[jj]) );
                  h_proton_CMA_lateral_charge_good_protons->Fill( event_store->fclustered_matched_charge[jj][ii],  event_store->flength[jj] * sqrt( 1 - pow( event_store->fcostheta_muon[jj] ,2 ) ) ); //lateral distance between muon and proton
                  h_proton_CMA_costheta_charge_good_protons->Fill( event_store->fclustered_matched_charge[jj][ii], event_store->fcostheta_muon[jj] ); //lateral distance between muon and proton
		  }
		  for ( unsigned ii=0; ii< event_store->fclustered_mismatched_charge[jj].size(); ii++ ) {
		  h_proton_clustered_mismatched_reco_charge_good_protons->Fill( event_store->fclustered_mismatched_charge[jj][ii], float(event_store->fclustered_mismatched[jj])/float(event_store->freco_mcp_collection_hits[jj]) );
                  h_proton_CMI_lateral_charge_good_protons->Fill( event_store->fclustered_mismatched_charge[jj][ii],  event_store->flength[jj] * sqrt( 1 - pow( event_store->fcostheta_muon[jj] ,2 ) ) ); //lateral distance between muon and proton
                  h_proton_CMI_costheta_charge_good_protons->Fill( event_store->fclustered_mismatched_charge[jj][ii], event_store->fcostheta_muon[jj] ); //lateral distance between muon and proton
		  }
	     }//proton

	  } //if all are tracked
	  
if ( count_not_tracked > 0 ) { //some are not tracked
	  
	if ( jj == unsigned(fmuon_pos) ) { //the muon is tracked by definition
		  if ( event_store->freco_mcp_collection_hits[jj] ) {
		  h_muon_clustering_prob_bad_protons->Fill(0., float(event_store->fnot_clustered[jj])/float(event_store->freco_mcp_collection_hits[jj]) );	
		  h_muon_clustering_prob_bad_protons->Fill(1., float(event_store->fclustered_matched[jj])/float(event_store->freco_mcp_collection_hits[jj]) );	
		  h_muon_clustering_prob_bad_protons->Fill(2., float(event_store->fclustered_mismatched[jj])/float(event_store->freco_mcp_collection_hits[jj]) );	
		  }

		  for (unsigned ii=0; ii<event_store->fhit_mismatch_pdg[jj].size(); ii++)
		  h_muon_clustering_mismatch_pdg_bad_protons->Fill(event_store->fhit_mismatch_pdg[jj][ii]);	
		 
		  if (event_store->fnot_clustered[jj])
		  h_muon_not_clustered_reco_hits_bad_protons->Fill( event_store->fnot_clustered[jj] , event_store->freco_mcp_collection_hits[jj] );
		  if (event_store->fclustered_matched[jj])
		  h_muon_clustered_matched_reco_hits_bad_protons->Fill( event_store->fclustered_matched[jj] , event_store->freco_mcp_collection_hits[jj] );
		  if (event_store->fclustered_mismatched[jj])
		  h_muon_clustered_mismatched_reco_hits_bad_protons->Fill( event_store->fclustered_mismatched[jj] , event_store->freco_mcp_collection_hits[jj] );
		  
		  for ( unsigned ii=0; ii< event_store->fnot_clustered_charge[jj].size(); ii++ ) 
		  h_muon_not_clustered_reco_charge_bad_protons->Fill( event_store->fnot_clustered_charge[jj][ii], float(event_store->fnot_clustered[jj])/float(event_store->freco_mcp_collection_hits[jj]) );
		  for ( unsigned ii=0; ii< event_store->fclustered_matched_charge[jj].size(); ii++ ) 
		  h_muon_clustered_matched_reco_charge_bad_protons->Fill( event_store->fclustered_matched_charge[jj][ii], float(event_store->fclustered_matched[jj])/float(event_store->freco_mcp_collection_hits[jj]) );
		  for ( unsigned ii=0; ii< event_store->fclustered_mismatched_charge[jj].size(); ii++ ) 
		  h_muon_clustered_mismatched_reco_charge_bad_protons->Fill( event_store->fclustered_mismatched_charge[jj][ii], float(event_store->fclustered_mismatched[jj])/float(event_store->freco_mcp_collection_hits[jj]) );

		  for ( unsigned zz=0; zz<event_store->fpdg.size(); zz++) {
			  if (event_store->fpdg[zz]!=2212) continue;
                  if (event_store->fnot_clustered[jj]) {
		  h_muon_NC_lateral_hits_bad_protons->Fill( float(event_store->fnot_clustered[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->flength[zz] * sqrt( 1 - pow( event_store->fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_NC_costheta_hits_bad_protons->Fill( float(event_store->fnot_clustered[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->fcostheta_muon[zz]  ); //lateral distance between muon and proton
		  }
		  if (event_store->fclustered_matched[jj]) {
                  h_muon_CMA_lateral_hits_bad_protons->Fill( float(event_store->fclustered_matched[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->flength[zz] * sqrt( 1 - pow( event_store->fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_CMA_costheta_hits_bad_protons->Fill( float(event_store->fclustered_matched[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->fcostheta_muon[zz]  ); //lateral distance between muon and proton
		  }
		  if (event_store->fclustered_mismatched[jj]){
                  h_muon_CMI_lateral_hits_bad_protons->Fill( float(event_store->fclustered_mismatched[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->flength[zz] * sqrt( 1 - pow( event_store->fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_CMI_costheta_hits_bad_protons->Fill( float(event_store->fclustered_mismatched[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->fcostheta_muon[zz]  ); //lateral distance between muon and proton
		  }
		 
		  for ( unsigned ii=0; ii< event_store->fnot_clustered_charge[jj].size(); ii++ ) {
                  h_muon_NC_lateral_charge_bad_protons->Fill( event_store->fnot_clustered_charge[jj][ii],  event_store->flength[zz] * sqrt( 1 - pow( event_store->fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_NC_costheta_charge_bad_protons->Fill( event_store->fnot_clustered_charge[jj][ii], event_store->fcostheta_muon[zz] ); //lateral distance between muon and proton
		  }
		  for ( unsigned ii=0; ii< event_store->fclustered_matched_charge[jj].size(); ii++ ) {
                  h_muon_CMA_lateral_charge_bad_protons->Fill( event_store->fclustered_matched_charge[jj][ii],  event_store->flength[zz] * sqrt( 1 - pow( event_store->fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_CMA_costheta_charge_bad_protons->Fill( event_store->fclustered_matched_charge[jj][ii], event_store->fcostheta_muon[zz] ); //lateral distance between muon and proton
		  }
		  for ( unsigned ii=0; ii< event_store->fclustered_mismatched_charge[jj].size(); ii++ ) {
                  h_muon_CMI_lateral_charge_bad_protons->Fill( event_store->fclustered_mismatched_charge[jj][ii],  event_store->flength[zz] * sqrt( 1 - pow( event_store->fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_CMI_costheta_charge_bad_protons->Fill( event_store->fclustered_mismatched_charge[jj][ii], event_store->fcostheta_muon[zz] ); //lateral distance between muon and proton
		  }
		  } //fpdg
		} //muon
	  
	if (event_store->fpdg[jj] ==2212 && event_store->fis_tracked[jj]) {
		  if ( event_store->freco_mcp_collection_hits[jj] ) {
		  h_proton_clustering_prob_bad_protons->Fill(0., float(event_store->fnot_clustered[jj])/float(event_store->freco_mcp_collection_hits[jj]) );	
		  h_proton_clustering_prob_bad_protons->Fill(1., float(event_store->fclustered_matched[jj])/float(event_store->freco_mcp_collection_hits[jj]) );	
		  h_proton_clustering_prob_bad_protons->Fill(2., float(event_store->fclustered_mismatched[jj])/float(event_store->freco_mcp_collection_hits[jj]) );	
		  } 
		  for (unsigned ii=0; ii<event_store->fhit_mismatch_pdg[jj].size(); ii++) 
		  h_proton_clustering_mismatch_pdg_bad_protons->Fill(event_store->fhit_mismatch_pdg[jj][ii]);	

		  if (event_store->fnot_clustered[jj]) {
		  h_proton_not_clustered_reco_hits_bad_protons->Fill( event_store->fnot_clustered[jj] , event_store->freco_mcp_collection_hits[jj] );
                  h_proton_NC_lateral_hits_bad_protons->Fill( float(event_store->fnot_clustered[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->flength[jj] * sqrt( 1 - pow( event_store->fcostheta_muon[jj] ,2 ) ) ); //lateral distance between proton and proton
                  h_proton_NC_costheta_hits_bad_protons->Fill( float(event_store->fnot_clustered[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->fcostheta_muon[jj]  ); //lateral distance between muon and proton
		  }
		  if (event_store->fclustered_matched[jj]){
		  h_proton_clustered_matched_reco_hits_bad_protons->Fill( event_store->fclustered_matched[jj] , event_store->freco_mcp_collection_hits[jj] );
                  h_proton_CMA_lateral_hits_bad_protons->Fill( float(event_store->fclustered_matched[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->flength[jj] * sqrt( 1 - pow( event_store->fcostheta_muon[jj] ,2 ) ) ); //lateral distance between proton and proton
                  h_proton_CMA_costheta_hits_bad_protons->Fill( float(event_store->fclustered_matched[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->fcostheta_muon[jj]  ); //lateral distance between muon and proton
		  }
		  if (event_store->fclustered_mismatched[jj]){
		  h_proton_clustered_mismatched_reco_hits_bad_protons->Fill( event_store->fclustered_mismatched[jj] , event_store->freco_mcp_collection_hits[jj] );
                  h_proton_CMI_lateral_hits_bad_protons->Fill( float(event_store->fclustered_mismatched[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->flength[jj] * sqrt( 1 - pow( event_store->fcostheta_muon[jj] ,2 ) ) ); //lateral distance between proton and proton
                  h_proton_CMI_costheta_hits_bad_protons->Fill( float(event_store->fclustered_mismatched[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->fcostheta_muon[jj]  ); //lateral distance between muon and proton
		  }

		  for ( unsigned ii=0; ii< event_store->fnot_clustered_charge[jj].size(); ii++ ) {
		  h_proton_not_clustered_reco_charge_bad_protons->Fill( event_store->fnot_clustered_charge[jj][ii], float(event_store->fnot_clustered[jj])/float(event_store->freco_mcp_collection_hits[jj]) );
                  h_proton_NC_lateral_charge_bad_protons->Fill( event_store->fnot_clustered_charge[jj][ii],  event_store->flength[jj] * sqrt( 1 - pow( event_store->fcostheta_muon[jj] ,2 ) ) ); //lateral distance between muon and proton
                  h_proton_NC_costheta_charge_bad_protons->Fill( event_store->fnot_clustered_charge[jj][ii], event_store->fcostheta_muon[jj] ); //lateral distance between muon and proton
		  }
		  for ( unsigned ii=0; ii< event_store->fclustered_matched_charge[jj].size(); ii++ ) {
		  h_proton_clustered_matched_reco_charge_bad_protons->Fill( event_store->fclustered_matched_charge[jj][ii], float(event_store->fclustered_matched[jj])/float(event_store->freco_mcp_collection_hits[jj]) );
                  h_proton_CMA_lateral_charge_bad_protons->Fill( event_store->fclustered_matched_charge[jj][ii],  event_store->flength[jj] * sqrt( 1 - pow( event_store->fcostheta_muon[jj] ,2 ) ) ); //lateral distance between muon and proton
                  h_proton_CMA_costheta_charge_bad_protons->Fill( event_store->fclustered_matched_charge[jj][ii], event_store->fcostheta_muon[jj] ); //lateral distance between muon and proton
		  }
		  for ( unsigned ii=0; ii< event_store->fclustered_mismatched_charge[jj].size(); ii++ ) {
		  h_proton_clustered_mismatched_reco_charge_bad_protons->Fill( event_store->fclustered_mismatched_charge[jj][ii], float(event_store->fclustered_mismatched[jj])/float(event_store->freco_mcp_collection_hits[jj]) );
                  h_proton_CMI_lateral_charge_bad_protons->Fill( event_store->fclustered_mismatched_charge[jj][ii],  event_store->flength[jj] * sqrt( 1 - pow( event_store->fcostheta_muon[jj] ,2 ) ) ); //lateral distance between muon and proton
                  h_proton_CMI_costheta_charge_bad_protons->Fill( event_store->fclustered_mismatched_charge[jj][ii], event_store->fcostheta_muon[jj] ); //lateral distance between muon and proton
		  }
	  } //protons && is tracked 
	
	if ( event_store->fpdg[jj] == 2212 && !event_store->fis_tracked[jj] ) {
		  if ( event_store->freco_mcp_collection_hits[jj] ) {
		  h_proton_clustering_prob_bad_protons_not_tracked->Fill(0., float(event_store->fnot_clustered[jj])/float(event_store->freco_mcp_collection_hits[jj]) );	
		  h_proton_clustering_prob_bad_protons_not_tracked->Fill(1., float(event_store->fclustered_matched[jj])/float(event_store->freco_mcp_collection_hits[jj]) );	
		  h_proton_clustering_prob_bad_protons_not_tracked->Fill(2., float(event_store->fclustered_mismatched[jj])/float(event_store->freco_mcp_collection_hits[jj]) );	
		  }
		  for (unsigned ii=0; ii<event_store->fhit_mismatch_pdg[jj].size(); ii++)
		  h_proton_clustering_mismatch_pdg_bad_protons_not_tracked->Fill(event_store->fhit_mismatch_pdg[jj][ii]);	
		 
		  if (event_store->fnot_clustered[jj]){
		  h_proton_not_clustered_reco_hits_bad_protons_not_tracked->Fill( event_store->fnot_clustered[jj] , event_store->freco_mcp_collection_hits[jj] );
                  h_proton_NC_lateral_hits_bad_protons_not_tracked->Fill( float(event_store->fnot_clustered[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->flength[jj] * sqrt( 1 - pow( event_store->fcostheta_muon[jj] ,2 ) ) ); //lateral distance between proton and proton
                  h_proton_NC_costheta_hits_bad_protons_not_tracked->Fill( float(event_store->fnot_clustered[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->fcostheta_muon[jj]  ); //lateral distance between muon and proton
		  }
		  if (event_store->fclustered_matched[jj]){
		  h_proton_clustered_matched_reco_hits_bad_protons_not_tracked->Fill( event_store->fclustered_matched[jj] , event_store->freco_mcp_collection_hits[jj] );
                  h_proton_CMA_lateral_hits_bad_protons_not_tracked->Fill( float(event_store->fclustered_matched[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->flength[jj] * sqrt( 1 - pow( event_store->fcostheta_muon[jj] ,2 ) ) ); //lateral distance between proton and proton
                  h_proton_CMA_costheta_hits_bad_protons_not_tracked->Fill( float(event_store->fclustered_matched[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->fcostheta_muon[jj]  ); //lateral distance between muon and proton
		  }
		  if (event_store->fclustered_mismatched[jj]){
		  h_proton_clustered_mismatched_reco_hits_bad_protons_not_tracked->Fill( event_store->fclustered_mismatched[jj] , event_store->freco_mcp_collection_hits[jj] );
                  h_proton_CMI_lateral_hits_bad_protons_not_tracked->Fill( float(event_store->fclustered_mismatched[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->flength[jj] * sqrt( 1 - pow( event_store->fcostheta_muon[jj] ,2 ) ) ); //lateral distance between proton and proton
                  h_proton_CMI_costheta_hits_bad_protons_not_tracked->Fill( float(event_store->fclustered_mismatched[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->fcostheta_muon[jj]  ); //lateral distance between muon and proton
		  }

		  for ( unsigned ii=0; ii< event_store->fnot_clustered_charge[jj].size(); ii++ ) {
		  h_proton_not_clustered_reco_charge_bad_protons_not_tracked->Fill( event_store->fnot_clustered_charge[jj][ii], float(event_store->fnot_clustered[jj])/float(event_store->freco_mcp_collection_hits[jj]) );
                  h_proton_NC_lateral_charge_bad_protons_not_tracked->Fill( event_store->fnot_clustered_charge[jj][ii],  event_store->flength[jj] * sqrt( 1 - pow( event_store->fcostheta_muon[jj] ,2 ) ) ); //lateral distance between muon and proton
                  h_proton_NC_costheta_charge_bad_protons_not_tracked->Fill( event_store->fnot_clustered_charge[jj][ii], event_store->fcostheta_muon[jj] ); //lateral distance between muon and proton
		  }
		  for ( unsigned ii=0; ii< event_store->fclustered_matched_charge[jj].size(); ii++ ) {
		  h_proton_clustered_matched_reco_charge_bad_protons_not_tracked->Fill( event_store->fclustered_matched_charge[jj][ii], float(event_store->fclustered_matched[jj])/float(event_store->freco_mcp_collection_hits[jj]) );
                  h_proton_CMA_lateral_charge_bad_protons_not_tracked->Fill( event_store->fclustered_matched_charge[jj][ii],  event_store->flength[jj] * sqrt( 1 - pow( event_store->fcostheta_muon[jj] ,2 ) ) ); //lateral distance between muon and proton
                  h_proton_CMA_costheta_charge_bad_protons_not_tracked->Fill( event_store->fclustered_matched_charge[jj][ii], event_store->fcostheta_muon[jj] ); //lateral distance between muon and proton
		  }
		  for ( unsigned ii=0; ii< event_store->fclustered_mismatched_charge[jj].size(); ii++ ) {
		  h_proton_clustered_mismatched_reco_charge_bad_protons_not_tracked->Fill( event_store->fclustered_mismatched_charge[jj][ii], float(event_store->fclustered_mismatched[jj])/float(event_store->freco_mcp_collection_hits[jj]) );
                  h_proton_CMI_lateral_charge_bad_protons_not_tracked->Fill( event_store->fclustered_mismatched_charge[jj][ii],  event_store->flength[jj] * sqrt( 1 - pow( event_store->fcostheta_muon[jj] ,2 ) ) ); //lateral distance between muon and proton
                  h_proton_CMI_costheta_charge_bad_protons_not_tracked->Fill( event_store->fclustered_mismatched_charge[jj][ii], event_store->fcostheta_muon[jj] ); //lateral distance between muon and proton
		  }
		} //proton && ! is tracked

}  //some are not tracked
	  
  if ( is_lowmomentum_p ) { //only protons below 20MeV, none is tracked
	if ( jj == unsigned(fmuon_pos) ) { //the muon is tracked by definition
		  if ( event_store->freco_mcp_collection_hits[jj] ) {
		  h_muon_clustering_prob_low_protons->Fill(0., float(event_store->fnot_clustered[jj])/float(event_store->freco_mcp_collection_hits[jj]) );	
		  h_muon_clustering_prob_low_protons->Fill(1., float(event_store->fclustered_matched[jj])/float(event_store->freco_mcp_collection_hits[jj]) );	
		  h_muon_clustering_prob_low_protons->Fill(2., float(event_store->fclustered_mismatched[jj])/float(event_store->freco_mcp_collection_hits[jj]) );	
		  }
		  for (unsigned ii=0; ii<event_store->fhit_mismatch_pdg[jj].size(); ii++)
		  h_muon_clustering_mismatch_pdg_low_protons->Fill(event_store->fhit_mismatch_pdg[jj][ii]);	
		  
		  for ( unsigned ii=0; ii< event_store->fnot_clustered_charge[jj].size(); ii++ )
		  h_muon_not_clustered_reco_charge_low_protons->Fill( event_store->fnot_clustered_charge[jj][ii], float(event_store->fnot_clustered[jj])/float(event_store->freco_mcp_collection_hits[jj]) );
		  for ( unsigned ii=0; ii< event_store->fclustered_matched_charge[jj].size(); ii++ )
		  h_muon_clustered_matched_reco_charge_low_protons->Fill( event_store->fclustered_matched_charge[jj][ii], float(event_store->fclustered_matched[jj])/float(event_store->freco_mcp_collection_hits[jj]) );
		  for ( unsigned ii=0; ii< event_store->fclustered_mismatched_charge[jj].size(); ii++ )
		  h_muon_clustered_mismatched_reco_charge_low_protons->Fill( event_store->fclustered_mismatched_charge[jj][ii], float(event_store->fclustered_mismatched[jj])/float(event_store->freco_mcp_collection_hits[jj]) );
		  
		  for ( unsigned zz=0; zz<event_store->fpdg.size(); zz++) {
			  if (event_store->fpdg[zz]!=2212) continue;
		  if (event_store->fnot_clustered[jj]) {
		  h_muon_not_clustered_reco_hits_low_protons->Fill( event_store->fnot_clustered[jj] , event_store->freco_mcp_collection_hits[jj] );
                  h_muon_NC_lateral_hits_low_protons->Fill( float(event_store->fnot_clustered[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->flength[zz] * sqrt( 1 - pow( event_store->fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_NC_costheta_hits_low_protons->Fill( float(event_store->fnot_clustered[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->fcostheta_muon[zz]  ); //lateral distance between muon and proton
		  } 
		  if (event_store->fclustered_matched[jj]) {
		  h_muon_clustered_matched_reco_hits_low_protons->Fill( event_store->fclustered_matched[jj] , event_store->freco_mcp_collection_hits[jj] );
                  h_muon_CMA_lateral_hits_low_protons->Fill( float(event_store->fclustered_matched[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->flength[zz] * sqrt( 1 - pow( event_store->fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_CMA_costheta_hits_low_protons->Fill( float(event_store->fclustered_matched[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->fcostheta_muon[zz]  ); //lateral distance between muon and proton
		  }
		  if (event_store->fclustered_mismatched[jj]) {
		  h_muon_clustered_mismatched_reco_hits_low_protons->Fill( event_store->fclustered_mismatched[jj] , event_store->freco_mcp_collection_hits[jj] );
                  h_muon_CMI_lateral_hits_low_protons->Fill( float(event_store->fclustered_mismatched[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->flength[zz] * sqrt( 1 - pow( event_store->fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_CMI_costheta_hits_low_protons->Fill( float(event_store->fclustered_mismatched[jj])/float(event_store->freco_mcp_collection_hits[jj]),  event_store->fcostheta_muon[zz]  ); //lateral distance between muon and proton
		  }

		  for ( unsigned ii=0; ii< event_store->fnot_clustered_charge[jj].size(); ii++ ) {
                  h_muon_NC_lateral_charge_low_protons->Fill( event_store->fnot_clustered_charge[jj][ii],  event_store->flength[zz] * sqrt( 1 - pow( event_store->fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_NC_costheta_charge_low_protons->Fill( event_store->fnot_clustered_charge[jj][ii], event_store->fcostheta_muon[zz] ); //lateral distance between muon and proton
		  }
		  for ( unsigned ii=0; ii< event_store->fclustered_matched_charge[jj].size(); ii++ ) {
                  h_muon_CMA_lateral_charge_low_protons->Fill( event_store->fclustered_matched_charge[jj][ii],  event_store->flength[zz] * sqrt( 1 - pow( event_store->fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_CMA_costheta_charge_low_protons->Fill( event_store->fclustered_matched_charge[jj][ii], event_store->fcostheta_muon[zz] ); //lateral distance between muon and proton
		  }
		  for ( unsigned ii=0; ii< event_store->fclustered_mismatched_charge[jj].size(); ii++ ) {
                  h_muon_CMI_lateral_charge_low_protons->Fill( event_store->fclustered_mismatched_charge[jj][ii],  event_store->flength[zz] * sqrt( 1 - pow( event_store->fcostheta_muon[zz] ,2 ) ) ); //lateral distance between muon and proton
                  h_muon_CMI_costheta_charge_low_protons->Fill( event_store->fclustered_mismatched_charge[jj][ii], event_store->fcostheta_muon[zz] ); //lateral distance between muon and proton
		  }
		    } //fpdg
		  } //muon
	  	} //low p

  } //loop on MCP for hit merging histograms
}//end fill hit

void HistoMaker::ScalePlots( int n_events ) {
  //scale some histos	
  if (n_events) {
  hproton_multi_all->Scale(1./n_events);
  hproton_multi_above20MeV->Scale(1./n_events);
  hproton_multi_below20MeV->Scale(1./n_events);
  }
  if (tot_n_protons) {
  hproton_merged_not_merged->Scale(1./tot_n_protons);
  h_fraction_pdgs_not_tracked_not_clustered->SetBinContent( 3 , h_fraction_pdgs_not_tracked_not_clustered->GetBinContent(3)/float(tot_n_protons) ); //>20MeV
  }

  if (n_all_protons)
  h_fraction_pdgs_not_tracked_not_clustered->SetBinContent( 2 , h_fraction_pdgs_not_tracked_not_clustered->GetBinContent(2)/float(n_all_protons) ); //all protons

  if (low_protons)
  h_fraction_pdgs_not_tracked_not_clustered->SetBinContent( 4 , h_fraction_pdgs_not_tracked_not_clustered->GetBinContent(4)/float(low_protons) ); //<20MeV

  h_fraction_pdgs_not_tracked_not_clustered->SetBinError( 2 , h_fraction_pdgs_not_tracked_not_clustered->GetBinContent(2)/sqrt(h_fraction_pdgs_not_tracked_not_clustered->GetEntries()) ); //all protons
  h_fraction_pdgs_not_tracked_not_clustered->SetBinError( 3 , h_fraction_pdgs_not_tracked_not_clustered->GetBinContent(3)/sqrt(h_fraction_pdgs_not_tracked_not_clustered->GetEntries()) ); //>20MeV
  h_fraction_pdgs_not_tracked_not_clustered->SetBinError( 4 , h_fraction_pdgs_not_tracked_not_clustered->GetBinContent(4)/sqrt(h_fraction_pdgs_not_tracked_not_clustered->GetEntries()) ); //<20MeV
}//end scale plots
   

}//end namespace

#endif
