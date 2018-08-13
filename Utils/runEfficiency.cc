#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "TEfficiency.h"
#include <string>
#include <vector>
#include <map>
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TH2.h"
#include "TApplication.h"
#include "TKey.h"

//TApplication app("app",0, NULL);

class qualifier
{   
	public:
	  qualifier(){
		  version="";
		  dir="";
		  histo="";
		  counter=0;
	  }
	  qualifier(std::string v, std::string d, std::string h, int c = 0){
		  version=v;
		  dir=d;
		  histo=h;
		  counter=c;
	  }
	  ~qualifier() {};
		

	  void Set(std::string a, std::string b, std::string c){ version=a; dir=b; histo=c; counter=0; }      
	  void setVersion(std::string xin){ version=xin; }      
	  void setHisto(std::string xin){ histo=xin; }      
	  void setDir(std::string xin){ dir=xin; }      
	  std::string getVersion() const { return version; }      
	  std::string getHisto() const { return histo; }      
	  std::string getDir() const { return dir; }      
	  int getCounter() const { return counter; }      

	  bool operator==(const qualifier& right)const{
		return (version == right.getVersion() && dir == right.getDir() && histo == right.getHisto() ) ;
	  }
	  bool operator<(const qualifier& right)const{
		return ( counter < right.getCounter() ) ;
	  }
	  bool operator>(const qualifier& right)const{
		return ( counter > right.getCounter() ) ;
	  }
	  friend std::ostream& operator<< (std::ostream& stream, const qualifier& matrix){
		  stream << "Version=" << matrix.getVersion() << "; Dir=" << matrix.getDir() << "; Histo=" << matrix.getHisto() << "; Counter=" << matrix.getCounter();
		  return stream;
	  }
		  
	private:
		std::string version;
		std::string dir;
		std::string histo;
		int counter;
};

TGraphAsymmErrors* find_efficiency(  std::map<qualifier, TEfficiency> saved_tot, qualifier qq) {
	for (auto const &saved : saved_tot) {
		//std::cout << saved.first << std::endl;
		//std::cout << qq << std::endl;
		if (saved.first==qq)
			return saved.second.CreateGraph();
	}
	 std::cout << "Error! map element not found for qualifier " << qq << std::endl;
	 exit(-1);
}

TH1D* find_plot (  std::map<qualifier, TH1D> saved_tot, qualifier qq) {
	for (auto const &saved : saved_tot) {
		if (saved.first==qq){
			TH1D* htemp = (TH1D*) saved.second.Clone();
			return htemp;
		}
	}
	 std::cout << "Error! map element not found for qualifier " << qq << std::endl;
	 exit(-1);
}
			
TH1D* get_tail_to_tot_2d( TFile* inFile, std::string file, std::string dir, bool is_merged ) {
	std::string name;
	if (is_merged)
		name = "dqdx_tailtotot_length_merged";
	else
		name = "dqdx_tailtotot_length_not_merged";
		
	TH2D* tail = (TH2D*) inFile->Get( std::string(dir+"/"+name).c_str() )->Clone();	
	TH1D* service = (TH1D*) tail->ProjectionX()->Clone(std::string( file + "_" + dir+"_"+name).c_str());
	service->Reset();
	for (int i = 0; i<=tail->GetNbinsX(); i++) {
	service->SetBinContent(i, tail->ProjectionY("tmp", i, i)->GetMean());
	//std::cout << tail->ProjectionY("tmp", i, i)->GetBinContent(1) << std::endl;
	}
	//TCanvas c1;
	//service->Draw();
	//app.Run(true);
	return service;
}




int main () {

TFile* outFile = new TFile("eff_test.root","recreate");
TDirectory* dir_eff = outFile->mkdir("dir_efficiency");
TDirectory* dir_plot = outFile->mkdir("dir_plots");
TEfficiency* pEff = 0;
std::vector<std::string> numerator;
std::vector<std::string> denominator;
std::vector<std::string> plots;
std::vector<std::string> dirs; std::vector<std::string> versions;
std::map<std::string, std::string> files;
std::map<qualifier, TEfficiency> saved_efficiency;
std::map<qualifier, TH1D> saved_plot;
files["develop"]="../v06_74_01/build_slf6.x86_64/numu_recobench_output.root";
files["mcc8"]="../new_truth_reco/build_slf6.x86_64/numu_recobench_mcc8.root";

dirs.push_back("pandoraNu"); versions.push_back("mcc8");
//dirs.push_back("pmtrack"); versions.push_back("mcc8");
dirs.push_back("pandoraNu"); versions.push_back("develop");
dirs.push_back("trajcluster"); versions.push_back("develop");

numerator.push_back("proton_kinE"); denominator.push_back("proton_kinE_all");
numerator.push_back("proton_nhits"); denominator.push_back("proton_nhits_all");
numerator.push_back("proton_l"); denominator.push_back("proton_l_all");

numerator.push_back("proton_l_tracked_angle1"); denominator.push_back("proton_l_all_angle1");
numerator.push_back("proton_l_tracked_angle2"); denominator.push_back("proton_l_all_angle2");
numerator.push_back("proton_l_tracked_angle3"); denominator.push_back("proton_l_all_angle3");
numerator.push_back("proton_kinE_tracked_angle1"); denominator.push_back("proton_kinE_all_angle1");
numerator.push_back("proton_kinE_tracked_angle2"); denominator.push_back("proton_kinE_all_angle2");
numerator.push_back("proton_kinE_tracked_angle3"); denominator.push_back("proton_kinE_all_angle3");
numerator.push_back("proton_nhits_tracked_angle1"); denominator.push_back("proton_nhits_all_angle1");
numerator.push_back("proton_nhits_tracked_angle2"); denominator.push_back("proton_nhits_all_angle2");
numerator.push_back("proton_nhits_tracked_angle3"); denominator.push_back("proton_nhits_all_angle3");

numerator.push_back("theta_mu_tracked"); denominator.push_back("proton_theta_mu");
numerator.push_back("pmu_end_tracked"); denominator.push_back("pmu_end_not_tracked");

numerator.push_back("shower_proton_kinE"); denominator.push_back("proton_kinE_all");
numerator.push_back("shower_proton_nhits"); denominator.push_back("proton_nhits_all");
numerator.push_back("shower_proton_costheta_muon"); denominator.push_back("proton_theta_mu");

plots.push_back("pmu_end_tracked");
plots.push_back("pmu_end_not_tracked");
plots.push_back("theta_mu_tracked");
plots.push_back("theta_mu_not_tracked");
plots.push_back("n_proton_showers");
plots.push_back("shower_proton_kinE");
plots.push_back("shower_proton_nhits");
plots.push_back("shower_proton_l");
plots.push_back("shower_proton_costheta_muon");
plots.push_back("proton_merged_not_merged");
plots.push_back("muon_pos_res_goodprotons");
plots.push_back("muon_pos_res_badprotons");
plots.push_back("proton_multi_above20MeV");
plots.push_back("proton_multi_below20MeV");
plots.push_back("tail_to_tot_2d_merged");
plots.push_back("tail_to_tot_2d_not_merged");

int cont = 0;
for (auto const &file : files ) {
TFile* inFile = new TFile(file.second.c_str(),"read");
for (auto const &dir : dirs ) {
	std::string vv = versions.at(&dir - &dirs[0]);
	if ( vv != file.first ) continue; //select only matching
	for (auto const &histo : numerator ) { //efficiency plots
		
		TH1* h_pass = (TH1D*) inFile->Get( std::string(dir+"/"+histo.c_str() ).c_str() )->Clone();
		TH1* h_total = (TH1D*) inFile->Get( std::string(dir+"/"+denominator[ &histo - &numerator[0] ]).c_str() )->Clone();
		if (!h_pass || !h_total) {
			std::cout << histo << " not found in file " << file.first << " and folder " << dir <<". Skipping." << std::endl;
			continue;
		}
		if ( histo.find("theta_mu") != std::string::npos ) { //hack
			h_pass->Rebin(20);
			h_total->Rebin(20);
		}
		if ( denominator[ &histo - &numerator[0] ].find("not_tracked") != std::string::npos ) { //hack2
			h_total->Add( h_pass, 1 );
		}
		if ( histo.find("shower") != std::string::npos ) { //hack3
			if  ( histo.find("costheta") != std::string::npos ) {
				h_pass->Rebin(2);
				h_total->Rebin(2);
			} else {
			h_pass->Rebin(20);
			h_total->Rebin(20);
			}
		}
	//h_pass and h_total are valid and consistent histograms
		if( TEfficiency::CheckConsistency( *h_pass, *h_total ) ) {
  		pEff = new TEfficiency( *h_pass, *h_total );
		pEff->SetName( std::string( file.first + "_" + dir+"_"+histo).c_str() );
		pEff->SetTitle( std::string( file.first + "_" + dir+"_"+histo).c_str() );
		// this will write the TEfficiency object to "myfile.root"
		// AND pEff will be attached to the current directory
		qualifier q(file.first, dir, histo, cont++);
		if (saved_efficiency.find(q) != saved_efficiency.end()) { std::cout << "Error! Map element ALREADY EXISTING for qualifier " << q << std::endl; return -1;}
		saved_efficiency[q] = * (TEfficiency*) pEff->Clone(std::string( file.first + "_" + dir+"_"+histo).c_str());
		//std::cout << "Creating object with qualifier " << q << std::endl;
  		dir_eff->cd();
		pEff->Write();
 		} else {
			std::cout << "non consistent " << histo << " " << denominator[ &histo - &numerator[0] ] << std::endl;
		}
	}

	for ( auto const &histo : plots ) { //std plots
		TH1D* h_pass;
		if ( histo.find("tail_to_tot_2d") != std::string::npos ) {
			if ( histo.find("not_merged") != std::string::npos )
				h_pass = get_tail_to_tot_2d( inFile, file.first, dir, 0 );
			else
				h_pass = get_tail_to_tot_2d( inFile, file.first, dir, 1 );
		//TCanvas c1;
	        //h_pass->Draw();
		//app.Run(true);
		}
		else
			h_pass = (TH1D*) inFile->Get( std::string(dir+"/"+histo.c_str() ).c_str() )->Clone();
		
		if (!h_pass) {
			std::cout << histo << " not found in file " << file.first << " and folder " << dir <<". Skipping." << std::endl;
			continue;
		}
		qualifier q(file.first, dir, histo, cont++);
		if (saved_plot.find(q) != saved_plot.end()) { std::cout << "Error! Map element ALREADY EXISTING for qualifier " << q << std::endl; return -1;}
		saved_plot[q] = * (TH1D*) h_pass->Clone(std::string( file.first + "_" + dir+"_"+histo).c_str());
  		dir_plot->cd();
		h_pass->Write();
	}
}
inFile->Close();

}

gStyle->SetOptStat(0);
//now make comparison canvas
//
//FIRST CANVAS
TCanvas c1; c1.SetName("mcc8_develop_pandoraNu_kinE");
qualifier q("develop", "pandoraNu", "proton_kinE");
TGraphAsymmErrors* eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency (PandoraNu); Proton Kinetic Energy (GeV); Efficiency");
eff1->SetName("develop");
eff1->GetXaxis()->SetRangeUser(0,0.4);
eff1->Draw();
qualifier q2("mcc8", "pandoraNu", "proton_kinE");
TGraphAsymmErrors* eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("mcc8");
eff2->SetLineColor(kRed); eff2->Draw("same");
auto legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(eff1,"develop","pl");
legend->AddEntry(eff2,"mcc8","pl");
legend->SetLineWidth(0);
legend->Draw();
c1.SetGridx();
c1.SetGridy();
outFile->cd(); c1.Write();

//SECOND CANVAS
TCanvas c2; c2.SetName("mcc8_develop_pandoraNu_nhits");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_nhits");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency (PandoraNu); nhits (all planes); Efficiency");
eff1->SetName("develop");
eff1->GetXaxis()->SetRangeUser(0,150);
eff1->Draw();
q2.setVersion("mcc8"); q2.setDir("pandoraNu"); q2.setHisto("proton_nhits");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("mcc8");
eff2->SetLineColor(kRed); eff2->Draw("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(eff1,"develop","pl");
legend->AddEntry(eff2,"mcc8","pl");
legend->SetLineWidth(0);
legend->Draw();
c2.SetGridx();
c2.SetGridy();
outFile->cd(); c2.Write();

//THRID CANVAS
TCanvas c3; c3.SetName("mcc8_develop_pandoraNu_l");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_l");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency (PandoraNu); length (cm); Efficiency");
eff1->SetName("develop");
eff1->GetXaxis()->SetRangeUser(0,30);
eff1->Draw();
q2.setVersion("mcc8"); q2.setDir("pandoraNu"); q2.setHisto("proton_l");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("mcc8");
eff2->SetLineColor(kRed); eff2->Draw("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(eff1,"develop","pl");
legend->AddEntry(eff2,"mcc8","pl");
legend->SetLineWidth(0);
legend->Draw();
c3.SetGridx();
c3.SetGridy();
outFile->cd(); c3.Write();

//FOURTH CANVAS - compare different angles in the same version - nhits/develop
TCanvas c4; c4.SetName("develop_pandoraNu_angles_nhits");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_nhits_tracked_angle1");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency for different track angles (develop, pandoraNu); nhits (all planes); Efficiency");
eff1->SetName("angle1");
eff1->GetXaxis()->SetRangeUser(0,60);
eff1->Draw();
q2.setVersion("develop"); q2.setDir("pandoraNu"); q2.setHisto("proton_nhits_tracked_angle2");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("angle2");
eff2->SetLineColor(kRed); eff2->Draw("same");
qualifier q3;
q3.setVersion("develop"); q3.setDir("pandoraNu"); q3.setHisto("proton_nhits_tracked_angle3");
TGraphAsymmErrors* eff3 = find_efficiency( saved_efficiency, q3);
eff3->SetName("angle3");
eff3->SetLineColor(kBlue); eff3->Draw("same");
legend = new TLegend(0.5,0.1,0.8,0.3);
legend->AddEntry(eff1,"#theta<30","pl");
legend->AddEntry(eff2,"#theta>30 && #theta<60","pl");
legend->AddEntry(eff3,"#theta>60 && #theta<90","pl");
legend->SetLineWidth(0);
legend->Draw();
c4.SetGridx();
c4.SetGridy();
outFile->cd(); c4.Write();


//FIFTH CANVAS - compare different angles in the same version - length/develop
TCanvas c5; c5.SetName("develop_pandoraNu_angles_l");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_l_tracked_angle1");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency for different track angles (develop, pandoraNu); Length (cm); Efficiency");
eff1->SetName("angle1");
eff1->GetXaxis()->SetRangeUser(0,10);
eff1->Draw();
q2.setVersion("develop"); q2.setDir("pandoraNu"); q2.setHisto("proton_l_tracked_angle2");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("angle2");
eff2->SetLineColor(kRed); eff2->Draw("same");
q3.setVersion("develop"); q3.setDir("pandoraNu"); q3.setHisto("proton_l_tracked_angle3");
eff3 = find_efficiency( saved_efficiency, q3);
eff3->SetName("angle3");
eff3->SetLineColor(kBlue); eff3->Draw("same");
legend = new TLegend(0.5,0.1,0.8,0.3);
legend->AddEntry(eff1,"#theta<30","pl");
legend->AddEntry(eff2,"#theta>30 && #theta<60","pl");
legend->AddEntry(eff3,"#theta>60 && #theta<90","pl");
legend->SetLineWidth(0);
legend->Draw();
c5.SetGridx();
c5.SetGridy();
outFile->cd(); c5.Write();

//SIXTH CANVAS - compare different angles in the same version - kinE/develop
TCanvas c6; c6.SetName("develop_pandoraNu_angles_kinE");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_kinE_tracked_angle1");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency for different track angles (develop, pandoraNu); Kinetic Energy (GeV); Efficiency");
eff1->SetName("angle1");
eff1->GetXaxis()->SetRangeUser(0,0.15);
eff1->Draw();
q2.setVersion("develop"); q2.setDir("pandoraNu"); q2.setHisto("proton_kinE_tracked_angle2");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("angle2");
eff2->SetLineColor(kRed); eff2->Draw("same");
q3.setVersion("develop"); q3.setDir("pandoraNu"); q3.setHisto("proton_kinE_tracked_angle3");
eff3 = find_efficiency( saved_efficiency, q3);
eff3->SetName("angle3");
eff3->SetLineColor(kBlue); eff3->Draw("same");
legend = new TLegend(0.5,0.1,0.8,0.3);
legend->AddEntry(eff1,"#theta<30","pl");
legend->AddEntry(eff2,"#theta>30 && #theta<60","pl");
legend->AddEntry(eff3,"#theta>60 && #theta<90","pl");
legend->SetLineWidth(0);
legend->Draw();
c6.SetGridx();
c6.SetGridy();
outFile->cd(); c6.Write();


//SEVENTH CANVAS - compare different angles in the same version - nhits/mcc8
TCanvas c7; c7.SetName("mcc8_pandoraNu_angles_nhits");
q.setVersion("mcc8"); q.setDir("pandoraNu"); q.setHisto("proton_nhits_tracked_angle1");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency for different track angles (mcc8, pandoraNu); nhits (all planes); Efficiency");
eff1->SetName("angle1");
eff1->GetXaxis()->SetRangeUser(0,60);
eff1->Draw();
q2.setVersion("mcc8"); q2.setDir("pandoraNu"); q2.setHisto("proton_nhits_tracked_angle2");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("angle2");
eff2->SetLineColor(kRed); eff2->Draw("same");
q3.setVersion("mcc8"); q3.setDir("pandoraNu"); q3.setHisto("proton_nhits_tracked_angle3");
eff3 = find_efficiency( saved_efficiency, q3);
eff3->SetName("angle3");
eff3->SetLineColor(kBlue); eff3->Draw("same");
legend = new TLegend(0.5,0.1,0.8,0.3);
legend->AddEntry(eff1,"#theta<30","pl");
legend->AddEntry(eff2,"#theta>30 && #theta<60","pl");
legend->AddEntry(eff3,"#theta>60 && #theta<90","pl");
legend->SetLineWidth(0);
legend->Draw();
c7.SetGridx();
c7.SetGridy();
outFile->cd(); c7.Write();


//EIGHTH CANVAS - compare different angles in the same version - length/mcc8
TCanvas c8; c8.SetName("mcc8_pandoraNu_angles_l");
q.setVersion("mcc8"); q.setDir("pandoraNu"); q.setHisto("proton_l_tracked_angle1");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency for different track angles (mcc8, pandoraNu); Length (cm); Efficiency");
eff1->SetName("angle1");
eff1->GetXaxis()->SetRangeUser(0,10);
eff1->Draw();
q2.setVersion("mcc8"); q2.setDir("pandoraNu"); q2.setHisto("proton_l_tracked_angle2");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("angle2");
eff2->SetLineColor(kRed); eff2->Draw("same");
q3.setVersion("mcc8"); q3.setDir("pandoraNu"); q3.setHisto("proton_l_tracked_angle3");
eff3 = find_efficiency( saved_efficiency, q3);
eff3->SetName("angle3");
eff3->SetLineColor(kBlue); eff3->Draw("same");
legend = new TLegend(0.5,0.1,0.8,0.3);
legend->AddEntry(eff1,"#theta<30","pl");
legend->AddEntry(eff2,"#theta>30 && #theta<60","pl");
legend->AddEntry(eff3,"#theta>60 && #theta<90","pl");
legend->SetLineWidth(0);
legend->Draw();
c8.SetGridx();
c8.SetGridy();
outFile->cd(); c8.Write();

//NINETH CANVAS - compare different angles in the same version - kinE/mcc8
TCanvas c9; c9.SetName("mcc8_pandoraNu_angles_kinE");
q.setVersion("mcc8"); q.setDir("pandoraNu"); q.setHisto("proton_kinE_tracked_angle1");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency for different track angles (mcc8, pandoraNu); Kinetic Energy (GeV); Efficiency");
eff1->SetName("angle1");
eff1->GetXaxis()->SetRangeUser(0,0.15);
eff1->Draw();
q2.setVersion("mcc8"); q2.setDir("pandoraNu"); q2.setHisto("proton_kinE_tracked_angle2");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("angle2");
eff2->SetLineColor(kRed); eff2->Draw("same");
q3.setVersion("mcc8"); q3.setDir("pandoraNu"); q3.setHisto("proton_kinE_tracked_angle3");
eff3 = find_efficiency( saved_efficiency, q3);
eff3->SetName("angle3");
eff3->SetLineColor(kBlue); eff3->Draw("same");
legend = new TLegend(0.5,0.1,0.8,0.3);
legend->AddEntry(eff1,"#theta<30","pl");
legend->AddEntry(eff2,"#theta>30 && #theta<60","pl");
legend->AddEntry(eff3,"#theta>60 && #theta<90","pl");
legend->SetLineWidth(0);
legend->Draw();
c9.SetGridx();
c9.SetGridy();
outFile->cd(); c9.Write();

//TENTH CANVAS - compare mcc8 and develop in a fixed angle
TCanvas c10; c10.SetName("mcc8_develop_pandoraNu_kinE_angle1");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_kinE_tracked_angle1");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency for #theta<30; Kinetic Energy (GeV); Efficiency");
eff1->SetName("develop");
eff1->GetXaxis()->SetRangeUser(0,0.15);
eff1->Draw();
q2.setVersion("mcc8"); q2.setDir("pandoraNu"); q2.setHisto("proton_kinE_tracked_angle1");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("mcc8");
eff2->SetLineColor(kRed); eff2->Draw("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(eff1,"develop","pl");
legend->AddEntry(eff2,"mcc8","pl");
legend->SetLineWidth(0);
legend->Draw();
c10.SetGridx();
c10.SetGridy();
outFile->cd(); c10.Write();

//ELEVENTH CANVAS - compare mcc8 and develop in a fixed angle
TCanvas c11; c11.SetName("mcc8_develop_pandoraNu_kinE_angle2");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_kinE_tracked_angle2");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency for #theta>30 && #theta<60; Kinetic Energy (GeV); Efficiency");
eff1->SetName("develop");
eff1->GetXaxis()->SetRangeUser(0,0.15);
eff1->Draw();
q2.setVersion("mcc8"); q2.setDir("pandoraNu"); q2.setHisto("proton_kinE_tracked_angle2");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("mcc8");
eff2->SetLineColor(kRed); eff2->Draw("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(eff1,"develop","pl");
legend->AddEntry(eff2,"mcc8","pl");
legend->SetLineWidth(0);
legend->Draw();
c11.SetGridx();
c11.SetGridy();
outFile->cd(); c11.Write();

//TWELFTH CANVAS - compare mcc8 and develop in a fixed angle
TCanvas c12; c12.SetName("mcc8_develop_pandoraNu_kinE_angle3");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_kinE_tracked_angle3");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency for #theta>60 && #theta<90; Kinetic Energy (GeV); Efficiency");
eff1->SetName("develop");
eff1->GetXaxis()->SetRangeUser(0,0.15);
eff1->Draw();
q2.setVersion("mcc8"); q2.setDir("pandoraNu"); q2.setHisto("proton_kinE_tracked_angle3");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("mcc8");
eff2->SetLineColor(kRed); eff2->Draw("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(eff1,"develop","pl");
legend->AddEntry(eff2,"mcc8","pl");
legend->SetLineWidth(0);
legend->Draw();
c12.SetGridx();
c12.SetGridy();
outFile->cd(); c12.Write();


//13th CANVAS - compare mcc8 and develop in a fixed angle
TCanvas c13; c13.SetName("mcc8_develop_pandoraNu_l_angle1");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_l_tracked_angle1");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency for #theta<30; Length (cm); Efficiency");
eff1->SetName("develop");
eff1->GetXaxis()->SetRangeUser(0,10);
eff1->Draw();
q2.setVersion("mcc8"); q2.setDir("pandoraNu"); q2.setHisto("proton_l_tracked_angle1");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("mcc8");
eff2->SetLineColor(kRed); eff2->Draw("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(eff1,"develop","pl");
legend->AddEntry(eff2,"mcc8","pl");
legend->SetLineWidth(0);
legend->Draw();
c13.SetGridx();
c13.SetGridy();
outFile->cd(); c13.Write();

//14th CANVAS - compare mcc8 and develop in a fixed angle
TCanvas c14; c14.SetName("mcc8_develop_pandoraNu_l_angle2");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_l_tracked_angle2");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency for #theta>30 && #theta<60; Length (cm); Efficiency");
eff1->SetName("develop");
eff1->GetXaxis()->SetRangeUser(0,10);
eff1->Draw();
q2.setVersion("mcc8"); q2.setDir("pandoraNu"); q2.setHisto("proton_l_tracked_angle2");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("mcc8");
eff2->SetLineColor(kRed); eff2->Draw("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(eff1,"develop","pl");
legend->AddEntry(eff2,"mcc8","pl");
legend->SetLineWidth(0);
legend->Draw();
c14.SetGridx();
c14.SetGridy();
outFile->cd(); c14.Write();

//15th CANVAS - compare mcc8 and develop in a fixed angle
TCanvas c15; c15.SetName("mcc8_develop_pandoraNu_l_angle3");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_l_tracked_angle3");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency for #theta>60 && #theta<90; Length (cm); Efficiency");
eff1->SetName("develop");
eff1->GetXaxis()->SetRangeUser(0,10);
eff1->Draw();
q2.setVersion("mcc8"); q2.setDir("pandoraNu"); q2.setHisto("proton_l_tracked_angle3");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("mcc8");
eff2->SetLineColor(kRed); eff2->Draw("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(eff1,"develop","pl");
legend->AddEntry(eff2,"mcc8","pl");
legend->SetLineWidth(0);
legend->Draw();
c15.SetGridx();
c15.SetGridy();
outFile->cd(); c15.Write();

//16th CANVAS - compare mcc8 and develop in a fixed angle
TCanvas c16; c16.SetName("mcc8_develop_pandoraNu_nhits_angle1");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_nhits_tracked_angle1");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency for #theta<30; nhits (all planes); Efficiency");
eff1->SetName("develop");
eff1->GetXaxis()->SetRangeUser(0,60);
eff1->Draw();
q2.setVersion("mcc8"); q2.setDir("pandoraNu"); q2.setHisto("proton_nhits_tracked_angle1");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("mcc8");
eff2->SetLineColor(kRed); eff2->Draw("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(eff1,"develop","pl");
legend->AddEntry(eff2,"mcc8","pl");
legend->SetLineWidth(0);
legend->Draw();
c16.SetGridx();
c16.SetGridy();
outFile->cd(); c16.Write();

//17th CANVAS - compare mcc8 and develop in a fixed angle
TCanvas c17; c17.SetName("mcc8_develop_pandoraNu_nhits_angle2");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_nhits_tracked_angle2");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency for #theta>30 && #theta<60; nhits (all planes); Efficiency");
eff1->SetName("develop");
eff1->GetXaxis()->SetRangeUser(0,60);
eff1->Draw();
q2.setVersion("mcc8"); q2.setDir("pandoraNu"); q2.setHisto("proton_nhits_tracked_angle2");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("mcc8");
eff2->SetLineColor(kRed); eff2->Draw("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(eff1,"develop","pl");
legend->AddEntry(eff2,"mcc8","pl");
legend->SetLineWidth(0);
legend->Draw();
c17.SetGridx();
c17.SetGridy();
outFile->cd(); c17.Write();

//18th CANVAS - compare mcc8 and develop in a fixed angle
TCanvas c18; c18.SetName("mcc8_develop_pandoraNu_nhits_angle3");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_nhits_tracked_angle3");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency for #theta>60 && #theta<90; nhits (all planes); Efficiency");
eff1->SetName("develop");
eff1->GetXaxis()->SetRangeUser(0,60);
eff1->Draw();
q2.setVersion("mcc8"); q2.setDir("pandoraNu"); q2.setHisto("proton_nhits_tracked_angle3");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("mcc8");
eff2->SetLineColor(kRed); eff2->Draw("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(eff1,"develop","pl");
legend->AddEntry(eff2,"mcc8","pl");
legend->SetLineWidth(0);
legend->Draw();
c18.SetGridx();
c18.SetGridy();
outFile->cd(); c18.Write();

//19th CANVAS - compare different angles in the same version - nhits/develop
TCanvas c19; c19.SetName("develop_trajcluster_angles_nhits");
q.setVersion("develop"); q.setDir("trajcluster"); q.setHisto("proton_nhits_tracked_angle1");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency for different track angles (develop, trajcluster); nhits (all planes); Efficiency");
eff1->SetName("angle1");
eff1->GetXaxis()->SetRangeUser(0,60);
eff1->Draw();
q2.setVersion("develop"); q2.setDir("trajcluster"); q2.setHisto("proton_nhits_tracked_angle2");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("angle2");
eff2->SetLineColor(kRed); eff2->Draw("same");
q3.setVersion("develop"); q3.setDir("trajcluster"); q3.setHisto("proton_nhits_tracked_angle3");
eff3 = find_efficiency( saved_efficiency, q3);
eff3->SetName("angle3");
eff3->SetLineColor(kBlue); eff3->Draw("same");
legend = new TLegend(0.5,0.1,0.8,0.3);
legend->AddEntry(eff1,"#theta<30","pl");
legend->AddEntry(eff2,"#theta>30 && #theta<60","pl");
legend->AddEntry(eff3,"#theta>60 && #theta<90","pl");
legend->SetLineWidth(0);
legend->Draw();
c19.SetGridx();
c19.SetGridy();
outFile->cd(); c19.Write();


//20th CANVAS - compare different angles in the same version - length/develop
TCanvas c20; c20.SetName("develop_trajcluster_angles_l");
q.setVersion("develop"); q.setDir("trajcluster"); q.setHisto("proton_l_tracked_angle1");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency for different track angles (develop, trajcluster); Length (cm); Efficiency");
eff1->SetName("angle1");
eff1->GetXaxis()->SetRangeUser(0,10);
eff1->Draw();
q2.setVersion("develop"); q2.setDir("trajcluster"); q2.setHisto("proton_l_tracked_angle2");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("angle2");
eff2->SetLineColor(kRed); eff2->Draw("same");
q3.setVersion("develop"); q3.setDir("trajcluster"); q3.setHisto("proton_l_tracked_angle3");
eff3 = find_efficiency( saved_efficiency, q3);
eff3->SetName("angle3");
eff3->SetLineColor(kBlue); eff3->Draw("same");
legend = new TLegend(0.5,0.1,0.8,0.3);
legend->AddEntry(eff1,"#theta<30","pl");
legend->AddEntry(eff2,"#theta>30 && #theta<60","pl");
legend->AddEntry(eff3,"#theta>60 && #theta<90","pl");
legend->SetLineWidth(0);
legend->Draw();
c20.SetGridx();
c20.SetGridy();
outFile->cd(); c20.Write();

//21st CANVAS - compare different angles in the same version - kinE/develop
TCanvas c21; c21.SetName("develop_trajcluster_angles_kinE");
q.setVersion("develop"); q.setDir("trajcluster"); q.setHisto("proton_kinE_tracked_angle1");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency for different track angles (develop); Kinetic Energy (GeV); Efficiency");
eff1->SetName("angle1");
eff1->GetXaxis()->SetRangeUser(0,0.15);
eff1->Draw();
q2.setVersion("develop"); q2.setDir("trajcluster"); q2.setHisto("proton_kinE_tracked_angle2");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("angle2");
eff2->SetLineColor(kRed); eff2->Draw("same");
q3.setVersion("develop"); q3.setDir("trajcluster"); q3.setHisto("proton_kinE_tracked_angle3");
eff3 = find_efficiency( saved_efficiency, q3);
eff3->SetName("angle3");
eff3->SetLineColor(kBlue); eff3->Draw("same");
legend = new TLegend(0.5,0.1,0.8,0.3);
legend->AddEntry(eff1,"#theta<30","pl");
legend->AddEntry(eff2,"#theta>30 && #theta<60","pl");
legend->AddEntry(eff3,"#theta>60 && #theta<90","pl");
legend->SetLineWidth(0);
legend->Draw();
c21.SetGridx();
c21.SetGridy();
outFile->cd(); c21.Write();

//22nd CANVAS - compare trajcluster and pandoraNu in develop in a fixed angle
TCanvas c22; c22.SetName("develop_pandoraNu_trajcluster_nhits_angle1");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_nhits_tracked_angle1");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency for #theta<30 (develop); nhits (all planes); Efficiency");
eff1->SetName("pandoraNu");
eff1->GetXaxis()->SetRangeUser(0,60);
eff1->Draw();
q2.setVersion("develop"); q2.setDir("trajcluster"); q2.setHisto("proton_nhits_tracked_angle1");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("trajcluster");
eff2->SetLineColor(kRed); eff2->Draw("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(eff1,"pandoraNu","pl");
legend->AddEntry(eff2,"trajcluster","pl");
legend->SetLineWidth(0);
legend->Draw();
c22.SetGridx();
c22.SetGridy();
outFile->cd(); c22.Write();

//23rd CANVAS - compare trajcluster and pandoraNu in develop in a fixed angle
TCanvas c23; c23.SetName("develop_pandoraNu_trajcluster_nhits_angle2");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_nhits_tracked_angle2");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency for #theta>30 && #theta<60 (develop); nhits (all planes); Efficiency");
eff1->SetName("pandoraNu");
eff1->GetXaxis()->SetRangeUser(0,60);
eff1->Draw();
q2.setVersion("develop"); q2.setDir("trajcluster"); q2.setHisto("proton_nhits_tracked_angle2");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("trajcluster");
eff2->SetLineColor(kRed); eff2->Draw("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(eff1,"pandoraNu","pl");
legend->AddEntry(eff2,"trajcluster","pl");
legend->SetLineWidth(0);
legend->Draw();
c23.SetGridx();
c23.SetGridy();
outFile->cd(); c23.Write();

//24th CANVAS - compare trajcluster and pandoraNu in develop in a fixed angle
TCanvas c24; c24.SetName("develop_pandoraNu_trajcluster_nhits_angle3");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_nhits_tracked_angle3");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency for #theta>60 && #theta<90 (develop); nhits (all planes); Efficiency");
eff1->SetName("pandoraNu");
eff1->GetXaxis()->SetRangeUser(0,60);
eff1->Draw();
q2.setVersion("develop"); q2.setDir("trajcluster"); q2.setHisto("proton_nhits_tracked_angle3");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("trajcluster");
eff2->SetLineColor(kRed); eff2->Draw("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(eff1,"pandoraNu","pl");
legend->AddEntry(eff2,"trajcluster","pl");
legend->SetLineWidth(0);
legend->Draw();
c24.SetGridx();
c24.SetGridy();
outFile->cd(); c24.Write();

//25th CANVAS - compare trajcluster and pandoraNu in develop in a fixed angle
TCanvas c25; c25.SetName("develop_pandoraNu_trajcluster_l_angle1");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_l_tracked_angle1");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency for #theta<30 (develop); Length (cm); Efficiency");
eff1->SetName("pandoraNu");
eff1->GetXaxis()->SetRangeUser(0,10);
eff1->Draw();
q2.setVersion("develop"); q2.setDir("trajcluster"); q2.setHisto("proton_l_tracked_angle1");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("trajcluster");
eff2->SetLineColor(kRed); eff2->Draw("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(eff1,"pandoraNu","pl");
legend->AddEntry(eff2,"trajcluster","pl");
legend->SetLineWidth(0);
legend->Draw();
c25.SetGridx();
c25.SetGridy();
outFile->cd(); c25.Write();

//26th CANVAS - compare trajcluster and pandoraNu in develop in a fixed angle
TCanvas c26; c26.SetName("develop_pandoraNu_trajcluster_l_angle2");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_l_tracked_angle2");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency for #theta>30 && #theta<60 (develop); Length (cm); Efficiency");
eff1->SetName("pandoraNu");
eff1->GetXaxis()->SetRangeUser(0,10);
eff1->Draw();
q2.setVersion("develop"); q2.setDir("trajcluster"); q2.setHisto("proton_l_tracked_angle2");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("trajcluster");
eff2->SetLineColor(kRed); eff2->Draw("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(eff1,"pandoraNu","pl");
legend->AddEntry(eff2,"trajcluster","pl");
legend->SetLineWidth(0);
legend->Draw();
c26.SetGridx();
c26.SetGridy();
outFile->cd(); c26.Write();

//27th CANVAS - compare trajcluster and pandoraNu in develop in a fixed angle
TCanvas c27; c27.SetName("develop_pandoraNu_trajcluster_l_angle3");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_l_tracked_angle3");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency for #theta>60 && #theta<90 (develop); Length (cm); Efficiency");
eff1->SetName("pandoraNu");
eff1->GetXaxis()->SetRangeUser(0,10);
eff1->Draw();
q2.setVersion("develop"); q2.setDir("trajcluster"); q2.setHisto("proton_l_tracked_angle3");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("trajcluster");
eff2->SetLineColor(kRed); eff2->Draw("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(eff1,"pandoraNu","pl");
legend->AddEntry(eff2,"trajcluster","pl");
legend->SetLineWidth(0);
legend->Draw();
c27.SetGridx();
c27.SetGridy();
outFile->cd(); c27.Write();

//28th CANVAS - compare trajcluster and pandoraNu in develop in a fixed angle
TCanvas c28; c28.SetName("develop_pandoraNu_trajcluster_kinE_angle1");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_kinE_tracked_angle1");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency for #theta<30 (develop); Kinetic Energy (GeV); Efficiency");
eff1->SetName("pandoraNu");
eff1->GetXaxis()->SetRangeUser(0,0.15);
eff1->Draw();
q2.setVersion("develop"); q2.setDir("trajcluster"); q2.setHisto("proton_kinE_tracked_angle1");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("trajcluster");
eff2->SetLineColor(kRed); eff2->Draw("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(eff1,"pandoraNu","pl");
legend->AddEntry(eff2,"trajcluster","pl");
legend->SetLineWidth(0);
legend->Draw();
c28.SetGridx();
c28.SetGridy();
outFile->cd(); c28.Write();

//29th CANVAS - compare trajcluster and pandoraNu in develop in a fixed angle
TCanvas c29; c29.SetName("develop_pandoraNu_trajcluster_kinE_angle2");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_kinE_tracked_angle2");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency for #theta>30 && #theta<60 (develop); Kinetic Energy (GeV); Efficiency");
eff1->SetName("pandoraNu");
eff1->GetXaxis()->SetRangeUser(0,0.15);
eff1->Draw();
q2.setVersion("develop"); q2.setDir("trajcluster"); q2.setHisto("proton_kinE_tracked_angle2");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("trajcluster");
eff2->SetLineColor(kRed); eff2->Draw("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(eff1,"pandoraNu","pl");
legend->AddEntry(eff2,"trajcluster","pl");
legend->SetLineWidth(0);
legend->Draw();
c29.SetGridx();
c29.SetGridy();
outFile->cd(); c29.Write();

//30th CANVAS - compare trajcluster and pandoraNu in develop in a fixed angle
TCanvas c30; c30.SetName("develop_pandoraNu_trajcluster_kinE_angle3");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_kinE_tracked_angle3");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency for #theta>60 && #theta<90 (develop); Kinetic Energy (GeV); Efficiency");
eff1->SetName("pandoraNu");
eff1->GetXaxis()->SetRangeUser(0,0.15);
eff1->Draw();
q2.setVersion("develop"); q2.setDir("trajcluster"); q2.setHisto("proton_kinE_tracked_angle3");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("trajcluster");
eff2->SetLineColor(kRed); eff2->Draw("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(eff1,"pandoraNu","pl");
legend->AddEntry(eff2,"trajcluster","pl");
legend->SetLineWidth(0);
legend->Draw();
c30.SetGridx();
c30.SetGridy();
outFile->cd(); c30.Write();

//31st CANVAS - compare trajcluster/pandoraNu in develop - length
TCanvas c31; c31.SetName("develop_pandoraNu_trajcluster_l");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_l");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency (develop); Length (cm); Efficiency");
eff1->SetName("pandora");
eff1->GetXaxis()->SetRangeUser(0,10);
eff1->Draw();
q2.setVersion("develop"); q2.setDir("trajcluster"); q2.setHisto("proton_l");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("trajcluster");
eff2->SetLineColor(kRed); eff2->Draw("same");
legend = new TLegend(0.5,0.1,0.8,0.3);
legend->AddEntry(eff1,"pandoraNu","pl");
legend->AddEntry(eff2,"trajcluster","pl");
legend->SetLineWidth(0);
legend->Draw();
c31.SetGridx();
c31.SetGridy();
outFile->cd(); c31.Write();

//32nd CANVAS - compare trajcluster/pandoraNu in develop - kinE
TCanvas c32; c32.SetName("develop_pandoraNu_trajcluster_kinE");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_kinE");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency (develop); Kinetic Energy (GeV); Efficiency");
eff1->SetName("pandora");
eff1->GetXaxis()->SetRangeUser(0,0.15);
eff1->Draw();
q2.setVersion("develop"); q2.setDir("trajcluster"); q2.setHisto("proton_kinE");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("trajcluster");
eff2->SetLineColor(kRed); eff2->Draw("same");
legend = new TLegend(0.5,0.1,0.8,0.3);
legend->AddEntry(eff1,"pandoraNu","pl");
legend->AddEntry(eff2,"trajcluster","pl");
legend->SetLineWidth(0);
legend->Draw();
c32.SetGridx();
c32.SetGridy();
outFile->cd(); c32.Write();

//33rd CANVAS - compare trajcluster/pandoraNu in develop - nhits
TCanvas c33; c33.SetName("develop_pandoraNu_trajcluster_nhits");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_nhits");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency (develop); nhits (all planes); Efficiency");
eff1->SetName("pandora");
eff1->GetXaxis()->SetRangeUser(0,60);
eff1->Draw();
q2.setVersion("develop"); q2.setDir("trajcluster"); q2.setHisto("proton_nhits");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("trajcluster");
eff2->SetLineColor(kRed); eff2->Draw("same");
legend = new TLegend(0.5,0.1,0.8,0.3);
legend->AddEntry(eff1,"pandoraNu","pl");
legend->AddEntry(eff2,"trajcluster","pl");
legend->SetLineWidth(0);
legend->Draw();
c33.SetGridx();
c33.SetGridy();
outFile->cd(); c33.Write();

//34th CANVAS - compare trajcluster/pandoraNu in develop - theta mu
TCanvas c34; c34.SetName("develop_pandoraNu_trajcluster_theta_mu");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("theta_mu_tracked");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency (develop); cos #theta_{#mu p}; Efficiency");
eff1->SetName("pandora");
eff1->GetXaxis()->SetRangeUser(-1,1);
eff1->Draw();
q2.setVersion("develop"); q2.setDir("trajcluster"); q2.setHisto("theta_mu_tracked");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("trajcluster");
eff2->SetLineColor(kRed); eff2->Draw("same");
legend = new TLegend(0.5,0.1,0.8,0.3);
legend->AddEntry(eff1,"pandoraNu","pl");
legend->AddEntry(eff2,"trajcluster","pl");
legend->SetLineWidth(0);
legend->Draw();
c34.SetGridx();
c34.SetGridy();
outFile->cd(); c34.Write();

//35th CANVAS - compare mcc8 and develop - theta mu
TCanvas c35; c35.SetName("mcc8_develop_pandoraNu_theta_mu");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("theta_mu_tracked");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency (develop); cos #theta_{#mu p}; Efficiency");
eff1->SetName("develop");
eff1->GetXaxis()->SetRangeUser(-1,1);
eff1->Draw();
q2.setVersion("mcc8"); q2.setDir("pandoraNu"); q2.setHisto("theta_mu_tracked");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("mcc8");
eff2->SetLineColor(kRed); eff2->Draw("same");
legend = new TLegend(0.5,0.1,0.8,0.3);
legend->AddEntry(eff1,"develop","pl");
legend->AddEntry(eff2,"mcc8","pl");
legend->SetLineWidth(0);
legend->Draw();
c35.SetGridx();
c35.SetGridy();
outFile->cd(); c35.Write();

//36th CANVAS - compare trajcluster/pandoraNu in develop - pmu
TCanvas c36; c36.SetName("develop_pandoraNu_trajcluster_pmu");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("pmu_end_tracked");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency (develop); d (cm); Efficiency");
eff1->SetName("pandora");
eff1->GetXaxis()->SetRangeUser(0,15);
eff1->Draw();
q2.setVersion("develop"); q2.setDir("trajcluster"); q2.setHisto("pmu_end_tracked");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("trajcluster");
eff2->SetLineColor(kRed); eff2->Draw("same");
legend = new TLegend(0.5,0.1,0.8,0.3);
legend->AddEntry(eff1,"pandoraNu","pl");
legend->AddEntry(eff2,"trajcluster","pl");
legend->SetLineWidth(0);
legend->Draw();
c36.SetGridx();
c36.SetGridy();
outFile->cd(); c36.Write();

//37th CANVAS - compare mcc8 and develop - theta mu
TCanvas c37; c37.SetName("mcc8_develop_pandoraNu_pmu");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("pmu_end_tracked");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton tracking efficiency (develop); d (cm); Efficiency");
eff1->SetName("develop");
eff1->GetXaxis()->SetRangeUser(0,30);
eff1->Draw();
q2.setVersion("mcc8"); q2.setDir("pandoraNu"); q2.setHisto("pmu_end_tracked");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("mcc8");
eff2->SetLineColor(kRed); eff2->Draw("same");
legend = new TLegend(0.5,0.1,0.8,0.3);
legend->AddEntry(eff1,"develop","pl");
legend->AddEntry(eff2,"mcc8","pl");
legend->SetLineWidth(0);
legend->Draw();
c37.SetGridx();
c37.SetGridy();
outFile->cd(); c37.Write();

//38th CANVAS
TCanvas c38; c38.SetName("mcc8_develop_pandoraNu_shower_proton_kinE");
q.Set("develop", "pandoraNu", "shower_proton_kinE");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton \"showering\" efficiency; Proton Kinetic Energy (GeV); Efficiency");
eff1->SetName("develop");
eff1->GetXaxis()->SetRangeUser(0,0.4);
eff1->GetYaxis()->SetRangeUser(0,0.3);
eff1->Draw();
q2.Set("mcc8", "pandoraNu", "shower_proton_kinE");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("mcc8");
eff2->SetLineColor(kRed); eff2->Draw("same");
q3.Set("develop", "trajcluster", "shower_proton_kinE");
eff3 = find_efficiency( saved_efficiency, q3);
eff3->SetName("trajcluster");
eff3->SetLineColor(kBlue); eff3->Draw("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(eff1,"develop - pandoraNu","pl");
legend->AddEntry(eff2,"mcc8 - pandoraNu","pl");
legend->AddEntry(eff3,"develop - trajcluster","pl");
legend->SetLineWidth(0);
legend->Draw();
c38.SetGridx();
c38.SetGridy();
outFile->cd(); c38.Write();

//39th CANVAS
TCanvas c39; c39.SetName("mcc8_develop_pandoraNu_shower_proton_nhits");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("shower_proton_nhits");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton \"showering\" efficiency; nhits (all planes); Efficiency");
eff1->SetName("develop");
eff1->GetXaxis()->SetRangeUser(0,150);
eff1->GetYaxis()->SetRangeUser(0,0.3);
eff1->Draw();
q2.setVersion("mcc8"); q2.setDir("pandoraNu"); q2.setHisto("shower_proton_nhits");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetName("mcc8");
eff2->SetLineColor(kRed); eff2->Draw("same");
q3.setVersion("develop"); q3.setDir("trajcluster"); q3.setHisto("shower_proton_nhits");
eff3 = find_efficiency( saved_efficiency, q3);
eff3->SetName("trajcluster");
eff3->SetLineColor(kBlue); eff3->Draw("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(eff1,"develop - pandoraNu","pl");
legend->AddEntry(eff2,"mcc8 - pandoraNu","pl");
legend->AddEntry(eff3,"develop - trajcluster","pl");
legend->SetLineWidth(0);
legend->Draw();
c39.SetGridx();
c39.SetGridy();
outFile->cd(); c39.Write();

//40th CANVAS - starting with the plots
TCanvas c40; c40.SetName("develop_pandoraNu_theta_mu");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("theta_mu_tracked");
TH1D* h1 = find_plot( saved_plot, q);
h1->SetTitle("Proton muon angle distribution (develop, pandoraNu); cos #theta_{#mu p}; Events");
h1->SetName("tracked");
h1->Rebin(10);
h1->GetXaxis()->SetRangeUser(-1,1);
//h1->GetYaxis()->SetRangeUser(0,900);
h1->Draw();
q2.setVersion("develop"); q2.setDir("pandoraNu"); q2.setHisto("theta_mu_not_tracked");
TH1D* h2 = find_plot( saved_plot, q2);
h2->SetName("not_tracked");
h2->Rebin(10);
h2->SetLineColor(kRed); h2->Draw("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(h1,"tracked protons","pl");
legend->AddEntry(h2,"not tracked protons","pl");
legend->SetLineWidth(0);
legend->Draw();
c40.SetGridx();
c40.SetGridy();
outFile->cd(); c40.Write();

//41st CANVAS - starting with the plots
TCanvas c41; c41.SetName("develop_pandoraNu_pmu");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("pmu_end_tracked");
h1 = find_plot( saved_plot, q);
h1->SetTitle("Proton \"d\" distance distribution (develop, pandoraNu); d (cm); Events");
h1->SetName("tracked");
q2.setVersion("develop"); q2.setDir("pandoraNu"); q2.setHisto("pmu_end_not_tracked");
h2 = find_plot( saved_plot, q2);
h2->SetName("not_tracked");
h2->SetTitle("Proton \"d\" distance distribution (develop, pandoraNu); d (cm); Events");
h2->GetXaxis()->SetRangeUser(0,50);
h2->SetLineColor(kRed); h2->Draw("");
h1->Draw("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(h1,"tracked protons","pl");
legend->AddEntry(h2,"not tracked protons","pl");
legend->SetLineWidth(0);
legend->Draw();
c41.SetGridx();
c41.SetGridy();
outFile->cd(); c41.Write();

//42th CANVAS - still shower efficiency
TCanvas c42; c42.SetName("mcc8_develop_pandoraNu_shower_proton_costheta_mu");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("shower_proton_costheta_muon");
eff1 = find_efficiency( saved_efficiency, q);
eff1->SetTitle("Proton \"showering\" efficiency; cos #theta_{#mu p}; Efficiency");
eff1->SetName("develop");
eff1->GetXaxis()->SetRangeUser(-1,1);
q2.setVersion("mcc8"); q2.setDir("pandoraNu"); q2.setHisto("shower_proton_costheta_muon");
eff2 = find_efficiency( saved_efficiency, q2);
eff2->SetTitle("Proton \"showering\" efficiency; cos #theta_{#mu p}; Efficiency");
eff2->SetName("mcc8");
eff2->SetLineColor(kRed); 
//q3.setVersion("develop"); q3.setDir("trajcluster"); q3.setHisto("shower_proton_costheta_muon");
//eff3 = find_efficiency( saved_efficiency, q3);
//eff3->SetName("trajcluster");
//eff3->SetLineColor(kBlue); 
eff2->GetXaxis()->SetRangeUser(-1,1);
eff2->GetYaxis()->SetRangeUser(0,0.22);
eff2->Draw();
//eff3->Draw("same");
eff1->Draw("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(eff1,"develop - pandoraNu","pl");
legend->AddEntry(eff2,"mcc8 - pandoraNu","pl");
//legend->AddEntry(eff3,"develop - trajcluster","pl");
legend->SetLineWidth(0);
legend->Draw();
c42.SetGridx();
c42.SetGridy();
outFile->cd(); c42.Write();

//43rd CANVAS - starting with the plots
TCanvas c43; c43.SetName("develop_pandoraNu_n_proton_showers");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("n_proton_showers");
h1 = find_plot( saved_plot, q);
h1->SetTitle("Probability of N proton showers in an event (pandoraNu); N; Probability");
h1->SetName("develop");
h1->GetXaxis()->SetRangeUser(0,10);
h1->DrawNormalized();
q2.setVersion("mcc8"); q2.setDir("pandoraNu"); q2.setHisto("n_proton_showers");
h2 = find_plot( saved_plot, q2);
h2->SetName("mcc8");
h2->SetLineColor(kRed); h2->DrawNormalized("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(h1,"develop - pandoraNu","pl");
legend->AddEntry(h2,"mcc8 - pandoraNu","pl");
legend->SetLineWidth(0);
legend->Draw();
c43.SetGridx();
c43.SetGridy();
outFile->cd(); c43.Write();

//44th CANVAS - starting with the plots
TCanvas c44; c44.SetName("develop_pandoraNu_trajcluster_proton_total_efficiency");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_merged_not_merged");
h1 = find_plot( saved_plot, q);
h1->SetTitle("Global tracking efficiency; 0=tracked, 1=not tracked; Probability");
h1->SetName("develop");
h1->GetXaxis()->SetRangeUser(0,1.);
h1->GetYaxis()->SetRangeUser(0.35,0.65);
h1->DrawNormalized();
q2.setVersion("mcc8"); q2.setDir("pandoraNu"); q2.setHisto("proton_merged_not_merged");
h2 = find_plot( saved_plot, q2);
h2->SetName("mcc8");
h2->SetLineColor(kRed); h2->DrawNormalized("same");
q3.Set("develop","trajcluster","proton_merged_not_merged");
TH1D* h3 = find_plot(saved_plot, q3);
h3->SetName("trajcluster");
h3->SetLineColor(kBlack); h3->DrawNormalized("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(h1,"develop - pandoraNu","pl");
legend->AddEntry(h2,"mcc8 - pandoraNu","pl");
legend->AddEntry(h3,"develop - trajcluster","pl");
legend->SetLineWidth(0);
legend->Draw();
c44.SetGridx();
c44.SetGridy();
outFile->cd(); c44.Write();

//45th CANVAS - starting with the plots
TCanvas c45; c45.SetName("develop_pandoraNu_trajcluster_muon_posres_good");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("muon_pos_res_goodprotons");
h1 = find_plot( saved_plot, q);
h1->SetTitle("Muon position resolution for events w/o lost protons > 20MeV; Distance from true vertex (cm); Events");
h1->SetName("develop");
h1->GetXaxis()->SetRangeUser(0,10);
h1->DrawNormalized();
q2.setVersion("mcc8"); q2.setDir("pandoraNu"); q2.setHisto("muon_pos_res_goodprotons");
h2 = find_plot( saved_plot, q2);
h2->SetName("mcc8");
h2->SetLineColor(kRed); h2->DrawNormalized("same");
q3.Set("develop","trajcluster","muon_pos_res_goodprotons");
h3 = find_plot(saved_plot, q3);
h3->SetName("trajcluster");
h3->SetLineColor(kBlack); h3->DrawNormalized("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(h1,"develop - pandoraNu","pl");
legend->AddEntry(h2,"mcc8 - pandoraNu","pl");
legend->AddEntry(h3,"develop - trajcluster","pl");
legend->SetLineWidth(0);
legend->Draw();
c45.SetGridx();
c45.SetGridy();
outFile->cd(); c45.Write();

//46th CANVAS - starting with the plots
TCanvas c46; c46.SetName("develop_pandoraNu_trajcluster_muon_posres_bad");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("muon_pos_res_badprotons");
h1 = find_plot( saved_plot, q);
h1->SetTitle("Muon position resolution for events w/o lost protons > 20MeV; Distance from true vertex (cm); Events");
h1->SetName("develop");
h1->GetXaxis()->SetRangeUser(0,10);
h1->DrawNormalized();
q2.setVersion("mcc8"); q2.setDir("pandoraNu"); q2.setHisto("muon_pos_res_badprotons");
h2 = find_plot( saved_plot, q2);
h2->SetName("mcc8");
h2->SetLineColor(kRed); h2->DrawNormalized("same");
q3.Set("develop","trajcluster","muon_pos_res_badprotons");
h3 = find_plot(saved_plot, q3);
h3->SetName("trajcluster");
h3->SetLineColor(kBlack); h3->DrawNormalized("same");
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(h1,"develop - pandoraNu","pl");
legend->AddEntry(h2,"mcc8 - pandoraNu","pl");
legend->AddEntry(h3,"develop - trajcluster","pl");
legend->SetLineWidth(0);
legend->Draw();
c46.SetGridx();
c46.SetGridy();
outFile->cd(); c46.Write();

//47th CANVAS - starting with the plots
TCanvas c47; c47.SetName("proton_multiplicity");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("proton_multi_above20MeV");
h1 = find_plot( saved_plot, q);
h1->SetTitle("Proton multiplicity; N_{p}; Area Normalized");
h1->SetName("develop");
h1->GetXaxis()->SetRangeUser(0,20);
q2.setVersion("mcc8"); q2.setDir("pandoraNu"); q2.setHisto("proton_multi_below20MeV");
h2 = find_plot( saved_plot, q2);
h2->SetTitle("Proton multiplicity; N_{p}; Area Normalized");
h2->SetName("mcc8");
h2->SetLineColor(kRed); h2->SetLineWidth(2); h2->DrawNormalized("hist");
h1->SetLineWidth(2);
h1->DrawNormalized("hist same");
h1->SetLineWidth(1);
h2->SetLineWidth(1);
legend = new TLegend(0.6,0.2,0.8,0.3);
legend->AddEntry(h1,"protons w/ kinE>20MeV","pl");
legend->AddEntry(h2,"protons w/ kinE<20MeV","pl");
legend->SetLineWidth(0);
legend->Draw();
c47.SetGridx();
c47.SetGridy();
outFile->cd(); c47.Write();

//48th CANVAS - starting with the plots
TCanvas c48; c48.SetName("develop_pandoraNu_trajcluster_tail_to_tot");
q.setVersion("develop"); q.setDir("pandoraNu"); q.setHisto("tail_to_tot_2d_merged");
h1 = find_plot( saved_plot, q);
h1->SetTitle("Tail to total (develop); Length (mm); A. U.");
h1->SetName("pandora_merged");
h1->GetXaxis()->SetRangeUser(0,400);
h1->SetLineColor(kBlack);
h1->Draw();
q2.setVersion("develop"); q2.setDir("pandoraNu"); q2.setHisto("tail_to_tot_2d_not_merged");
h2 = find_plot( saved_plot, q2);
h2->SetName("pandora_notmerged");
h2->SetLineColor(kRed); h2->Draw("same");
q3.setVersion("develop"); q3.setDir("trajcluster"); q3.setHisto("tail_to_tot_2d_merged");
h3 = find_plot( saved_plot, q3);
h3->SetName("trajcluster_merged");
h3->SetLineColor(kBlue); h3->Draw("same");
qualifier q4("develop", "trajcluster", "tail_to_tot_2d_not_merged");
TH1D* h4 = find_plot( saved_plot, q4);
h4->SetName("trajcluster_notmerged");
h4->SetLineColor(kGreen); h4->Draw("same");
qualifier q5("mcc8", "pandoraNu", "tail_to_tot_2d_merged");
TH1D* h5 = find_plot( saved_plot, q5);
h5->SetName("mcc8_merged");
h5->SetLineColor(kViolet); h5->Draw("same");
qualifier q6("mcc8", "pandoraNu", "tail_to_tot_2d_not_merged");
TH1D* h6 = find_plot( saved_plot, q6);
h6->SetName("mcc8_notmerged");
h6->SetLineColor(kOrange); h6->Draw("same");
legend = new TLegend(0.4,0.1,0.8,0.3);
legend->AddEntry(h1,"merged, pandoraNu","pl");
legend->AddEntry(h2,"not merged, pandoraNu","pl");
legend->AddEntry(h3,"merged, trajcluster","pl");
legend->AddEntry(h4,"not merged, trajcluster","pl");
legend->AddEntry(h5,"merged, pandoraNu mcc8","pl");
legend->AddEntry(h6,"not merged, pandoraNu mcc8","pl");
legend->SetLineWidth(0);
legend->SetLineWidth(0);
legend->Draw();
c48.SetGridx();
c48.SetGridy();
outFile->cd(); c48.Write();



//save files
TIter next(outFile->GetListOfKeys());
TKey *key;
while ((key = (TKey*)next())) {
	  TClass *clsPtr = gROOT->GetClass(key->GetClassName());
	  TString name = key->GetClassName();
	  if (! name.Contains("TCanvas")) continue;
	  // do something with the key
	  TCanvas* save_canvas;
	  save_canvas = (TCanvas*) key->ReadObj();
	  save_canvas->SaveAs(("plots/"+std::string(key->GetName())+".pdf").c_str());
}















































































































outFile->Close();

}
