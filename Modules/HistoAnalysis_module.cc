////////////////////////////////////////////////////////////////////////
// Class:       HistoAnalysis
// Plugin Type: analyzer (art v2_05_01)
// File:        HistoAnalysis_module.cc
//
// Generated at Mon Aug 13 10:38:01 2018 by Simone Marcocci using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"


//ROOT
#include "TFile.h"
#include "TTree.h"

//local
#include "ubana/ProtonBenchmarker/Algos/HistoMaker.h"
#include "ubana/ProtonBenchmarker/Datatypes/StoredEvent.h"

namespace recohelper {
class HistoAnalysis;
}


class HistoAnalysis : public art::EDAnalyzer {
public:
  explicit HistoAnalysis(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  HistoAnalysis(HistoAnalysis const &) = delete;
  HistoAnalysis(HistoAnalysis &&) = delete;
  HistoAnalysis & operator = (HistoAnalysis const &) = delete;
  HistoAnalysis & operator = (HistoAnalysis &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:

  // Declare member data here.
  std::string fFileName;
  std::string fDirName;
  std::string fTreeName;

  art::ServiceHandle< art::TFileService > tfs;
  StoredEvent* event_store;
  reco_histo::HistoMaker* histo_maker;



};


HistoAnalysis::HistoAnalysis(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  fFileName = p.get<std::string> ("FileName");
  fDirName = p.get<std::string> ("DirName");
  fTreeName = p.get<std::string> ("TreeName");

}

void HistoAnalysis::analyze(art::Event const & e)
{
  // Implementation of required member function here.
  
	std::string fullname;
	cet::search_path sp("FW_SEARCH_PATH");
	sp.find_file(fFileName, fullname);
	
	if ( fullname.empty() )
		throw cet::exception("FileOpenError") << "Can't find " << fFileName;
	
	TFile* infile = new TFile(fullname.c_str(),"read");
	if (infile->IsZombie() )
		throw cet::exception("FileReadError") << "Can't open " << fullname;

	TTree* tree;
	tree = (TTree*) infile->Get( std::string( fDirName + "/" + fTreeName).c_str() );
	if (!tree)
		throw cet::exception("FileReadError") << "Can't open tree " << fDirName << "/" << fTreeName;

	event_store = new StoredEvent();
	tree->SetBranchAddress("stored", &event_store);

  	histo_maker = new reco_histo::HistoMaker();
  	histo_maker->Init( tfs );
  	art::TFileDirectory hits_dir = tfs->mkdir("hits_dir");
  	histo_maker->Init_Hit(  hits_dir );

	int n_events = 0;
	std::cout << "Processing " << tree->GetEntries() << " entries" << std:: endl;
	for (long i=0; i<tree->GetEntries(); i++) {
		if ( i!=0 && i%5000==0) std::cout << i << "/" << tree->GetEntries() << " entries completed. "<< std::endl;
		n_events++;
		tree->GetEntry(i);
				
   		histo_maker->Fill_Truth_Histos( event_store );
  		int muon_pos;
  		histo_maker->Fill_Analysis_Histos( event_store, muon_pos, true ); //plots tracking and vertexing information
  		histo_maker->Fill_Hit_Histos( event_store );

	}

  	histo_maker->ScalePlots( n_events );
  	histo_maker->FillCumulativeHistograms();


}

DEFINE_ART_MODULE(HistoAnalysis)
