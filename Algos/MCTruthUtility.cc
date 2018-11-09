#include "MCTruthUtility.h"

namespace mctruthutil {


void MCTruthUtility::FillTrackMCPVectors( std::unique_ptr<truth::IMCTruthMatching> &fMCTruthMatching, art::FindManyP<recob::Hit> &hits_from_tracks, std::vector< art::Ptr<recob::Track> > &tracklist, 
		art::Handle< std::vector<simb::MCParticle> > &mcpartHandle, std::vector< art::Ptr<simb::MCParticle> > &mcpar_v, std::vector< art::Ptr<recob::Track> > &track_v, std::vector< anab::BackTrackerMatchingData > &btdata_v) {
	
    double maxe = -1;
    double tote = 0;
    double totenergy = 0;
    anab::BackTrackerMatchingData btdata;

    size_t NTracks = tracklist.size();
    std::cout << "NTRACKS " << NTracks << std::endl;
      
      // Now to access MCTruth for each track... 
      for(size_t iTrk=0; iTrk < NTracks; ++iTrk) { 
	int TrueTrackID = 0;
	int TrackID     = 0;
	std::vector< art::Ptr<recob::Hit> > allHits = hits_from_tracks.at(iTrk);
	
	std::map<int,double> trkide;
	std::map<int,double> trkides;
	for(size_t h = 0; h < allHits.size(); ++h){
	  art::Ptr<recob::Hit> hit = allHits[h];
	  std::vector<sim::TrackIDE> TrackIDs = fMCTruthMatching->HitToTrackID(hit);
      	  std::vector<sim::TrackIDE> eveIDs = fMCTruthMatching->HitToEveID(hit);
	  
	  for(size_t e = 0; e < TrackIDs.size(); ++e){
  	    if ( e==0 ) 	trkide[TrackIDs[e].trackID] = TrackIDs[e].energy;
	    else		trkide[TrackIDs[e].trackID] += TrackIDs[e].energy;
	  }
	  for(size_t e = 0; e < eveIDs.size(); ++e){
	  if (e==0)	trkides[eveIDs[e].trackID] = eveIDs[e].energy;
	  else		trkides[eveIDs[e].trackID] += eveIDs[e].energy;
	  }

	TrackIDs.clear();
	eveIDs.clear();
	  
	}
	// Work out which IDE despoited the most charge in the hit if there was more than one.
	maxe = -1;
	tote = 0;
      for (std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
	tote += ii->second;
	if ((ii->second)>maxe){
	  maxe = ii->second;
	  TrackID = ii->first;
	}
      }
      
      btdata.cleanliness = maxe/tote;//purity
      
      const simb::MCParticle *tmpParticle = fMCTruthMatching->TrackIDToParticle(TrackID);
      if (!tmpParticle) continue; // Retain this check that the BackTracker can find the right particle
      // Now, loop through the MCParticle's myself to find the correct match
      int mcpart_i(-1);
      for (auto const particle : *mcpartHandle){
        mcpart_i++;
        if (TrackID == particle.TrackId()){
          break;
        }
      }
	
      const simb::MCParticle particle = mcpartHandle.product()->at(mcpart_i);
      TrueTrackID = particle.TrackId();
      
      auto diff = mcpart_i; // check I have a sensible value for this counter
      if (diff >= (int)mcpartHandle->size()){
        std::cout << "Error, the backtracker is doing weird things to your pointers!" << std::endl;
        throw std::exception();
      }
     
      totenergy=trkides[ TrueTrackID ];
      btdata.completeness = maxe/totenergy;

      art::Ptr<simb::MCParticle> mcpartPtr(mcpartHandle, mcpart_i);
	
      std::cout << "pushing back...." << std::endl;
      mcpar_v.push_back( mcpartPtr );
      track_v.push_back( tracklist[iTrk] );
      btdata_v.push_back( btdata );
      std::cout << "Size after pushing back " << mcpar_v.size() << std::endl;
    
      trkide.clear();
      trkides.clear();
    } // Loop over tracks   
}
  
void MCTruthUtility::FillShowerMCPVectors( std::unique_ptr<truth::IMCTruthMatching> &fMCTruthMatching, art::FindManyP<recob::Hit> &hits_from_showers, std::vector< art::Ptr<recob::Shower> > &showerlist, 
		art::Handle< std::vector<simb::MCParticle> > &mcpartHandle, std::vector< art::Ptr<simb::MCParticle> > &mcpar_v, std::vector< art::Ptr<recob::Shower> > &shower_v, std::vector< anab::BackTrackerMatchingData > &btdata_v) {
  
    double maxe = -1;
    double tote = 0;
    double totenergy = 0;
    anab::BackTrackerMatchingData btdata;
	
    // Now Loop over showers....
    size_t NShowers = showerlist.size();
    for (size_t Shower = 0; Shower < NShowers; ++Shower) {
      int ShowerID          = 0;
      std::vector< art::Ptr<recob::Hit> > allHits = hits_from_showers.at(Shower);
      
      std::map<int,double> showeride;
      std::map<int,double> showerides;
      for(size_t h = 0; h < allHits.size(); ++h){
	art::Ptr<recob::Hit> hit = allHits[h];
	std::vector<sim::IDE> ides;
	std::vector<sim::TrackIDE> TrackIDs = fMCTruthMatching->HitToTrackID(hit);
      	std::vector<sim::TrackIDE> eveIDs = fMCTruthMatching->HitToEveID(hit);
	
	for(size_t e = 0; e < TrackIDs.size(); ++e){
	if (e==0)	  showeride[TrackIDs[e].trackID] = TrackIDs[e].energy;
	else		  showeride[TrackIDs[e].trackID] += TrackIDs[e].energy;
	}
	for(size_t e = 0; e < eveIDs.size(); ++e){
	if (e==0)	showerides[eveIDs[e].trackID] = eveIDs[e].energy;
	else		showerides[eveIDs[e].trackID] += eveIDs[e].energy;
	}

	TrackIDs.clear();
	eveIDs.clear();
      }
      // Work out which IDE despoited the most charge in the hit if there was more than one.
      maxe = -1;
      tote = 0;
      for (std::map<int,double>::iterator ii = showeride.begin(); ii!=showeride.end(); ++ii){
	tote += ii->second;
	if ((ii->second)>maxe){
	  maxe = ii->second;
	  ShowerID = ii->first;
	}
      }
      btdata.cleanliness = maxe/tote;

      // Now have MCParticle trackID corresponding to shower, so get PdG code and T0 etc.
      const simb::MCParticle *tmpParticle = fMCTruthMatching->TrackIDToParticle(ShowerID);
      if (!tmpParticle) continue; // Retain this check that the BackTracker can find the right particle
      // Now, loop through the MCParticle's myself to find the correct match
      int mcpart_i(-1);
      for (auto const particle : *mcpartHandle){
        mcpart_i++;
        if (ShowerID == particle.TrackId()){
          break;
        }
      }
      const simb::MCParticle particle = mcpartHandle.product()->at(mcpart_i);
      ShowerID = particle.TrackId();
      
      auto diff = mcpart_i; // check I have a sensible value for this counter
      if (diff >= (int)mcpartHandle->size()){
        std::cout << "Error, the backtracker is doing weird things to your pointers!" << std::endl;
        throw std::exception();
      }
      
      totenergy= showerides[ ShowerID ];
      btdata.completeness = maxe/totenergy;

      art::Ptr<simb::MCParticle> mcpartPtr(mcpartHandle, mcpart_i);
	
      mcpar_v.push_back( mcpartPtr );
      shower_v.push_back( showerlist[Shower] );
      btdata_v.push_back( btdata );

      showeride.clear();
      showerides.clear();
    
    }// Loop over showers
}
  

void MCTruthUtility::FillPFParticleMCPVectors( std::unique_ptr<truth::IMCTruthMatching> &fMCTruthMatching, art::FindManyP<recob::Cluster> &clusters_from_pfps, art::FindManyP<recob::Hit> &hits_from_clusters,
		std::vector< art::Ptr<recob::PFParticle> > &pfplist, 
		art::Handle< std::vector<simb::MCParticle> > &mcpartHandle, std::vector< art::Ptr<simb::MCParticle> > &mcpar_v, std::vector< art::Ptr<recob::PFParticle> > &pfp_v, std::vector< anab::BackTrackerMatchingData > &btdata_v) {
  
    
    size_t NPfparticles = pfplist.size();
    
    // Now to access MCTruth for each pfparticle... 
    for(size_t iPfp=0; iPfp < NPfparticles; ++iPfp) { 
      int TrackID     = 0;
      anab::BackTrackerMatchingData btdata;
      
      std::vector< art::Ptr<recob::Hit> > allHits;
      //Get all hits through associated clusters
      std::vector< art::Ptr<recob::Cluster>> allClusters = clusters_from_pfps.at(iPfp);
      for (size_t iclu = 0; iclu<allClusters.size(); ++iclu){
        std::vector< art::Ptr<recob::Hit>> hits = hits_from_clusters.at(iclu);
        allHits.insert(allHits.end(), hits.begin(), hits.end());
      }

      std::map<int,double> trkide;
      std::map<int,double> trkides;
      for(size_t h = 0; h < allHits.size(); ++h){
	art::Ptr<recob::Hit> hit = allHits[h];
	std::vector<sim::TrackIDE> TrackIDs = fMCTruthMatching->HitToTrackID(hit);
      	std::vector<sim::TrackIDE> eveIDs = fMCTruthMatching->HitToEveID(hit);
	
	for(size_t e = 0; e < TrackIDs.size(); ++e){
	if (e==0)  trkide[TrackIDs[e].trackID] = TrackIDs[e].energy;
	else  trkide[TrackIDs[e].trackID] += TrackIDs[e].energy;
	}
	for(size_t e = 0; e < eveIDs.size(); ++e){
	if (e==0)	  trkides[eveIDs[e].trackID] = eveIDs[e].energy;
	else		  trkides[eveIDs[e].trackID] += eveIDs[e].energy;
	}

	TrackIDs.clear();
	eveIDs.clear();
      }
      // Work out which IDE despoited the most charge in the hit if there was more than one.
      double maxe = -1;
      double tote = 0;
      double totenergy = 0;
      for (std::map<int,double>::iterator ii = trkide.begin(); ii!=trkide.end(); ++ii){
	tote += ii->second;
	if ((ii->second)>maxe){
	  maxe = ii->second;
	  TrackID = ii->first;
	}
      }
      btdata.cleanliness = maxe/tote;
      
      const simb::MCParticle *tmpParticle = fMCTruthMatching->TrackIDToParticle(TrackID);
      if (!tmpParticle) continue; // Retain this check that the BackTracker can find the right particle
      // Now, loop through the MCParticle's myself to find the correct match
      int mcpart_i(-1);
      for (auto const particle : *mcpartHandle){
        mcpart_i++;
        if (TrackID == particle.TrackId()){
          break;
        }
      }
	
      const simb::MCParticle particle = mcpartHandle.product()->at(mcpart_i);
      int TrueTrackID = particle.TrackId();
      
      auto diff = mcpart_i; // check I have a sensible value for this counter
      if (diff >= (int)mcpartHandle->size()){
        std::cout << "Error, the backtracker is doing weird things to your pointers!" << std::endl;
        throw std::exception();
      }
     
      totenergy = trkides[ TrueTrackID ];
      btdata.completeness = maxe/totenergy;

      art::Ptr<simb::MCParticle> mcpartPtr(mcpartHandle, mcpart_i);
	
      mcpar_v.push_back( mcpartPtr );
      pfp_v.push_back( pfplist[iPfp] );
      btdata_v.push_back( btdata );
	
      trkide.clear();
      trkides.clear();

    } //loop over PFParticles

}

}//namespace
