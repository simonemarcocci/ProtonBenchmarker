#ifndef RECOBENCHMARKERUTILITY_H
#define RECOBENCHMARKERUTILITY_H

// larsoft includes
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/GeometryCore.h"
#include "lardata/Utilities/LArFFT.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// ROOT includes
#include "TLorentzVector.h"

// cc includes
#include <numeric>
#include <vector>
#include <functional>
#include <stdexcept>
#include <cmath>

namespace rbutil{

  class recoBenchmarkerUtility{

    public:

      bool isInTPC(const art::Ptr< simb::MCParticle >& );
      
      bool isInTPC(const simb::MCParticle & );
            
      std::vector<double> getMomentumVector(const art::Ptr< simb::MCParticle >& a);

      std::vector<double> getMomentumVector(const recob::Track& a);

      float getAngle(const std::vector<double> a, const std::vector<double> b, recoBenchmarkerUtility rbutil, std::string proj);
      
      std::vector<float> getUnitVector(std::vector<double> a); 

      float getDotProduct(std::vector<float> a, std::vector<float> b);

      std::vector<float> getHitXZPosition(const recob::Hit& thisHit, recoBenchmarkerUtility rbutil);

      float convertTicksToX(const recob::Hit& thisHit);

      bool isHitNearVertex(std::vector<float> v, std::vector<float> h);
  };

}

#endif
