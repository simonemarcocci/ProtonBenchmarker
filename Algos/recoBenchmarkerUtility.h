#ifndef RECOBENCHMARKERUTILITY_H
#define RECOBENCHMARKERUTILITY_H

// larsoft includes
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "canvas/Persistency/Common/Ptr.h"

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

      std::vector<double> getMomentumVector(const art::Ptr< simb::MCParticle >& a);

      std::vector<double> getMomentumVector(const recob::Track& a);

      float getAngle(const std::vector<double> a, const std::vector<double> b, recoBenchmarkerUtility rbutil, std::string proj);
      
      std::vector<float> getUnitVector(std::vector<double> a); 

      float getDotProduct(std::vector<float> a, std::vector<float> b);

  };

}

#endif
