#ifndef RECOBENCHMARKERUTILITY_H
#define RECOBENCHMARKERUTILITY_H

// larsoft includes
#include "nusimdata/SimulationBase/MCParticle.h"
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

      float getAngle(const art::Ptr<simb::MCParticle>& a, const art::Ptr<simb::MCParticle>& b, recoBenchmarkerUtility rbutil, std::string proj);
      
      std::vector<float> getUnitVector(std::vector<double> a); 

      float getDotProduct(std::vector<float> a, std::vector<float> b);

  };

}

#endif
