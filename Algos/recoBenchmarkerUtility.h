#ifndef RECOBENCHMARKERUTILITY_H
#define RECOBENCHMARKERUTILITY_H

// larsoft includes
#include "nusimdata/SimulationBase/MCParticle.h"
#include "canvas/Persistency/Common/fwd.h"

// ROOT includes
#include "TLorentzVector.h"

// cc includes
#include <numeric>
#include <vector>
#include <functional>

namespace rbutil{

  class recoBenchmarkerUtility{

    public:

      float getAngle(const art::Ptr<simb::MCParticle>& a, const art::Ptr<simb::MCParticle>& b);

  };

}

#endif
