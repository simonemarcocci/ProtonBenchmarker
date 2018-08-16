#include "canvas/Persistency/Common/PtrVector.h" 
#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Common/Assns.h"

#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"

#include "ubana/ProtonBenchmarker/Datatypes/StoredEvent.h"

#include "nusimdata/SimulationBase/MCParticle.h"



#include <vector>

template class std::vector<double>;
template class std::vector<float>;
template class std::vector<int>;
template class std::vector<bool>;
template class std::vector< std::vector<int>>;
template class std::vector< std::vector<float>>;
template class std::map<int, int>;

class StoredEvent;

