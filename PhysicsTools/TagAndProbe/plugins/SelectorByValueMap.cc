#include "PhysicsTools/TagAndProbe/plugins/SelectorByValueMap.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

typedef SelectorByValueMap<pat::Electron, bool> PatElectronSelectorByValueMap;
DEFINE_FWK_MODULE(PatElectronSelectorByValueMap);

typedef SelectorByValueMap<pat::Photon, bool> PatPhotonSelectorByValueMap;
DEFINE_FWK_MODULE(PatPhotonSelectorByValueMap);
