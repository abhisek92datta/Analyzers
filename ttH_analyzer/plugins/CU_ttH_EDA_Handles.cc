#ifndef CU_ttH_EDA_Handles_CC
#define CU_ttH_EDA_Handles_CC

#include "CU_ttH_EDA_Handles.h"

/// Set up handles with getByToken from edm::Event
void Set_up_handles(const Event &iEvent, const EventSetup &iSetup, edm_Handles &handle, edm_Tokens &token,
                    int isdata)
{
    iEvent.getByToken(token.triggerResults, handle.triggerResults);
    iEvent.getByToken(token.filterResults, handle.filterResults);

    iEvent.getByToken(token.vertices, handle.vertices);
    iEvent.getByToken(token.sec_vertices, handle.sec_vertices);
    iEvent.getByToken(token.PU_info, handle.PU_info);
    iEvent.getByToken(token.srcRho, handle.srcRho);

    iEvent.getByToken(token.electrons, handle.electrons);
    iEvent.getByToken(token.muons, handle.muons);
    iEvent.getByToken(token.jets, handle.jets);
    iEvent.getByToken(token.METs, handle.METs);	
    iEvent.getByToken(token.genjets, handle.genjets);

    iEvent.getByToken(token.PF_candidates, handle.PF_candidates);

    iEvent.getByToken(token.BS, handle.BS);

    iEvent.getByToken(token.mvaValuesMapToken_, handle.mvaValues);
    iEvent.getByToken(token.mvaCategoriesMapToken_, handle.mvaCategories);
    iEvent.getByToken(token.electrons_for_mva_token, handle.electrons_for_mva);
    iEvent.getByToken(token.muon_h_token, handle.muon_h);

    iEvent.getByToken(token.genTtbarIdToken_, handle.genTtbarId);
    iEvent.getByToken(token.puInfoToken, handle.PupInfo);
    
    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", handle.ttrkbuilder);
}

#endif
