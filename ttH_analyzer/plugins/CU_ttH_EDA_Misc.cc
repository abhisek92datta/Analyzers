#ifndef CU_ttH_EDA_Misc_cc
#define CU_ttH_EDA_Misc_cc

/// Includes
#include "CU_ttH_EDA.h"

/*
int CU_ttH_EDA::End_Run_hist_fill_triggers()
{
    TAxis *axis = h_hlt->GetXaxis();
    if (!axis)
        return 1;

    for (std::map<std::string, unsigned long>::const_iterator iter =
             n_trigger_fired.begin();
         iter != n_trigger_fired.end(); ++iter) {

        int bin_num = axis->FindBin((iter->first).c_str());
        h_hlt->Fill(bin_num - 1, (iter->second));
    }

    axis = h_flt->GetXaxis();
    if (!axis)
        return 1;

    for (std::map<std::string, unsigned long>::const_iterator iter =
             n_filter_fired.begin();
         iter != n_filter_fired.end(); ++iter) {

        int bin_num = axis->FindBin((iter->first).c_str());
        h_flt->Fill(bin_num - 1, (iter->second));
    }

    return 0;
}
*/

void CU_ttH_EDA::Update_common_vars(const edm::Event &iEvent,
                                    CU_ttH_EDA_event_vars &local)
{
    local.run_nr = iEvent.id().run();
    local.event_nr = iEvent.id().event();
    local.lumisection_nr = iEvent.id().luminosityBlock();
}

/*
int CU_ttH_EDA::Check_beam_spot(edm::Handle<reco::BeamSpot> BS)
{
    if (!BS.isValid())
        return 1;

    if (verbose_)
        printf("\t BeamSpot: x = %.2f,\t y = %.2f,\t z = %.2f \n", BS->x0(),
               BS->y0(), BS->z0());

    return 0;
}
*/

int CU_ttH_EDA::Check_triggers(edm::Handle<edm::TriggerResults> triggerResults,
                               CU_ttH_EDA_event_vars &local)
{
    if (!triggerResults.isValid()) {
        std::cerr << "Trigger results not valid for tag " << hltTag
                  << std::endl;
        return 1;
    }

    /// scan: trigger_on_HLT_<type>
    local.pass_single_e =
        Check_triggers_iterator(trigger_on_HLT_e, triggerResults);
    local.pass_single_mu =
        Check_triggers_iterator(trigger_on_HLT_mu, triggerResults);
    local.pass_double_e =
        Check_triggers_iterator(trigger_on_HLT_ee, triggerResults);
    local.pass_elemu =
        Check_triggers_iterator(trigger_on_HLT_emu, triggerResults);
    local.pass_double_mu =
        Check_triggers_iterator(trigger_on_HLT_mumu, triggerResults);
	/*
    if (trigger_stats) {
        for (std::vector<std::string>::const_iterator trigger =
                 trigger_names.begin();
             trigger != trigger_names.end(); ++trigger) {
            std::string pathName = *trigger;
            unsigned int hltIndex = hlt_config.triggerIndex(pathName);

            if (hltIndex >= triggerResults->size())
                continue;

            bool trigger_accept = triggerResults->accept(hltIndex);
            int prescale = -1; // hlt_config.prescaleValue(iEvent, iSetup,
                               // pathName);

            if (verbose_ && dumpHLT_)
                std::cout << " =====>  HLT: path name = " << pathName
                          << ",\t prescale = " << prescale
                          << ",\t pass = " << trigger_accept << std::endl;

            std::string pathNameNoVer = hlt_config.removeVersion(pathName);

            if (trigger_accept)
                ++n_trigger_fired[pathNameNoVer];
        }
    }
    */

    return 0;
}

inline bool CU_ttH_EDA::Check_triggers_iterator(
    const vector<string> &triggers,
    edm::Handle<edm::TriggerResults> triggerResults)
{
    for (std::vector<std::string>::const_iterator trigger = triggers.begin();
         trigger != triggers.end(); ++trigger) {

        std::string trigger_it;
        char s[100];
        unsigned int hltIndex;
        for (int i = 10; i >= 1; i--) {
            trigger_it.assign(*trigger);
            sprintf(s, "%d", i);
            trigger_it.append(s);
            hltIndex = hlt_config.triggerIndex(trigger_it);
            if (hltIndex >= triggerResults->size())
                continue;
            if (triggerResults->accept(hltIndex))
                return true;
        }
    }

    return false;
}


int CU_ttH_EDA::Check_filters(edm::Handle<edm::TriggerResults> filterResults, CU_ttH_EDA_event_vars &local)
{
    if (!filterResults.isValid()) {
        std::cerr << "Trigger results not valid for tag " << filterTag
                  << std::endl;
        return 1;
    }
     
    bool pass = 1;
    for (std::vector<std::string>::const_iterator filter = MET_filter_names.begin();
         filter != MET_filter_names.end(); ++filter) {
	    
        unsigned int filterIndex;  
	    std::string pathName = *filter;
            filterIndex = filter_config.triggerIndex(pathName);
            if (filterIndex >= filterResults->size()){
                pass = pass*0;
		break;
	    }
            if (filterResults->accept(filterIndex))
                pass=pass*1;
	    else {
	    	pass = pass*0;
		break;
	    }		   
    }
	 
    local.MET_filters = pass;
	
    /*	
    if (trigger_stats) {
        for (std::vector<std::string>::const_iterator trigger =
                 filter_names.begin();
             trigger != filter_names.end(); ++trigger) {
            std::string pathName = *trigger;
            unsigned int hltIndex = filter_config.triggerIndex(pathName);

            if (hltIndex >= filterResults->size())
                continue;

            bool filter_accept = filterResults->accept(hltIndex);
            int prescale = -1; // filter_config.prescaleValue(iEvent, iSetup,
                               // pathName);

            if (verbose_ && dumpHLT_)
                std::cout << " =====>  Filter: path name = " << pathName
                          << ",\t prescale = " << prescale
                          << ",\t pass = " << filter_accept << std::endl;

            std::string pathNameNoVer = filter_config.removeVersion(pathName);

            if (filter_accept)
                ++n_filter_fired[pathNameNoVer];
        }
    }
    */
    return 0;
}


int CU_ttH_EDA::Check_vertices_set_MAODhelper(
    edm::Handle<reco::VertexCollection> vertices)
{
    /// Primary vertex handling
    if (!vertices.isValid())
        return 1;

    reco::Vertex vertex;
    int n_PVs = 0;

    for (reco::VertexCollection::const_iterator vtx = vertices->begin();
         vtx != vertices->end(); ++vtx) {

        if (vtx->isFake() || vtx->ndof() < 4.0 || abs(vtx->z()) > 24.0 ||
            abs(vtx->position().Rho()) > 2.0)
            continue;

        if (n_PVs == 0)
            vertex = *vtx;

        ++n_PVs;
    }

    if (verbose_)
        printf("\t Event PV: x = %.3f,\t y = %.3f,\t z = %.3f \n", vertex.x(),
               vertex.y(), vertex.z());

    if (n_PVs > 0)
        miniAODhelper.SetVertex(
            vertex); // FIXME?: overload miniAODhelper::SetVertex(reco::Vertex&)
    else
        miniAODhelper.SetVertex(*(vertices->begin()));

    return 0;
}

int CU_ttH_EDA::Check_PV(const edm::Handle<reco::VertexCollection> &vertices)
{
    reco::VertexCollection::const_iterator vtx = vertices->begin();
    if (vtx->isFake() || vtx->ndof() < 4.0 || abs(vtx->z()) > 24.0 ||
        abs(vtx->position().Rho()) > 2.0)
        return 0;
    else
        return 1;
}

/// Other functions

void CU_ttH_EDA::Check_Fill_Print_single_lepton(
    const CU_ttH_EDA_event_vars &local)
{
    fprintf(events_combined, "%d,%d,%d,", local.run_nr, local.lumisection_nr,
            local.event_nr);
    fprintf(events_combined, "%d,%d,%d,%d,%d,", local.is_e, local.is_mu,
            local.is_ee, local.is_emu, local.is_mumu);
    fprintf(events_combined, "%d,%d,", local.n_sl_jets, local.n_sl_btags);
    if (local.n_electrons == 1) {
        fprintf(events_combined, "%.4f,%.4f,%d,", local.e_selected[0].pt(),
                miniAODhelper.GetElectronRelIso(local.e_selected[0],
                                                coneSize::R03, corrType::rhoEA,
                                                effAreaType::spring15),
                local.e_selected[0].pdgId());
    } else if (local.n_muons == 1) {
        fprintf(events_combined, "%.4f,%.4f,%d,", local.mu_selected[0].pt(),
                miniAODhelper.GetMuonRelIso(local.mu_selected[0], coneSize::R04,
                                            corrType::deltaBeta),
                local.mu_selected[0].pdgId());
    }
    fprintf(events_combined, "-1,-1,-1,");
    fprintf(
        events_combined, "%.4f,%.4f,%.4f,",
        local.jets_sl_selected_sorted[0].pt(),
        local.jets_sl_selected_sorted[0].eta(),
        local.jets_sl_selected_sorted[0].phi());
    fprintf(events_combined, "%.4f,%.4f,%.4f,%.4f,", local.jet1SF_sl,
            local.jet1SF_up_sl, local.jet1SF_down_sl,
            miniAODhelper.GetJetCSV(local.jets_sl_selected_sorted[0],
                                "pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    fprintf(
        events_combined, "%.4f,%.4f,%.4f,",
        local.jets_sl_selected_sorted[1].pt(),
        local.jets_sl_selected_sorted[1].eta(),
        local.jets_sl_selected_sorted[1].phi());
    fprintf(events_combined, "%.4f,%.4f,%.4f,%.4f,", local.jet2SF_sl,
            local.jet2SF_up_sl, local.jet2SF_down_sl,
            miniAODhelper.GetJetCSV(local.jets_sl_selected_sorted[1],
                                "pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    fprintf(events_combined, "%.4f,%.4f,-1,", local.met_pt, local.met_phi);
    fprintf(events_combined, "%d,", local.ttHf_cat);
    fprintf(events_combined, "%d,", local.truenpv);
    if (!isdata)
        fprintf(events_combined, "%.4f,", local.PU_weight);
    else
        fprintf(events_combined, "-1,");
    if (!isdata) {
        fprintf(events_combined, "%.4f,%.4f,%.4f,%.4f,", local.b_weight_sl,
                local.b_weight_sl_lfup, local.b_weight_sl_hfdown, 
                local.b_weight_sl_cErr1_down);
        //fprintf(events_combined, "%.4f,%.4f,%.4f,", local.lep_sf_trig_sl, local.lep_sf_id_sl,
        //        local.lep_sf_iso_sl);
    } else
        //fprintf(events_combined, "-1,-1,-1,-1,-1,-1,-1,");
        fprintf(events_combined, "-1,-1,-1,-1,");
    if (!isdata)
        fprintf(events_combined, "%.4f,%.4f,", local.pdf_weight_up,
                local.pdf_weight_down);
    else
        fprintf(events_combined, "-1,-1,");
    if (!isdata)
        fprintf(events_combined, "%.4f,%.4f,", local.q2_weight_up,
                local.q2_weight_down);
    else
        fprintf(events_combined, "-1,-1,");      
    //fprintf(events_combined, "-1,-1,-1,-1\n");
    //fprintf(events_combined, "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",local.bjetnessFV_num_leps, local.bjetnessFV_npvTrkOVcollTrk, local.bjetnessFV_avip3d_val, local.bjetnessFV_avip3d_sig, local.bjetnessFV_avsip3d_sig, local.bjetnessFV_avip1d_sig);
}

void CU_ttH_EDA::Check_Fill_Print_di_lepton(const CU_ttH_EDA_event_vars &local)
{
    fprintf(events_combined, "%d,%d,%d,", local.run_nr, local.lumisection_nr,
            local.event_nr);
    fprintf(events_combined, "%d,%d,%d,%d,%d,", local.is_e, local.is_mu,
            local.is_ee, local.is_emu, local.is_mumu);
    fprintf(events_combined, "%d,%d,", local.n_di_jets, local.n_di_btags);
    if (local.n_di_electrons == 2) {
        fprintf(events_combined, "%.4f,%.4f,%d,",
                local.e_di_selected_sorted[0].pt(),
                miniAODhelper.GetElectronRelIso(local.e_di_selected_sorted[0],
                                                coneSize::R03, corrType::rhoEA,
                                                effAreaType::spring15),
                local.e_di_selected_sorted[0].pdgId());
        fprintf(events_combined, "%.4f,%.4f,%d,",
                local.e_di_selected_sorted[1].pt(),
                miniAODhelper.GetElectronRelIso(local.e_di_selected_sorted[1],
                                                coneSize::R03, corrType::rhoEA,
                                                effAreaType::spring15),
                local.e_di_selected_sorted[1].pdgId());
    } else if (local.n_di_muons == 2) {
        fprintf(events_combined, "%.4f,%.4f,%d,",
                local.mu_di_selected_sorted[0].pt(),
                miniAODhelper.GetMuonRelIso(local.mu_di_selected_sorted[0],
                                            coneSize::R04, corrType::deltaBeta),
                local.mu_di_selected_sorted[0].pdgId());
        fprintf(events_combined, "%.4f,%.4f,%d,",
                local.mu_di_selected_sorted[1].pt(),
                miniAODhelper.GetMuonRelIso(local.mu_di_selected_sorted[1],
                                            coneSize::R04, corrType::deltaBeta),
                local.mu_di_selected_sorted[1].pdgId());
    } else {
        if (local.mu_di_selected_sorted[0].pt() >
            local.e_di_selected_sorted[0].pt()) {
            fprintf(events_combined, "%.4f,%.4f,%d,",
                    local.mu_di_selected_sorted[0].pt(),
                    miniAODhelper.GetMuonRelIso(local.mu_di_selected_sorted[0],
                                                coneSize::R04,
                                                corrType::deltaBeta),
                    local.mu_di_selected_sorted[0].pdgId());
            fprintf(events_combined, "%.4f,%.4f,%d,",
                    local.e_di_selected_sorted[0].pt(),
                    miniAODhelper.GetElectronRelIso(
                        local.e_di_selected_sorted[0], coneSize::R03,
                        corrType::rhoEA, effAreaType::spring15),
                    local.e_di_selected_sorted[0].pdgId());
        } else {
            fprintf(events_combined, "%.4f,%.4f,%d,",
                    local.e_di_selected_sorted[0].pt(),
                    miniAODhelper.GetElectronRelIso(
                        local.e_di_selected_sorted[0], coneSize::R03,
                        corrType::rhoEA, effAreaType::spring15),
                    local.e_di_selected_sorted[0].pdgId());
            fprintf(events_combined, "%.4f,%.4f,%d,",
                    local.mu_di_selected_sorted[0].pt(),
                    miniAODhelper.GetMuonRelIso(local.mu_di_selected_sorted[0],
                                                coneSize::R04,
                                                corrType::deltaBeta),
                    local.mu_di_selected_sorted[0].pdgId());
        }
    }
    fprintf(
        events_combined, "%.4f,%.4f,%.4f,",
        local.jets_di_selected_sorted[0].pt(),
        local.jets_di_selected_sorted[0].eta(),
        local.jets_di_selected_sorted[0].phi());
    fprintf(events_combined, "%.4f,%.4f,%.4f,%.4f,", local.jet1SF_di,
            local.jet1SF_up_di, local.jet1SF_down_di,
            miniAODhelper.GetJetCSV(local.jets_di_selected_sorted[0],
                                "pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    fprintf(
        events_combined, "%.4f,%.4f,%.4f,",
        local.jets_di_selected_sorted[1].pt(),
        local.jets_di_selected_sorted[1].eta(),
        local.jets_di_selected_sorted[1].phi());
    fprintf(events_combined, "%.4f,%.4f,%.4f,%.4f,", local.jet2SF_di,
            local.jet2SF_up_di, local.jet2SF_down_di,
            miniAODhelper.GetJetCSV(local.jets_di_selected_sorted[1],
                                "pfCombinedInclusiveSecondaryVertexV2BJetTags"));   
    fprintf(events_combined, "%.4f,%.4f,%.4f,", local.met_pt, local.met_phi,
            local.mll);  
    fprintf(events_combined, "%d,", local.ttHf_cat);
    fprintf(events_combined, "%d,", local.truenpv);
    if (!isdata)
        fprintf(events_combined, "%.4f,", local.PU_weight);
    else
        fprintf(events_combined, "-1,");
    if (!isdata) {
        fprintf(events_combined, "%.4f,%.4f,%.4f,%.4f,", local.b_weight_di,
                local.b_weight_di_lfup, local.b_weight_di_hfdown, 
                local.b_weight_di_cErr1_down);
        fprintf(events_combined, "%.4f,%.4f,%.4f,", local.lep_sf_trig_di, local.lep_sf_id_di,
                local.lep_sf_iso_di);
    } else
        //fprintf(events_combined, "-1,-1,-1,-1,-1,-1,-1,");
        fprintf(events_combined, "-1,-1,-1,-1,");
    if (!isdata)
        fprintf(events_combined, "%.4f,%.4f,", local.pdf_weight_up,
                local.pdf_weight_down);
    else
        fprintf(events_combined, "-1,-1,");
    if (!isdata)
        fprintf(events_combined, "%.4f,%.4f,", local.q2_weight_up,
                local.q2_weight_down);
    else
        fprintf(events_combined, "-1,-1,");      
    //fprintf(events_combined, "-1,-1,-1,-1\n");
    
    //fprintf(events_combined, "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",local.bjetnessFV_num_leps, local.bjetnessFV_npvTrkOVcollTrk, local.bjetnessFV_avip3d_val, local.bjetnessFV_avip3d_sig, local.bjetnessFV_avsip3d_sig, local.bjetnessFV_avip1d_sig);
}

void CU_ttH_EDA::Select_Leptons(CU_ttH_EDA_event_vars &local,
                                const edm_Handles &handle)
{

    local.e_with_id = miniAODhelper.GetElectronsWithMVAid(
        handle.electrons_for_mva, handle.mvaValues, handle.mvaCategories);
    // Single Lepton
    local.mu_selected = miniAODhelper.GetSelectedMuons(
        *(handle.muons), min_mu_pT, muonID::muonTight, coneSize::R04,
        corrType::deltaBeta, max_mu_eta);
    local.mu_veto_selected = miniAODhelper.GetSelectedMuons(
        *(handle.muons), min_veto_mu_pT, muonID::muonTightDL, coneSize::R04,
        corrType::deltaBeta, max_veto_mu_eta);
    local.e_selected = miniAODhelper.GetSelectedElectrons(
        local.e_with_id, min_ele_pT, electronID::electron80XCutBasedM,
        max_ele_eta);
    local.e_veto_selected = miniAODhelper.GetSelectedElectrons(
        local.e_with_id, min_veto_ele_pT,
        electronID::electron80XCutBasedM,
	max_veto_ele_eta);
    //local.e_selected = GetSelectedElectrons(
    //    *(handle.electrons_for_mva), min_ele_pT, electronID::electronGenPurposeMVAid80,
    //    handle.medium_id_decisions, max_ele_eta);
    //local.e_veto_selected = GetSelectedElectrons(
    //    *(handle.electrons_for_mva), min_veto_ele_pT,
    //    electronID::electronGenPurposeMVAid80, handle.medium_id_decisions,
	//max_veto_ele_eta);
    local.n_electrons = static_cast<int>(local.e_selected.size());
    local.n_veto_electrons = static_cast<int>(local.e_veto_selected.size());
    local.n_muons = static_cast<int>(local.mu_selected.size());
    local.n_veto_muons = static_cast<int>(local.mu_veto_selected.size());
    local.n_leptons = local.n_electrons + local.n_muons;

    // Dilepton
    local.mu_di_selected = miniAODhelper.GetSelectedMuons(
        *(handle.muons), min_di_mu2_pT, muonID::muonTightDL, coneSize::R04,
        corrType::deltaBeta, max_di_mu2_eta);
    local.e_di_selected = miniAODhelper.GetSelectedElectrons(
        local.e_with_id, min_di_ele2_pT,
        electronID::electron80XCutBasedM,
	max_di_ele2_eta);
    //local.e_di_selected = GetSelectedElectrons(
    //    *(handle.electrons_for_mva), min_di_ele2_pT,
    //    electronID::electronGenPurposeMVAid80, handle.medium_id_decisions,
	//max_di_ele2_eta);
    local.n_di_electrons = static_cast<int>(local.e_di_selected.size());
    local.n_di_muons = static_cast<int>(local.mu_di_selected.size());
    /// Sort leptons by pT
    local.mu_di_selected_sorted =
        miniAODhelper.GetSortedByPt(local.mu_di_selected);
    local.e_di_selected_sorted =
        miniAODhelper.GetSortedByPt(local.e_di_selected);
    local.n_di_leptons = local.n_di_electrons + local.n_di_muons;
}

inline std::vector<pat::Electron>
CU_ttH_EDA::GetSelectedElectrons(const edm::View<pat::Electron>& inputElectrons, const float iMinPt, const electronID::electronID iElectronID, const edm::Handle<edm::ValueMap<bool>>& medium_id_decisions, const float iMaxEta){
  std::vector<pat::Electron> selectedElectrons;
  bool passesID;
  
  for (size_t i = 0; i < inputElectrons.size(); ++i){
	const auto el = inputElectrons.ptrAt(i);
    passesID = false;    
    passesID = (*medium_id_decisions)[el];
   
    if(miniAODhelper.isGoodElectron(*el,iMinPt,iMaxEta,iElectronID) ) {
  		if( passesID == true )
  			selectedElectrons.push_back(*el); 
  	}
  }
  return selectedElectrons;
}

void CU_ttH_EDA::Select_Jets(CU_ttH_EDA_event_vars &local,
                             const edm::Event &iEvent,
                             const edm::EventSetup &iSetup,
                             const edm_Handles &handle, const double &rho,
                             const JME::JetResolution &resolution)
{

    // ID Check
    local.jets_raw = miniAODhelper.GetSelectedJets(*(handle.jets), 0., 999,
                                                   jetID::jetTight, '-');

    // Overlap removal
    local.jets_sl_raw = miniAODhelper.GetDeltaRCleanedJets(
        local.jets_raw, local.mu_veto_selected, local.e_veto_selected, 0.4);
    local.jets_di_raw = miniAODhelper.GetDeltaRCleanedJets(
        local.jets_raw, local.mu_di_selected, local.e_di_selected, 0.4);

    // Uncorrected jets
    local.jets_sl_raw = miniAODhelper.GetUncorrectedJets(local.jets_sl_raw);
    local.jets_di_raw = miniAODhelper.GetUncorrectedJets(local.jets_di_raw);

    // Jet Energy Correction

    bool doJER;
    if (isdata)
        doJER = 0;
    else
        doJER = 1;

    // using my jet correction function
    local.jets_sl_corrected =
        GetCorrectedJets(local.jets_sl_raw, handle.genjets, rho, resolution,
                         sysType::NA, 1, doJER);
    local.jets_di_corrected =
        GetCorrectedJets(local.jets_di_raw, handle.genjets, rho, resolution,
                         sysType::NA, 1, doJER);

    // using MiniAODHelper's jet correction function
    // local.jets_sl_corrected =
    // miniAODhelper.GetCorrectedJets(local.jets_sl_raw, iEvent, iSetup, handle.genjets,r,
    // sysType::NA, 1, doJER);
    // local.jets_di_corrected =
    // miniAODhelper.GetCorrectedJets(local.jets_di_raw, iEvent, iSetup, handle.genjets,r,
    // sysType::NA, 1, doJER);

    // for b-weight
    local.iSys.push_back(0); // none 
    local.iSys.push_back(9); // LF up
    local.iSys.push_back(12);   // HF down
    local.iSys.push_back(22);   // CErr1 down
	
    // Jet Selection
    local.jets_sl_selected = miniAODhelper.GetSelectedJets(
        local.jets_sl_corrected, min_jet_pT, max_jet_eta, jetID::none, '-');
    local.jets_di_selected = miniAODhelper.GetSelectedJets(
        local.jets_di_corrected, min_jet2_pT, max_jet_eta, jetID::none, '-');

    // Storing selected raw jets for JEC SF calculation

    for (unsigned int j = 0; j < local.jets_sl_corrected.size(); ++j) {
        if ((local.jets_sl_corrected[j].pt() > min_jet_pT) &&
            (fabs(local.jets_sl_corrected[j].eta()) < max_jet_eta)) {
            local.jets_sl_selected_raw.push_back(local.jets_sl_raw[j]);
        }
    }

    for (unsigned int j = 0; j < local.jets_di_corrected.size(); ++j) {
        if ((local.jets_di_corrected[j].pt() > min_jet2_pT) &&
            (fabs(local.jets_di_corrected[j].eta()) < max_jet_eta)) {
            local.jets_di_selected_raw.push_back(local.jets_di_raw[j]);
        }
    }

    // b-tagged jet selection
    local.jets_sl_selected_tag = miniAODhelper.GetSelectedJets(
        local.jets_sl_selected, min_bjet_pT, max_bjet_eta, jetID::none,
        MAODHelper_b_tag_strength);
    local.jets_di_selected_tag = miniAODhelper.GetSelectedJets(
        local.jets_di_selected, min_bjet_pT, max_bjet_eta, jetID::none,
        MAODHelper_b_tag_strength);

    local.n_sl_jets = static_cast<int>(local.jets_sl_selected.size());
    local.n_sl_btags = static_cast<int>(local.jets_sl_selected_tag.size());
    local.n_di_jets = static_cast<int>(local.jets_di_selected.size());
    local.n_di_btags = static_cast<int>(local.jets_di_selected_tag.size());

    /// Sort jets by pT
    local.jets_sl_selected_sorted =
        miniAODhelper.GetSortedByPt(local.jets_sl_selected);
    local.jets_sl_selected_tag_sorted =
        miniAODhelper.GetSortedByPt(local.jets_sl_selected_tag);

    local.jets_di_selected_sorted =
        miniAODhelper.GetSortedByPt(local.jets_di_selected);
    local.jets_di_selected_tag_sorted =
        miniAODhelper.GetSortedByPt(local.jets_di_selected_tag);
}

void CU_ttH_EDA::Init_Mets(CU_ttH_EDA_event_vars &local,
                           const edm_Handles &handle)
{

    local.pfMET = handle.METs->front();
    local.met_pt = local.pfMET.pt();
    local.met_phi = local.pfMET.phi();
    local.met_passed = 0;
    local.mll_passed = 0;
}

inline std::vector<pat::Jet>
CU_ttH_EDA::CheckJetID(const std::vector<pat::Jet> &inputJets,
                       const std::vector<pat::Jet> &inputJets_old)
{
    std::vector<pat::Jet> outputJets;
    bool loose = false;
    int N = static_cast<int>(inputJets.size());
    double scale;
    for (int i = 0; i < N; ++i) {
        scale = inputJets[i].pt() / inputJets_old[i].pt();
        loose = (inputJets[i].neutralHadronEnergyFraction() * scale < 0.99 &&
                 inputJets[i].chargedEmEnergyFraction() * scale < 0.99 &&
                 inputJets[i].neutralEmEnergyFraction() * scale < 0.99 &&
                 (inputJets[i].neutralMultiplicity() +
                  inputJets[i].chargedMultiplicity()) > 1);

        if (fabs(inputJets[i].eta()) < 2.4) {
            loose = (loose &&
                     inputJets[i].chargedHadronEnergyFraction() * scale > 0.0 &&
                     inputJets[i].chargedMultiplicity() > 0);
        }
        if (loose == true)
            outputJets.push_back(inputJets[i]);
    }
    return outputJets;
}

void CU_ttH_EDA::SetFactorizedJetCorrector(const sysType::sysType iSysType, CU_ttH_EDA_event_vars &local)
{

    std::vector<JetCorrectorParameters> corrParams;
    if (isdata) {

        if( local.run_nr>=272007 && local.run_nr<=276811 ) {
            JetCorrectorParameters *L3JetPar = new JetCorrectorParameters(
                "data/JEC/Summer16_23Sep2016BCDV3_DATA_L3Absolute_AK4PFchs.txt");
            JetCorrectorParameters *L2JetPar = new JetCorrectorParameters(
                "data/JEC/Summer16_23Sep2016BCDV3_DATA_L2Relative_AK4PFchs.txt");
            JetCorrectorParameters *L1JetPar = new JetCorrectorParameters(
                "data/JEC/Summer16_23Sep2016BCDV3_DATA_L1FastJet_AK4PFchs.txt");
            JetCorrectorParameters *L2L3JetPar = new JetCorrectorParameters(
                "data/JEC/Summer16_23Sep2016BCDV3_DATA_L2L3Residual_AK4PFchs.txt");
        }
        else if( local.run_nr>=276831 && local.run_nr<=278801 ) {
            JetCorrectorParameters *L3JetPar = new JetCorrectorParameters(
               "data/JEC/Summer16_23Sep2016EFV3_DATA_L3Absolute_AK4PFchs.txt");
            JetCorrectorParameters *L2JetPar = new JetCorrectorParameters(
               "data/JEC/Summer16_23Sep2016EFV3_DATA_L2Relative_AK4PFchs.txt");
            JetCorrectorParameters *L1JetPar = new JetCorrectorParameters(
               "data/JEC/Summer16_23Sep2016EFV3_DATA_L1FastJet_AK4PFchs.txt");
            JetCorrectorParameters *L2L3JetPar = new JetCorrectorParameters(
               "data/JEC/Summer16_23Sep2016EFV3_DATA_L2L3Residual_AK4PFchs.txt");
        }
        else if( local.run_nr>=278802 && local.run_nr<=280385 ) {
            JetCorrectorParameters *L3JetPar = new JetCorrectorParameters(
                "data/JEC/Summer16_23Sep2016GV3_DATA_L3Absolute_AK4PFchs.txt");
            JetCorrectorParameters *L2JetPar = new JetCorrectorParameters(
                "data/JEC/Summer16_23Sep2016GV3_DATA_L2Relative_AK4PFchs.txt");
            JetCorrectorParameters *L1JetPar = new JetCorrectorParameters(
                "data/JEC/Summer16_23Sep2016GV3_DATA_L1FastJet_AK4PFchs.txt");
            JetCorrectorParameters *L2L3JetPar = new JetCorrectorParameters(
                 "data/JEC/Summer16_23Sep2016GV3_DATA_L2L3Residual_AK4PFchs.txt");
        }
        else if( local.run_nr>=280919 && local.run_nr<=284044 ) {
            JetCorrectorParameters *L3JetPar = new JetCorrectorParameters(
                "data/JEC/Summer16_23Sep2016HV3_DATA_L3Absolute_AK4PFchs.txt");
            JetCorrectorParameters *L2JetPar = new JetCorrectorParameters(
                "data/JEC/Summer16_23Sep2016HV3_DATA_L2Relative_AK4PFchs.txt");
            JetCorrectorParameters *L1JetPar = new JetCorrectorParameters(
                "data/JEC/Summer16_23Sep2016HV3_DATA_L1FastJet_AK4PFchs.txt");
            JetCorrectorParameters *L2L3JetPar = new JetCorrectorParameters(
                "data/JEC/Summer16_23Sep2016HV3_DATA_L2L3Residual_AK4PFchs.txt");
        }

        corrParams.push_back(*L1JetPar);
        corrParams.push_back(*L2JetPar);
        corrParams.push_back(*L3JetPar);
        corrParams.push_back(*L2L3JetPar);
        _jetCorrector = new FactorizedJetCorrector(corrParams);

        std::string _JESUncFile =
            "data/JEC/Spring16_25nsV6_DATA_Uncertainty_AK4PFchs.txt";
        _jetCorrectorUnc = new JetCorrectionUncertainty(_JESUncFile);

        delete L3JetPar;
        delete L2JetPar;
        delete L1JetPar;
        delete L2L3JetPar;
    }

    else {
        JetCorrectorParameters *L3JetPar = new JetCorrectorParameters(
            "data/JEC/Summer16_23Sep2016V3_MC_L3Absolute_AK4PFchs.txt");
        JetCorrectorParameters *L2JetPar = new JetCorrectorParameters(
            "data/JEC/Summer16_23Sep2016V3_MC_L2Relative_AK4PFchs.txt");
        JetCorrectorParameters *L1JetPar = new JetCorrectorParameters(
            "data/JEC/Summer16_23Sep2016V3_MC_L1FastJet_AK4PFchs.txt");

        corrParams.push_back(*L1JetPar);
        corrParams.push_back(*L2JetPar);
        corrParams.push_back(*L3JetPar);
        _jetCorrector = new FactorizedJetCorrector(corrParams);

        std::string _JESUncFile =
            "data/JEC/Summer16_23Sep2016V3_MC_Uncertainty_AK4PFchs.txt";
        _jetCorrectorUnc = new JetCorrectionUncertainty(_JESUncFile);

        delete L3JetPar;
        delete L2JetPar;
        delete L1JetPar;
    }
}

inline double CU_ttH_EDA::GetJetSF(pat::Jet jet,
                                   const sysType::sysType iSysType,
                                   const double &rho)
{
    double scale = 1;
    _jetCorrector->setJetPt(jet.pt());
    _jetCorrector->setJetEta(jet.eta());
    _jetCorrector->setJetA(jet.jetArea());
    _jetCorrector->setRho(rho); //=fixedGridRhoFastjetAll

    scale = _jetCorrector->getCorrection();
    jet.scaleEnergy(scale);

    if (iSysType == sysType::JESup || iSysType == sysType::JESdown) {
        _jetCorrectorUnc->setJetPt(jet.pt());
        _jetCorrectorUnc->setJetEta(jet.eta());
        double unc = 1;
        double jes = 1;
        if (iSysType == sysType::JESup) {
            unc = _jetCorrectorUnc->getUncertainty(true);
            jes = 1 + unc;
        } else if (iSysType == sysType::JESdown) {
            unc = _jetCorrectorUnc->getUncertainty(false);
            jes = 1 - unc;
        }
        jet.scaleEnergy(jes);
        return jes;
    } else
        return scale;
}

inline std::vector<pat::Jet> CU_ttH_EDA::GetCorrectedJets(
    const std::vector<pat::Jet> &inputJets,
    const edm::Handle<reco::GenJetCollection> &genjets, const double &rho,
    const JME::JetResolution &resolution, const sysType::sysType iSysType,
    const bool &doJES, const bool &doJER, const float &corrFactor,
    const float &uncFactor)
{

    std::vector<pat::Jet> outputJets;

    for (std::vector<pat::Jet>::const_iterator it = inputJets.begin(),
                                               ed = inputJets.end();
         it != ed; ++it) {

        pat::Jet jet = (*it);
        double scale = 1.;

        // JEC
        if (doJES == 1) {

            _jetCorrector->setJetPt(jet.pt());
            _jetCorrector->setJetEta(jet.eta());
            _jetCorrector->setJetA(jet.jetArea());
            _jetCorrector->setRho(rho); //=fixedGridRhoFastjetAll

            scale = _jetCorrector->getCorrection();
            jet.scaleEnergy(scale);

            if (iSysType == sysType::JESup || iSysType == sysType::JESdown) {
                _jetCorrectorUnc->setJetPt(jet.pt());
                _jetCorrectorUnc->setJetEta(
                    jet.eta()); // here you must use the CORRECTED jet pt
                double unc = 1;
                double jes = 1;
                if (iSysType == sysType::JESup) {
                    unc = _jetCorrectorUnc->getUncertainty(true);
                    jes = 1 + unc;
                } else if (iSysType == sysType::JESdown) {
                    unc = _jetCorrectorUnc->getUncertainty(false);
                    jes = 1 - unc;
                }

                jet.scaleEnergy(jes);
            }
        }
        // JER
        if (doJER == 1) {

            double jerSF = 1.;
            bool genjet_match = 0;
            double dpt_min = 99999;
            double dpt;
            double dR;

            JME::JetParameters parameters_1;
            parameters_1.setJetPt(jet.pt());
            parameters_1.setJetEta(jet.eta());
            parameters_1.setRho(rho);
            float res = resolution.getResolution(parameters_1) * jet.pt();
            reco::GenJet matched_genjet;

            for (reco::GenJetCollection::const_iterator iter = genjets->begin();
                 iter != genjets->end(); ++iter) {
                dpt = fabs(jet.pt() - iter->pt());
                dR = miniAODhelper.DeltaR(&jet, iter);
                if (dR < (0.4 / 2)) {
                    if (dpt < (3 * fabs(res))) {
                        genjet_match = 1;
                        if (dpt <= dpt_min) {
                            matched_genjet = *(iter);
                            dpt_min = dpt;
                        }
                    }
                }
            }

            if (genjet_match == 1) {
                if (iSysType == sysType::JERup) {
                    jerSF = getJERfactor(uncFactor, fabs(jet.eta()),
                                         matched_genjet.pt(), jet.pt());
                } else if (iSysType == sysType::JERdown) {
                    jerSF = getJERfactor(-uncFactor, fabs(jet.eta()),
                                         matched_genjet.pt(), jet.pt());
                } else {
                    jerSF = getJERfactor(0, fabs(jet.eta()),
                                         matched_genjet.pt(), jet.pt());
                }
            } else if (genjet_match == 0) {
                 jerSF = r->Gaus();
                 r->SetSeed(1);
                 //jerSF = 1;
            }
            jet.scaleEnergy(jerSF * corrFactor);
        }
        outputJets.push_back(jet);
    }

    return outputJets;
}

// JER factors for 80X
inline double CU_ttH_EDA::getJERfactor(const int returnType,
                                       const double jetAbsETA,
                                       const double genjetPT,
                                       const double recojetPT)
{

    double factor = 1.;

    double scale_JER = 1., scale_JERup = 1., scale_JERdown = 1.;
    double extrauncertainty = 1.5;

    if (jetAbsETA < 0.5) {
        scale_JER = 1.109;
        scale_JERup = 1.109 + 0.008 * extrauncertainty;
        scale_JERdown = 1.109 - 0.008 * extrauncertainty;
    } else if (jetAbsETA < 0.8) {
        scale_JER = 1.138;
        scale_JERup = 1.138 + 0.013 * extrauncertainty;
        scale_JERdown = 1.138 - 0.013 * extrauncertainty;
    } else if (jetAbsETA < 1.1) {
        scale_JER = 1.114;
        scale_JERup = 1.114 + 0.013 * extrauncertainty;
        scale_JERdown = 1.114 - 0.013 * extrauncertainty;
    } else if (jetAbsETA < 1.3) {
        scale_JER = 1.123;
        scale_JERup = 1.123 + 0.024 * extrauncertainty;
        scale_JERdown = 1.123 - 0.024 * extrauncertainty;
    } else if (jetAbsETA < 1.7) {
        scale_JER = 1.084;
        scale_JERup = 1.084 + 0.011 * extrauncertainty;
        scale_JERdown = 1.084 - 0.011 * extrauncertainty;
    } else if (jetAbsETA < 1.9) {
        scale_JER = 1.082;
        scale_JERup = 1.082 + 0.035 * extrauncertainty;
        scale_JERdown = 1.082 - 0.035 * extrauncertainty;
    } else if (jetAbsETA < 2.1) {
        scale_JER = 1.140;
        scale_JERup = 1.140 + 0.047 * extrauncertainty;
        scale_JERdown = 1.140 - 0.047 * extrauncertainty;
    } else if (jetAbsETA < 2.3) {
        scale_JER = 1.067;
        scale_JERup = 1.067 + 0.053 * extrauncertainty;
        scale_JERdown = 1.067 - 0.053 * extrauncertainty;
    } else if (jetAbsETA < 2.5) {
        scale_JER = 1.177;
        scale_JERup = 1.177 + 0.041 * extrauncertainty;
        scale_JERdown = 1.177 - 0.041 * extrauncertainty;
    } else if (jetAbsETA < 2.8) {
        scale_JER = 1.364;
        scale_JERup = 1.364 + 0.039 * extrauncertainty;
        scale_JERdown = 1.364 - 0.039 * extrauncertainty;
    } else if (jetAbsETA < 3.0) {
        scale_JER = 1.857;
        scale_JERup = 1.857 + 0.071 * extrauncertainty;
        scale_JERdown = 1.857 - 0.071 * extrauncertainty;
    } else if (jetAbsETA < 3.2) {
        scale_JER = 1.328;
        scale_JERup = 1.328 + 0.022 * extrauncertainty;
        scale_JERdown = 1.328 - 0.022 * extrauncertainty;
    } else if (jetAbsETA < 4.7) {
        scale_JER = 1.160;
        scale_JERup = 1.160 + 0.029 * extrauncertainty;
        scale_JERdown = 1.160 - 0.029 * extrauncertainty;
    }

    double jetPt_JER = recojetPT;
    double jetPt_JERup = recojetPT;
    double jetPt_JERdown = recojetPT;

    double diff_recojet_genjet = recojetPT - genjetPT;

    jetPt_JER = std::max(0., genjetPT + scale_JER * (diff_recojet_genjet));
    jetPt_JERup = std::max(0., genjetPT + scale_JERup * (diff_recojet_genjet));
    jetPt_JERdown =
        std::max(0., genjetPT + scale_JERdown * (diff_recojet_genjet));

    if (returnType == 1)
        factor = jetPt_JERup / recojetPT;
    else if (returnType == -1)
        factor = jetPt_JERdown / recojetPT;
    else
        factor = jetPt_JER / recojetPT;

    return factor;
}

void CU_ttH_EDA::Check_SL_Event_Selection(CU_ttH_EDA_event_vars &local)
{
   
    if (local.pass_single_e == 1 || local.pass_single_mu == 1) {
        
        if (local.MET_filters ==0 || local.filterbadChCandidate == 0 || local.filterbadPFMuon == 0)
        	return;
        if (local.n_prim_V <= 0)
            return;
        if (local.n_leptons != 1)
            return;

        if (local.n_electrons == 1) {
            if (local.n_veto_electrons != 1 || local.n_veto_muons != 0 ||
                local.pass_single_e != 1)
                return;
            if (local.n_sl_jets < min_njets || local.n_sl_btags < min_nbtags)
                return;
            local.event_selection_SL = true;
            local.is_e = true;
            ++sl_e;
        }

        else if (local.n_muons == 1) {
            if (local.n_veto_muons != 1 || local.n_veto_electrons != 0 ||
                local.pass_single_mu != 1)
                return;
            if (local.n_sl_jets < min_njets || local.n_sl_btags < min_nbtags)
                return;
            local.event_selection_SL = true;
            local.is_mu = true;
            ++sl_mu;
        }
    }
}

void CU_ttH_EDA::Check_DL_Event_Selection(CU_ttH_EDA_event_vars &local)
{
    if (local.pass_double_e == 1 || local.pass_double_mu == 1 ||
        local.pass_elemu == 1) {
        
        if (local.MET_filters ==0 || local.filterbadChCandidate == 0 || local.filterbadPFMuon == 0)
        	return;
        if (local.n_prim_V <= 0)
            return;
        if (local.n_di_leptons != 2)
            return;

        if (local.n_di_electrons == 2) { // di e
            if ((local.e_di_selected_sorted[0].pdgId() *
                     local.e_di_selected_sorted[1].pdgId() >
                 0) ||
                local.pass_double_e != 1)
                return;
            if (local.e_di_selected_sorted[0].pt() <= min_di_ele1_pT)
                return;
            if (local.n_di_jets < min_di_njets ||
                local.n_di_btags < min_di_nbtags)
                return;
            if (local.jets_di_selected_sorted[0].pt() <= min_jet_pT ||
                local.jets_di_selected_sorted[1].pt() <= min_jet_pT)
                return;
            E = local.e_di_selected_sorted[0].energy() +
                local.e_di_selected_sorted[1].energy();
            px = local.e_di_selected_sorted[0].px() +
                 local.e_di_selected_sorted[1].px();
            py = local.e_di_selected_sorted[0].py() +
                 local.e_di_selected_sorted[1].py();
            pz = local.e_di_selected_sorted[0].pz() +
                 local.e_di_selected_sorted[1].pz();
            p = sqrt(px * px + py * py + pz * pz);
            local.mll = sqrt(E * E - p * p);
            if (local.met_pt >= min_di_met)
                local.met_passed = 1;
            if (local.mll > min_di_mll) {
                if (local.mll < 76 || local.mll > 106)
                    local.mll_passed = 1;
            }
            if (local.met_passed == 1 && local.mll_passed == 1) {
                local.event_selection_DL = true;
                local.is_ee = true;
                ++dl_ee;
            }
        }

        else if (local.n_di_muons == 2) { // di mu
            if ((local.mu_di_selected_sorted[0].pdgId() *
                     local.mu_di_selected_sorted[1].pdgId() >
                 0) ||
                local.pass_double_mu != 1)
                return;
            if (local.mu_di_selected_sorted[0].pt() <= min_di_mu1_pT)
                return;
            if (local.n_di_jets < min_di_njets ||
                local.n_di_btags < min_di_nbtags)
                return;
            if (local.jets_di_selected_sorted[0].pt() <= min_jet_pT ||
                local.jets_di_selected_sorted[1].pt() <= min_jet_pT)
                return;
            E = local.mu_di_selected_sorted[0].energy() +
                local.mu_di_selected_sorted[1].energy();
            px = local.mu_di_selected_sorted[0].px() +
                 local.mu_di_selected_sorted[1].px();
            py = local.mu_di_selected_sorted[0].py() +
                 local.mu_di_selected_sorted[1].py();
            pz = local.mu_di_selected_sorted[0].pz() +
                 local.mu_di_selected_sorted[1].pz();
            p = sqrt(px * px + py * py + pz * pz);
            local.mll = sqrt(E * E - p * p);
            if (local.met_pt >= min_di_met)
                local.met_passed = 1;
            if (local.mll > min_di_mll) {
                if (local.mll < 76 || local.mll > 106)
                    local.mll_passed = 1;
            }
            if (local.met_passed == 1 && local.mll_passed == 1) {
                local.event_selection_DL = true;
                local.is_mumu = true;
                ++dl_mumu;
            }
        }

        else if (local.n_di_muons == 1 &&
                 local.n_di_electrons == 1) { // 1 e 1 mu
            if ((local.mu_di_selected_sorted[0].pdgId() *
                     local.e_di_selected_sorted[0].pdgId() >
                 0) ||
                local.pass_elemu != 1)
                return;
            if ((local.mu_di_selected_sorted[0].pt() <= min_di_mu1_pT) &&
                (local.e_di_selected_sorted[0].pt() <= min_di_ele1_pT))
                return;
            if (local.n_di_jets < min_di_njets ||
                local.n_di_btags < min_di_nbtags)
                return;
            if (local.jets_di_selected_sorted[0].pt() <= min_jet_pT ||
                local.jets_di_selected_sorted[1].pt() <= min_jet_pT)
                return;
            E = local.e_di_selected_sorted[0].energy() +
                local.mu_di_selected_sorted[0].energy();
            px = local.e_di_selected_sorted[0].px() +
                 local.mu_di_selected_sorted[0].px();
            py = local.e_di_selected_sorted[0].py() +
                 local.mu_di_selected_sorted[0].py();
            pz = local.e_di_selected_sorted[0].pz() +
                 local.mu_di_selected_sorted[0].pz();
            p = sqrt(px * px + py * py + pz * pz);
            local.mll = sqrt(E * E - p * p);
            local.met_passed = 1;
            if (local.mll > min_di_mll) {
                local.mll_passed = 1;
            }
            if (local.met_passed == 1 && local.mll_passed == 1) {
                local.event_selection_DL = true;
                local.is_emu = true;
                ++dl_emu;
            }
        }
    }
}

void CU_ttH_EDA::Fill_addn_quant(CU_ttH_EDA_event_vars &local,
                                 const edm::Event &iEvent,
                                 const edm::EventSetup &iSetup,
                                 const double &rho, const edm_Handles &handle)
{

    getJECSF(local, rho, handle); // to get JEC scale factors

    getbweight(local); // to get b-weight

    getLeptonSF(local); // to get Lepton ID, Iso and Trigger SFs

    // ttHf Category
    local.ttHf_cat = -1;
    if (handle.genTtbarId.isValid())
        local.ttHf_cat = *handle.genTtbarId % 100;
   
    // Generator Weight
    local.gen_weight = -1;
    if (!isdata)
    	local.gen_weight = handle.event_gen_info->weight();

    // PDF Weight
    if (!isdata)
        getPDFweight(local, handle.event_gen_info);

    // PU Weight
    if (!isdata)
        local.PU_weight = getPUweight(handle.PupInfo, local);

    // Q2 Weight
    if (!isdata) {
        local.q2_weight_up =
            getQ2weight(handle.event_gen_info, handle.EvtHandle, "1005");
        local.q2_weight_down =
            getQ2weight(handle.event_gen_info, handle.EvtHandle, "1009");
    }
    
    if (local.n_prim_V == 1) {
    	set_bjetness_input(local, handle.vertices);
    	get_bjetness_vars(local.jets_inp_bjetness, local.PV, *handle.ttrkbuilder, handle.electrons_for_mva, handle.muon_h, local.bjetnessFV_num_leps, local.bjetnessFV_npvTrkOVcollTrk, local.bjetnessFV_avip3d_val, local.bjetnessFV_avip3d_sig, local.bjetnessFV_avsip3d_sig, local.bjetnessFV_avip1d_sig);
     }
}

inline void CU_ttH_EDA::getJECSF(CU_ttH_EDA_event_vars &local, const double &rho, const edm_Handles &handle)
{

    int jet_index = 0;
    int jet_index2 = 0;
    if (local.event_selection_SL != 0) {

        double temp_pt = local.jets_sl_selected[0].pt();
        for (int i = 0; i < local.n_sl_jets; ++i) {
            if (local.jets_sl_selected[i].pt() > temp_pt) {
                jet_index = i;
                temp_pt = local.jets_sl_selected[i].pt();
            }
        }
	double temp2_pt = local.jets_sl_selected[0].pt();
        for (int i = 0; i < local.n_sl_jets; ++i) {
	    if (i==jet_index)
		continue;		
            if (local.jets_sl_selected[i].pt() > temp2_pt ) {
                jet_index2 = i;
                temp2_pt = local.jets_sl_selected[i].pt();
            }
        }
        // Jet SF
        pat::Jet jet1 = local.jets_sl_selected_raw[jet_index];
	pat::Jet jet2 = local.jets_sl_selected_raw[jet_index2];
        /*
        local.jet1SF_sl = miniAODhelper.GetJetCorrectionFactor(jet, iEvent,
        iSetup, handle.genjets, r, sysType::NA, 1, 0);
        local.jet1SF_up_sl = miniAODhelper.GetJetCorrectionFactor(jet,iEvent,
        iSetup, handle.genjets, r, sysType::JESup, 1, 0)/local.jet1SF_sl;
        local.jet1SF_down_sl = miniAODhelper.GetJetCorrectionFactor(jet,iEvent,
        iSetup, handle.genjets, r, sysType::JESdown, 1, 0)/local.jet1SF_sl;
        */

        local.jet1SF_sl = GetJetSF(jet1, sysType::NA, rho);
        local.jet2SF_sl = GetJetSF(jet2, sysType::NA, rho);

        if (!isdata){
            local.jet1SF_up_sl = GetJetSF(jet1, sysType::JESup, rho);
            local.jet1SF_down_sl = GetJetSF(jet1, sysType::JESdown, rho);
            local.jet2SF_up_sl = GetJetSF(jet2, sysType::JESup, rho);
            local.jet2SF_down_sl = GetJetSF(jet2, sysType::JESdown, rho);
        }
        else {
            local.jet1SF_up_sl = -1;
            local.jet1SF_down_sl = -1;
            local.jet2SF_up_sl = -1;
            local.jet2SF_down_sl = -1;
        }

    }

    else if (local.event_selection_DL != 0) {

        double temp_pt = local.jets_di_selected[0].pt();
        for (int i = 0; i < local.n_di_jets; ++i) {
            if (local.jets_di_selected[i].pt() > temp_pt) {
                jet_index = i;
                temp_pt = local.jets_di_selected[i].pt();
            }
        }
	double temp2_pt = local.jets_di_selected[0].pt();
        for (int i = 0; i < local.n_di_jets; ++i) {
	    if (i==jet_index)
		continue;		
            if (local.jets_di_selected[i].pt() > temp2_pt ) {
                jet_index2 = i;
                temp2_pt = local.jets_di_selected[i].pt();
            }
        }
        // Jet SF
        pat::Jet jet1 = local.jets_di_selected_raw[jet_index];
        pat::Jet jet2 = local.jets_di_selected_raw[jet_index2];
        /*
        local.jet1SF_di = miniAODhelper.GetJetCorrectionFactor(jet, iEvent,
        iSetup, handle.genjets, r, sysType::NA, 1, 0);
        local.jet1SF_up_di = miniAODhelper.GetJetCorrectionFactor(jet,iEvent,
        iSetup, handle.genjets, r, sysType::JESup, 1, 0)/local.jet1SF_di;
        local.jet1SF_down_di = miniAODhelper.GetJetCorrectionFactor(jet,iEvent,
        iSetup, handle.genjets, r, sysType::JESdown, 1, 0)/local.jet1SF_di;
        */

        local.jet1SF_di = GetJetSF(jet1, sysType::NA, rho);
        local.jet2SF_di = GetJetSF(jet2, sysType::NA, rho);

        if (!isdata){
            local.jet1SF_up_di = GetJetSF(jet1, sysType::JESup, rho);
            local.jet1SF_down_di = GetJetSF(jet1, sysType::JESdown, rho);
            local.jet2SF_up_di = GetJetSF(jet2, sysType::JESup, rho);
            local.jet2SF_down_di = GetJetSF(jet2, sysType::JESdown, rho);
        }
        else {
            local.jet1SF_up_di = -1;
            local.jet1SF_down_di = -1;
            local.jet2SF_up_di = -1;
            local.jet2SF_down_di = -1;
        }
    }
}

inline void CU_ttH_EDA::getLeptonSF(CU_ttH_EDA_event_vars &local)
{

    if (local.event_selection_SL != 0) {

        if (local.is_e) {
            local.lep_sf_id_sl = leptonSFhelper.GetElectronSF(
                local.e_selected[0].pt(),
                local.e_selected[0].superCluster()->position().eta(), 0, "ID");
            local.lep_sf_iso_sl = leptonSFhelper.GetElectronSF(
                local.e_selected[0].pt(),
                local.e_selected[0].superCluster()->position().eta(), 0, "Iso");
            local.lep_sf_gsf_sl = leptonSFhelper.GetElectronSF(
                local.e_selected[0].pt(),
                local.e_selected[0].superCluster()->position().eta(), 0, "Gsf");
            local.lep_sf_hip_sl = -1;
            local.lep_sf_trig_sl = leptonSFhelper.GetElectronSF(
                local.e_selected[0].pt(),
                local.e_selected[0].superCluster()->position().eta(), 0, "Trigger");
        } else if (local.is_mu) {
            local.lep_sf_id_sl = leptonSFhelper.GetMuonSF(
                local.mu_selected[0].pt(), local.mu_selected[0].eta(), 0, "ID");
            local.lep_sf_iso_sl =
                leptonSFhelper.GetMuonSF(local.mu_selected[0].pt(),
                                         local.mu_selected[0].eta(), 0, "Iso");
            local.lep_sf_gsf_sl = -1;
            local.lep_sf_hip_sl = leptonSFhelper.GetMuonSF(local.mu_selected[0].pt(),
                                         local.mu_selected[0].eta(), 0, "HIP");                                                            
            local.lep_sf_trig_sl = leptonSFhelper.GetMuonSF(
                local.mu_selected[0].pt(), local.mu_selected[0].eta(), 0,
                "Trigger");
        }
    }

    if (local.event_selection_DL != 0) {

        if (local.is_ee) {
            local.lep_sf_id_di =
                leptonSFhelper.GetElectronSF(
                    local.e_di_selected[0].pt(),
                    local.e_di_selected[0].superCluster()->position().eta(), 0,
                    "ID") *
                leptonSFhelper.GetElectronSF(
                    local.e_di_selected[1].pt(),
                    local.e_di_selected[1].superCluster()->position().eta(), 0,
                    "ID");
            local.lep_sf_iso_di =
                leptonSFhelper.GetElectronSF(
                    local.e_di_selected[0].pt(),
                    local.e_di_selected[0].superCluster()->position().eta(), 0,
                    "Iso") *
                leptonSFhelper.GetElectronSF(
                    local.e_di_selected[1].pt(),
                    local.e_di_selected[1].superCluster()->position().eta(), 0,
                    "Iso");
            local.lep_sf_gsf_di =
                leptonSFhelper.GetElectronSF(
                    local.e_di_selected[0].pt(),
                    local.e_di_selected[0].superCluster()->position().eta(), 0,
                    "Gsf") *
                leptonSFhelper.GetElectronSF(
                    local.e_di_selected[1].pt(),
                    local.e_di_selected[1].superCluster()->position().eta(), 0,
                    "Gsf");
            local.lep_sf_hip_di = -1;
            local.lep_sf_trig_di = leptonSFhelper.GetElectronElectronSF(
                local.e_di_selected[0].eta(), local.e_di_selected[1].eta(), 0,
                "Trigger");
        } else if (local.is_mumu) {
            local.lep_sf_id_di =
                leptonSFhelper.GetMuonSF(local.mu_di_selected[0].pt(),
                                         local.mu_di_selected[0].eta(), 0,
                                         "ID") *
                leptonSFhelper.GetMuonSF(local.mu_di_selected[1].pt(),
                                         local.mu_di_selected[1].eta(), 0,
                                         "ID");
            local.lep_sf_iso_di =
                leptonSFhelper.GetMuonSF(local.mu_di_selected[0].pt(),
                                         local.mu_di_selected[0].eta(), 0,
                                         "Iso") *
                leptonSFhelper.GetMuonSF(local.mu_di_selected[1].pt(),
                                         local.mu_di_selected[1].eta(), 0,
                                         "Iso");
            local.lep_sf_gsf_di = -1;
            local.lep_sf_hip_di =
                leptonSFhelper.GetMuonSF(local.mu_di_selected[0].pt(),
                                         local.mu_di_selected[0].eta(), 0,
                                         "HIP") *
                leptonSFhelper.GetMuonSF(local.mu_di_selected[1].pt(),
                                         local.mu_di_selected[1].eta(), 0,
                                         "HIP");
            local.lep_sf_trig_di = leptonSFhelper.GetMuonMuonSF(
                local.mu_di_selected[0].eta(), local.mu_di_selected[1].eta(), 0,
                "Trigger");
        } else if (local.is_emu) {
            local.lep_sf_id_di =
                leptonSFhelper.GetElectronSF(
                    local.e_di_selected[0].pt(),
                    local.e_di_selected[0].superCluster()->position().eta(), 0,
                    "ID") *
                leptonSFhelper.GetMuonSF(local.mu_di_selected[0].pt(),
                                         local.mu_di_selected[0].eta(), 0,
                                         "ID");
            local.lep_sf_iso_di =
                leptonSFhelper.GetElectronSF(
                    local.e_di_selected[0].pt(),
                    local.e_di_selected[0].superCluster()->position().eta(), 0,
                    "Iso") *
                leptonSFhelper.GetMuonSF(local.mu_di_selected[0].pt(),
                                         local.mu_di_selected[0].eta(), 0,
                                         "Iso");
            local.lep_sf_gsf_di = -1;
            local.lep_sf_hip_di = -1;
            local.lep_sf_trig_di = leptonSFhelper.GetElectronMuonSF(
                local.e_di_selected[0].eta(), local.mu_di_selected[0].eta(), 0,
                "Trigger");
        }
    }
}

inline double
CU_ttH_EDA::getPUweight(edm::Handle<std::vector<PileupSummaryInfo>> PupInfo, CU_ttH_EDA_event_vars &local)
{
    double pu_weight = -1;
    double numTruePV = -1;
    if ((PupInfo.isValid())) {
        for (std::vector<PileupSummaryInfo>::const_iterator PVI =
                 PupInfo->begin();
             PVI != PupInfo->end(); ++PVI) {
            int BX = PVI->getBunchCrossing();
            if (BX == 0) {
                numTruePV = PVI->getTrueNumInteractions();
            }
        }
    }
    local.truenpv =  numTruePV;
    for (int i = 0; i < 75; ++i) {
        if (numTruePV < (PU_x[i] + 1)) {
            pu_weight = PU_y[i];
            break;
        }
    }
    return pu_weight;
}

inline double
CU_ttH_EDA::getQ2weight(const edm::Handle<GenEventInfoProduct> &event_gen_info,
                        const edm::Handle<LHEEventProduct> &EvtHandle,
                        const string &ud)
{
    double theWeight;
    theWeight = event_gen_info->weight();
    unsigned int i;
    for (i = 0; i < EvtHandle->weights().size(); ++i) {
        if (!(ud.compare(EvtHandle->weights()[i].id)))
            theWeight *=
                EvtHandle->weights()[i].wgt / EvtHandle->originalXWGTUP();
    }
    return theWeight;
}

void CU_ttH_EDA::fillCSVHistos(TFile *fileHF, TFile *fileLF)
{
    for (int iSys = 0; iSys < 9; ++iSys) {
        for (int iPt = 0; iPt < 5; ++iPt)
            h_csv_wgt_hf[iSys][iPt] = NULL;
        for (int iPt = 0; iPt < 3; ++iPt) {
            for (int iEta = 0; iEta < 3; ++iEta)
                h_csv_wgt_lf[iSys][iPt][iEta] = NULL;
        }
    }
    for (int iSys = 0; iSys < 5; ++iSys) {
        for (int iPt = 0; iPt < 5; ++iPt)
            c_csv_wgt_hf[iSys][iPt] = NULL;
    }

    // CSV reweighting /// only care about the nominal ones
    for (int iSys = 0; iSys < 9; ++iSys) {
        TString syst_csv_suffix_hf = "final";
        TString syst_csv_suffix_c = "final";
        TString syst_csv_suffix_lf = "final";

        switch (iSys) {
        case 0:
            // this is the nominal case
            break;
        case 1:
            // JESUp
            syst_csv_suffix_hf = "final_JESUp";
            syst_csv_suffix_lf = "final_JESUp";
            syst_csv_suffix_c = "final_cErr1Up";
            break;
        case 2:
            // JESDown
            syst_csv_suffix_hf = "final_JESDown";
            syst_csv_suffix_lf = "final_JESDown";
            syst_csv_suffix_c = "final_cErr1Down";
            break;
        case 3:
            // purity up
            syst_csv_suffix_hf = "final_LFUp";
            syst_csv_suffix_lf = "final_HFUp";
            syst_csv_suffix_c = "final_cErr2Up";
            break;
        case 4:
            // purity down
            syst_csv_suffix_hf = "final_LFDown";
            syst_csv_suffix_lf = "final_HFDown";
            syst_csv_suffix_c = "final_cErr2Down";
            break;
        case 5:
            // stats1 up
            syst_csv_suffix_hf = "final_Stats1Up";
            syst_csv_suffix_lf = "final_Stats1Up";
            break;
        case 6:
            // stats1 down
            syst_csv_suffix_hf = "final_Stats1Down";
            syst_csv_suffix_lf = "final_Stats1Down";
            break;
        case 7:
            // stats2 up
            syst_csv_suffix_hf = "final_Stats2Up";
            syst_csv_suffix_lf = "final_Stats2Up";
            break;
        case 8:
            // stats2 down
            syst_csv_suffix_hf = "final_Stats2Down";
            syst_csv_suffix_lf = "final_Stats2Down";
            break;
        }

        for (int iPt = 0; iPt < 5; ++iPt)
            h_csv_wgt_hf[iSys][iPt] = (TH1D *)fileHF->Get(
                Form("csv_ratio_Pt%i_Eta0_%s", iPt, syst_csv_suffix_hf.Data()));

        if (iSys < 5) {
            for (int iPt = 0; iPt < 5; ++iPt)
                c_csv_wgt_hf[iSys][iPt] = (TH1D *)fileHF->Get(Form(
                    "c_csv_ratio_Pt%i_Eta0_%s", iPt, syst_csv_suffix_c.Data()));
        }

        for (int iPt = 0; iPt < 4; ++iPt) {
            for (int iEta = 0; iEta < 3; ++iEta)
                h_csv_wgt_lf[iSys][iPt][iEta] =
                    (TH1D *)fileLF->Get(Form("csv_ratio_Pt%i_Eta%i_%s", iPt,
                                             iEta, syst_csv_suffix_lf.Data()));
        }
    }

    return;
}

double CU_ttH_EDA::getCSVWeight(std::vector<double> jetPts,
                                std::vector<double> jetEtas,
                                std::vector<double> jetCSVs,
                                std::vector<int> jetFlavors, int iSys,
                                double &csvWgtHF, double &csvWgtLF,
                                double &csvWgtCF)
{
    int iSysHF = 0;
    switch (iSys) {
    case 7:
        iSysHF = 1;
        break; // JESUp
    case 8:
        iSysHF = 2;
        break; // JESDown
    case 9:
        iSysHF = 3;
        break; // LFUp
    case 10:
        iSysHF = 4;
        break; // LFDown
    case 13:
        iSysHF = 5;
        break; // Stats1Up
    case 14:
        iSysHF = 6;
        break; // Stats1Down
    case 15:
        iSysHF = 7;
        break; // Stats2Up
    case 16:
        iSysHF = 8;
        break; // Stats2Down
    default:
        iSysHF = 0;
        break; // NoSys
    }

    int iSysC = 0;
    switch (iSys) {
    case 21:
        iSysC = 1;
        break;
    case 22:
        iSysC = 2;
        break;
    case 23:
        iSysC = 3;
        break;
    case 24:
        iSysC = 4;
        break;
    default:
        iSysC = 0;
        break;
    }

    int iSysLF = 0;
    switch (iSys) {
    case 7:
        iSysLF = 1;
        break; // JESUp
    case 8:
        iSysLF = 2;
        break; // JESDown
    case 11:
        iSysLF = 3;
        break; // HFUp
    case 12:
        iSysLF = 4;
        break; // HFDown
    case 17:
        iSysLF = 5;
        break; // Stats1Up
    case 18:
        iSysLF = 6;
        break; // Stats1Down
    case 19:
        iSysLF = 7;
        break; // Stats2Up
    case 20:
        iSysLF = 8;
        break; // Stats2Down
    default:
        iSysLF = 0;
        break; // NoSys
    }

    double csvWgthf = 1.;
    double csvWgtC = 1.;
    double csvWgtlf = 1.;

    for (int iJet = 0; iJet < int(jetPts.size()); ++iJet) {

        double csv = jetCSVs[iJet];
        double jetPt = jetPts[iJet];
        double jetAbsEta = fabs(jetEtas[iJet]);
        int flavor = jetFlavors[iJet];

        int iPt = -1;
        int iEta = -1;
        if (jetPt >= 19.99 && jetPt < 30)
            iPt = 0;
        else if (jetPt >= 30 && jetPt < 40)
            iPt = 1;
        else if (jetPt >= 40 && jetPt < 60)
            iPt = 2;
        else if (jetPt >= 60 && jetPt < 100)
            iPt = 3;
        else if (jetPt >= 100)
            iPt = 4;

        if (jetAbsEta >= 0 && jetAbsEta < 0.8)
            iEta = 0;
        else if (jetAbsEta >= 0.8 && jetAbsEta < 1.6)
            iEta = 1;
        else if (jetAbsEta >= 1.6 && jetAbsEta < 2.41)
            iEta = 2;

        if (iPt < 0 || iEta < 0)
            std::cout << "Error, couldn't find Pt, Eta bins for this b-flavor "
                         "jet, jetPt = "
                      << jetPt << ", jetAbsEta = " << jetAbsEta << std::endl;

        if (abs(flavor) == 5) {
            int useCSVBin =
                (csv >= 0.) ? h_csv_wgt_hf[iSysHF][iPt]->FindBin(csv) : 1;
            double iCSVWgtHF =
                h_csv_wgt_hf[iSysHF][iPt]->GetBinContent(useCSVBin);
            if (iCSVWgtHF != 0)
                csvWgthf *= iCSVWgtHF;
        } else if (abs(flavor) == 4) {
            int useCSVBin =
                (csv >= 0.) ? c_csv_wgt_hf[iSysC][iPt]->FindBin(csv) : 1;
            double iCSVWgtC =
                c_csv_wgt_hf[iSysC][iPt]->GetBinContent(useCSVBin);
            if (iCSVWgtC != 0)
                csvWgtC *= iCSVWgtC;
        } else {
            if (iPt >= 3)
                iPt =
                    3; /// [30-40], [40-60] and [60-10000] only 3 Pt bins for lf
            int useCSVBin =
                (csv >= 0.) ? h_csv_wgt_lf[iSysLF][iPt][iEta]->FindBin(csv) : 1;
            double iCSVWgtLF =
                h_csv_wgt_lf[iSysLF][iPt][iEta]->GetBinContent(useCSVBin);
            if (iCSVWgtLF != 0)
                csvWgtlf *= iCSVWgtLF;
        }
    }

    double csvWgtTotal = csvWgthf * csvWgtC * csvWgtlf;

    csvWgtHF = csvWgthf;
    csvWgtLF = csvWgtlf;
    csvWgtCF = csvWgtC;

    return csvWgtTotal;
}

inline void CU_ttH_EDA::getbweight(CU_ttH_EDA_event_vars &local)
{
    double csvWgtHF, csvWgtLF, csvWgtCF;

    if (local.event_selection_SL == 1) {
        for (std::vector<pat::Jet>::const_iterator iJet =
                 local.jets_sl_selected.begin();
             iJet != local.jets_sl_selected.end(); ++iJet) {
            local.vec_jet_pt.push_back(iJet->pt());
            local.vec_jet_eta.push_back(iJet->eta());
            local.vec_jet_csv.push_back(miniAODhelper.GetJetCSV(
                *iJet, "pfCombinedInclusiveSecondaryVertexV2BJetTags"));
            local.vec_jet_hadronFlavour.push_back(iJet->hadronFlavour());
        }

	csvWgtHF = csvWgtLF = csvWgtCF = 0;
        local.b_weight_sl =
            getCSVWeight(local.vec_jet_pt, local.vec_jet_eta, local.vec_jet_csv,
                         local.vec_jet_hadronFlavour, local.iSys[0], csvWgtHF,
                         csvWgtLF, csvWgtCF);
	csvWgtHF = csvWgtLF = csvWgtCF = 0;
	local.b_weight_sl_lfup =
            getCSVWeight(local.vec_jet_pt, local.vec_jet_eta, local.vec_jet_csv,
                         local.vec_jet_hadronFlavour, local.iSys[1], csvWgtHF,
                         csvWgtLF, csvWgtCF);
	csvWgtHF = csvWgtLF = csvWgtCF = 0;
	local.b_weight_sl_hfdown =
            getCSVWeight(local.vec_jet_pt, local.vec_jet_eta, local.vec_jet_csv,
                         local.vec_jet_hadronFlavour, local.iSys[2], csvWgtHF,
                         csvWgtLF, csvWgtCF);
	csvWgtHF = csvWgtLF = csvWgtCF = 0;
	local.b_weight_sl_cErr1_down =
            getCSVWeight(local.vec_jet_pt, local.vec_jet_eta, local.vec_jet_csv,
                         local.vec_jet_hadronFlavour, local.iSys[3], csvWgtHF,
                         csvWgtLF, csvWgtCF);
        local.b_weight_di = local.b_weight_di_lfup = local.b_weight_di_hfdown = local.b_weight_di_cErr1_down = -1;
    }

    else if (local.event_selection_DL == 1) {
        for (std::vector<pat::Jet>::const_iterator iJet =
                 local.jets_di_selected.begin();
             iJet != local.jets_di_selected.end(); ++iJet) {
            local.vec_jet_pt.push_back(iJet->pt());
            local.vec_jet_eta.push_back(iJet->eta());
            local.vec_jet_csv.push_back(miniAODhelper.GetJetCSV(
                *iJet, "pfCombinedInclusiveSecondaryVertexV2BJetTags"));
            local.vec_jet_hadronFlavour.push_back(iJet->hadronFlavour());
        }

        csvWgtHF = csvWgtLF = csvWgtCF = 0;
        local.b_weight_di =
            getCSVWeight(local.vec_jet_pt, local.vec_jet_eta, local.vec_jet_csv,
                         local.vec_jet_hadronFlavour, local.iSys[0], csvWgtHF,
                         csvWgtLF, csvWgtCF);
	csvWgtHF = csvWgtLF = csvWgtCF = 0;
	local.b_weight_di_lfup =
            getCSVWeight(local.vec_jet_pt, local.vec_jet_eta, local.vec_jet_csv,
                         local.vec_jet_hadronFlavour, local.iSys[1], csvWgtHF,
                         csvWgtLF, csvWgtCF);
	csvWgtHF = csvWgtLF = csvWgtCF = 0;
	local.b_weight_di_hfdown =
            getCSVWeight(local.vec_jet_pt, local.vec_jet_eta, local.vec_jet_csv,
                         local.vec_jet_hadronFlavour, local.iSys[2], csvWgtHF,
                         csvWgtLF, csvWgtCF);
	csvWgtHF = csvWgtLF = csvWgtCF = 0;
	local.b_weight_di_cErr1_down =
            getCSVWeight(local.vec_jet_pt, local.vec_jet_eta, local.vec_jet_csv,
                         local.vec_jet_hadronFlavour, local.iSys[3], csvWgtHF,
                         csvWgtLF, csvWgtCF);
        local.b_weight_sl = local.b_weight_sl_lfup = local.b_weight_sl_hfdown = local.b_weight_sl_cErr1_down = -1;
    }
}

inline void
CU_ttH_EDA::getPDFweight(CU_ttH_EDA_event_vars &local,
                         const edm::Handle<GenEventInfoProduct> &genInfos)
{

    auto pdfInfos = genInfos->pdf();
    double pdfNominal = pdfInfos->xPDF.first * pdfInfos->xPDF.second;

    std::vector<double> pdfs;
    for (size_t j = 0; j < NNPDF30_nlo_as_0118_PDFSet->size(); ++j) {
        double xpdf1 = _systPDFs[j]->xfxQ(pdfInfos->id.first, pdfInfos->x.first,
                                          pdfInfos->scalePDF);
        double xpdf2 = _systPDFs[j]->xfxQ(
            pdfInfos->id.second, pdfInfos->x.second, pdfInfos->scalePDF);
        pdfs.push_back(xpdf1 * xpdf2);
    }

    const LHAPDF::PDFUncertainty pdfUnc =
        NNPDF30_nlo_as_0118_PDFSet->uncertainty(pdfs, 68.);

    double weight_up = 1.0;
    double weight_down = 1.0;
    if (std::isfinite(1. / pdfNominal)) {
        weight_up = (pdfUnc.central + pdfUnc.errplus) / pdfNominal;
        weight_down = (pdfUnc.central - pdfUnc.errminus) / pdfNominal;
    }
    local.pdf_weight_up = weight_up;
    local.pdf_weight_down = weight_down;
}

inline void 
CU_ttH_EDA::set_bjetness_input(CU_ttH_EDA_event_vars &local, const edm::Handle<reco::VertexCollection> &vertices)
{
	local.PV = *vertices->begin();	
	std::vector<pat::Jet> jets;
	pat::Jet temp_jet;
	double n_jets=0;
	double n_max;
	if(local.event_selection_SL == true) {
		 jets = local.jets_sl_selected;
		 n_jets = local.n_sl_jets;
	}
	else if (local.event_selection_DL == true) {	
		jets = local.jets_di_selected;
		n_jets = local.n_di_jets;
	}
	
	
	for(int i=0; i<n_jets; i++) {
		for (int j=(n_jets-1); j>i; j--) {
			if ( miniAODhelper.GetJetCSV(jets[j],"pfCombinedInclusiveSecondaryVertexV2BJetTags") > miniAODhelper.GetJetCSV(jets[j-1],"pfCombinedInclusiveSecondaryVertexV2BJetTags") ) {
			
			temp_jet = jets[j-1];
			jets[j-1] = jets[j];
			jets[j] = temp_jet;
                               
			}
		}
	}
	
	if(n_jets > 6) n_max = 6;
	else n_max = n_jets;
	
	for (int k=1; k<n_max; k++)
		local.jets_inp_bjetness.push_back(jets[k]);

}

inline vector<reco::TransientTrack> 
CU_ttH_EDA::get_ttrks(vector<reco::Track> trks, const TransientTrackBuilder& ttrkbuilder)
{
  vector<reco::TransientTrack> ttrks;
  for(uint tr=0; tr<trks.size(); tr++){
   reco::TransientTrack ttrk = ttrkbuilder.build(&trks[tr]);
   ttrks.push_back(ttrk);
  }
 return ttrks;
}

inline bool 
CU_ttH_EDA::is_goodtrk(reco::Track trk,const reco::Vertex& vtx)
{
	bool isgoodtrk = false;
 	if(trk.pt()>1 &&
   		trk.hitPattern().numberOfValidHits()>=8 &&
   		trk.hitPattern().numberOfValidPixelHits()>=2 &&
   		trk.normalizedChi2()<5 &&
   		std::abs(trk.dxy(vtx.position()))<0.2 &&
   		std::abs(trk.dz(vtx.position()))<17) 
   			isgoodtrk = true;
	
	return isgoodtrk;
}

inline bool 
CU_ttH_EDA::is_loosePOG_jetmuon(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Muon> > muon_h)
{
  bool ismu = false;
  for(const pat::Muon &mu : *muon_h){
    if(deltaR(jcand.p4(),mu.p4())<0.1 && fabs(jcand.pt()-mu.pt())/mu.pt()<0.05){
     if(mu.isLooseMuon()) ismu = true;
     if(ismu) break;
    }
  }  
  return ismu;    
}   
    
inline bool 
CU_ttH_EDA::is_softLep_jetelectron(const pat::PackedCandidate &jcand, edm::Handle<edm::View<pat::Electron> > electron_pat, const reco::Vertex& vtx)
{
  bool isele = false;
  for(edm::View<pat::Electron>::const_iterator ele = electron_pat->begin(); ele != electron_pat->end(); ele++){
    const pat::Electron &lele = *ele;
    if(deltaR(jcand.p4(),lele.p4())<0.1 && fabs(jcand.pt()-lele.pt())/lele.pt()<0.05 ){
      const reco::HitPattern &hitPattern = lele.gsfTrack().get()->hitPattern();
      uint32_t hit = hitPattern.getHitPattern(reco::HitPattern::TRACK_HITS, 0);
      bool hitCondition = !(reco::HitPattern::validHitFilter(hit) && ((reco::HitPattern::pixelBarrelHitFilter(hit) && reco::HitPattern::getLayer(hit) < 3) || reco::HitPattern::pixelEndcapHitFilter(hit)));
      if(!hitCondition && lele.passConversionVeto()) isele = true;
      if(isele) break;
    }
  }
  return isele;
}

inline void 
CU_ttH_EDA::get_bjetness_trkinfos(vector<pat::Jet> evtjets, const reco::Vertex& vtx, vector<reco::Track>& jetchtrks, double& bjetness_num_pvtrks, double& bjetness_num_npvtrks, edm::Handle<edm::View<pat::Electron> > electron_pat, edm::Handle<edm::View<pat::Muon> > muon_h, double& bjetness_num_eles, double& bjetness_num_mus, vector<tuple<double, double, double> >& jetsdir)
{
  //Loop over evt jet
  for(uint j=0; j<evtjets.size(); j++){
    pat::Jet jet = evtjets[j];
    //Access jet daughters
    vector<reco::CandidatePtr> jdaus(jet.daughterPtrVector());
    sort(jdaus.begin(), jdaus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) {return p1->pt() > p2->pt();});
    for(uint jd=0; jd<jdaus.size(); jd++){
      const pat::PackedCandidate &jcand = dynamic_cast<const pat::PackedCandidate &>(*jdaus[jd]);
      //dR requirement
      if(deltaR(jcand.p4(),jet.p4())>0.4) continue;
      reco::Track trk = reco::Track(jcand.pseudoTrack());
      bool isgoodtrk = is_goodtrk(trk,vtx);
      //Minimal conditions for a BJetness jet constituent 
      if(isgoodtrk && jcand.charge()!=0 && jcand.fromPV()>1){
        jetchtrks.push_back(trk);
        if(jcand.fromPV()==3) bjetness_num_pvtrks++;
        if(jcand.fromPV()==2) bjetness_num_npvtrks++;
        jetsdir.push_back(make_tuple(jet.px(),jet.py(),jet.pz()));
        if(fabs(jcand.pdgId())==13 && is_loosePOG_jetmuon(jcand,muon_h)) bjetness_num_mus++;
        if(fabs(jcand.pdgId())==11 && is_softLep_jetelectron(jcand,electron_pat,vtx)) bjetness_num_eles++;       
        //if(fabs(jcand.pdgId())==11 && is_loosePOGNoIPNoIso_jetelectron(jcand,electron_pat,vtx)) bjetness_num_eles++;
      }//Ch trks 
    }//Loop on jet daus 
  }//Loop on evt jet
}

inline void 
CU_ttH_EDA::get_avip3d(vector<reco::Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& jetsdir, double& jetchtrks_avip3d_val, double& jetchtrks_avip3d_sig, double& jetchtrks_avsip3d_sig)
{
  double valtemp = 0;
  vector<reco::TransientTrack> ttrks = get_ttrks(trks,ttrkbuilder);
  for(uint t=0; t<ttrks.size(); t++){
    valtemp = IPTools::absoluteImpactParameter3D(ttrks[t],vtx).second.value();
    if(valtemp==valtemp) jetchtrks_avip3d_val  += valtemp;
    valtemp = IPTools::absoluteImpactParameter3D(ttrks[t],vtx).second.significance();
    if(valtemp==valtemp) jetchtrks_avip3d_sig  += valtemp;
    GlobalVector jetsdirgv(get<0>(jetsdir[t]),get<1>(jetsdir[t]),get<2>(jetsdir[t]));
    valtemp = IPTools::signedImpactParameter3D(ttrks[t],jetsdirgv,vtx).second.significance();
    if(valtemp==valtemp) jetchtrks_avsip3d_sig += valtemp;
  }
}

inline void 
CU_ttH_EDA::get_avip1d(vector<reco::Track> trks, const TransientTrackBuilder& ttrkbuilder, reco::Vertex vtx, vector<tuple<double, double, double> >& jetsdir, double& jetchtrks_avip1d_sig)
{
  double valtemp = 0;
  vector<reco::TransientTrack> ttrks = get_ttrks(trks,ttrkbuilder);
  SignedTransverseImpactParameter stip;
  for(uint t=0; t<ttrks.size(); t++){
    GlobalVector jetsdirgv(get<0>(jetsdir[t]),get<1>(jetsdir[t]),get<2>(jetsdir[t]));
    valtemp = fabs(stip.zImpactParameter(ttrks[t],jetsdirgv,vtx).second.significance());
    if(valtemp==valtemp) jetchtrks_avip1d_sig  += valtemp;
  }
}

inline void 
CU_ttH_EDA::get_bjetness_vars( vector<pat::Jet> evtjets, const reco::Vertex& vtx, const TransientTrackBuilder& ttrkbuilder, edm::Handle<edm::View<pat::Electron> > electron_pat, edm::Handle<edm::View<pat::Muon> > muon_h, double& bjetnessFV_num_leps, double& bjetnessFV_npvTrkOVcollTrk, double& bjetnessFV_avip3d_val, double& bjetnessFV_avip3d_sig, double& bjetnessFV_avsip3d_sig, double& bjetnessFV_avip1d_sig )
{
	//Get BJetness trk info
  vector<reco::Track> jetschtrks; 
  jetschtrks.clear(); 
  double num_pvtrks  = 0;
  double num_npvtrks = 0;
  double num_eles    = 0;  
  double num_mus     = 0;           
  vector<tuple<double, double, double> > jetsdir; 
  jetsdir.clear(); 
  
  get_bjetness_trkinfos(evtjets, vtx, jetschtrks, num_pvtrks, num_npvtrks, electron_pat, muon_h, num_eles, num_mus, jetsdir);
  
  bjetnessFV_num_leps = num_eles+num_mus;
  if(jetschtrks.size()!=0){
    bjetnessFV_npvTrkOVcollTrk       = num_npvtrks/double(jetschtrks.size()); 
    //Get BJetness Impact Parameters
    double ip_valtemp = 0;
    //3D
    double jetchtrks_avip3d_val  = 0;
    double jetchtrks_avip3d_sig  = 0;
    double jetchtrks_avsip3d_sig = 0;
    
    get_avip3d(jetschtrks, ttrkbuilder, vtx, jetsdir, jetchtrks_avip3d_val,jetchtrks_avip3d_sig,jetchtrks_avsip3d_sig);
    ip_valtemp = jetchtrks_avip3d_val/jetschtrks.size();
    if(ip_valtemp==ip_valtemp) bjetnessFV_avip3d_val = ip_valtemp;
    else                       bjetnessFV_avip3d_val = -996;
    ip_valtemp = jetchtrks_avip3d_sig/jetschtrks.size();
    if(ip_valtemp==ip_valtemp) bjetnessFV_avip3d_sig = ip_valtemp;
    else                       bjetnessFV_avip3d_sig = -996; 
    ip_valtemp = jetchtrks_avsip3d_sig/jetschtrks.size();
    if(ip_valtemp==ip_valtemp) bjetnessFV_avsip3d_sig = ip_valtemp;
    else                       bjetnessFV_avsip3d_sig = -996;
    //1D
    double jetchtrks_avip1d_sig  = 0;
    get_avip1d(jetschtrks, ttrkbuilder, vtx, jetsdir, jetchtrks_avip1d_sig);
    ip_valtemp = jetchtrks_avip1d_sig/jetschtrks.size();
    if(ip_valtemp==ip_valtemp) bjetnessFV_avip1d_sig = ip_valtemp;
    else                       bjetnessFV_avip1d_sig = -996;    
  }else{
    bjetnessFV_npvTrkOVcollTrk       = -998;
    bjetnessFV_avip3d_val            = -998;
    bjetnessFV_avip3d_sig            = -998;
    bjetnessFV_avsip3d_sig           = -998;
    bjetnessFV_avip1d_sig            = -998;
  }
}


#endif
