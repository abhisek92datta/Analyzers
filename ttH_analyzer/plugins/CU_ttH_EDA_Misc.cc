#ifndef CU_ttH_EDA_Misc_cc
#define CU_ttH_EDA_Misc_cc

/// Includes
#include "CU_ttH_EDA.h"

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

void CU_ttH_EDA::Update_common_vars(const edm::Event &iEvent,
									CU_ttH_EDA_event_vars &local)
{
	local.run_nr = iEvent.id().run();
	local.event_nr = iEvent.id().event();
	local.lumisection_nr = iEvent.id().luminosityBlock();
}

int CU_ttH_EDA::Check_beam_spot(edm::Handle<reco::BeamSpot> BS)
{
	if (!BS.isValid())
		return 1;

	// 	math::XYZPoint beamSpotPosition;
	// 	beamSpotPosition.SetCoordinates(0, 0, 0);
	// //		double BSx = 0, BSy = 0, BSz = 0;
	//
	// 	double BSx = BS->x0();
	// 	double BSy = BS->y0();
	// 	double BSz = BS->z0();
	// 	beamSpotPosition = BS->position();
	//
	//		if (verbose_)
	//			printf("\t BeamSpot: x = %.2f,\t y = %.2f,\t z = %.2f \n",
	//				BSx, BSy, BSz );
	if (verbose_)
		printf("\t BeamSpot: x = %.2f,\t y = %.2f,\t z = %.2f \n", BS->x0(),
			   BS->y0(), BS->z0());

	return 0;
}

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

	return 0;
}

bool CU_ttH_EDA::Check_triggers_iterator(
	const vector<string> &triggers,
	edm::Handle<edm::TriggerResults> triggerResults)
{
	for (std::vector<std::string>::const_iterator trigger = triggers.begin();
		 trigger != triggers.end(); ++trigger) {
		 
		 std::string trigger_it;
		 char s[100];
		 unsigned int hltIndex;
		 for (int i=6; i>=1; i--) {
		 	trigger_it.assign(*trigger);
		 	sprintf(s,"%d",i);
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

int CU_ttH_EDA::Check_filters(edm::Handle<edm::TriggerResults> filterResults)
{
	if (!filterResults.isValid()) {
		std::cerr << "Trigger results not valid for tag " << filterTag
				  << std::endl;
		return 1;
	}

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

	return 0;
}

int CU_ttH_EDA::Check_vertices_set_MAODhelper(edm::Handle<reco::VertexCollection> vertices)
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

int CU_ttH_EDA::Check_PV(edm::Handle<reco::VertexCollection> vertices){
	reco::VertexCollection::const_iterator vtx = vertices->begin();
	if (vtx->isFake() || vtx->ndof() < 4.0 || abs(vtx->z()) > 24.0 || abs(vtx->position().Rho()) > 2.0)
		return 0;
	else 
		return 1;
}



/*
/// Taggers
int CU_ttH_EDA::Higgs_tagger(
	Handle<boosted::SubFilterJetCollection> subfilter_jets,
	CU_ttH_EDA_event_vars &local)
{
	local.n_Htags = 0;

	if (!subfilter_jets.isValid())
		return 1;

	boosted::SubFilterJetCollection subfilterjets =
		BoostedUtils::GetSortedByPt(*subfilter_jets);

	for (boosted::SubFilterJetCollection::iterator higgsJet =
			 subfilterjets.begin();
		 higgsJet != subfilterjets.end(); ++higgsJet) {
		// pt and eta requirements for top jet
		if (higgsJet->fatjet.pt() <= 250. || abs(higgsJet->fatjet.eta()) >= 1.8)
			continue;

		int numBtagFiltJets = 0;
		std::vector<pat::Jet> filterjets = higgsJet->filterjets;
		int numFiltJets = filterjets.size();
		for (int ijet = 0; ijet < numFiltJets; ++ijet) {
			if (verbose_) {
				printf("\t\t filt jet %2d:\t pT = %.1f,\t eta = %.2f,\t phi = "
					   "%.2f,\t CSVv2 = %+5.3f,\t CSVv1 = %+5.3f \n",
					   ijet, filterjets[ijet].pt(), filterjets[ijet].eta(),
					   filterjets[ijet].phi(),
					   filterjets[ijet].bDiscriminator(
						   "combinedInclusiveSecondaryVertexV2BJetTags"),
					   filterjets[ijet].bDiscriminator(
						   "combinedSecondaryVertexBJetTags"));
			}

			if (filterjets[ijet].pt() <= 20. ||
				abs(filterjets[ijet].eta()) >= 2.5)
				continue;

			// b-tag medium WP
			if (filterjets[ijet].bDiscriminator(
					"combinedInclusiveSecondaryVertexV2BJetTags") < 0.814)
				continue;

			++numBtagFiltJets;
		}

		if (verbose_) {
			printf("\t Higgs jet %2d:\t pT = %.1f,\t eta = %.2f,\t phi = "
				   "%.2f,\t numFiltJets = %2d,\t numBtagFiltJets = %2d\n",
				   int(higgsJet - subfilterjets.begin()), higgsJet->fatjet.pt(),
				   higgsJet->fatjet.eta(), higgsJet->fatjet.phi(), numFiltJets,
				   numBtagFiltJets);
		}

		if (numBtagFiltJets >= 2)
			++local.n_Htags;
	}

	return 0;
}

int CU_ttH_EDA::Top_tagger(Handle<boosted::HTTTopJetCollection> top_jets,
						   CU_ttH_EDA_event_vars &local)
{
	local.n_ttags = 0;

	if (!top_jets.isValid())
		return 1;

	boosted::HTTTopJetCollection heptopjets =
		BoostedUtils::GetSortedByPt(*top_jets);

	for (boosted::HTTTopJetCollection::iterator topJet = heptopjets.begin();
		 topJet != heptopjets.end(); ++topJet) {
		// pt and eta requirements on top jet
		if (topJet->fatjet.pt() <= 250. || abs(topJet->fatjet.eta()) >= 1.8)
			continue;

		// pt and eta requirements on subjets
		if (topJet->nonW.pt() <= 20 || abs(topJet->nonW.eta()) >= 2.5 ||
			topJet->W1.pt() <= 20 || abs(topJet->W1.eta()) >= 2.5 ||
			topJet->W2.pt() <= 20 || abs(topJet->W2.eta()) >= 2.5)
			continue;

		// must be top-tagged
		if (toptagger.GetTopTaggerOutput(*topJet)<=-1) 
			continue;

		++local.n_ttags;
	}

	return 0;
}
*/

/// Other functions

void CU_ttH_EDA::Check_Fill_Print_single_lepton(CU_ttH_EDA_event_vars &local)
{
	fprintf(events_combined, "%d,%d,%d,", local.run_nr, local.lumisection_nr, local.event_nr);	
	fprintf(events_combined, "%d,%d,%d,%d,%d,", local.is_e, local.is_mu, local.is_ee, local.is_emu, local.is_mumu);
	fprintf(events_combined, "%d,%d,", local.n_sl_jets, local.n_sl_btags );
	if (local.n_electrons == 1) {
		fprintf(events_combined, "%.4f,%.4f,%d,", local.e_selected[0].pt(), miniAODhelper.GetElectronRelIso(local.e_selected[0], coneSize::R03, corrType::rhoEA,effAreaType::spring15), local.e_selected[0].pdgId());
	}
	else if (local.n_muons == 1) {
		fprintf(events_combined, "%.4f,%.4f,%d,", local.mu_selected[0].pt(), miniAODhelper.GetMuonRelIso(local.mu_selected[0], coneSize::R04, corrType::deltaBeta), local.mu_selected[0].pdgId());
	}
	fprintf(events_combined, "-1,-1,-1,");
	fprintf(events_combined, "%.4f,%.4f,%.4f,%.4f,", local.jets_sl_selected_sorted[0].pt(), local.jets_sl_selected_sorted[1].pt(), miniAODhelper.GetJetCSV(local.jets_sl_selected_sorted[0],"pfCombinedInclusiveSecondaryVertexV2BJetTags"), miniAODhelper.GetJetCSV(local.jets_sl_selected_sorted[1],"pfCombinedInclusiveSecondaryVertexV2BJetTags"));
	fprintf(events_combined, "%.4f,%.4f,%.4f,", local.jet1SF_sl, local.jet1SF_up_sl, local.jet1SF_down_sl);
	fprintf(events_combined, "%.4f,%.4f,-1,", local.met_pt, local.met_phi);
	fprintf(events_combined, "%d,", local.ttHf_cat);
	if(!isdata)
		fprintf(events_combined, "%.4f,",local.PU_weight);
	else
		fprintf(events_combined, "-1,");
	if(!isdata){
		fprintf(events_combined, "%.4f,%.4f,", local.b_weight_sl, local.lep_sf_trig_sl);
		fprintf(events_combined, "%.4f,%.4f,", local.lep_sf_id_sl, local.lep_sf_iso_sl);
	}
	else
		fprintf(events_combined, "-1,-1,-1,-1,");
	if(!isdata)
		fprintf(events_combined, "%.4f,%.4f\n", local.q2_weight_up, local.q2_weight_down);
	else
		fprintf(events_combined, "-1,-1\n");	
	if(!isdata)
		fprintf(events_combined, "%.4f,%.4f\n", local.pdf_weight_up, local.pdf_weight_down);
	else
		fprintf(events_combined, "-1,-1\n");
}

void CU_ttH_EDA::Check_Fill_Print_di_lepton(CU_ttH_EDA_event_vars &local)
{
	fprintf(events_combined, "%d,%d,%d,", local.run_nr, local.lumisection_nr, local.event_nr);
	fprintf(events_combined, "%d,%d,%d,%d,%d,", local.is_e, local.is_mu, local.is_ee, local.is_emu, local.is_mumu);
	fprintf(events_combined, "%d,%d,", local.n_di_jets, local.n_di_btags );
	if (local.n_di_electrons == 2) {
		fprintf(events_combined, "%.4f,%.4f,%d,", local.e_di_selected_sorted[0].pt(), miniAODhelper.GetElectronRelIso(local.e_di_selected_sorted[0], coneSize::R03, corrType::rhoEA,effAreaType::spring15), local.e_di_selected_sorted[0].pdgId());
		fprintf(events_combined, "%.4f,%.4f,%d,", local.e_di_selected_sorted[1].pt(), miniAODhelper.GetElectronRelIso(local.e_di_selected_sorted[1], coneSize::R03, corrType::rhoEA,effAreaType::spring15), local.e_di_selected_sorted[1].pdgId());
	}
	else if (local.n_di_muons == 2) {
		fprintf(events_combined, "%.4f,%.4f,%d,", local.mu_di_selected_sorted[0].pt(), miniAODhelper.GetMuonRelIso(local.mu_di_selected_sorted[0], coneSize::R04, corrType::deltaBeta), local.mu_di_selected_sorted[0].pdgId());
		fprintf(events_combined, "%.4f,%.4f,%d,", local.mu_di_selected_sorted[1].pt(), miniAODhelper.GetMuonRelIso(local.mu_di_selected_sorted[1], coneSize::R04, corrType::deltaBeta), local.mu_di_selected_sorted[1].pdgId());
	}
	else {
		if (local.mu_di_selected_sorted[0].pt() > local.e_di_selected_sorted[0].pt()) {
			fprintf(events_combined, "%.4f,%.4f,%d,", local.mu_di_selected_sorted[0].pt(), miniAODhelper.GetMuonRelIso(local.mu_di_selected_sorted[0], coneSize::R04, corrType::deltaBeta), local.mu_di_selected_sorted[0].pdgId());
			fprintf(events_combined, "%.4f,%.4f,%d,", local.e_di_selected_sorted[0].pt(), miniAODhelper.GetElectronRelIso(local.e_di_selected_sorted[0], coneSize::R03, corrType::rhoEA,effAreaType::spring15), local.e_di_selected_sorted[0].pdgId());
		}
		else {
			fprintf(events_combined, "%.4f,%.4f,%d,", local.e_di_selected_sorted[0].pt(), miniAODhelper.GetElectronRelIso(local.e_di_selected_sorted[0], coneSize::R03, corrType::rhoEA,effAreaType::spring15), local.e_di_selected_sorted[0].pdgId());
			fprintf(events_combined, "%.4f,%.4f,%d,", local.mu_di_selected_sorted[0].pt(), miniAODhelper.GetMuonRelIso(local.mu_di_selected_sorted[0], coneSize::R04, corrType::deltaBeta), local.mu_di_selected_sorted[0].pdgId());
		}	
	}
	fprintf(events_combined, "%.4f,%.4f,%.4f,%.4f,", local.jets_di_selected_sorted[0].pt(), local.jets_di_selected_sorted[1].pt(), miniAODhelper.GetJetCSV(local.jets_di_selected_sorted[0],"pfCombinedInclusiveSecondaryVertexV2BJetTags"), miniAODhelper.GetJetCSV(local.jets_di_selected_sorted[1],"pfCombinedInclusiveSecondaryVertexV2BJetTags"));
	fprintf(events_combined, "%.4f,%.4f,%.4f,", local.jet1SF_di, local.jet1SF_up_di, local.jet1SF_down_di);
	fprintf(events_combined, "%.4f,%.4f,%.4f,", local.met_pt, local.met_phi, local.mll);
	fprintf(events_combined, "%d,", local.ttHf_cat);
	if(!isdata)
		fprintf(events_combined, "%.4f,",local.PU_weight);
	else
		fprintf(events_combined, "-1,");
	if(!isdata){
		fprintf(events_combined, "%.4f,%.4f,", local.b_weight_di, local.lep_sf_trig_di);
		fprintf(events_combined, "%.4f,%.4f,", local.lep_sf_id_di, local.lep_sf_iso_di);
	}
	else
		fprintf(events_combined, "-1,-1,-1,-1,");
	if(!isdata)
		fprintf(events_combined, "%.4f,%.4f\n", local.q2_weight_up, local.q2_weight_down);
	else
		fprintf(events_combined, "-1,-1\n");	
	if(!isdata)
		fprintf(events_combined, "%.4f,%.4f\n", local.pdf_weight_up, local.pdf_weight_down);
	else
		fprintf(events_combined, "-1,-1\n");
}

/*
void CU_ttH_EDA::printDecayChain(const reco::Candidate &p, int &index,
								 int mother_index, bool details)
{
	int ndaug = p.numberOfDaughters();

	for (int j = 0; j < ndaug; ++j) {
		const reco::Candidate *daug = p.daughter(j);
		index++;
		cout << index << "\t" << daug->pdgId() << "\t" << daug->status() << "\t"
			 << mother_index << "\t" << daug->numberOfDaughters();
		if (details) {
			cout << "\t" << daug->pt() << "\t" << daug->eta() << "\t"
				 << daug->phi() << endl;
		} else
			cout << endl;

		if (daug->status() != 1)
			printDecayChain(*daug, index, index, details);
	}
}
*/
/*
float CU_ttH_EDA::getMHT(CU_ttH_EDA_event_vars &local)
{
	float MHT_x = 0;
	float MHT_y = 0;

	for (auto & mu : local.mu_selected_sorted) {
		MHT_x -= mu.px();
		MHT_y -= mu.py();
	}

	for (auto & ele : local.e_selected_sorted) {
		MHT_x -= ele.px();
		MHT_y -= ele.py();
	}
	for (auto & jet : local.jets_selected_sorted) {
		MHT_x -= jet.px();
		MHT_y -= jet.py();
	}

	//local.met_pt = sqrt(MHT_x * MHT_x + MHT_y * MHT_y);
	//local.met_phi = atan(MHT_y/MHT_x);
	return sqrt(MHT_x * MHT_x + MHT_y * MHT_y);
}
*/
/*
std::vector<pat::Jet> 
CU_ttH_EDA::CheckJetID (const std::vector<pat::Jet>& inputJets) {
    std::vector<pat::Jet> outputJets;
    bool loose = false;
    for (const auto& iJet : inputJets) {
    	 loose = (
		  iJet.neutralHadronEnergyFraction() < 0.99 &&
		  iJet.chargedEmEnergyFraction() < 0.99 &&
		  iJet.neutralEmEnergyFraction() < 0.99 &&
		  (iJet.neutralMultiplicity() + iJet.chargedMultiplicity() )  > 1
		  );
      
    	if( fabs(iJet.eta())<2.4 ){
		 loose = ( loose &&
		 iJet.chargedHadronEnergyFraction() > 0.0 &&
		 iJet.chargedMultiplicity() > 0
	      	);
    	}
    if (loose == true)
    	outputJets.push_back(iJet);
    }
    return outputJets;
}
*/
std::vector<pat::Jet> 
CU_ttH_EDA::CheckJetID (const std::vector<pat::Jet>& inputJets, const std::vector<pat::Jet>& inputJets_old) {
    std::vector<pat::Jet> outputJets;
    bool loose = false;
    int N = static_cast<int>(inputJets.size());
    double scale;
    for(int i=0; i<N; i++) {
    	scale = inputJets[i].pt()/inputJets_old[i].pt();
    	loose = ( inputJets[i].neutralHadronEnergyFraction()*scale < 0.99 &&
		  inputJets[i].chargedEmEnergyFraction()*scale < 0.99 &&
		  inputJets[i].neutralEmEnergyFraction()*scale < 0.99 &&
		  (inputJets[i].neutralMultiplicity() + inputJets[i].chargedMultiplicity() )  > 1
		  );
      
    	if( fabs(inputJets[i].eta())<2.4 ){
		 loose = ( loose &&
		 inputJets[i].chargedHadronEnergyFraction()*scale > 0.0 &&
		 inputJets[i].chargedMultiplicity() > 0
	      	);
    	}
    	if (loose == true)
    	outputJets.push_back(inputJets[i]);
    	}
    return outputJets;
}

void CU_ttH_EDA::SetFactorizedJetCorrector(const sysType::sysType iSysType){

    std::vector<JetCorrectorParameters> corrParams;	
    if (isdata) {
    
    JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters(string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/Spring16_25nsV3_DATA_L3Absolute_AK4PFchs.txt");
    JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters(string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/Spring16_25nsV3_DATA_L2Relative_AK4PFchs.txt");
    JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters(string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/Spring16_25nsV3_DATA_L1FastJet_AK4PFchs.txt");
    JetCorrectorParameters *L2L3JetPar  = new JetCorrectorParameters(string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/Spring16_25nsV3_DATA_L2L3Residual_AK4PFchs.txt");

    corrParams.push_back(*L1JetPar);
    corrParams.push_back(*L2JetPar);
    corrParams.push_back(*L3JetPar);
    corrParams.push_back(*L2L3JetPar);
    _jetCorrector = new FactorizedJetCorrector(corrParams);

    std::string _JESUncFile = string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/Spring16_25nsV3_DATA_Uncertainty_AK4PFchs.txt";	
    _jetCorrectorUnc = new JetCorrectionUncertainty(_JESUncFile);
    	
    }
    
    else {	
    //setting up the JetCorrector
    JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters(string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/Spring16_25nsV3_MC_L3Absolute_AK4PFchs.txt");
    JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters(string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/Spring16_25nsV3_MC_L2Relative_AK4PFchs.txt");
    JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters(string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/Spring16_25nsV3_MC_L1FastJet_AK4PFchs.txt");

    //JetCorrectorParameters *L3JetPar  = new JetCorrectorParameters(string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/Summer15_25nsV6_MC_L3Absolute_AK4PFchs.txt");
    //JetCorrectorParameters *L2JetPar  = new JetCorrectorParameters(string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/Summer15_25nsV6_MC_L2Relative_AK4PFchs.txt");
    //JetCorrectorParameters *L1JetPar  = new JetCorrectorParameters(string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/Summer15_25nsV6_MC_L1FastJet_AK4PFchs.txt");

    corrParams.push_back(*L1JetPar);
    corrParams.push_back(*L2JetPar);
    corrParams.push_back(*L3JetPar);
    _jetCorrector = new FactorizedJetCorrector(corrParams);

    std::string _JESUncFile = string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/Spring16_25nsV3_MC_Uncertainty_AK4PFchs.txt";	
    //std::string _JESUncFile = string(getenv("CMSSW_BASE")) + "/src/Analyzers/ttH_analyzer/data/Summer15_25nsV6_MC_Uncertainty_AK4PFchs.txt";	
    _jetCorrectorUnc = new JetCorrectionUncertainty(_JESUncFile);
    }
    
}

double
CU_ttH_EDA::GetJetSF( pat::Jet jet, const sysType::sysType iSysType, double rho) {
	double scale = 1;
	_jetCorrector->setJetPt(jet.pt());
    	_jetCorrector->setJetEta(jet.eta());
    	_jetCorrector->setJetA(jet.jetArea());
    	_jetCorrector->setRho(rho); //=fixedGridRhoFastjetAll

	 scale = _jetCorrector->getCorrection();
	 jet.scaleEnergy( scale );
	
	if( iSysType == sysType::JESup || iSysType == sysType::JESdown ){
      		_jetCorrectorUnc->setJetPt(jet.pt());
      		_jetCorrectorUnc->setJetEta(jet.eta());
      		double unc = 1;
      		double jes = 1;
      		if( iSysType==sysType::JESup ){
			unc = _jetCorrectorUnc->getUncertainty(true);
			jes = 1 + unc;
		}
      		else if( iSysType==sysType::JESdown ){
			unc = _jetCorrectorUnc->getUncertainty(false);
			jes = 1 - unc;
      		}
      		//scale = scale*jes;
      		jet.scaleEnergy( jes );
      		return jes;
    	}
    	else
    		return scale;
}


//JEC
std::vector<pat::Jet> 
CU_ttH_EDA::GetCorrectedJets_JEC(const std::vector<pat::Jet>& inputJets, double rho, const sysType::sysType iSysType, const float& corrFactor , const float& uncFactor ){
	
  std::vector<pat::Jet> outputJets;

  for( std::vector<pat::Jet>::const_iterator it = inputJets.begin(), ed = inputJets.end(); it != ed; ++it ){
    
    pat::Jet jet = (*it);
    double scale = 1.;

    _jetCorrector->setJetPt(jet.pt());
    _jetCorrector->setJetEta(jet.eta());
    _jetCorrector->setJetA(jet.jetArea());
    _jetCorrector->setRho(rho); //=fixedGridRhoFastjetAll

    scale = _jetCorrector->getCorrection();
    jet.scaleEnergy( scale );

    if( iSysType == sysType::JESup || iSysType == sysType::JESdown ){
      _jetCorrectorUnc->setJetPt(jet.pt());
      _jetCorrectorUnc->setJetEta(jet.eta()); // here you must use the CORRECTED jet pt
      double unc = 1;
      double jes = 1;
      if( iSysType==sysType::JESup ){
	unc = _jetCorrectorUnc->getUncertainty(true);
	jes = 1 + unc;
      }
      else if( iSysType==sysType::JESdown ){
	unc = _jetCorrectorUnc->getUncertainty(false);
	jes = 1 - unc;
      }

      jet.scaleEnergy( jes );
    }
    
    outputJets.push_back(jet);
  }

  return outputJets;
}

// JER
std::vector<pat::Jet> 
CU_ttH_EDA::GetCorrectedJets_JER(const std::vector<pat::Jet>& inputJets, double rho, const sysType::sysType iSysType, const float& corrFactor , const float& uncFactor ){
	
  std::vector<pat::Jet> outputJets;

  for( std::vector<pat::Jet>::const_iterator it = inputJets.begin(), ed = inputJets.end(); it != ed; ++it ){
    
    pat::Jet jet = (*it);

      double jerSF = 1.;
      
      if( jet.genJet() ){
        if( iSysType == sysType::JERup ){
	      jerSF = getJERfactor(uncFactor, fabs(jet.eta()), jet.genJet()->pt(), jet.pt());
        }
        else if( iSysType == sysType::JERdown ){
	      jerSF = getJERfactor(-uncFactor, fabs(jet.eta()), jet.genJet()->pt(), jet.pt());
        }
        else {
  	      jerSF = getJERfactor(0, fabs(jet.eta()), jet.genJet()->pt(), jet.pt());
        }
      }
     
      jet.scaleEnergy( jerSF*corrFactor );
    
    outputJets.push_back(jet);
  }

  return outputJets;
}

// JER factors for 76X
// Use miniAODHelper's function for 74X factors
double CU_ttH_EDA::getJERfactor( const int returnType, const double jetAbsETA, const double genjetPT, const double recojetPT){

  double factor = 1.;
    
  double scale_JER = 1., scale_JERup = 1., scale_JERdown = 1.;
  double extrauncertainty=1.5;
  
  if( jetAbsETA<0.5 ){ 
    scale_JER = 1.095; scale_JERup = 1.095 + 0.018*extrauncertainty; scale_JERdown = 1.095 - 0.018*extrauncertainty;
  }
  else if( jetAbsETA<0.8 ){ 
    scale_JER = 1.120; scale_JERup = 1.120 + 0.028*extrauncertainty; scale_JERdown = 1.120 - 0.028*extrauncertainty;
  }
  else if( jetAbsETA<1.1 ){ 
    scale_JER = 1.097; scale_JERup = 1.097 + 0.017*extrauncertainty; scale_JERdown = 1.097 - 0.017*extrauncertainty;
  }
  else if( jetAbsETA<1.3 ){ 
    scale_JER = 1.103; scale_JERup = 1.103 + 0.033*extrauncertainty; scale_JERdown = 1.103 - 0.033*extrauncertainty;
  }
  else if( jetAbsETA<1.7 ){ 
    scale_JER = 1.118; scale_JERup = 1.118 + 0.014*extrauncertainty; scale_JERdown = 1.118 - 0.014*extrauncertainty;
  }
  else if( jetAbsETA<1.9 ){ 
    scale_JER = 1.100; scale_JERup = 1.100 + 0.033*extrauncertainty; scale_JERdown = 1.100 - 0.033*extrauncertainty;
  }
  else if( jetAbsETA<2.1 ){ 
    scale_JER = 1.162; scale_JERup = 1.162 + 0.044*extrauncertainty; scale_JERdown = 1.162 - 0.044*extrauncertainty;
  }
  else if( jetAbsETA<2.3 ){ 
    scale_JER = 1.160; scale_JERup = 1.160 + 0.048*extrauncertainty; scale_JERdown = 1.160 - 0.048*extrauncertainty;
  }
  else if( jetAbsETA<2.5 ){ 
    scale_JER = 1.161; scale_JERup = 1.161 + 0.060*extrauncertainty; scale_JERdown = 1.161 - 0.060*extrauncertainty;
  }
  else if( jetAbsETA<2.8 ){ 
    scale_JER = 1.209; scale_JERup = 1.209 + 0.059*extrauncertainty; scale_JERdown = 1.209 - 0.059*extrauncertainty;
  }
  else if( jetAbsETA<3.0 ){ 
    scale_JER = 1.564; scale_JERup = 1.564 + 0.321*extrauncertainty; scale_JERdown = 1.564 - 0.321*extrauncertainty;
  }
  else if( jetAbsETA<3.2 ){ 
    scale_JER = 1.384; scale_JERup = 1.384 + 0.033*extrauncertainty; scale_JERdown = 1.384 - 0.033*extrauncertainty;
  }
  else if( jetAbsETA<5.0 ){ 
    scale_JER = 1.216; scale_JERup = 1.216 + 0.050*extrauncertainty; scale_JERdown = 1.216 - 0.050*extrauncertainty;
  }

  double jetPt_JER = recojetPT;
  double jetPt_JERup = recojetPT;
  double jetPt_JERdown = recojetPT;

  double diff_recojet_genjet = recojetPT - genjetPT;

  //if( genjetPT>5 ){
    jetPt_JER = std::max( 0., genjetPT + scale_JER * ( diff_recojet_genjet ) );
    jetPt_JERup = std::max( 0., genjetPT + scale_JERup * ( diff_recojet_genjet ) );
    jetPt_JERdown = std::max( 0., genjetPT + scale_JERdown * ( diff_recojet_genjet ) );
  //}

  if( returnType==1 )       factor = jetPt_JERup/recojetPT;
  else if( returnType==-1 ) factor = jetPt_JERdown/recojetPT;
  else                      factor = jetPt_JER/recojetPT;

  //if( !(genjetPT>5) ) factor = 1.;

  return factor;
}

void CU_ttH_EDA::Check_SL_Event_Selection(CU_ttH_EDA_event_vars &local){
	if ( local.pass_single_e == 1 || local.pass_single_mu == 1 ) {
		if (local.n_prim_V > 0) {
			if (local.n_leptons == 1) {
				
				if (local.n_electrons == 1) {
					if (local.n_veto_electrons == 1 && local.n_veto_muons == 0 && local.pass_single_e == 1) {
						if (local.n_sl_jets >= min_njets && local.n_sl_btags >= min_nbtags) {
							local.event_selection_SL = true;
							local.is_e = true;
							sl_e++;
						}
					}
				}
				
				else if (local.n_muons == 1) {
					if (local.n_veto_muons == 1 && local.n_veto_electrons == 0 && local.pass_single_mu == 1) {
						if (local.n_sl_jets >= min_njets && local.n_sl_btags >= min_nbtags) {
							local.event_selection_SL = true;
							local.is_mu = true;
							sl_mu++;
						}
					}
				}
				
			}
		}
	}
}

void CU_ttH_EDA::Check_DL_Event_Selection(CU_ttH_EDA_event_vars &local){
	if ( local.pass_double_e == 1 || local.pass_double_mu == 1 || local.pass_elemu == 1 ) {
		if (local.n_prim_V > 0) {
			if (local.n_di_leptons == 2) {   
				if (local.n_di_electrons == 2) {   // di e
					if ( (local.e_di_selected_sorted[0].pdgId()*local.e_di_selected_sorted[1].pdgId() < 0)  &&  local.pass_double_e == 1) {
						if ( local.e_di_selected_sorted[0].pt() > min_di_ele1_pT ) {
							if (local.n_di_jets >= min_di_njets && local.n_di_btags >= min_di_nbtags) {
								//if (local.jets_di_selected_sorted[0].pt() > min_jet_pT && local.jets_di_selected_sorted[1].pt() > min_jet_pT && local.jets_di_selected_tag_sorted[0].pt() > min_jet_pT) {
								if (local.jets_di_selected_sorted[0].pt() > min_jet_pT && local.jets_di_selected_sorted[1].pt() > min_jet_pT) {
									E = local.e_di_selected_sorted[0].energy() + local.e_di_selected_sorted[1].energy();
									px = local.e_di_selected_sorted[0].px() + local.e_di_selected_sorted[1].px();
									py = local.e_di_selected_sorted[0].py() + local.e_di_selected_sorted[1].py();
									pz = local.e_di_selected_sorted[0].pz() + local.e_di_selected_sorted[1].pz();
									p = sqrt(px*px + py*py + pz*pz);
									local.mll = sqrt( E*E - p*p );
									if (local.met_pt >= min_di_met)
										local.met_passed = 1;
									if (local.mll > min_di_mll) {	
										if (local.mll < 76 || local.mll > 106)
											local.mll_passed = 1;
									}	
									if (local.met_passed == 1 && local.mll_passed == 1) {
										local.event_selection_DL = true;
										local.is_ee = true;
										dl_ee++;
									}
								}	
							}	
						}
					}
				}	
				else if (local.n_di_muons == 2) {     // di mu
					if ((local.mu_di_selected_sorted[0].pdgId()*local.mu_di_selected_sorted[1].pdgId() < 0)  &&  local.pass_double_mu == 1) {
						if ( local.mu_di_selected_sorted[0].pt() > min_di_mu1_pT ) {
							if (local.n_di_jets >= min_di_njets && local.n_di_btags >= min_di_nbtags) {
								//if (local.jets_di_selected_sorted[0].pt() > min_jet_pT && local.jets_di_selected_sorted[1].pt() > min_jet_pT && local.jets_di_selected_tag_sorted[0].pt() > min_jet_pT) {
								if (local.jets_di_selected_sorted[0].pt() > min_jet_pT && local.jets_di_selected_sorted[1].pt() > min_jet_pT) {
									E = local.mu_di_selected_sorted[0].energy() + local.mu_di_selected_sorted[1].energy();
									px = local.mu_di_selected_sorted[0].px() + local.mu_di_selected_sorted[1].px();
									py = local.mu_di_selected_sorted[0].py() + local.mu_di_selected_sorted[1].py();
									pz = local.mu_di_selected_sorted[0].pz() + local.mu_di_selected_sorted[1].pz();
									p = sqrt(px*px + py*py + pz*pz);
									local.mll = sqrt( E*E - p*p );
									if (local.met_pt >= min_di_met)
										local.met_passed = 1;
									if (local.mll > min_di_mll) {	
										if (local.mll < 76 || local.mll > 106)
											local.mll_passed = 1;
									}
									if (local.met_passed == 1 && local.mll_passed == 1) {
										local.event_selection_DL = true;
										local.is_mumu = true;
										dl_mumu++;
									}
								}
							}	
						}
					}
				}  
				else {     // 1 e 1 mu
					if ((local.mu_di_selected_sorted[0].pdgId()*local.e_di_selected_sorted[0].pdgId() < 0)  &&  local.pass_elemu == 1) {
						if ( (local.mu_di_selected_sorted[0].pt() > min_di_mu1_pT) || (local.e_di_selected_sorted[0].pt() > min_di_ele1_pT) ) {
							if (local.n_di_jets >= min_di_njets && local.n_di_btags >= min_di_nbtags) {
								//if (local.jets_di_selected_sorted[0].pt() > min_jet_pT && local.jets_di_selected_sorted[1].pt() > min_jet_pT && local.jets_di_selected_tag_sorted[0].pt() > min_jet_pT) {
								if (local.jets_di_selected_sorted[0].pt() > min_jet_pT && local.jets_di_selected_sorted[1].pt() > min_jet_pT) {
									E = local.e_di_selected_sorted[0].energy() + local.mu_di_selected_sorted[0].energy();
									px = local.e_di_selected_sorted[0].px() + local.mu_di_selected_sorted[0].px();
									py = local.e_di_selected_sorted[0].py() + local.mu_di_selected_sorted[0].py();
									pz = local.e_di_selected_sorted[0].pz() + local.mu_di_selected_sorted[0].pz();
									p = sqrt(px*px + py*py + pz*pz);
									local.mll = sqrt( E*E - p*p );
									local.met_passed = 1;
									if (local.mll > min_di_mll) {	
										local.mll_passed = 1;
									}
									if (local.met_passed == 1 && local.mll_passed == 1) {
										local.event_selection_DL = true;
										local.is_emu = true;
										dl_emu++;
									}
								}
							}	
						}
					}
				}
			}
		}
	}
}

void CU_ttH_EDA::Fill_addn_quant(CU_ttH_EDA_event_vars &local, double rho, edm_Handles handle) {
	
	//local.lep_sf_id_sl = local.lep_sf_iso_sl = local.lep_sf_id_di = local.lep_sf_iso_di = 0;
	
	if (local.event_selection_SL!=0) {
		
		// Jet SF
		pat::Jet jet = local.jets_sl_selected_JEC[0];
		jet.setP4(jet.correctedJet(0).p4());
		local.jet1SF_sl = GetJetSF(jet,sysType::NA, rho);
		local.jet1SF_up_sl = GetJetSF(jet,sysType::JESup, rho);
		local.jet1SF_down_sl = GetJetSF(jet,sysType::JESdown, rho);
		
		// to get b-weight
		getbweight(local);
		
		// Lepton SF
		if (local.is_e) {
			local.lep_sf_id_sl = leptonSFhelper.GetElectronSF(  local.e_selected[0].pt() , local.e_selected[0].eta() , 0 , "ID" );
			local.lep_sf_iso_sl = leptonSFhelper.GetElectronSF(  local.e_selected[0].pt() , local.e_selected[0].eta() , 0 , "Iso" );
			local.lep_sf_trig_sl = leptonSFhelper.GetElectronSF(  local.e_selected[0].pt() , local.e_selected[0].eta() , 0 , "Trigger" );
		}
		else if (local.is_mu) {
			local.lep_sf_id_sl = leptonSFhelper.GetMuonSF(  local.mu_selected[0].pt() , local.mu_selected[0].eta() , 0 , "ID" );
			local.lep_sf_iso_sl = leptonSFhelper.GetMuonSF(  local.mu_selected[0].pt() , local.mu_selected[0].eta() , 0 , "Iso" );
			local.lep_sf_trig_sl = leptonSFhelper.GetMuonSF(  local.mu_selected[0].pt() , local.mu_selected[0].eta() , 0 , "Trigger" );
		}
		
		// ttHf Category
		local.ttHf_cat = -1;
  		if( handle.genTtbarId.isValid() ) local.ttHf_cat = *handle.genTtbarId%100;
  		
  		//PDF Weight
  		//if(!isdata)
  		//	getPDFweight(local,handle.event_gen_info);
  		
  		// PU Weight
		if(!isdata)
			local.PU_weight = PU_weight(handle.PupInfo);
	
		// Q2 Weight
		if(!isdata) {
			local.q2_weight_up = getQ2weight( handle.EvtHandle, "1005");
			local.q2_weight_down = getQ2weight( handle.EvtHandle, "1009");
		}
	}

	if (local.event_selection_DL!=0) {
	
		// Jet SF
		pat::Jet jet = local.jets_di_selected_JEC[0];
		jet.setP4(jet.correctedJet(0).p4());
		local.jet1SF_di = GetJetSF(jet,sysType::NA, rho);
		local.jet1SF_up_di = GetJetSF(jet,sysType::JESup, rho);
		local.jet1SF_down_di = GetJetSF(jet,sysType::JESdown, rho);

		// to get b-weight
		getbweight(local);
		
		
		// Lepton SF
		if (local.is_ee) {
			local.lep_sf_id_di = leptonSFhelper.GetElectronSF(  local.e_di_selected[0].pt() , local.e_di_selected[0].eta() , 0 , "ID" )*leptonSFhelper.GetElectronSF(  local.e_di_selected[1].pt() , local.e_di_selected[1].eta() , 0 , "ID" );
			local.lep_sf_iso_di = leptonSFhelper.GetElectronSF(  local.e_di_selected[0].pt() , local.e_di_selected[0].eta() , 0 , "Iso" )*leptonSFhelper.GetElectronSF(  local.e_di_selected[1].pt() , local.e_di_selected[1].eta() , 0 , "Iso" );
			local.lep_sf_trig_di = leptonSFhelper.GetElectronSF(  local.e_di_selected[0].pt() , local.e_di_selected[0].eta() , 0 , "Trigger" )*leptonSFhelper.GetElectronSF(  local.e_di_selected[1].pt() , local.e_di_selected[1].eta() , 0 , "Trigger" );
		}
		else if (local.is_mumu) {
			local.lep_sf_id_di = leptonSFhelper.GetMuonSF(  local.mu_di_selected[0].pt() , local.mu_di_selected[0].eta() , 0 , "ID" )*leptonSFhelper.GetMuonSF(  local.mu_di_selected[1].pt() , local.mu_di_selected[1].eta() , 0 , "ID" );
			local.lep_sf_iso_di = leptonSFhelper.GetMuonSF(  local.mu_di_selected[0].pt() , local.mu_di_selected[0].eta() , 0 , "Iso" )*leptonSFhelper.GetMuonSF(  local.mu_di_selected[1].pt() , local.mu_di_selected[1].eta() , 0 , "Iso" );
			local.lep_sf_trig_di = leptonSFhelper.GetMuonSF(  local.mu_di_selected[0].pt() , local.mu_di_selected[0].eta() , 0 , "Trigger" )*leptonSFhelper.GetMuonSF(  local.mu_di_selected[1].pt() , local.mu_di_selected[1].eta() , 0 , "Trigger" );
		}
		else if (local.is_emu) {
			local.lep_sf_id_di = leptonSFhelper.GetElectronSF(  local.e_di_selected[0].pt() , local.e_di_selected[0].eta() , 0 , "ID" )*leptonSFhelper.GetMuonSF(  local.mu_di_selected[0].pt() , local.mu_di_selected[0].eta() , 0 , "ID" );
			local.lep_sf_iso_di = leptonSFhelper.GetElectronSF(  local.e_di_selected[0].pt() , local.e_di_selected[0].eta() , 0 , "Iso" )*leptonSFhelper.GetMuonSF(  local.mu_di_selected[0].pt() , local.mu_di_selected[0].eta() , 0 , "Iso" );
			local.lep_sf_trig_di = leptonSFhelper.GetElectronSF(  local.e_di_selected[0].pt() , local.e_di_selected[0].eta() , 0 , "Trigger" )*leptonSFhelper.GetMuonSF(  local.mu_di_selected[0].pt() , local.mu_di_selected[0].eta() , 0 , "Trigger" );
		}
		
		// ttHf Category
		local.ttHf_cat = -1;
  		if( handle.genTtbarId.isValid() ) local.ttHf_cat = *handle.genTtbarId%100;
  		
  		//PDF Weight
  		//if(!isdata)
  		//	getPDFweight(local,handle.event_gen_info);
  			
  		// PU Weight
		if(!isdata)
			local.PU_weight = PU_weight(handle.PupInfo);
	
		// Q2 Weight
		if(!isdata) {
			local.q2_weight_up = getQ2weight( handle.EvtHandle, "1005");
			local.q2_weight_down = getQ2weight( handle.EvtHandle, "1009");
		}
		
	}
}

double CU_ttH_EDA::PU_weight ( edm::Handle<std::vector< PileupSummaryInfo > > PupInfo  )
{
	double pu_weight = -1;
	double numTruePV = -1;
  	if( (PupInfo.isValid()) ){
    		for( std::vector<PileupSummaryInfo>::const_iterator PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI ) {
      			int BX = PVI->getBunchCrossing();
      			if( BX==0 ){
				numTruePV = PVI->getTrueNumInteractions();
			}
    		}
  	}
  	
  	for(int i=0; i<50; i++){
  		if(numTruePV < PU_x[i]) {
  			pu_weight = PU_y[i];
  			break;
  		}	
  	}
  	
	return pu_weight;
}

double getQ2weight( edm::Handle<LHEEventProduct> EvtHandle, string ud) {
	double theWeight;
	theWeight = handle.event_gen_info->weight();
	unsigned int i;
	for (i=0; i<handle.EvtHandle->weights().size(); i++) {
   		if ( !(ud.compare(handle.EvtHandle->weights()[i].id))) 
   			theWeight *= handle.EvtHandle->weights()[i].wgt/handle.EvtHandle->originalXWGTUP(); 
	}
	return theWeight;
}

void CU_ttH_EDA::fillCSVHistos(TFile *fileHF, TFile *fileLF)
{
    for( int iSys=0; iSys<9; iSys++ ){
    for( int iPt=0; iPt<5; iPt++ ) h_csv_wgt_hf[iSys][iPt] = NULL;
    for( int iPt=0; iPt<3; iPt++ ){
      for( int iEta=0; iEta<3; iEta++ )h_csv_wgt_lf[iSys][iPt][iEta] = NULL;
    }
  }
  for( int iSys=0; iSys<5; iSys++ ){
    for( int iPt=0; iPt<5; iPt++ ) c_csv_wgt_hf[iSys][iPt] = NULL;
  }

  // CSV reweighting /// only care about the nominal ones
  for( int iSys=0; iSys<9; iSys++ ){
    TString syst_csv_suffix_hf = "final";
    TString syst_csv_suffix_c = "final";
    TString syst_csv_suffix_lf = "final";
    
    switch( iSys ){
    case 0:
      // this is the nominal case
      break;
    case 1:
      // JESUp
      syst_csv_suffix_hf = "final_JESUp"; syst_csv_suffix_lf = "final_JESUp";
      syst_csv_suffix_c  = "final_cErr1Up";
      break;
    case 2:
      // JESDown
      syst_csv_suffix_hf = "final_JESDown"; syst_csv_suffix_lf = "final_JESDown";
      syst_csv_suffix_c  = "final_cErr1Down";
      break;
    case 3:
      // purity up
      syst_csv_suffix_hf = "final_LFUp"; syst_csv_suffix_lf = "final_HFUp";
      syst_csv_suffix_c  = "final_cErr2Up";
      break;
    case 4:
      // purity down
      syst_csv_suffix_hf = "final_LFDown"; syst_csv_suffix_lf = "final_HFDown";
      syst_csv_suffix_c  = "final_cErr2Down";
      break;
    case 5:
      // stats1 up
      syst_csv_suffix_hf = "final_Stats1Up"; syst_csv_suffix_lf = "final_Stats1Up";
      break;
    case 6:
      // stats1 down
      syst_csv_suffix_hf = "final_Stats1Down"; syst_csv_suffix_lf = "final_Stats1Down";
      break;
    case 7:
      // stats2 up
      syst_csv_suffix_hf = "final_Stats2Up"; syst_csv_suffix_lf = "final_Stats2Up";
      break;
    case 8:
      // stats2 down
      syst_csv_suffix_hf = "final_Stats2Down"; syst_csv_suffix_lf = "final_Stats2Down";
      break;
    }

    for( int iPt=0; iPt<5; iPt++ ) h_csv_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_hf.Data()) );

    if( iSys<5 ){
      for( int iPt=0; iPt<5; iPt++ ) c_csv_wgt_hf[iSys][iPt] = (TH1D*)fileHF->Get( Form("c_csv_ratio_Pt%i_Eta0_%s",iPt,syst_csv_suffix_c.Data()) );
    }
    
    for( int iPt=0; iPt<4; iPt++ ){
      for( int iEta=0; iEta<3; iEta++ )h_csv_wgt_lf[iSys][iPt][iEta] = (TH1D*)fileLF->Get( Form("csv_ratio_Pt%i_Eta%i_%s",iPt,iEta,syst_csv_suffix_lf.Data()) );
    }
  }

  return;
}

double CU_ttH_EDA::getCSVWeight(std::vector<double> jetPts, std::vector<double> jetEtas, std::vector<double> jetCSVs,
                       std::vector<int> jetFlavors, int iSys, double &csvWgtHF, double &csvWgtLF, double &csvWgtCF)
{
    int iSysHF = 0;
  switch(iSys){
  case 7:  iSysHF=1; break; //JESUp
  case 8:  iSysHF=2; break; //JESDown
  case 9:  iSysHF=3; break; //LFUp
  case 10: iSysHF=4; break; //LFDown
  case 13: iSysHF=5; break; //Stats1Up
  case 14: iSysHF=6; break; //Stats1Down
  case 15: iSysHF=7; break; //Stats2Up
  case 16: iSysHF=8; break; //Stats2Down
  default : iSysHF = 0; break; //NoSys
  }

  int iSysC = 0;
  switch(iSys){
  case 21: iSysC=1; break;
  case 22: iSysC=2; break;
  case 23: iSysC=3; break;
  case 24: iSysC=4; break;
  default : iSysC = 0; break;
  }

  int iSysLF = 0;
  switch(iSys){
  case 7:  iSysLF=1; break; //JESUp
  case 8:  iSysLF=2; break; //JESDown
  case 11: iSysLF=3; break; //HFUp
  case 12: iSysLF=4; break; //HFDown
  case 17: iSysLF=5; break; //Stats1Up
  case 18: iSysLF=6; break; //Stats1Down
  case 19: iSysLF=7; break; //Stats2Up
  case 20: iSysLF=8; break; //Stats2Down
  default : iSysLF = 0; break; //NoSys
  }

  double csvWgthf = 1.;
  double csvWgtC  = 1.;
  double csvWgtlf = 1.;

  for( int iJet=0; iJet<int(jetPts.size()); iJet++ ){

    double csv = jetCSVs[iJet];
    double jetPt = jetPts[iJet];
    double jetAbsEta = fabs(jetEtas[iJet]);
    int flavor = jetFlavors[iJet];

    int iPt = -1; int iEta = -1;
    if (jetPt >=19.99 && jetPt<30) iPt = 0;
    else if (jetPt >=30 && jetPt<40) iPt = 1;
    else if (jetPt >=40 && jetPt<60) iPt = 2;
    else if (jetPt >=60 && jetPt<100) iPt = 3;
    else if (jetPt >=100) iPt = 4;

    if (jetAbsEta >=0 &&  jetAbsEta<0.8 ) iEta = 0;
    else if ( jetAbsEta>=0.8 && jetAbsEta<1.6 )  iEta = 1;
    else if ( jetAbsEta>=1.6 && jetAbsEta<2.41 ) iEta = 2;

    if (iPt < 0 || iEta < 0) std::cout << "Error, couldn't find Pt, Eta bins for this b-flavor jet, jetPt = " << jetPt << ", jetAbsEta = " << jetAbsEta << std::endl;

    if (abs(flavor) == 5 ){
      int useCSVBin = (csv>=0.) ? h_csv_wgt_hf[iSysHF][iPt]->FindBin(csv) : 1;
      double iCSVWgtHF = h_csv_wgt_hf[iSysHF][iPt]->GetBinContent(useCSVBin);
      if( iCSVWgtHF!=0 ) csvWgthf *= iCSVWgtHF;
    }
    else if( abs(flavor) == 4 ){
      int useCSVBin = (csv>=0.) ? c_csv_wgt_hf[iSysC][iPt]->FindBin(csv) : 1;
      double iCSVWgtC = c_csv_wgt_hf[iSysC][iPt]->GetBinContent(useCSVBin);
      if( iCSVWgtC!=0 ) csvWgtC *= iCSVWgtC;
    }
    else {
      if (iPt >=3) iPt=3;       /// [30-40], [40-60] and [60-10000] only 3 Pt bins for lf
      int useCSVBin = (csv>=0.) ? h_csv_wgt_lf[iSysLF][iPt][iEta]->FindBin(csv) : 1;
      double iCSVWgtLF = h_csv_wgt_lf[iSysLF][iPt][iEta]->GetBinContent(useCSVBin);
      if( iCSVWgtLF!=0 ) csvWgtlf *= iCSVWgtLF;
    }
  }

  double csvWgtTotal = csvWgthf * csvWgtC * csvWgtlf;

  csvWgtHF = csvWgthf;
  csvWgtLF = csvWgtlf;
  csvWgtCF = csvWgtC;

  return csvWgtTotal;
}

void CU_ttH_EDA::getbweight (CU_ttH_EDA_event_vars &local) {
	double csvWgtHF, csvWgtLF, csvWgtCF;
	csvWgtHF = csvWgtLF = csvWgtCF = 0;
	
	if (local.event_selection_SL == 1) {
  		for( std::vector<pat::Jet>::const_iterator iJet = local.jets_sl_selected.begin(); iJet != local.jets_sl_selected.end(); iJet++ ){ 
			 local.vec_jet_pt.push_back(iJet->pt());
			 local.vec_jet_eta.push_back(iJet->eta());
    	 	  	 local.vec_jet_csv.push_back(miniAODhelper.GetJetCSV(*iJet,"pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    			 local.vec_jet_hadronFlavour.push_back(iJet->hadronFlavour());
	 	} 
	
		local.b_weight_sl = getCSVWeight(local.vec_jet_pt, local.vec_jet_eta, local.vec_jet_csv, local.vec_jet_hadronFlavour, local.iSys, csvWgtHF, csvWgtLF, csvWgtCF);
		local.b_weight_di = 0;
	}
	
	else if (local.event_selection_DL ==1) {
		for( std::vector<pat::Jet>::const_iterator iJet = local.jets_di_selected.begin(); iJet != local.jets_di_selected.end(); iJet++ ){ 
			 local.vec_jet_pt.push_back(iJet->pt());
			 local.vec_jet_eta.push_back(iJet->eta());
    	 		 local.vec_jet_csv.push_back(miniAODhelper.GetJetCSV(*iJet,"pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    			 local.vec_jet_hadronFlavour.push_back(iJet->hadronFlavour());
		 } 
	
		local.b_weight_di = getCSVWeight(local.vec_jet_pt, local.vec_jet_eta, local.vec_jet_csv, local.vec_jet_hadronFlavour, local.iSys, csvWgtHF, csvWgtLF, csvWgtCF);
		local.b_weight_sl = 0;
	}
}

/*
void CU_ttH_EDA::getPDFweight(CU_ttH_EDA_event_vars &local, edm::Handle<GenEventInfoProduct> genInfos) {

LHAPDF::PDFSet CT14nlo_PDFSet("CT14nlo");
_systPDFs = CT14nlo_PDFSet.mkPDFs();

auto pdfInfos = genInfos->pdf();
double pdfNominal = pdfInfos->xPDF.first * pdfInfos->xPDF.second;

std::vector<double> pdfs;
for (size_t j = 0; j < CT14nlo_PDFSet.size(); ++j) {
    double xpdf1 = CT14nlo_PDFSet[j]->xfxQ(pdfInfos->id.first, pdfInfos->x.first, pdfInfos->scalePDF);
    double xpdf2 = CT14nlo_PDFSet[j]->xfxQ(pdfInfos->id.second, pdfInfos->x.second, pdfInfos->scalePDF);
    pdfs.push_back(xpdf1 * xpdf2);
}


//const LHAPDF::PDFUncertainty pdfUnc = _systPDFSets[i].uncertainty(pdfs, 68.);
const LHAPDF::PDFUncertainty pdfUnc = _systPDFs[i].uncertainty(pdfs, 68.);

double weight_up = 1.0;
double weight_down = 1.0;
if (std::isfinite(1./pdfNominal)) {
  weight_up = (pdfUnc.central + pdfUnc.errplus) / pdfNominal;
  weight_down = (pdfUnc.central - pdfUnc.errminus) / pdfNominal;
}
local.pdf_weight_up = weight_up;
local.pdf_weight_down = weight_down;
}

*/



#endif
