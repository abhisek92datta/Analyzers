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
	fprintf(events_combined, "%d, %d, %d, ", local.run_nr, local.lumisection_nr, local.event_nr);	
	fprintf(events_combined, "%d, %d, %d, %d, %d, ", local.is_e, local.is_mu, local.is_ee, local.is_emu, local.is_mumu);
	fprintf(events_combined, "%d, %d, ", local.n_sl_jets, local.n_sl_btags );
	if (local.n_electrons == 1) {
		fprintf(events_combined, "%.4f, %.4f, %d, ", local.e_selected[0].pt(), miniAODhelper.GetElectronRelIso(local.e_selected[0], coneSize::R03, corrType::rhoEA,effAreaType::spring15), local.e_selected[0].pdgId());
	}
	else if (local.n_muons == 1) {
		fprintf(events_combined, "%.4f, %.4f, %d, ", local.mu_selected[0].pt(), miniAODhelper.GetMuonRelIso(local.mu_selected[0], coneSize::R04, corrType::deltaBeta), local.mu_selected[0].pdgId());
	}
	fprintf(events_combined, "0, 0, 0, ");
	fprintf(events_combined, "%.4f, %.4f, %.4f, %.4f, ", local.jets_sl_selected_sorted[0].pt(), local.jets_sl_selected_sorted[1].pt(), miniAODhelper.GetJetCSV(local.jets_sl_selected_sorted[0],"pfCombinedInclusiveSecondaryVertexV2BJetTags"), miniAODhelper.GetJetCSV(local.jets_sl_selected_sorted[1],"pfCombinedInclusiveSecondaryVertexV2BJetTags"));
	fprintf(events_combined, "%.4f, %.4f, %.4f, ", local.jet1SF_sl, local.jet1SF_up_sl, local.jet1SF_down_sl);
	fprintf(events_combined, "%.4f, %.4f, %.4f, ", local.met_pt, local.met_phi, local.mll);
	fprintf(events_combined, "-1, -1, ");
	fprintf(events_combined, "%.4f, ", local.b_weight_sl);
	fprintf(events_combined, "-1, -1, -1, -1, -1, -1, \n");
}

void CU_ttH_EDA::Check_Fill_Print_di_lepton(CU_ttH_EDA_event_vars &local)
{
	fprintf(events_combined, "%d, %d, %d, ", local.run_nr, local.lumisection_nr, local.event_nr);
	fprintf(events_combined, "%d, %d, %d, %d, %d, ", local.is_e, local.is_mu, local.is_ee, local.is_emu, local.is_mumu);
	fprintf(events_combined, "%d, %d, ", local.n_di_jets, local.n_di_btags );
	if (local.n_di_electrons == 2) {
		fprintf(events_combined, "%.4f, %.4f, %d, ", local.e_di_selected_sorted[0].pt(), miniAODhelper.GetElectronRelIso(local.e_di_selected_sorted[0], coneSize::R03, corrType::rhoEA,effAreaType::spring15), local.e_di_selected_sorted[0].pdgId());
		fprintf(events_combined, "%.4f, %.4f, %d, ", local.e_di_selected_sorted[1].pt(), miniAODhelper.GetElectronRelIso(local.e_di_selected_sorted[1], coneSize::R03, corrType::rhoEA,effAreaType::spring15), local.e_di_selected_sorted[1].pdgId());
	}
	else if (local.n_di_muons == 2) {
		fprintf(events_combined, "%.4f, %.4f, %d, ", local.mu_di_selected_sorted[0].pt(), miniAODhelper.GetMuonRelIso(local.mu_di_selected_sorted[0], coneSize::R04, corrType::deltaBeta), local.mu_di_selected_sorted[0].pdgId());
		fprintf(events_combined, "%.4f, %.4f, %d, ", local.mu_di_selected_sorted[1].pt(), miniAODhelper.GetMuonRelIso(local.mu_di_selected_sorted[1], coneSize::R04, corrType::deltaBeta), local.mu_di_selected_sorted[1].pdgId());
	}
	else {
		if (local.mu_di_selected_sorted[0].pt() > local.e_di_selected_sorted[0].pt()) {
			fprintf(events_combined, "%.4f, %.4f, %d, ", local.mu_di_selected_sorted[0].pt(), miniAODhelper.GetMuonRelIso(local.mu_di_selected_sorted[0], coneSize::R04, corrType::deltaBeta), local.mu_di_selected_sorted[0].pdgId());
			fprintf(events_combined, "%.4f, %.4f, %d, ", local.e_di_selected_sorted[0].pt(), miniAODhelper.GetElectronRelIso(local.e_di_selected_sorted[0], coneSize::R03, corrType::rhoEA,effAreaType::spring15), local.e_di_selected_sorted[0].pdgId());
		}
		else {
			fprintf(events_combined, "%.4f, %.4f, %d, ", local.e_di_selected_sorted[0].pt(), miniAODhelper.GetElectronRelIso(local.e_di_selected_sorted[0], coneSize::R03, corrType::rhoEA,effAreaType::spring15), local.e_di_selected_sorted[0].pdgId());
			fprintf(events_combined, "%.4f, %.4f, %d, ", local.mu_di_selected_sorted[0].pt(), miniAODhelper.GetMuonRelIso(local.mu_di_selected_sorted[0], coneSize::R04, corrType::deltaBeta), local.mu_di_selected_sorted[0].pdgId());
		}	
	}
	fprintf(events_combined, "%.4f, %.4f, %.4f, %.4f, ", local.jets_di_selected_sorted[0].pt(), local.jets_di_selected_sorted[1].pt(), miniAODhelper.GetJetCSV(local.jets_di_selected_sorted[0],"pfCombinedInclusiveSecondaryVertexV2BJetTags"), miniAODhelper.GetJetCSV(local.jets_di_selected_sorted[1],"pfCombinedInclusiveSecondaryVertexV2BJetTags"));
	fprintf(events_combined, "%.4f, %.4f, %.4f, ", local.jet1SF_di, local.jet1SF_up_di, local.jet1SF_down_di);
	fprintf(events_combined, "%.4f, %.4f, %.4f, ", local.met_pt, local.met_phi, local.mll);
	fprintf(events_combined, "-1, -1, ");
	fprintf(events_combined, "%.4f, ", local.b_weight_di);
	fprintf(events_combined, "-1, -1, -1, -1, -1, -1, \n");
}

/*
void CU_ttH_EDA::Check_Fill_Print_muj(CU_ttH_EDA_event_vars &local)
{
	int fill_itr = 0;

	h_tth_syncex1_mu->Fill(0.5 + fill_itr++); // fills 0.5 first
	if (local.pass_single_mu) {
		h_tth_syncex1_mu->Fill(0.5 + fill_itr++);
		Print_event_in_file1(events_mu_cut1, local.mu_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_muons == 1) {
		h_tth_syncex1_mu->Fill(0.5 + fill_itr++);
		Print_event_in_file1(events_mu_cut2, local.mu_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_electrons == 0) {
		h_tth_syncex1_mu->Fill(0.5 + fill_itr++);
		Print_event_in_file1(events_mu_cut3, local.mu_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_jets >= 4) {
		h_tth_syncex1_mu->Fill(0.5 + fill_itr++);
		Print_event_in_file1(events_mu_cut4, local.mu_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_btags >= 2) {
		h_tth_syncex1_mu->Fill(0.5 + fill_itr++);
		Print_event_in_file1(events_mu_cut5, local.mu_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_ttags >= 1) {
		h_tth_syncex1_mu->Fill(0.5 + fill_itr++);
		Print_event_in_file1(events_mu_cut6, local.mu_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_Htags >= 1) {
		h_tth_syncex1_mu->Fill(0.5 + fill_itr);
		Print_event_in_file1(events_mu_cut7, local.mu_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;
}

void CU_ttH_EDA::Check_Fill_Print_ej(CU_ttH_EDA_event_vars &local)
{
	int fill_itr = 0;

	h_tth_syncex1_ele->Fill(0.5 + fill_itr++); // fills 0.5 first

	if (local.pass_single_e) {
		h_tth_syncex1_ele->Fill(0.5 + fill_itr++);
		Print_event_in_file1(events_e_cut1, local.e_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_electrons == 1) {
		h_tth_syncex1_ele->Fill(0.5 + fill_itr++);
		Print_event_in_file1(events_e_cut2, local.e_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_muons == 0) {
		h_tth_syncex1_ele->Fill(0.5 + fill_itr++);
		Print_event_in_file1(events_e_cut3, local.e_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_jets >= 4) {
		h_tth_syncex1_ele->Fill(0.5 + fill_itr++);
		Print_event_in_file1(events_e_cut4, local.e_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_btags >= 2) {
		h_tth_syncex1_ele->Fill(0.5 + fill_itr++);
		Print_event_in_file1(events_e_cut5, local.e_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_ttags >= 1) {
		h_tth_syncex1_ele->Fill(0.5 + fill_itr++);
		Print_event_in_file1(events_e_cut6, local.e_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_Htags >= 1) {
		h_tth_syncex1_ele->Fill(0.5 + fill_itr);
		Print_event_in_file1(events_e_cut7, local.e_selected_sorted,
							 local.jets_selected_sorted, local);
	} else
		return;
}

// template<class lepton>
void CU_ttH_EDA::Check_Fill_Print_dimuj(CU_ttH_EDA_event_vars &local)
{
	std::vector<pat::Muon> muon1, muon2;
	if (local.n_muons >= 2) {
		muon1.push_back(local.mu_selected_sorted[0]);
		muon2.push_back(local.mu_selected_sorted[1]);

		TLorentzVector mu1, mu2;
		mu1.SetPtEtaPhiM(local.mu_selected_sorted[0].pt(),
						 local.mu_selected_sorted[0].eta(),
						 local.mu_selected_sorted[0].phi(),
						 local.mu_selected_sorted[0].mass());
		mu2.SetPtEtaPhiM(local.mu_selected_sorted[1].pt(),
						 local.mu_selected_sorted[1].eta(),
						 local.mu_selected_sorted[1].phi(),
						 local.mu_selected_sorted[1].mass());
		local.dimuon_mass = (mu1 + mu2).M();
	} else if (local.n_muons == 1) {
		muon1.push_back(local.mu_selected_sorted[0]);
	}

	int fill_itr = 0;

	h_tth_syncex1_dimu->Fill(0.5 + fill_itr++); // fills 0.5 first

	if (local.pass_double_mu) {
		h_tth_syncex1_dimu->Fill(0.5 + fill_itr++);
		Print_event_in_file1_dilepton(events_dimu_cut1, muon1, muon2,
									  local.dimuon_mass,
									  local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_muons >= 2) {
		if (local.mu_selected_sorted[0].charge() !=
			local.mu_selected_sorted[1].charge()) {
			h_tth_syncex1_dimu->Fill(0.5 + fill_itr++);
			Print_event_in_file1_dilepton(events_dimu_cut2, muon1, muon2,
										  local.dimuon_mass,
										  local.jets_selected_sorted, local);
		} else {
			std::cout << "Found event with two same-charge muons" << std::endl;
		}
		if (local.n_muons > 2) {
			std::cout << "Found event with more than 2 muons" << std::endl;
		}
	} else
		return;

	if (local.dimuon_mass >= 20) {
		h_tth_syncex1_dimu->Fill(0.5 + fill_itr++);
		Print_event_in_file1_dilepton(events_dimu_cut3, muon1, muon2,
									  local.dimuon_mass,
									  local.jets_selected_sorted, local);
	} else
		return;

	if (local.dimuon_mass <= 76 or local.dimuon_mass >= 106) {
		h_tth_syncex1_dimu->Fill(0.5 + fill_itr++);
		Print_event_in_file1_dilepton(events_dimu_cut4, muon1, muon2,
									  local.dimuon_mass,
									  local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_jets >= 2) {
		h_tth_syncex1_dimu->Fill(0.5 + fill_itr++);
		Print_event_in_file1_dilepton(events_dimu_cut5, muon1, muon2,
									  local.dimuon_mass,
									  local.jets_selected_sorted, local);
	} else
		return;

	if (local.MET_corrected.pt() > 40) {
		h_tth_syncex1_dimu->Fill(0.5 + fill_itr++);
		Print_event_in_file1_dilepton(events_dimu_cut6, muon1, muon2,
									  local.dimuon_mass,
									  local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_btags >= 1) {
		h_tth_syncex1_dimu->Fill(0.5 + fill_itr);
		Print_event_in_file1_dilepton(events_dimu_cut7, muon1, muon2,
									  local.dimuon_mass,
									  local.jets_selected_sorted, local);
	} else
		return;
}

void CU_ttH_EDA::Check_Fill_Print_dielej(CU_ttH_EDA_event_vars &local)
{
	std::vector<pat::Electron> electron1, electron2;
	if (local.n_electrons >= 2) {
		electron1.push_back(local.e_selected_sorted[0]);
		electron2.push_back(local.e_selected_sorted[1]);

		TLorentzVector ele1, ele2;
		ele1.SetPtEtaPhiM(local.e_selected_sorted[0].pt(),
						  local.e_selected_sorted[0].eta(),
						  local.e_selected_sorted[0].phi(),
						  local.e_selected_sorted[0].mass());
		ele2.SetPtEtaPhiM(local.e_selected_sorted[1].pt(),
						  local.e_selected_sorted[1].eta(),
						  local.e_selected_sorted[1].phi(),
						  local.e_selected_sorted[1].mass());
		local.dielectron_mass = (ele1 + ele2).M();
	} else if (local.n_electrons == 1) {
		electron1.push_back(local.e_selected_sorted[0]);
	}

	int fill_itr = 0;

	h_tth_syncex1_diele->Fill(0.5 + fill_itr++); // fills 0.5 first

	if (local.pass_double_e) {
		h_tth_syncex1_diele->Fill(0.5 + fill_itr++);
		Print_event_in_file1_dilepton(events_diele_cut1, electron1, electron2,
									  local.dielectron_mass,
									  local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_electrons >= 2) {
		if (local.e_selected_sorted[0].charge() !=
			local.e_selected_sorted[1].charge()) {
			h_tth_syncex1_diele->Fill(0.5 + fill_itr++);
			Print_event_in_file1_dilepton(events_diele_cut2, electron1,
										  electron2, local.dielectron_mass,
										  local.jets_selected_sorted, local);
		} else {
			std::cout << "Found event with two same-charge electrons"
					  << std::endl;
		}
		if (local.n_electrons > 2) {
			std::cout << "Found event with more than 2 electrons" << std::endl;
		}
	} else
		return;

	if (local.dielectron_mass >= 20) {
		h_tth_syncex1_diele->Fill(0.5 + fill_itr++);
		Print_event_in_file1_dilepton(events_diele_cut3, electron1, electron2,
									  local.dielectron_mass,
									  local.jets_selected_sorted, local);
	} else
		return;

	if (local.dielectron_mass <= 76 or local.dielectron_mass >= 106) {
		h_tth_syncex1_diele->Fill(0.5 + fill_itr++);
		Print_event_in_file1_dilepton(events_diele_cut4, electron1, electron2,
									  local.dielectron_mass,
									  local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_jets >= 2) {
		h_tth_syncex1_diele->Fill(0.5 + fill_itr++);
		Print_event_in_file1_dilepton(events_diele_cut5, electron1, electron2,
									  local.dielectron_mass,
									  local.jets_selected_sorted, local);
	} else
		return;

	if (local.MET_corrected.pt() > 40) {
		h_tth_syncex1_diele->Fill(0.5 + fill_itr++);
		Print_event_in_file1_dilepton(events_diele_cut6, electron1, electron2,
									  local.dielectron_mass,
									  local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_btags >= 1) {
		h_tth_syncex1_diele->Fill(0.5 + fill_itr);
		Print_event_in_file1_dilepton(events_diele_cut7, electron1, electron2,
									  local.dielectron_mass,
									  local.jets_selected_sorted, local);
	} else
		return;
}

void CU_ttH_EDA::Check_Fill_Print_elemuj(CU_ttH_EDA_event_vars &local)
{
	std::vector<pat::Electron> electron1;
	std::vector<pat::Muon> muon1;
	if (local.n_electrons >= 1) {
		electron1.push_back(local.e_selected_sorted[0]);
	}
	if (local.n_muons >= 1) {
		muon1.push_back(local.mu_selected_sorted[0]);
	}
	if ((local.n_electrons >= 1) and (local.n_muons >= 1)) {
		// electron1.push_back( local.e_selected_sorted[0] );
		// muon1.push_back( local.mu_selected_sorted[0] );

		TLorentzVector ele1, mu1;
		ele1.SetPtEtaPhiM(local.e_selected_sorted[0].pt(),
						  local.e_selected_sorted[0].eta(),
						  local.e_selected_sorted[0].phi(),
						  local.e_selected_sorted[0].mass());
		mu1.SetPtEtaPhiM(local.mu_selected_sorted[0].pt(),
						 local.mu_selected_sorted[0].eta(),
						 local.mu_selected_sorted[0].phi(),
						 local.mu_selected_sorted[0].mass());
		local.dilepton_mass = (ele1 + mu1).M();
	}

	int fill_itr = 0;

	h_tth_syncex1_elemu->Fill(0.5 + fill_itr++); // fills 0.5 first

	if (local.pass_elemu) {
		h_tth_syncex1_elemu->Fill(0.5 + fill_itr++);
		Print_event_in_file1_dilepton(events_elemu_cut1, muon1, electron1,
									  local.dilepton_mass,
									  local.jets_selected_sorted, local);
	} else
		return;

	if ((local.n_electrons >= 1) and (local.n_muons >= 1)) {
		if (local.e_selected_sorted[0].charge() !=
			local.mu_selected_sorted[0].charge()) {
			h_tth_syncex1_elemu->Fill(0.5 + fill_itr++);
			Print_event_in_file1_dilepton(events_elemu_cut2, muon1, electron1,
										  local.dilepton_mass,
										  local.jets_selected_sorted, local);
		} else {
			std::cout << "Found event with two same-charge leptons"
					  << std::endl;
		}
		if (local.n_electrons > 2) {
			std::cout << "Found event with more than 2 leptons" << std::endl;
		}
	} else
		return;

	if (local.dilepton_mass >= 20) {
		h_tth_syncex1_elemu->Fill(0.5 + fill_itr++);
		Print_event_in_file1_dilepton(events_elemu_cut3, muon1, electron1,
									  local.dilepton_mass,
									  local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_jets >= 2) {
		h_tth_syncex1_elemu->Fill(0.5 + fill_itr++);
		Print_event_in_file1_dilepton(events_elemu_cut4, muon1, electron1,
									  local.dilepton_mass,
									  local.jets_selected_sorted, local);
	} else
		return;

	if (local.n_btags >= 1) {
		h_tth_syncex1_elemu->Fill(0.5 + fill_itr);
		Print_event_in_file1_dilepton(events_elemu_cut5, muon1, electron1,
									  local.dilepton_mass,
									  local.jets_selected_sorted, local);
	} else
		return;
}


template <class lepton>
int CU_ttH_EDA::Print_event_in_file1(FILE *file, lepton &lpt,
									 std::vector<pat::Jet> &jets,
									 CU_ttH_EDA_event_vars &local)
{
	// print generic event info
	fprintf(file, "%6d %8d %10d   ", local.run_nr, local.lumisection_nr,
			local.event_nr);

	// print lepton info
	if (lpt.size() != 0)
		fprintf(file, "%6.2f %+4.2f %+4.2f   ", lpt[0].pt(), lpt[0].eta(),
				lpt[0].phi());
	else
		fprintf(file, "%6.2f %+4.2f %+4.2f   ", -99., -99., -99.);

	// print jet(s) info
	if (jets.size() >= 4) {
		fprintf(file,
				"%6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   ",
				jets[0].pt(), jets[1].pt(), jets[2].pt(), jets[3].pt(),
				jets[0].bDiscriminator(
					"combinedInclusiveSecondaryVertexV2BJetTags"),
				jets[1].bDiscriminator(
					"combinedInclusiveSecondaryVertexV2BJetTags"),
				jets[2].bDiscriminator(
					"combinedInclusiveSecondaryVertexV2BJetTags"),
				jets[3].bDiscriminator(
					"combinedInclusiveSecondaryVertexV2BJetTags"));
	} else {
		switch (jets.size()) {
		case 3:
			fprintf(file,
					"%6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   ",
					jets[0].pt(), jets[1].pt(), jets[2].pt(), -99.,
					jets[0].bDiscriminator(
						"combinedInclusiveSecondaryVertexV2BJetTags"),
					jets[1].bDiscriminator(
						"combinedInclusiveSecondaryVertexV2BJetTags"),
					jets[2].bDiscriminator(
						"combinedInclusiveSecondaryVertexV2BJetTags"),
					-99.);
			break;

		case 2:
			fprintf(file,
					"%6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   ",
					jets[0].pt(), jets[1].pt(), -99., -99.,
					jets[0].bDiscriminator(
						"combinedInclusiveSecondaryVertexV2BJetTags"),
					jets[1].bDiscriminator(
						"combinedInclusiveSecondaryVertexV2BJetTags"),
					-99., -99.);
			break;

		case 1:
			fprintf(file,
					"%6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   ",
					jets[0].pt(), -99., -99., -99.,
					jets[0].bDiscriminator(
						"combinedInclusiveSecondaryVertexV2BJetTags"),
					-99., -99., -99.);
			break;

		default:
			fprintf(file,
					"%6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   ",
					-99., -99., -99., -99., -99., -99., -99., -99.);
		}
	}

	// print number of tags
	fprintf(file, "%2d  %2d   %2d  %2d\n", local.n_jets, local.n_btags,
			local.n_ttags, local.n_Htags);

	return ferror(file);
}

template <class lep1, class lep2>
// template<class lepton2>
int CU_ttH_EDA::Print_event_in_file1_dilepton(FILE *file, lep1 &lepton1,
											  lep2 &lepton2,
											  double dilepton_mass,
											  std::vector<pat::Jet> &jets,
											  CU_ttH_EDA_event_vars &local)
{
	// std::vector<lepton> leptons;
	// switch (dilepton_type) {
	//        case 'dimuon':
	//        leptons

	// print generic event info
	fprintf(file, "%6d %8d %10d   ", local.run_nr, local.lumisection_nr,
			local.event_nr);

	// print lepton info
	if (lepton1.size() == 1) {
		fprintf(file, "%6.2f %+4.2f %+4.2f   ", lepton1[0].pt(),
				lepton1[0].eta(), lepton1[0].phi());
	} else {
		fprintf(file, "%6.2f %+4.2f %+4.2f   ", -99., -99., -99.);
	}

	if (lepton2.size() == 1) {
		fprintf(file, "%6.2f %+4.2f %+4.2f   ", lepton2[0].pt(),
				lepton2[0].eta(), lepton2[0].phi());
	} else {
		fprintf(file, "%6.2f %+4.2f %+4.2f   ", -99., -99., -99.);
	}

	// print MET and dilepton mass
	fprintf(file, "%6.2f  %6.2f   ", dilepton_mass, local.MET_corrected.pt());

	// print jet(s) info
	if (jets.size() >= 4) {
		fprintf(file,
				"%6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   ",
				jets[0].pt(), jets[1].pt(), jets[2].pt(), jets[3].pt(),
				jets[0].bDiscriminator(
					"combinedInclusiveSecondaryVertexV2BJetTags"),
				jets[1].bDiscriminator(
					"combinedInclusiveSecondaryVertexV2BJetTags"),
				jets[2].bDiscriminator(
					"combinedInclusiveSecondaryVertexV2BJetTags"),
				jets[3].bDiscriminator(
					"combinedInclusiveSecondaryVertexV2BJetTags"));
	} else {
		switch (jets.size()) {
		case 3:
			fprintf(file,
					"%6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   ",
					jets[0].pt(), jets[1].pt(), jets[2].pt(), -99.,
					jets[0].bDiscriminator(
						"combinedInclusiveSecondaryVertexV2BJetTags"),
					jets[1].bDiscriminator(
						"combinedInclusiveSecondaryVertexV2BJetTags"),
					jets[2].bDiscriminator(
						"combinedInclusiveSecondaryVertexV2BJetTags"),
					-99.);
			break;

		case 2:
			fprintf(file,
					"%6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   ",
					jets[0].pt(), jets[1].pt(), -99., -99.,
					jets[0].bDiscriminator(
						"combinedInclusiveSecondaryVertexV2BJetTags"),
					jets[1].bDiscriminator(
						"combinedInclusiveSecondaryVertexV2BJetTags"),
					-99., -99.);
			break;

		case 1:
			fprintf(file,
					"%6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   ",
					jets[0].pt(), -99., -99., -99.,
					jets[0].bDiscriminator(
						"combinedInclusiveSecondaryVertexV2BJetTags"),
					-99., -99., -99.);
			break;

		default:
			fprintf(file,
					"%6.2f %6.2f %6.2f %6.2f   %+7.3f %+7.3f %+7.3f %+7.3f   ",
					-99., -99., -99., -99., -99., -99., -99., -99.);
		}
	}

	// print number of tags
	fprintf(file, "%2d  %2d\n", local.n_jets, local.n_btags);

	return ferror(file);
}

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


void CU_ttH_EDA::SetFactorizedJetCorrector(const sysType::sysType iSysType){

    //setting up the JetCorrector
    std::vector<JetCorrectorParameters> corrParams;
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
      		scale = scale*jes;
      		jet.scaleEnergy( jes );
    	}
    	if( jet.genJet() )
    		scale = scale * miniAODhelper.getJERfactor(0, fabs(jet.eta()), jet.genJet()->pt(), jet.pt());
	return scale;
}



std::vector<pat::Jet> 
CU_ttH_EDA::GetCorrectedJets(const std::vector<pat::Jet>& inputJets, double rho, const sysType::sysType iSysType, const float& corrFactor , const float& uncFactor ){
	
  std::vector<pat::Jet> outputJets;

  for( std::vector<pat::Jet>::const_iterator it = inputJets.begin(), ed = inputJets.end(); it != ed; ++it ){
    
    pat::Jet jet = (*it);
    double scale = 1.;

    /// JES
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

    /// JER
      double jerSF = 1.;
      
      if( jet.genJet() ){
        if( iSysType == sysType::JERup ){
	      jerSF = miniAODhelper.getJERfactor(uncFactor, fabs(jet.eta()), jet.genJet()->pt(), jet.pt());
        }
        else if( iSysType == sysType::JERdown ){
	      jerSF = miniAODhelper.getJERfactor(-uncFactor, fabs(jet.eta()), jet.genJet()->pt(), jet.pt());
        }
        else {
  	      jerSF = miniAODhelper.getJERfactor(0, fabs(jet.eta()), jet.genJet()->pt(), jet.pt());
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

  if( genjetPT>10. ){
    jetPt_JER = std::max( 0., genjetPT + scale_JER * ( diff_recojet_genjet ) );
    jetPt_JERup = std::max( 0., genjetPT + scale_JERup * ( diff_recojet_genjet ) );
    jetPt_JERdown = std::max( 0., genjetPT + scale_JERdown * ( diff_recojet_genjet ) );
  }

  if( returnType==1 )       factor = jetPt_JERup/recojetPT;
  else if( returnType==-1 ) factor = jetPt_JERdown/recojetPT;
  else                      factor = jetPt_JER/recojetPT;

  if( !(genjetPT>10.) ) factor = 1.;

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
						}
					}
				}	
				else if (local.n_muons == 1) {
					if (local.n_veto_muons == 1 && local.n_veto_electrons == 0 && local.pass_single_mu == 1) {
						if (local.n_sl_jets >= min_njets && local.n_sl_btags >= min_nbtags) {
							local.event_selection_SL = true;
							local.is_mu = true;
						}
					}
				}
			}
		}
	}
	//if (local.met_pt > 30)
	//	local.met_passed = 1;
	//local.mll = 0;
}

void CU_ttH_EDA::Check_DL_Event_Selection(CU_ttH_EDA_event_vars &local){
	if ( local.pass_double_e == 1 || local.pass_double_mu == 1 || local.pass_elemu == 1 ) {
		if (local.n_prim_V > 0) {
			if (local.n_di_leptons == 2) {   
				if (local.n_di_electrons == 2) {   // di e
					if ( (local.e_di_selected_sorted[0].pdgId()*local.e_di_selected_sorted[1].pdgId() < 0)  &&  local.pass_double_e == 1) {
						if ( local.e_di_selected_sorted[0].pt() > min_di_ele1_pT ) {
							if (local.n_di_jets >= min_di_njets && local.n_di_btags >= min_di_nbtags) {
								if (local.jets_di_selected_sorted[0].pt() > min_jet_pT && local.jets_di_selected_sorted[1].pt() > min_jet_pT && local.jets_di_selected_tag_sorted[0].pt() > min_jet_pT) {
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
								if (local.jets_di_selected_sorted[0].pt() > min_jet_pT && local.jets_di_selected_sorted[1].pt() > min_jet_pT && local.jets_di_selected_tag_sorted[0].pt() > min_jet_pT) {
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
								if (local.jets_di_selected_sorted[0].pt() > min_jet_pT && local.jets_di_selected_sorted[1].pt() > min_jet_pT && local.jets_di_selected_tag_sorted[0].pt() > min_jet_pT) {
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
	
  	for( std::vector<pat::Jet>::const_iterator iJet = local.jets_sl_selected.begin(); iJet != local.jets_sl_selected.end(); iJet++ ){ 
		 local.vec_jet_pt.push_back(iJet->pt());
		 local.vec_jet_eta.push_back(iJet->eta());
    	 	 local.vec_jet_csv.push_back(miniAODhelper.GetJetCSV(*iJet,"pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    		 local.vec_jet_hadronFlavour.push_back(iJet->hadronFlavour());
	 } 
	
	local.b_weight_sl = getCSVWeight(local.vec_jet_pt, local.vec_jet_eta, local.vec_jet_csv, local.vec_jet_hadronFlavour, local.iSys, csvWgtHF, csvWgtLF, csvWgtCF);

	csvWgtHF = csvWgtLF = csvWgtCF = 0;
	
	for( std::vector<pat::Jet>::const_iterator iJet = local.jets_di_selected.begin(); iJet != local.jets_di_selected.end(); iJet++ ){ 
		 local.vec_jet_pt.push_back(iJet->pt());
		 local.vec_jet_eta.push_back(iJet->eta());
    	 	 local.vec_jet_csv.push_back(miniAODhelper.GetJetCSV(*iJet,"pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    		 local.vec_jet_hadronFlavour.push_back(iJet->hadronFlavour());
	 } 
	
	local.b_weight_di = getCSVWeight(local.vec_jet_pt, local.vec_jet_eta, local.vec_jet_csv, local.vec_jet_hadronFlavour, local.iSys, csvWgtHF, csvWgtLF, csvWgtCF);
}

#endif
