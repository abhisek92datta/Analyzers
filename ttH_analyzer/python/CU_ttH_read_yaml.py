import yaml

# path of yaml file with respect to the directory from where you are running the job eg. in this case from Analyzers/
with open("ttH_analyzer/CU_ttH_cfg_file.yml", 'r') as ymlfile:
    cfg = yaml.load(ymlfile)

analysis_type = cfg["Generic"]["analysis_type"]
verbosity = cfg["Generic"]["verbosity"]
print_HLT_event_path = cfg["Generic"]["print_HLT_event_path"]
HLT_config_tag = cfg["Generic"]["HLT_config_tag"]
filter_config_tag = cfg["Generic"]["filter_config_tag"]


collect_trigger_stats = cfg["Triggers"]["collect_trigger_stats"]
HLT_electron_triggers = cfg["Triggers"]["HLT_electron_triggers"]
HLT_muon_triggers = cfg["Triggers"]["HLT_muon_triggers"]
HLT_electron_electron_triggers = cfg["Triggers"]["HLT_electron_electron_triggers"]
HLT_electron_muon_triggers = cfg["Triggers"]["HLT_electron_muon_triggers"]
HLT_muon_muon_triggers = cfg["Triggers"]["HLT_muon_muon_triggers"]


#min_tight_lepton_pT = cfg["Cuts"]["min_tight_lepton_pT"]
min_ele_pT = cfg["Cuts"]["min_ele_pT"]
min_mu_pT = cfg["Cuts"]["min_mu_pT"]
min_veto_ele_pT = cfg["Cuts"]["min_veto_ele_pT"]
min_veto_mu_pT = cfg["Cuts"]["min_veto_mu_pT"]
min_jet_pT = cfg["Cuts"]["min_jet_pT"]
min_bjet_pT = cfg["Cuts"]["min_bjet_pT"]
max_ele_eta = cfg["Cuts"]["max_ele_eta"]
max_mu_eta = cfg["Cuts"]["max_mu_eta"]
max_veto_ele_eta = cfg["Cuts"]["max_veto_ele_eta"]
max_veto_mu_eta = cfg["Cuts"]["max_veto_mu_eta"]
max_jet_eta = cfg["Cuts"]["max_jet_eta"]
max_bjet_eta = cfg["Cuts"]["max_bjet_eta"]
min_njets = cfg["Cuts"]["min_njets"]
min_nbtags = cfg["Cuts"]["min_nbtags"]

#jet_corrector = cfg["Jets"]["jet_corrector"]

using_real_data = cfg["MiniAODhelper"]["using_real_data"]

b_tag_strength = cfg["B_tag"]["b_tag_strength"]






