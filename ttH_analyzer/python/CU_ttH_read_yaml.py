import yaml

# path of yaml file with respect to the directory from where you are running the job eg. in this case from Analyzers/
with open("CU_ttH_cfg_file.yml", 'r') as ymlfile:
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
filter_names = cfg["Triggers"]["filter_names"]

min_ele_pT = cfg["Cuts"]["min_ele_pT"]
min_mu_pT = cfg["Cuts"]["min_mu_pT"]
min_veto_ele_pT = cfg["Cuts"]["min_veto_ele_pT"]
min_veto_mu_pT = cfg["Cuts"]["min_veto_mu_pT"]
min_di_ele1_pT = cfg["Cuts"]["min_di_ele1_pT"]
min_di_ele2_pT = cfg["Cuts"]["min_di_ele2_pT"]
min_di_mu1_pT = cfg["Cuts"]["min_di_mu1_pT"]
min_di_mu2_pT = cfg["Cuts"]["min_di_mu2_pT"]
min_jet_pT = cfg["Cuts"]["min_jet_pT"]
min_jet2_pT = cfg["Cuts"]["min_jet2_pT"]
min_bjet_pT = cfg["Cuts"]["min_bjet_pT"]
max_ele_eta = cfg["Cuts"]["max_ele_eta"]
max_mu_eta = cfg["Cuts"]["max_mu_eta"]
max_veto_ele_eta = cfg["Cuts"]["max_veto_ele_eta"]
max_veto_mu_eta = cfg["Cuts"]["max_veto_mu_eta"]
max_di_ele1_eta = cfg["Cuts"]["max_di_ele1_eta"]
max_di_ele2_eta = cfg["Cuts"]["max_di_ele2_eta"]
max_di_mu1_eta = cfg["Cuts"]["max_di_mu1_eta"]
max_di_mu2_eta = cfg["Cuts"]["max_di_mu2_eta"]
max_jet_eta = cfg["Cuts"]["max_jet_eta"]
max_bjet_eta = cfg["Cuts"]["max_bjet_eta"]
min_njets = cfg["Cuts"]["min_njets"]
min_di_njets = cfg["Cuts"]["min_di_njets"]
min_nbtags = cfg["Cuts"]["min_nbtags"]
min_di_nbtags = cfg["Cuts"]["min_di_nbtags"]
min_di_mll = cfg["Cuts"]["min_di_mll"]
min_di_met = cfg["Cuts"]["min_di_met"]

using_real_data = cfg["MiniAODhelper"]["using_real_data"]
dataset = cfg["MiniAODhelper"]["dataset"]

b_tag_strength = cfg["B_tag"]["b_tag_strength"]






