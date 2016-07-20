# Analyzers

This is for "ttH " -  H to bb" analysis

#Installation :

cmsrel CMSSW_8_0_12

cd CMSSW_8_0_12/src/

cmsenv

cd ../../

git clone https://github.com/abhisek92datta/Analyzers.git

cp Analyzers/Analyzer_run_recipe.sh .

chmod 755 Analyzer_run_recipe.sh

./Analyzer_run_recipe.sh


Run :

cd CMSSW_8_0_12/src/Analyzers/ttH_analyzer/

For MC :

1. in yaml config file -- using_real_data : false
                          HLT_config_tag: 'HLT2'
2. put desired filename in : test/CU_ttH_mc_EDA_cfg.py
3. cmsRun test/CU_ttH_mc_EDA_cfg.py 

For Data :

1. in yaml config file -- using_real_data : true
                          HLT_config_tag: 'HLT'
2. put desired filename in : test/CU_ttH_data_EDA_cfg.py
3. cmsRun test/CU_ttH_data_EDA_cfg.py 

