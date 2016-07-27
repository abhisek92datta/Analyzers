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

cd CMSSW_8_0_12/src/Analyzers/ttH_analyzer/

scram b -j 32

RUN :

For MC :

in yaml config file -- using_real_data : false and HLT_config_tag: 'HLT2'

To Run Locally:

1. put desired MC filename in : test/CU_ttH_mc_EDA_cfg.py
2. cmsRun test/CU_ttH_mc_EDA_cfg.py 

To Run on CRAB :

1. remove filename from test/CU_ttH_mc_EDA_cfg.py
2. put desired MC Dataset in : crabConfig_analyzer.py
3. in crabConfig_analyzer.py : config.Data.splitting = 'FileBased'
                             : config.Data.unitsPerJob = 10  (put desired number)
                             : config.Data.outLFNDirBase = '/store/user/abdatta/' (your own storage area)
                             : config.Site.storageSite = 'T3_US_FNALLPC' (your own storage site)
4. crab submit -c test/CU_ttH_mc_EDA_cfg.py

For DATA :

in yaml config file -- using_real_data : true and HLT_config_tag: 'HLT'

To Run Locally:

1. put desired DATA filename in : test/CU_ttH_data_EDA_cfg.py
2. cmsRun test/CU_ttH_data_EDA_cfg.py 

To Run on CRAB :

1. remove filename from test/CU_ttH_data_EDA_cfg.py
2. put desired DATA filename in : crabConfig_analyzer.py
3. put LUMI file name in : crabConfig_analyzer.py
4. in crabConfig_analyzer.py : config.Data.splitting = 'LumiBased'
                             : config.Data.unitsPerJob = 10  (put desired number)
                             : config.Data.outLFNDirBase = '/store/user/abdatta/' (your own storage area)
                             : config.Site.storageSite = 'T3_US_FNALLPC' (your own storage site)
5. crab submit -c test/CU_ttH_data_EDA_cfg.py

To check CRAB status :
crab status -d \<crab_output_directory_name\>

To resubmit :
crab resubmit -d \<crab_output_directory_name\>



