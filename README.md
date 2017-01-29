# Analyzers

This is for "ttH " -  H to bb" analysis

#Installation :

git clone https://github.com/abhisek92datta/Analyzers.git

cp Analyzers/Analyzer_run_recipe.sh .

chmod 755 Analyzer_run_recipe.sh

( Put the latest LUMI File name in : $CMSSW_VERSION/src/Analyzers/ttH_analyzer/macros/Pile_Up_Calc/Pileup_calculation_script.sh)

./Analyzer_run_recipe.sh

cd $CMSSW_VERSION/src/Analyzers/ttH_analyzer/

RUN :

For MC :

in yaml config file -- using_real_data : false 

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
4. crab submit -c crabConfig_analyzer.py

For DATA :

in yaml config file -- using_real_data : true 

Select the appropriate event selection for the dataset
in CU_ttH_EDA.cc and CU_ttH_EDA_Misc.cc

To Run Locally:

1. put desired DATA filename in : test/CU_ttH_data_EDA_cfg.py
2. put latest LUMI filename in : test/CU_ttH_data_EDA_cfg.py
3. cmsRun test/CU_ttH_data_EDA_cfg.py 

To Run on CRAB :

1. remove filename from test/CU_ttH_data_EDA_cfg.py
2. put desired DATA dataset in : crabConfig_analyzer.py
3. put LUMI file name in : crabConfig_analyzer.py
4. in crabConfig_analyzer.py : config.Data.splitting = 'LumiBased'
 
                             : config.Data.unitsPerJob = 10  (put desired number)

                             : config.Data.outLFNDirBase = '/store/user/abdatta/' (your own storage area)
                             
                             : config.Site.storageSite = 'T3_US_FNALLPC' (your own storage site)
5. crab submit -c crabConfig_analyzer.py

To check CRAB status :
crab status -d \<crab_output_directory_name\>

To resubmit :
crab resubmit -d \<crab_output_directory_name\>


Formatting :

clang-format -style=file -i *.h
clang-format -style=file -i *.cc



