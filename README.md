# Analyzers

This is for "ttH " -  H to bb" analysis

#Installation :

Setup CMSSW environment and get Analyzer repository:

cmsrel CMSSW_8_0_8_patch1

cd CMSSW_8_0_8_patch1/src/

cmsenv

git cms-addpkg PhysicsTools/JetMCAlgos

cd PhysicsTools/

git checkout CMSSW_8_0_8_patchX

cd ../

git clone https://github.com/abhisek92datta/Analyzers.git

mv Analyzers/GenTtbarCategorizer_cfi.py PhysicsTools/JetMCAlgos/python/

mv Analyzers/GenTtbarCategorizer.cc PhysicsTools/JetMCAlgos/plugins/

Get dependencies:

git clone https://github.com/abhisek92datta/MiniAOD.git

cd MiniAOD

git checkout CMSSW_8_0_8

cd ..

Compile:

scram b -j 32

Run :

cd Analyzers/ttH_analyzer/

For MC :

1. using_real_data : false

2. cmsRun test/CU_ttH_mc_EDA_cfg.py 

For Data :

1. using_real_data : true

2. cmsRun test/CU_ttH_data_EDA_cfg.py
