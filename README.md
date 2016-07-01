# Analyzers

This is for "ttH " -  H to bb" analysis

#Installation :

Setup CMSSW environment and get Analyzer repository:

cmsrel CMSSW_8_0_8_patch1
cd CMSSW_8_0_8_patch1/src/
cmsenv

git cms-addpkg PhysicsTools/JetMCAlgos
cd PhysicsTools/
git checkout CMSSW_8_0_8
cd ../

mv GenTtbarCategorizer_cfi.py PhysicsTools/JetMCAlgos/python/
mv GenTtbarCategorizer.cc PhysicsTools/JetMCAlgos/plugins/

git clone https://github.com/abhisek92datta/Analyzers.git

Get dependencies:

git clone https://github.com/cms-ttH/MiniAOD.git
cd MiniAOD
git checkout CMSSW_8_0_8
cd ..

Compile:

scram b -j 32

Run :

cd Analyzers/

cmsRun ttH_analyzer/test/CU_ttH_mc_EDA_cfg.py
or 
cmsRun ttH_analyzer/test/CU_ttH_data_EDA_cfg.py
