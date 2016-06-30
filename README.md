# Analyzers

This is for "ttH " -  H to bb" analysis

#Installation :

Setup CMSSW environment and get Analyzer repository:

cmsrel CMSSW_8_0_8_patch1

cd CMSSW_8_0_8_patch1/src/

cmsenv

git cms-init

git cms-addpkg PhysicsTools/JetMCAlgos

git clone https://github.com/abhisek92datta/Analyzers.git

cd Analyzers

git checkout master

cd ..

Get dependencies:

git clone https://github.com/cms-ttH/MiniAOD.git

Then switch to CMSSW_8_0_8 branch:

cd MiniAOD

git checkout CMSSW_8_0_8

cd ..

Compile:

scram b -j 32

Run :

cd Analyzers/

cmsRun ttH_analyzer/test/CU_ttH_EDA_cfg.py
