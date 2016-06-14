# Analyzers

This is for "ttH " -  H to bb" analysis

#Installation :

Setup CMSSW environment and get Analyzer repository:

cmsrel CMSSW_7_6_3_patch2

cd CMSSW_7_6_3_patch2/src/

cmsenv

git cms-init

git clone https://github.com/abhisek92datta/Analyzers.git

cd Analyzers

git checkout master

cd ..

Get dependencies:

git clone https://github.com/cms-ttH/MiniAOD.git

Then switch to CMSSW_7_6_3 branch:

cd MiniAOD

git checkout CMSSW_7_6_3

cd ..

For ttH, get LeptonID package shared with ttH-Multilepton group:

git clone https://github.com/cms-ttH/ttH-LeptonID.git

mkdir ttH

mv ttH-LeptonID ttH/LeptonID

Compile:

scram b -j 16

Run :

cmsRun <config_file_name> 
