voms-proxy-init -voms cms
cd $CMSSW_VERSION/src
git cms-init
git cms-merge-topic riga:deterministicSeeds
git cms-merge-topic riga:badGlobalMuonTagger_fix
git cms-merge-topic cms-met:METRecipe_80X_part2
git clone https://github.com/cms-met/MetTools.git
git cms-merge-topic riga:deterministicEGMSmearer_v2
cd EgammaAnalysis/ElectronTools/data
git clone https://github.com/ECALELFS/ScalesSmearings.git -b Moriond17_gainSwitch_unc
cd $CMSSW_BASE/src
git cms-merge-topic riga:ttHFGenFilter_tagging
git remote add ahinzmann https://github.com/ahinzmann/cmssw.git
git fetch ahinzmann PUidMiniAODfix80
git cherry-pick ca33756e1747aec27d13971bcfd0874b16724e7f
#git clone https://gitlab.cern.ch/ttH/CommonClassifier.git TTH/CommonClassifier
#source TTH/CommonClassifier/setup/install_mem.sh
cd $CMSSW_BASE/src
mv ../../Analyzers .
git cms-addpkg PhysicsTools/JetMCAlgos
cp Analyzers/Extras/PhysicsTools_JetMCAlgos_python/GenTtbarCategorizer_cfi.py PhysicsTools/JetMCAlgos/python/
cp Analyzers/Extras/PhysicsTools_JetMCAlgos_plugins/GenTtbarCategorizer.cc PhysicsTools/JetMCAlgos/plugins/
git clone https://github.com/abhisek92datta/MiniAOD.git
cd MiniAOD
git checkout CMSSW_8_0_8_ICHEP_Leptons
cd ..
scram b
cd Analyzers/ttH_analyzer/macros/Pile_Up_Calc/
chmod 755 Pileup_calculation_script.sh
