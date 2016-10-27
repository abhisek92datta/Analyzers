voms-proxy-init -voms cms
cd CMSSW_8_0_20/src/
git cms-merge-topic cms-met:METRecipe_8020
git cms-addpkg PhysicsTools/JetMCAlgos
mv ../../Analyzers .
cp Analyzers/Extras/PhysicsTools_JetMCAlgos_python/GenTtbarCategorizer_cfi.py PhysicsTools/JetMCAlgos/python/
cp Analyzers/Extras/PhysicsTools_JetMCAlgos_plugins/GenTtbarCategorizer.cc PhysicsTools/JetMCAlgos/plugins/
git clone https://github.com/abhisek92datta/MiniAOD.git
cd MiniAOD
git checkout CMSSW_8_0_8_ICHEP_Leptons
cd ..
scram b -j 32
scram b -j 32
scram b -j 32
cd Analyzers/ttH_analyzer/macros/Pile_Up_Calc/
chmod 755 Pileup_calculation_script.sh
./Pileup_calculation_script.sh
