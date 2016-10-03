cd CMSSW_8_0_17/src/
git cms-merge-topic -u cms-btv-pog:bTagHIPMitigation-PR_from-CMSSW_8_0_17
git cms-addpkg PhysicsTools/JetMCAlgos
mv ../../Analyzers .
cp Analyzers/Extras/PhysicsTools_JetMCAlgos_python/GenTtbarCategorizer_cfi.py PhysicsTools/JetMCAlgos/python/
cp Analyzers/Extras/PhysicsTools_JetMCAlgos_plugins/GenTtbarCategorizer.cc PhysicsTools/JetMCAlgos/plugins/
#echo /PhysicsTools/PatUtils/ >> .git/info/sparse-checkout
#git cms-merge-topic cms-met:metTool80X
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
