cd CMSSW_8_0_12/src/
git cms-addpkg PhysicsTools/JetMCAlgos
cd PhysicsTools/
cd ../
mv ../../Analyzers .
mv Analyzers/GenTtbarCategorizer_cfi.py PhysicsTools/JetMCAlgos/python/
mv Analyzers/GenTtbarCategorizer.cc PhysicsTools/JetMCAlgos/plugins/
#echo /PhysicsTools/PatUtils/ >> .git/info/sparse-checkout
#git cms-merge-topic cms-met:metTool80X
git clone https://github.com/abhisek92datta/MiniAOD.git
cd MiniAOD
git checkout CMSSW_8_0_8_LeptonSF_updated
cd ..
scram b -j 32
scram b -j 32
cd Analyzers/ttH_analyzer/macros/Pile_Up_Calc/
chmod 755 Pileup_calculation_script.sh
./Pileup_calculation_script.sh
