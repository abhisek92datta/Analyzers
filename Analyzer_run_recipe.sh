cd CMSSW_8_0_12/src/
git cms-addpkg PhysicsTools/JetMCAlgos
cd PhysicsTools/
git checkout CMSSW_8_0_X
cd ../
mv ../../Analyzers .
mv Analyzers/GenTtbarCategorizer_cfi.py PhysicsTools/JetMCAlgos/python/
mv Analyzers/GenTtbarCategorizer.cc PhysicsTools/JetMCAlgos/plugins/
