cd CMSSW_8_0_12/src/
git cms-addpkg PhysicsTools/JetMCAlgos
cd PhysicsTools/
git checkout CMSSW_8_0_X
cd ../
mv ../../Analyzers .
mv Analyzers/GenTtbarCategorizer_cfi.py PhysicsTools/JetMCAlgos/python/
mv Analyzers/GenTtbarCategorizer.cc PhysicsTools/JetMCAlgos/plugins/
git clone https://github.com/abhisek92datta/MiniAOD.git
cd MiniAOD
git checkout CMSSW_8_0_8
cd ..
scram b -j 32
scram b -j 32
