scram b -j 32

python pileupCalc.py -i '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/Cert_271036-275783_13TeV_PromptReco_Collisions16_JSON.txt' --inputLumiJSON '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt' --calcMode 'true' --minBiasXsec 71300 --maxPileupBin 50 --numPileupBins 50  output.root

root -l -q PU_data_hist_prod.C

g++ PU_hist_calc.cc -o PU_hist_calc

./PU_hist_calc

#cp PU_weights.txt ../../data/PU_weight/


