for i in {1..57}
do
   echo "getting file $i "
   #xrdcp -v root://se.cis.gov.pl:11000//store/user/kbunkow/ZprimeToMuMu_M-6000_TuneCP5_14TeV-pythia8/CRAB_omtf_nn_MC_analysis_ZprimeToMuMu_PU140_v3_t98/201125_175704/0000/omtfAnalysis2_$i.root ./crab_omtf_nn_MC_analysis_ZprimeToMuMu_PU140_v3_t98/results/.
   xrdcp -v root://se.cis.gov.pl:11000//store/user/kbunkow/Nu_E10-pythia8-gun/CRAB3_omtf_nn_MC_analysis_SingleNeutrino_PU200_v3_t100/201126_215339/0000/omtfAnalysis2_$i.root ./crab_omtf_nn_MC_analysis_SingleNeutrino_PU200_v3_t100/results/.

done
