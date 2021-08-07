#import subprocess
#subprocess.call(["xrdcp", "root://cms-xrd-global.cern.ch//store/user/kbunkow/Mu_FlatPt2to100-pythia8-gun/CRAB3_omtf_nn_MC_analysis_MuFlatPt_PU200_v1_t26/200403_164429/0000/omtfAnalysis2_7.root ./crab_omtf_nn_MC_analysis_MuFlatPt_PU200_v1_t26/results/"])
   
import os

for i in range(1, 65, 1) :
    #os.system("xrdcp root://cms-xrd-global.cern.ch//store/user/kbunkow/Mu_FlatPt2to100-pythia8-gun/CRAB3_omtf_nn_MC_analysis_MuFlatPt_PU200_v1_t26/200403_164429/0000/omtfAnalysis2_"  + str(i) +".root ./crab_omtf_nn_MC_analysis_MuFlatPt_PU200_v1_t26/results/")
    print ("doing file " + str(i))
    os.system("xrdcp -v root://cms-xrd-global.cern.ch//store/user/kbunkow/Nu_E10-pythia8-gun/CRAB3_omtf_nn_MC_analysis_MuFlatPt_PU200_v2_t29/200408_202824/0000/omtfAnalysis2_"  + str(i) +".root ./crab_omtf_nn_MC_analysis_MuFlatPt_PU200_v2_t29/results/")

