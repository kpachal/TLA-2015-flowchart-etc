import sys
import os
import subprocess

seeds = [5]#[1,2,3,4,5]
functions = ["our3Par","our4Par","UA2","Multijet9","GammaGamma3ParThird","GammaGamma4ParThird"] 

masses = ["0p55"]#["0p45","0p55","0p65","0p75","0p85","0p95","1p05"]
couplings = ["0p20","0p30","0p40"]#["0p10","0p15","0p20","0p30","0p40"]

permitWindow = True

markLows = [442,459,477,495,514]
markHigh = 1237

commandTemplate = "SearchPhase --config {0} --file /cluster/warehouse/kpachal/TLA2016/samples/bkgPlusSig/MixedSamples/dataLikeHistos_bkgPlusSig_seed{4}.root --histName mjj_3.5_withSignal_mass{2}_gSM{3} --outputfile TestResults_dressRehearsal/withSignal_noWindow/SearchResultWithSig_{1}_mass{2}_gSM{3}_seed{4}_{5}.root 2>/dev/null\n"

batchScript = "scripts/Step2_BatchScript_Template.sh"

for seed in seeds :
 for function in functions :
  for markLow in markLows :

    extraString = "withSig_from{0}".format(markLow+1)
    if permitWindow :
      extraString = extraString+"_permitWindow"
      commandTemplate = commandTemplate.replace("noWindow","permitWindow")

    for mass in masses :
      for coupling in couplings :

        configInName = "configurations/CompareFunctionConfigs/Step1_SearchPhase_{0}.config".format(function)
        configOutName = "submitConfigs/Step1_SearchPhase_{0}_mass{1}_{2}_{3}_seed{4}.config".format(function,mass,coupling,seed,extraString)
        configOut = open(configOutName,'w')
        with open(configInName) as configInData :
          for line in configInData :
            if "minXForFit" in line :
              line = "minXForFit  {0}\n".format(markLow)
            elif "maxXForFit" in line :
              line = "maxXForFit  {0}\n".format(markHigh)
            if permitWindow and "permitWindow" in line :
              line = "permitWindow true\n" 
            if "nPseudoExp" in line :
              line = "nPseudoExp 100\n"
            configOut.write(line)
        configOut.close()

        modcommand = "chmod 744 {0}".format(configOutName)
        subprocess.call(modcommand,shell=True)
 
        thisCommand = commandTemplate.format(configOutName,function,mass,coupling,seed,extraString)
        batchtempname = "submitScripts/submitScript_{0}_{1}_{2}_seed{3}_{4}.sh".format(function,mass,coupling,seed,extraString)

        fbatchin = open(batchScript,'r')
        fbatchindata = fbatchin.read()
        fbatchin.close()
        fbatchout = open(batchtempname,'w')
        fbatchoutdata = fbatchindata.replace("ZZZ",thisCommand)
        fbatchout.write(fbatchoutdata)
        fbatchout.close()

        modcommand = "chmod 744 {0}".format(batchtempname)
        subprocess.call(modcommand,shell=True)
 
        submitcommand = "qsub {0}".format(batchtempname)
        print submitcommand
#        subprocess.call(submitcommand,shell=True)
