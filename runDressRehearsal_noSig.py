import sys
import os
import subprocess

seed = 3
functions = ["our3Par","our4Par","UA2","Multijet9","GammaGamma3ParThird","GammaGamma4ParThird"] 

permitWindow = False

#markLows = [442,459,477,495,514]
markLows = [442,459,477,495,514]
markHigh = 1237

commandTemplate = "SearchPhase --config {0} --file /cluster/warehouse/kpachal/TLA2016/samples/bkgPlusSig_yStar0p3/MixedSamples/dataLikeHistos_bkgPlusSig_seed{2}.root --histName mjj_lumi3500 --outputfile TestResults_dressRehearsal_yStar0p3/noSignal/SearchResultNoSig_{1}_seed{2}_{3}.root 2>/dev/null\n"

batchScript = "scripts/Step2_BatchScript_Template.sh"

if permitWindow :
  commandTemplate = commandTemplate.replace("noSignal/","noSignal_permitWindow/")
else :
  commandTemplate = commandTemplate.replace("noSignal/","noSignal_noWindow/")

for function in functions :

  for markLow in markLows :

    extraString = "noSig_from{0}".format(markLow+1)
    if permitWindow :
      extraString = extraString+"_permitWindow"

    configInName = "configurations/CompareFunctionConfigs/Step1_SearchPhase_{0}.config".format(function)
    configOutName = "submitConfigs/Step1_SearchPhase_{0}_seed{1}_{2}.config".format(function,seed,extraString)
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
 
    thisCommand = commandTemplate.format(configOutName,function,seed,extraString)
    batchtempname = "submitScripts/submitScript_{0}_seed{1}_{2}.sh".format(function,seed,extraString)

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
    subprocess.call(submitcommand,shell=True)
