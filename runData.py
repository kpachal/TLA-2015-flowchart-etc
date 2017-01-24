import sys
import os
import subprocess

functions = ["our3Par","our4Par","UA2","Multijet9","GammaGamma3ParThird","GammaGamma4ParThird"] 
#functions = ["our4Par"]

#markLows = [442,459,477,495,514,533,553,573,594]
#markLows = [363,378,393,409,425,442,459,477,495]
markLows = [363,378,393,409,425,442,459,477,495]
markHigh = 1237

permitWindow = True

#commandTemplate = "SearchPhase --config {0} --file /cluster/warehouse/kpachal/TLA2016/samples/data/mjj_fullDataset.root --histName mjj --outputfile Results_Data_shorter/windowStat/SearchResultData_{1}_{2}.root 2>/dev/null\n"
#commandTemplate = "SearchPhase --config {0} --file /cluster/warehouse/kpachal/TLA2016/samples/data/mjj_DEF.root --histName mjj --outputfile Results_Data_shorter/windowStat/SearchResultData_{1}_{2}.root 2>/dev/null\n"
#commandTemplate = "SearchPhase --config {0} --file /cluster/warehouse/kpachal/TLA2016/samples/data/mjj_fullDataset_03.root --histName mjj --outputfile Results_Data/windowStat/SearchResultData_{1}_{2}.root 2>/dev/null\n"
commandTemplate = "SearchPhase --config {0} --file /cluster/warehouse/kpachal/TLA2016/samples/data/mjj_DEF_03.root --histName mjj --outputfile Results_Data/windowStat/SearchResultData_{1}_{2}.root 2>/dev/null\n"

batchScript = "scripts/Step2_BatchScript_Template.sh"

if permitWindow :
  commandTemplate = commandTemplate.replace("windowStat/","permitWindow/")
else :
  commandTemplate = commandTemplate.replace("windowStat/","noWindow/")

for function in functions :

  for markLow in markLows :

#    extraString = "fullDataset_from{0}".format(markLow+1)
#    extraString = "periodsDEF_from{0}".format(markLow+1)
#    extraString = "fullDataset_yStar0p3_from{0}".format(markLow+1)
    extraString = "periodsDEF_yStar0p3_from{0}".format(markLow+1)
    if permitWindow :
      extraString = extraString+"_permitWindow"

    configInName = "configurations/CompareFunctionConfigs/Step1_SearchPhase_{0}.config".format(function)
    configOutName = "submitConfigs/Step1_SearchPhase_{0}_{1}.config".format(function,extraString)
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
          line = "nPseudoExp 10000\n"
        configOut.write(line)
    configOut.close()

    modcommand = "chmod 744 {0}".format(configOutName)
    subprocess.call(modcommand,shell=True)
 
    thisCommand = commandTemplate.format(configOutName,function,extraString)
    batchtempname = "submitScripts/submitScript_{0}_{1}.sh".format(function,extraString)

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
