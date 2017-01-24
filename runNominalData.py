import sys
import os
import subprocess

function = "UA2"
altfunction = "our4Par"

low =442
#low = 393

permitWindow = True

commandTemplate = "SearchPhase --config {0} --file /cluster/warehouse/kpachal/TLA2016/samples/data/mjj_fullDataset.root --histName mjj --outputfile Results_Data_Nominal/SearchResultData_{1}_{2}.root 2>/dev/null\n"
#commandTemplate = "SearchPhase --config {0} --file /cluster/warehouse/kpachal/TLA2016/samples/data/mjj_fullDataset_03.root --histName mjj --outputfile Results_Data_Nominal/SearchResultData_{1}_{2}.root 2>/dev/null\n"

batchScript = "scripts/Step2_BatchScript_Template.sh"

extraString = "fullDataset_from{0}".format(low+1)
#extraString = "fullDataset_yStar0p3_from{0}".format(low+1)

configInName = "configurations/CompareFunctionConfigs/Step1_SearchPhase_{0}.config".format(function)
altConfigName = "configurations/CompareFunctionConfigs/Step1_SearchPhase_{0}.config".format(altfunction)
configOutName = "submitConfigs/Step1_SearchPhase_nominal_{0}_{1}.config".format(function,extraString)
configOut = open(configOutName,'w')
with open(configInName) as configInData :
 with open(altConfigName) as altConfigData :
  for line in configInData :
    if "minXForFit" in line :
      line = "minXForFit  {0}\n".format(low)
    elif "maxXForFit" in line :
      line = "maxXForFit  {0}\n".format(1237)
    if permitWindow and "permitWindow" in line :
      line = "permitWindow true\n" 
    if "nPseudoExp" in line :
      line = "nPseudoExp 10000\n"

    if "doAlternateFunction" in line :
      line = "doAlternateFunction true\n"
      for nline in altConfigData :
        if "functionCode" in nline :
          line = line+nline.replace("function","alternateFunction")+"\n"
        elif "nParameters" in nline :
          line = line+nline.replace("nParameters","alternateNParameters")+"\n"
        elif "parameter" in nline and "altpar" not in nline :
          line = line+nline.replace("par","altpar")+"\n"
    elif "alt" in line :
      continue

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
