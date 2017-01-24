#!/usr/bin/env python

import os
import ROOT
from art.morisot import Morisot
from array import array
import sys
import numpy as np
import math


# Initialize painter
myPainter = Morisot()
myPainter.setColourPalette("Teals")
myPainter.setLabelType(2)
myPainter.setEPS(True)

luminosity = 3500
CME = 13000

doIndividualPlots = False
doSignalVersion = True

class searchFileData :

  def __init__(self,filename,permitWindow=False) :

    self.permitWindow = permitWindow

    searchInputFile = ROOT.TFile.Open(filename,"READ")

    # Retrieve search phase inputs
    self.basicData = searchInputFile.Get("basicData")
    self.basicData.SetDirectory(0)
    self.basicBkgFrom4ParamFit = searchInputFile.Get("basicBkgFrom4ParamFit")
    self.basicBkgFrom4ParamFit.SetDirectory(0)
    self.residualHist = searchInputFile.Get("residualHist")
    self.residualHist.SetDirectory(0)
    self.relativeDiffHist = searchInputFile.Get("relativeDiffHist")
    self.relativeDiffHist.SetDirectory(0)
    self.sigOfDiffHist = searchInputFile.Get("sigOfDiffHist")
    self.sigOfDiffHist.SetDirectory(0)
    self.logLikelihoodPseudoStatHist = searchInputFile.Get("logLikelihoodStatHistNullCase")
    self.logLikelihoodPseudoStatHist.SetDirectory(0)
    self.chi2PseudoStatHist = searchInputFile.Get("chi2StatHistNullCase")
    self.chi2PseudoStatHist.SetDirectory(0)
    self.bumpHunterStatHist = searchInputFile.Get("bumpHunterStatHistNullCase")
    self.bumpHunterStatHist.SetDirectory(0)
    self.bumpHunterTomographyPlot = searchInputFile.Get('bumpHunterTomographyFromPseudoexperiments')

    bumpHunterStatOfFitToData = searchInputFile.Get('bumpHunterStatOfFitToData')
    logLOfFitToDataVec = searchInputFile.Get('logLOfFitToData')
    chi2OfFitToDataVec = searchInputFile.Get('chi2OfFitToData')
    statOfFitToData = searchInputFile.Get('bumpHunterPLowHigh')
    self.logLOfFitToData = logLOfFitToDataVec[0]
    self.logLPVal = logLOfFitToDataVec[1]
    self.chi2OfFitToData = chi2OfFitToDataVec[0]
    self.chi2PVal = chi2OfFitToDataVec[1]
    self.bumpHunterStatFitToData = statOfFitToData[0]
    self.bumpHunterPVal = bumpHunterStatOfFitToData[1]
    self.bumpLowEdge = statOfFitToData[1]
    self.bumpHighEdge = statOfFitToData[2]

    self.NDF = searchInputFile.Get('NDF')[0]

    excludeWindowNums = searchInputFile.Get('excludeWindowNums')
    self.excludeWindow = int(excludeWindowNums[0]+0.5)
    self.bottomWindowEdge = excludeWindowNums[1]
    self.topWindowEdge = excludeWindowNums[2]

    if (self.excludeWindow and self.permitWindow) :
      statsOfRemainingSpectrum = searchInputFile.Get("BHLogLAndChi2OfRemainderAfterWindow")
      self.BHPValRemainder = statsOfRemainingSpectrum[0]
      self.LogLPValRemainder = statsOfRemainingSpectrum[1]
      self.Chi2PValRemainder = self.calculateRemainingChi2() #statsOfRemainingSpectrum[2]

    searchInputFile.Close()

  def getPValErrs(self) :

    # (DeltaX/X)^2 = (1/DeltaX)^2 = 1/X: set errors
    nRightBH = self.bumpHunterStatHist.Integral(self.bumpHunterStatHist.FindBin(self.bumpHunterStatFitToData),self.bumpHunterStatHist.GetNbinsX())
    nLeftBH = self.bumpHunterStatHist.Integral() - nRightBH
    if nRightBH > 0 and nLeftBH > 0 : deltaPvalBH = self.bumpHunterPVal * math.sqrt(1/nRightBH + 1/nLeftBH)
    else : deltaPvalBH = 0
    nRightChi2 = self.chi2PseudoStatHist.Integral(self.chi2PseudoStatHist.FindBin(self.chi2OfFitToData),self.chi2PseudoStatHist.GetNbinsX())
    nLeftChi2 = self.chi2PseudoStatHist.Integral() - nRightChi2
    if nRightChi2 > 0 and nLeftChi2 > 0 : deltaPvalChi2 = self.chi2PVal * math.sqrt(1/nRightChi2 + 1/nLeftChi2)
    else : deltaPvalChi2 = 0
    nRightLogL = self.logLikelihoodPseudoStatHist.Integral(self.logLikelihoodPseudoStatHist.FindBin(self.logLOfFitToData),self.logLikelihoodPseudoStatHist.GetNbinsX())
    nLeftLogL = self.logLikelihoodPseudoStatHist.Integral() - nRightLogL
    if nRightLogL > 0 and nLeftLogL > 0 : deltaPvalLogL = self.logLPVal * math.sqrt(1/nRightLogL + 1/nLeftLogL)
    else : deltaPvalLogL = 0

    return deltaPvalBH,deltaPvalChi2,deltaPvalLogL

  def calculateRemainingChi2(self) :

    firstBin = 0
    for bin in range(1,self.basicBkgFrom4ParamFit.GetNbinsX()+2) :
      firstBin = bin
      if self.basicBkgFrom4ParamFit.GetBinContent(bin) > 0 :
        break
    lastBin = 0
    for bin in range(self.basicBkgFrom4ParamFit.GetNbinsX()+1,0,-1) :
      lastBin = bin
      if self.basicBkgFrom4ParamFit.GetBinContent(bin) > 0 :
        break
    firstWindowBin = 0
    lastWindowBin = 0
    if self.excludeWindow :
      for bin in range(1,self.basicBkgFrom4ParamFit.GetNbinsX()+2) :
        if math.fabs(self.basicBkgFrom4ParamFit.GetBinLowEdge(bin) - self.bottomWindowEdge) < 0.1 :
          firstWindowBin = bin
        if math.fabs(self.basicBkgFrom4ParamFit.GetBinLowEdge(bin)+self.basicBkgFrom4ParamFit.GetBinWidth(bin) - self.topWindowEdge) < 0.1 :
          lastWindowBin = bin

    answer = 0
    for bin in range(firstBin,lastBin+1) :

      if self.excludeWindow and bin >= firstWindowBin and bin <= lastWindowBin : continue

      d = self.basicData.GetBinContent(bin)
      if (d==0) : continue
      b = self.basicBkgFrom4ParamFit.GetBinContent(bin)
      deltaB = self.basicBkgFrom4ParamFit.GetBinError(bin)

      term = (d - b) / math.sqrt(b+deltaB*deltaB)
      answer = answer + (term*term)

    nRightChi2 = self.chi2PseudoStatHist.Integral(self.chi2PseudoStatHist.FindBin(answer),self.chi2PseudoStatHist.GetNbinsX())
    nTotal = self.chi2PseudoStatHist.Integral()
    return float(nRightChi2)/float(nTotal)


  def makeSearchPhasePlots(self,low,folder,ext,funcName=[]) :
 
    firstBin = self.basicData.FindBin(low)
    lastBin = self.basicData.FindBin(1234)+1

    if self.excludeWindow and self.permitWindow :
      myPainter.drawDataAndFitOverSignificanceHist(self.basicData,self.basicBkgFrom4ParamFit,self.residualHist,\
         'm_{jj} [TeV]','Events','Significance','{0}/figure1'.format(folder)+ext,\
         luminosity,13,low,1234,firstBin,lastBin,True,self.bumpLowEdge,self.bumpHighEdge,doWindowLimits=self.excludeWindow,windowLow=self.bottomWindowEdge,windowHigh=self.topWindowEdge,extraLegendLines=funcName)
    else :
      myPainter.drawDataAndFitOverSignificanceHist(self.basicData,self.basicBkgFrom4ParamFit,self.residualHist,\
         'm_{jj} [TeV]','Events','Significance','{0}/figure1'.format(folder)+ext,\
         luminosity,13,low,1234,firstBin,lastBin,True,self.bumpLowEdge,self.bumpHighEdge,extraLegendLines=funcName)
    myPainter.drawPseudoExperimentsWithObservedStat(self.logLikelihoodPseudoStatHist,float(self.logLOfFitToData),self.logLPVal,0,luminosity,13,\
       'logL statistic','Pseudo-exeperiments',"{0}/logLStatPlot".format(folder)+ext)
    myPainter.drawPseudoExperimentsWithObservedStat(self.chi2PseudoStatHist,float(self.chi2OfFitToData),self.chi2PVal,0,luminosity,13,\
       "#chi^{2}",'Pseudo-exeperiments',"{0}/chi2StatPlot".format(folder)+ext)
    myPainter.drawPseudoExperimentsWithObservedStat(self.bumpHunterStatHist,float(self.bumpHunterStatFitToData),self.bumpHunterPVal,0,luminosity,13,\
       'BumpHunter','Pseudo-exeperiments',"{0}/bumpHunterStatPlot".format(folder)+ext)
    myPainter.drawBumpHunterTomographyPlot(self.bumpHunterTomographyPlot,"{0}/bumpHunterTomographyPlot".format(folder)+ext)


def GetZVal (p, excess) :
  #the function normal_quantile converts a p-value into a significance,
  #i.e. the number of standard deviations corresponding to the right-tail of 
  #a Gaussian
  if excess :
    zval = ROOT.Math.normal_quantile(1-p,1);
  else :
    zval = ROOT.Math.normal_quantile(p,1);

  return zval


def MakeHistoFromStats(statistics) :

  nentries = len(statistics)
  nBins = int(float(nentries)/10.0)

  maxVal = max(statistics)
  minVal = min(statistics)
  axisrange = maxVal - minVal;

  thismin = minVal-0.05*axisrange;
  thismax = maxVal+0.05*axisrange;

  statPlot = ROOT.TH1D("statPlot","",nBins,thismin,thismax)
  for val in range(len(statistics)) :
    statPlot.Fill(statistics[val])

  return statPlot

def rms(x):
  mean = np.mean(x)
  sumvec = [(val-mean)*(val-mean) for val in x]
  stddev = np.sqrt(np.sum(sumvec)/np.size(x))
  return stddev

massCouplingDict = {}
seedRange = []

if doSignalVersion :
  massCouplingDict["0p10"] = ["0p45","0p55","0p65","0p75","0p85"]
  massCouplingDict["0p15"] = ["0p45","0p55","0p65","0p75","0p85","0p95"]
  massCouplingDict["0p20"] = ["0p55","0p65","0p75","0p85","0p95"]
  massCouplingDict["0p30"] = ["0p45","0p55","0p65","0p75","0p85","0p95","1p05"]
  massCouplingDict["0p40"] = ["0p65","0p75","0p85","0p95","1p05"]

  seedRange = [1,5]

else :
  massCouplingDict["None"] = ["None"]

  seedRange = [1,99]
#  seedRange = [12,12]

#markLows = [378,393,409,425,442,459,477]
#markLows = [425,442,459,477]
markLows = [442,459,477,495,514,1000]
markHigh = 1250

pValCutoff = 0.05

functions = ["our3Par","our4Par","UA2","Multijet9","GammaGamma3ParThird","GammaGamma4ParThird"]
#functions = ["our4Par","UA2"]
specialext = ""

# Get input
filenametemp = "/cluster/warehouse/kpachal/TLA2016/StatisticalAnalysis/Bayesian/TestResults_dressRehearsal/{0}/{1}"
#filenametemp = "/cluster/warehouse/kpachal/TLA2016/StatisticalAnalysis/Bayesian/TestResults_dressRehearsal_yStar0p3/{0}/{1}"

countValid=0
countBumps=0
countProblems=0
countNeedMoreBins=0
countMissingFiles=0
needMoreBins=[]
countTotal = 0

BHPValDistribution = ROOT.TH1D("BHPValDist","BHPValDist",20,0.0,1.0)

resultDict = {}

for coupling in massCouplingDict.keys() :
  resultDict[coupling] = {}
  for mass in massCouplingDict[coupling] :
    resultDict[coupling][mass] = {}
    resultDict[coupling][mass]["bumps"] = {}
    resultDict[coupling][mass]["funcs"] = {}
    resultDict[coupling][mass]["nomFile"] = {}

# Loop over no signal and every signal we had.
for seed in range(seedRange[0],seedRange[1]+1) :

 if doSignalVersion : basefolder = "Plots_dressRehearsal_withSig/" 
 else : basefolder = "Plots_dressRehearsal_yStar0p3/"
# else : basefolder = "Plots_dressRehearsal/"
 folder = basefolder + "/seed{0}/".format(seed)

 CrossCheckFolder = folder+"/IndividualPlots"
 
 # make plots folder i.e. make folder extension
 if not os.path.exists(folder):
   os.makedirs(folder)
 if not os.path.exists(CrossCheckFolder) :
   os.makedirs(CrossCheckFolder)

 for coupling in massCouplingDict.keys() :
  for mass in massCouplingDict[coupling] :

    if (mass is "None" and coupling is not "None") or (mass is not "None" and coupling is "None") :
      continue
    anyFileReal = False
    countTotal = countTotal+1

    windowOK = False
    leadHasBump = False
    subleadHasBump = False

    functionGraphDicts = {}
    withWindowGraphDicts = {}
    for function in functions :

      chi2Graph =  ROOT.TGraphErrors()
      logLGraph =  ROOT.TGraphErrors()
      BHGraph =  ROOT.TGraphErrors()
      windowChi2Graph = ROOT.TGraphErrors()
      windowLogLGraph = ROOT.TGraphErrors()
      windowBHGraph = ROOT.TGraphErrors()
      functionGraphDicts[function] = {}
      functionGraphDicts[function]["Chi2"] = chi2Graph
      functionGraphDicts[function]["BH"] = BHGraph
      functionGraphDicts[function]["LogL"] = logLGraph
      withWindowGraphDicts[function] = {}
      withWindowGraphDicts[function]["Chi2"] = chi2Graph
      withWindowGraphDicts[function]["BH"] = BHGraph
      withWindowGraphDicts[function]["LogL"] = logLGraph

    # Increase left-hand side step by step until two functions at least pass cutoff,
    # under consideration of certain limit requirements.
    pickLow = -1
    index = -1
    for low in sorted(markLows) :
      low = low+1  

      print "\nBeginning range from",low
 
      index = index+1
      nBetter = 0
      bestPVal = -1
      bestFunc = ""
      secondBestPVal = -1
      secondBestFunc = ""

      # Step 1: check if I can find two functions with good enough chi2 p-values
      # without allowing a window exclusion.
      foundRealBump = False

      worksWithoutWindow = []
      worksTotal = []

      for function in functions :

        midfolder = "withSignal_noWindow"
        filename = "SearchResultWithSig_{0}_mass{1}_gSM{2}_seed{3}_withSig_from{4}.root".format(function,mass,coupling,seed,low)
        if mass is "None" and coupling is "None" :
          midfolder = "noSignal_noWindow"
          filename = "SearchResultNoSig_{0}_seed{1}_noSig_from{2}.root".format(function,seed,low)
        infilename = filenametemp.format(midfolder,filename)

        try :       
          theseData = searchFileData(infilename,False)
          anyFileReal = True
        except :
          functionGraphDicts[function]["Chi2"].SetPoint(index,low,0.0)
          functionGraphDicts[function]["LogL"].SetPoint(index,low,0.0)
          functionGraphDicts[function]["BH"].SetPoint(index,low,0.0)
          functionGraphDicts[function]["Chi2"].SetPointError(index,0.0,0.0)
          functionGraphDicts[function]["LogL"].SetPointError(index,0.0,0.0)
          functionGraphDicts[function]["BH"].SetPointError(index,0.0,0.0)
          continue

        functionGraphDicts[function]["Chi2"].SetPoint(index,low,theseData.chi2PVal)
        functionGraphDicts[function]["LogL"].SetPoint(index,low,theseData.logLPVal)
        functionGraphDicts[function]["BH"].SetPoint(index,low,theseData.bumpHunterPVal)

        deltaPvalBH,deltaPvalChi2,deltaPvalLogL = theseData.getPValErrs()
        functionGraphDicts[function]["BH"].SetPointError(index,0.0,deltaPvalBH)
        functionGraphDicts[function]["Chi2"].SetPointError(index,0.0,deltaPvalChi2)
        functionGraphDicts[function]["LogL"].SetPointError(index,0.0,deltaPvalLogL)

        ext = "{0}_mass{1}_gSM{2}_seed{3}_{4}up_noWindow".format(function,mass,coupling,seed,low)+specialext

        # Search phase plots
        if (doIndividualPlots) :
          theseData.makeSearchPhasePlots(low,CrossCheckFolder,ext)

        # How many successful functions do we have here? Which ones?
        if theseData.chi2PVal > pValCutoff :
          nBetter = nBetter+1
          worksWithoutWindow.append(function)
          worksTotal.append(function)
          if theseData.chi2PVal > bestPVal :
            # Bump best values down to second best
            secondBestPVal = bestPVal
            secondBestFunc = bestFunc
            # Store new vals
            bestPVal = theseData.chi2PVal
            bestFunc = function
            print "Selecting best function",bestFunc,"without window exclusion."
          elif theseData.chi2PVal > secondBestPVal :
            secondBestPVal = theseData.chi2PVal
            secondBestFunc = function
            print "Selecting second best function",secondBestFunc,"without window exclusion"

      # Step 2: Look at results for all functions for this start point. Do we have at least two good functions?
      if nBetter > 1 and pickLow < 0 :
        pickLow = low
        windowOK = True
        print "List of working functions without allowing windows:",worksTotal
        break

      # Step 3: If no, allow windows and then check again
      for function in functions :

        midfolder = "withSignal_permitWindow"
        filename = "SearchResultWithSig_{0}_mass{1}_gSM{2}_seed{3}_withSig_from{4}_permitWindow.root".format(function,mass,coupling,seed,low)
        if mass is "None" and coupling is "None" :
          midfolder = "noSignal_permitWindow"
          filename = "SearchResultNoSig_{0}_seed{1}_noSig_from{2}_permitWindow.root".format(function,seed,low)
        infilename = filenametemp.format(midfolder,filename)

        try :
          theseData = searchFileData(infilename,True)
          anyFileReal = True
        except :
          withWindowGraphDicts[function]["Chi2"].SetPoint(index,low,0.0)
          withWindowGraphDicts[function]["LogL"].SetPoint(index,low,0.0)
          withWindowGraphDicts[function]["BH"].SetPoint(index,low,0.0)
          withWindowGraphDicts[function]["Chi2"].SetPointError(index,0.0,0.0)
          withWindowGraphDicts[function]["LogL"].SetPointError(index,0.0,0.0)
          withWindowGraphDicts[function]["BH"].SetPointError(index,0.0,0.0)
          continue

        withWindowGraphDicts[function]["Chi2"].SetPoint(index,low,theseData.chi2PVal)
        withWindowGraphDicts[function]["LogL"].SetPoint(index,low,theseData.logLPVal)
        withWindowGraphDicts[function]["BH"].SetPoint(index,low,theseData.bumpHunterPVal)

        deltaPvalBH,deltaPvalChi2,deltaPvalLogL = theseData.getPValErrs()
        withWindowGraphDicts[function]["BH"].SetPointError(index,0.0,deltaPvalBH)
        withWindowGraphDicts[function]["Chi2"].SetPointError(index,0.0,deltaPvalChi2)
        withWindowGraphDicts[function]["LogL"].SetPointError(index,0.0,deltaPvalLogL)

        ext = "{0}_mass{1}_gSM{2}_seed{3}_{4}up_permitWindow".format(function,mass,coupling,seed,low)+specialext

        # Search phase plots
        if (doIndividualPlots) :
          theseData.makeSearchPhasePlots(low,CrossCheckFolder,ext)

        # How many successful functions do we have here? Which ones?
        # Function is not counted if the window it chooses is > 20% width/central mass ratio!
        # Skip to the next function.
        foundBumpThisFunc = False
        if theseData.excludeWindow :
          width = (theseData.topWindowEdge - theseData.bottomWindowEdge)/2.0
          testmass = theseData.bottomWindowEdge + width
          if width/testmass > 0.20 :
            continue

          # Step 4: what is p-value of remaining spectrum? 
          # If still bad after 20% window exclusion we do not
          # accept this fit.
          if theseData.BHPValRemainder > 0.01 and theseData.Chi2PValRemainder > 0.05 : 
            foundBumpThisFunc = True

        # Step 5: If we know that either no window needed
        # to be removed or a window was removed that was narrower than 20%
        # and left a remaining spectrum with p-vale > 0.01, then
        # the window is OK and this is our final background estimate!
        # If no window was excluded, we get the same thing as before.
        # If a window was excluded, need to compare p-value of remainder of 
        # spectrum when deciding if this should be nominal.
        if (not theseData.excludeWindow) or foundBumpThisFunc :
          windowOK = True

        if foundBumpThisFunc :

          hasCurrentLeadWindow = False
          hasCurrentSubleadWindow = False

          # Only allowed to add it where no solution was found without a window.
          # Careful not to double count functions we already chose: always prefer
          # to use version without a window exclusion if we can
          if function in worksWithoutWindow :
            continue

          # Haven't already counted a success for this function, so we can add it 
          if theseData.Chi2PValRemainder > pValCutoff :
            nBetter = nBetter+1
            worksTotal.append(function)
            if theseData.Chi2PValRemainder > bestPVal and (hasCurrentLeadWindow or bestPVal < 0) :
              # Bump best values down to second best
              secondBestPVal = bestPVal
              secondBestFunc = bestFunc
              hasCurrentSubleadWindow = hasCurrentLeadWindow
              # Store new bests
              bestPVal = theseData.Chi2PValRemainder
              bestFunc = function
              hasCurrentLeadWindow = True
              foundRealBump = True
              print "Replaced lead with function",function,"with window. Remaining p-val is",theseData.Chi2PValRemainder
            elif theseData.Chi2PValRemainder > secondBestPVal and (hasCurrentSubleadWindow or secondBestPVal < 0) :
              secondBestPVal = theseData.Chi2PValRemainder
              secondBestFunc = function
              hasCurrentSubleadWindow = True
              print "Replaced sublead with function",function,"with window. Remaining p-val is",theseData.Chi2PValRemainder

      # Step 2: Look at results for all functions for this start point. Do we have at least two good functions?
      # If so don't look at other start points.
      if windowOK and nBetter > 1 and pickLow < 0:
        pickLow = low
        break

      print "List of working functions after allowing windows:",worksTotal

    # Have now looped over all of the starting edges and we can tell if we have anything that's good enough.
    # Plot the graphs of p-values.
    Chi2GraphList = []
    LogLGraphList = []
    BHGraphList = []
    for function in functions :
      Chi2GraphList.append(functionGraphDicts[function]["Chi2"])
      BHGraphList.append(functionGraphDicts[function]["BH"])
      LogLGraphList.append(functionGraphDicts[function]["LogL"])

    if not anyFileReal :
      countMissingFiles = countMissingFiles+1
      resultDict[coupling][mass]["bumps"][seed] = 0
      continue

    ext = "_withSignal_mass{0}_gSM{1}".format(mass,coupling)
    if mass is "None" and coupling is "None" :
      ext = "_noSignal"
    myPainter.drawSeveralObservedLimits(Chi2GraphList,functions,folder+"/Chi2Graph"+ext+specialext,"Low edge of fit [GeV]","#chi^{2} p-value",luminosity,13,430,\
                        550,0.0,1.0,[],False,False,True,"Right","BottomL",True,[0.05])
    myPainter.drawSeveralObservedLimits(LogLGraphList,functions,folder+"/LogLGraph"+ext+specialext,"Low edge of fit [GeV]","LogL p-value",luminosity,13,430,\
                        550,0.0,1.0,[],False,False,True,"Right","BottomL",True,[0.05])
    myPainter.drawSeveralObservedLimits(BHGraphList,functions,folder+"/BHGraph"+ext+specialext,"Low edge of fit [GeV]","BH p-value",luminosity,13,430,\
                        550,0.0,1.0,[],False,False,True,"Right","BottomL",True,[0.05])

    # If we had a bump: we want to test n more bins for ONLY cases with no window. 
    # If we find one that works, use it instead.
    print "About to peek in bins"
    if foundRealBump :
      nExtraBins = 2
      extraBinList = markLows[markLows.index(pickLow-1)+1:markLows.index(pickLow-1)+1+nExtraBins]
      print "Found a bump with start point",pickLow
      print "So now look for better fits with no window in start vals:",extraBinList
      for low in sorted(extraBinList) :
        low = low+1  
        print "Checking start",low

        nBonusBetter = 0
        bestPValB = -1
        bestFuncB = ""
        secondBestPValB = -1
        secondBestFuncB = ""

        for function in functions :

          midfolder = "withSignal_noWindow"
          filename = "SearchResultWithSig_{0}_mass{1}_gSM{2}_seed{3}_withSig_from{4}.root".format(function,mass,coupling,seed,low)
          if mass is "None" and coupling is "None" :
            midfolder = "noSignal_noWindow"
            filename = "SearchResultNoSig_{0}_seed{1}_noSig_from{2}.root".format(function,seed,low)
          infilename = filenametemp.format(midfolder,filename)

          try :       
            theseData = searchFileData(infilename,False)
            anyFileReal = True
          except :
            continue

          # How many successful functions do we have here? Which ones?
          if theseData.chi2PVal > pValCutoff :
            nBonusBetter = nBonusBetter+1
            if theseData.chi2PVal > bestPValB :
              # Bump best values down to second best
              secondBestPValB = bestPValB
              secondBestFuncB = bestFuncB
              # Store new vals
              bestPValB = theseData.chi2PVal
              bestFuncB = function
              print "Found best function",bestFunc,"in bonus bin without window exclusion."
            elif theseData.chi2PVal > secondBestPValB :
              secondBestPValB = theseData.chi2PVal
              secondBestFuncB = function
              print "Found second best function",secondBestFunc,"in bonus bin without window exclusion"

        print "nBonusBetter is",nBonusBetter
        if nBonusBetter > 0 :
          # we keep these fit results instead
          pickLow = low
          secondBestFunc = bestFunc
          bestFunc = bestFuncB
          bestPVal = bestPValB
          if nBonusBetter > 1 :
            secondBestFunc = secondBestFuncB
            secondBestPVal = secondBestPValB
          foundRealBump = False

          # and if there are two, we don't need the other bonus bins, if any
          if nBonusBetter > 1 :
            break  

    # Now we know what was the last mass for which our fit was good enough.
    # Our nominal search phase result for this signal is the one which has
    # a window permitted and corresponds to this point.

    foundGoodStart = True
    if pickLow < 0 :
      foundGoodStart = False
      pickLow = markLows[-1]+1

    nomfolder = "withSignal_permitWindow"
    nomfilename = "SearchResultWithSig_{0}_mass{1}_gSM{2}_seed{3}_withSig_from{4}_permitWindow.root".format(bestFunc,mass,coupling,seed,pickLow)
    if mass is "None" and coupling is "None" :
      nomfolder = "noSignal_permitWindow"
      nomfilename = "SearchResultNoSig_{0}_seed{1}_noSig_from{2}_permitWindow.root".format(bestFunc,seed,pickLow)
    nomfilename = filenametemp.format(nomfolder,nomfilename)

    if os.path.isfile(nomfilename) :
      theseData = searchFileData(nomfilename,True)
    else :
      countMissingFiles = countMissingFiles+1
      resultDict[coupling][mass]["bumps"][seed] = 0
      continue

    ext = "nominalChoice_mass{0}_gSM{1}_seed{2}".format(mass,coupling,seed)+specialext

    theseData.makeSearchPhasePlots(pickLow,CrossCheckFolder,ext,[bestFunc])

    # Count
    if foundRealBump: 
      countBumps = countBumps+1
    elif foundGoodStart :
      countValid = countValid+1
    else :
      countNeedMoreBins = countNeedMoreBins+1
      needMoreBins.append(seed)

    resultDict[coupling][mass]["bumps"][seed] = foundRealBump
    resultDict[coupling][mass]["funcs"][seed] = [bestFunc,secondBestFunc]
    resultDict[coupling][mass]["nomFile"][seed] = nomfilename

  print "Cases with valid bumps removed:",countBumps
  print "Cases with good start values found:",countValid
  print "Cases where more data bins required:",countNeedMoreBins,". Relevant seeds:",needMoreBins
  print "Cases where relevant files are missing:",countMissingFiles
  print "Total number of cases was:",countTotal

  BHPValDistribution.Fill(bestPVal)

counter = 0
print "\nNominal and alternate functions chosen were:"
for coupling in sorted(resultDict.keys()) :
  for mass in sorted(resultDict[coupling].keys()) :
    for index in resultDict[coupling][mass]["funcs"].keys() :
      list = resultDict[coupling][mass]["funcs"][index]
      print list[0],list[1]
      if (list[0] == "our4Par" and list[1] == "UA2") or (list[0] == "UA2" and list[1] == "our4Par") : counter = counter+1
print "\n Percent of cases where 4Par and UA2 are selected:",float(counter)/99.0


if doSignalVersion :
  print "\n"
  print "Summary of discoverability:\n"
  print "Coupling\tMass\tFraction of cases discovered",
  for coupling in sorted(resultDict.keys()) :
    print "\n"
    for mass in sorted(resultDict[coupling].keys()) :
      bumpDict = resultDict[coupling][mass]["bumps"].values()
      print coupling,"\t",mass,"\t",float(sum(bumpDict))/float(len(bumpDict))


else : 
  print "\n"
  print "Bumps found in seeds",
  for coupling in sorted(resultDict.keys()) :
    for mass in sorted(resultDict[coupling].keys()) :
      for seed in range(seedRange[0],seedRange[1]+1) :
        if resultDict[coupling][mass]["bumps"][seed] == 1 :
          print seed,",",

  myPainter.drawPseudoExperimentsWithObservedStat(BHPValDistribution,-1,0.5,0.5,luminosity,CME,"BH p-value","Test cases",basefolder+"BHPValDistAcrossSeeds")

# Print a file with info on what the nominal selection was
filename = "chosenNominalFiles_{0}.txt"
if doSignalVersion : filename = filename.format("withSignal")
else : filename = filename.format("noSignal")
with open(filename,"w") as fout :
  for coupling in sorted(resultDict.keys()) :
    for mass in sorted(resultDict[coupling].keys()) :
      for seed in range(seedRange[0],seedRange[1]+1) :
        try :
          line = "{0} {1} {2} {3}".format(seed,coupling,mass,resultDict[coupling][mass]["nomFile"][seed])
          fout.write(line+"\n")
        except :
          continue

# For table
if doSignalVersion :

  print "\n"
  print "Table for note: fraction of cases discovered\n"
  print "Mass & gSM=0.1 & gSM=0.15 & gSM=0.2 & gSM=0.3 & gSM=0.4 \\\\"
  newdict = {}
  for coupling in sorted(resultDict.keys()) :
    for mass in sorted(resultDict[coupling].keys()) :
      bumpDict = resultDict[coupling][mass]["bumps"].values()
      if not mass in newdict.keys() :
        newdict[mass] = {}
      newdict[mass][coupling] = float(sum(bumpDict))/float(len(bumpDict))

  for mass in sorted(newdict.keys()) :
    string = "{0} &".format(mass)
    if "0p10" in newdict[mass].keys() : string = string + " {0} &".format(newdict[mass]["0p10"])
    else : string = string+ " - &"
    if "0p15" in newdict[mass].keys() : string = string + " {0} &".format(newdict[mass]["0p15"])
    else : string = string+ " - &"
    if "0p20" in newdict[mass].keys() : string = string + " {0} &".format(newdict[mass]["0p20"])
    else : string = string+ " - &"
    if "0p30" in newdict[mass].keys() : string = string + " {0} &".format(newdict[mass]["0p30"])
    else : string = string+ " - &"
    if "0p40" in newdict[mass].keys() : string = string + " {0} &".format(newdict[mass]["0p40"])
    else : string = string+ " - &"
    print string+"\\\\"

print "Done."
