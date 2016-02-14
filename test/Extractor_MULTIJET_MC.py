import sys, os
sys.path.insert(1, os.getcwd())

from Extractor_MULTIJET_common import *

#with open('../../../PatTopProduction/multijet_QCD_HT_fewEvents.list') as f:
#  files = f.readlines()


process = createExtractorProcess(True, True, useShiftCorrectedMET = True, globalTag = "76X_mcRun2_asymptotic_v12")

#process.source.fileNames = cms.untracked.vstring(files)
process.source.fileNames = cms.untracked.vstring(
        #'file:patTuple_1.root' # local file
        #'/store/mc/RunIIFall15MiniAODv2/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/PU25nsData2015v1_magnetOn_76X_mcRun2_asymptotic_v12-v1/00000/02B9457D-47BC-E511-B18F-E41D2D08DE30.root' # QCD Pt Flat test
        '/store/mc/RunIIFall15MiniAODv2/QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/10000/04A66FE4-47B9-E511-85C6-002590DB9216.root' # QCD HT test

    )
process.maxEvents.input = 1000
