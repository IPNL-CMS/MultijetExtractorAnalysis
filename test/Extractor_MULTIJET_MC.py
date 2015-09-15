import sys, os
sys.path.insert(1, os.getcwd())

from Extractor_MULTIJET_common import *

#with open('../../../PatTopProduction/multijet_QCD_HT_fewEvents.list') as f:
#  files = f.readlines()

#process = createExtractorProcess(True, True, useShiftCorrectedMET = True, globalTag = "PHYS14_25_V2")
process = createExtractorProcess(True, True, useShiftCorrectedMET = True, globalTag = "74X_mcRun2_asymptotic_v2")

#process.source.fileNames = cms.untracked.vstring(files)
process.source.fileNames = cms.untracked.vstring(
        #'file:patTuple_1.root'
        #'/store/user/apequegn/RelValTTbar_13/RelValTTbar_13TeV_18Mar15-v1/150318_170818/0000/patTuple_1.root' # relVal
        #'/store/user/apequegn/QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6/QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6_25Jun15-v1/150625_204806/0000/patTuple_10.root' 
        #'/store/mc/RunIISpring15DR74/QCD_Pt-15TTo7000_TuneZ2star-Flat_13TeV_pythia6/MINIAODSIM/Asympt50ns_MCRUN2_74_V9A-v1/00000/0AAC75BF-39FE-E411-BD62-0025905A60A0.root' 
        #'/store/mc/RunIISpring15DR74/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/Asympt50nsRaw_MCRUN2_74_V9A-v3/50000/CC91FB27-F308-E511-A212-00266CFBE43C.root' 
        '/store/mc/RunIISpring15DR74/QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/MINIAODSIM/Asympt25nsRaw_MCRUN2_74_V9-v3/00000/08547C6F-CD08-E511-81F2-0025905A60AA.root' 
    )
process.maxEvents.input = 1000
