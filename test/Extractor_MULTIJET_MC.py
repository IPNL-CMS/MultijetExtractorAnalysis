from Extractor_MULTIJET_common import *

with open('../../../PatTopProduction/multijet_QCD_HT_fewEvents.list') as f:
  files = f.readlines()

process = createExtractorProcess(True, True, useShiftCorrectedMET = True, globalTag = "START53_V27")

process.source.fileNames = cms.untracked.vstring(files)
#process.source.fileNames = cms.untracked.vstring( 'file:/gridgroup/cms/pequegnot/CMSSW/CMSSW_5_3_9_patch2/src/PatTopProduction/patTuple.root'
    #)
process.maxEvents.input = -1
