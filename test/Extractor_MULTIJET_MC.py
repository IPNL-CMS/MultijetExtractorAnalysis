from Extractor_MULTIJET_common import *

process = createExtractorProcess(True, True, useShiftCorrectedMET = True, globalTag = "START53_V23")

process.source.fileNames = cms.untracked.vstring( 'file:/gridgroup/cms/pequegnot/CMSSW/CMSSW_5_3_9_patch2/src/PatTopProduction/patTuple.root'
    )
process.maxEvents.input = 5000
