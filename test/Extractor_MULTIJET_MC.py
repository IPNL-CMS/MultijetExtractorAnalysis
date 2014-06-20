from Extractor_MULTIJET_common import *

#with open('../../../PatTopProduction/multijet_QCD_HT_fewEvents.list') as f:
#  files = f.readlines()

process = createExtractorProcess(True, True, useShiftCorrectedMET = True, globalTag = "START53_V27")

#process.source.fileNames = cms.untracked.vstring(files)
process.source.fileNames = cms.untracked.vstring( '/store/user/apequegn/QCD_HT-100To250_TuneZ2star_8TeV-madgraph-pythia/QCD_HT-100To250_START53_V7A_13Mar14-v1/c5f9c59e100f883a59cec8d8908af608/patTuple_976_1_4Lv.root'
    )
process.maxEvents.input = -1
