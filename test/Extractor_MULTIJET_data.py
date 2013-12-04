from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing()

options.register ('globalTag',
    'START53_V23',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "The globaltag to use")

options.parseArguments()
if len(options.globalTag) == 0:
  raise Exception("You _must_ pass a globalTag options to this script. Use --help for more informations")

from Extractor_MULTIJET_common import *

process = createExtractorProcess(False, True, useShiftCorrectedMET = True, globalTag = options.globalTag)


process.source.fileNames = cms.untracked.vstring( 'file:/gridgroup/cms/pequegnot/CMSSW/CMSSW_5_3_9_patch2/src/PatTopProduction/patTupleData.root'
)

process.maxEvents.input = -1

