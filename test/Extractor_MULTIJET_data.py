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


process.source.fileNames = cms.untracked.vstring( '/store/user/apequegn/Jet/Jet_Run2012A-22Jan2013_04Dec13-v1/2b71cc75519c435218af9c37f5843477/patTuple_92_1_Syy.root'
)

process.maxEvents.input = -1

