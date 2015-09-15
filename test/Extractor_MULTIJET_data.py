from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing()

options.register ('globalTag',
    '74X_dataRun2_Prompt_v0',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "The globaltag to use")

options.parseArguments()
if len(options.globalTag) == 0:
  raise Exception("You _must_ pass a globalTag options to this script. Use --help for more informations")

from Extractor_MULTIJET_common import *

process = createExtractorProcess(False, True, useShiftCorrectedMET = True, globalTag = options.globalTag)


process.source.fileNames = cms.untracked.vstring( 
        #'file:/gridgroup/cms/pequegnot/CMSSW/CMSSW_5_3_9_patch2/src/PatTopProduction/patTuple.root'
        #'/store/user/apequegn/JetHT/JetHT_7_4_X_RunD_25Jun15-v1/339c8c37e62a7df1890df9cb21e1059e/patTuple_11_3_QvH.root' # GT: PHYS14_25_V2
        #'/store/data/Run2015B/JetHT/MINIAOD/PromptReco-v1/000/251/164/00000/503874CF-A826-E511-A22B-02163E011A34.root'
        '/store/data/Run2015C/JetHT/MINIAOD/PromptReco-v1/000/253/620/00000/D4BD3FF3-1D40-E511-8375-02163E0142EE.root'
)

process.maxEvents.input = 10

