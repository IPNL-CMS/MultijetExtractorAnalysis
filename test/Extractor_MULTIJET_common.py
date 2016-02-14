#########################################
#
# Base macro for launching the PatExtractor
#
# The macro is for tests
#
#########################################


import FWCore.ParameterSet.Config as cms

def readFile(file):
  return cms.untracked.string(open(file).read())

def createExtractorProcess(isMC, isSemiMu, useShiftCorrectedMET, globalTag):
  process = cms.Process("PATextractor2")


  #########################################
  #
  # Main configuration statements
  #
  #########################################

  process.load('Configuration/StandardSequences/Services_cff')
  #process.load('Configuration/StandardSequences/GeometryIdeal_cff')
  #process.load('Configuration/StandardSequences/MagneticField_38T_cff')
  process.load('Configuration/StandardSequences/EndOfProcess_cff')
  #process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
  process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
  process.load('Configuration.StandardSequences.MagneticField_38T_cff')
  process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
  #process.load("Configuration.EventContent.EventContent_cff")
  process.load("FWCore.MessageLogger.MessageLogger_cfi")
  #process.load("Extractors.PatExtractor.PAT_extractor_lyonPatTuples_cff")
  process.load("Extractors.PatExtractor.PAT_extractor_miniAOD_cff")

  process.maxEvents = cms.untracked.PSet(
      input = cms.untracked.int32(10) #
      )

  #Global tag and data type choice
  #process.GlobalTag.globaltag = '%s::All' % globalTag
  process.GlobalTag.globaltag = '%s' % globalTag
  process.PATextraction.isMC  = isMC
  process.PATextraction.extractors.MC.enable  = isMC

  #Input PAT file to extract
  process.source = cms.Source("PoolSource",
      fileNames = cms.untracked.vstring(
        ),                           
      duplicateCheckMode = cms.untracked.string( 'noDuplicateCheck' )
      )

  #Output extracted file name
  if isMC:
    process.PATextraction.extractedRootFile = cms.string('extracted_mc.root')
  else:
    process.PATextraction.extractedRootFile = cms.string('extracted.root')

  #########################################
  #
  # PAT extractor main options statements
  #
  #########################################

  #
  # Adapt it to your needs
  #
  # If you are lost, see the example here (PART 3.2):
  # http://sviret.web.cern.ch/sviret/Welcome.php?n=CMS.PHYTuto
  #
  # Here we just extract, and don't perform any analysis

  process.PATextraction.extractors.photon.enable = True
  process.PATextraction.extractors.HLT.enable = True

  if not isMC:
    process.PATextraction.extractors.HLT.parameters.triggers = readFile("triggers.xml")

  # Jets correction : needs a valid global tags, or an external DB where JEC are stored
  process.PATextraction.extractors.jetmet.parameters.redoJetCorrection = True

  if isMC:
    process.PATextraction.extractors.jetmet.parameters.jetCorrectorLabel = "ak4PFCHSL1FastL2L3"
  else:
    # For closure:
    #process.PATextraction.extractors.jetmet.parameters.jetCorrectorLabel = "ak4PFCHSL1FastL2L3Residual"
    # For L3Res extraction:
    process.PATextraction.extractors.jetmet.parameters.jetCorrectorLabel = "ak4PFCHSL1FastL2L3"

  process.PATextraction.extractors.jetmet.parameters.doJER = True # Disabled automatically on data
  process.PATextraction.extractors.jetmet.parameters.doLooseJetID = True
  #process.PATextraction.extractors.jetmet.parameters.useGlobalTagForJEC = True
  process.PATextraction.extractors.jetmet.parameters.useGlobalTagForJEC = False

  # only L2Res on data (exclude L3Res)
  # Fall15_25nsV1_DATA_L2Residual_AK4PFchs.txt is a coy of Fall15_25ns_COMB_LOGLIN_L2Residual_v2_AK4PFchs_nokFSR.txt Fall15_25nsV1_DATA_L2Residual_AK4PFchs.txt
  process.PATextraction.extractors.jetmet.parameters.jecPayload = "Extractors/PatExtractor/data/jec_payloads_76X_Fall15_25nsV1_L2ResOnlyOnData.xml"
  process.PATextraction.extractors.jetmet.parameters.jecJetAlgo = "AK4PFchs"
  process.PATextraction.extractors.jetmet.parameters.useType1Fix = True
  #process.PATextraction.extractors.jetmet.parameters.useType1Fix = False
  #process.PATextraction.extractors.jetmet.parameters.useGlobalTagForType1Fix = True
  process.PATextraction.extractors.jetmet.parameters.useGlobalTagForType1Fix = False
  process.PATextraction.extractors.jetmet.parameters.jetCorrectorLabelForType1Fix = "ak4PFCHSL1RC"
  process.PATextraction.extractors.jetmet.parameters.jecPayload_L1ForType1Fix = "Extractors/PatExtractor/data/jec_payloads_L1ForType1Fix_Fall15_25nsV1.xml"

  # JER systematics:
  # Use -1 for 1-sigma down, 0 for nominal correction, and 1 for 1-sigma up
  process.PATextraction.extractors.jetmet.parameters.jerSign = 0

  # JES systematics:
  # Use -1 for 1-sigma down, 0 for nominal correction, and 1 for 1-sigma up
  process.PATextraction.extractors.jetmet.parameters.jesSign = 0
  process.PATextraction.extractors.jetmet.parameters.jes_uncertainties_file = cms.untracked.string("")

  process.PATextraction.extractors.jetmet.parameters.redoMetPhiCorrection   = True
  process.PATextraction.extractors.jetmet.parameters.redoMetTypeICorrection = True # Automatically true if redoJetCorrection is True


  # Multijet analysis configuration
  process.PATextraction.plugins = cms.PSet(
    multijetExtractorAnalysis = cms.PSet(
      PUJets = cms.PSet(
	removePUJets = cms.bool(False),
    # ID flag for PU jet identification
    # 4: jet is not PU with a Loose CL
    # 6: jet is not PU with a Medium CL 
    # 7: jet is not PU with a Tight CL
    # 0: jet is PU (?)
    id_min = cms.double(7)
	),

      recoilJets = cms.PSet(
	pt_min = cms.double(30.),
	eta_max = cms.double(3.0),
	),
	
      vertex = cms.PSet(
	#number_min = cms.double(1),
	tracks_min = cms.double(5)
	),

      muons = cms.PSet( # see https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonIsolationForRun2
	isolation_max = cms.double(0.40)
	)
	
      ) 

      )



  #########################################
  #
  # Launch the job
  #
  #########################################
  from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
  # turn on VID producer, indicate data format  to be
  # DataFormat.AOD or DataFormat.MiniAOD, as appropriate

  # Set up input/output depending on the format
  # You can list here either AOD or miniAOD files, but not both types mixed
  #
  useAOD = False

  if useAOD == True :
      dataFormat = DataFormat.AOD
  else :
      dataFormat = DataFormat.MiniAOD
 
  switchOnVIDPhotonIdProducer(process, dataFormat)
 
  # define which IDs we want to produce
  my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring15_25ns_V1_cff']

  # running on pattuple 
#  process.egmPhotonIDs.physicsObjectSrc = cms.InputTag("slimmedPhotonsPFlow")
  #process.photonIDValueMapProducer.ebReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma:reducedEBRecHits")
  #process.photonIDValueMapProducer.eeReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma:reducedEERecHits")
  #process.photonIDValueMapProducer.esReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma:reducedESRecHits")
  #process.photonIDValueMapProducer.verticesMiniAOD = cms.InputTag("offlineSlimmedPrimaryVertices")
  #process.photonIDValueMapProducer.pfCandidatesMiniAOD = cms.InputTag("packedPFCandidatesPFlow")
  #process.photonIDValueMapProducer.srcMiniAOD = cms.InputTag("slimmedPhotonsPFlow")

  # running on miniAOD
  process.egmPhotonIDs.physicsObjectSrc = cms.InputTag("slimmedPhotons")
  process.photonIDValueMapProducer.ebReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma:reducedEBRecHits")
  process.photonIDValueMapProducer.eeReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma:reducedEERecHits")
  process.photonIDValueMapProducer.esReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma:reducedESRecHits")
  process.photonIDValueMapProducer.verticesMiniAOD = cms.InputTag("offlineSlimmedPrimaryVertices")
  process.photonIDValueMapProducer.pfCandidatesMiniAOD = cms.InputTag("packedPFCandidates")
  process.photonIDValueMapProducer.srcMiniAOD = cms.InputTag("slimmedPhotons")

  #add them to the VID producer
  for idmod in my_id_modules:
      setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)
  

  switchOnVIDElectronIdProducer(process, dataFormat)

  id_modules = [
          'RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff'
          ]

  for mod in id_modules:
      setupAllVIDIdsInModule(process, mod, setupVIDElectronSelection)

  process.p = cms.Path(process.egmPhotonIDs+process.egmGsfElectronIDSequence+process.PATextraction)
  process.MessageLogger.cerr.FwkReport.reportEvery = 1000

  process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False), allowUnscheduled = cms.untracked.bool(True))
  return process
