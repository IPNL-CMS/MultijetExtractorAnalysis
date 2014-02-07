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
  process.load('Configuration/StandardSequences/GeometryIdeal_cff')
  process.load('Configuration/StandardSequences/MagneticField_38T_cff')
  process.load('Configuration/StandardSequences/EndOfProcess_cff')
  process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
  process.load("FWCore.MessageLogger.MessageLogger_cfi")
  process.load("Extractors.PatExtractor.PAT_extractor_cff")

  process.maxEvents = cms.untracked.PSet(
      input = cms.untracked.int32(10) #
      )

  #Global tag and data type choice
  process.GlobalTag.globaltag = '%s::All' % globalTag
  process.PATextraction.isMC  = isMC
  process.PATextraction.doMC  = isMC

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

  process.PATextraction.doMuon     = True
  process.PATextraction.doElectron = True
  process.PATextraction.doJet      = True
  process.PATextraction.doPhoton   = True

  process.PATextraction.doMET      = True
  if useShiftCorrectedMET:
    process.PATextraction.MET_PF.input = cms.InputTag("patMETsShiftCorrPFlow")
  else:
    process.PATextraction.MET_PF.input  = cms.InputTag("patMETsPFlow")

  process.PATextraction.doVertex   = True
  process.PATextraction.vtx_tag    = cms.InputTag( "goodOfflinePrimaryVertices" )
  process.PATextraction.doHLT      = True

  if not isMC:
    process.PATextraction.triggersXML = readFile("triggers.xml")

  # Jets correction : needs a valid global tags, or an external DB where JEC are stored
  process.PATextraction.jet_PF.redoJetCorrection = True

  if isMC:
    process.PATextraction.jet_PF.jetCorrectorLabel = "ak5PFchsL1FastL2L3"
  else:
    process.PATextraction.jet_PF.jetCorrectorLabel = "ak5PFchsL1FastL2L3Residual"

  process.PATextraction.jet_PF.doJER = True # Disable automatically on data
  process.PATextraction.jet_PF.doLooseJetID = False 

  # JER systematics:
  # Use -1 for 1-sigma down, 0 for nominal correction, and 1 for 1-sigma up
  process.PATextraction.jet_PF.jerSign = 0

  # JES systematics:
  # Use -1 for 1-sigma down, 0 for nominal correction, and 1 for 1-sigma up
  process.PATextraction.jet_PF.jesSign = 0
  process.PATextraction.jet_PF.jes_uncertainties_file = cms.untracked.string("")

  process.PATextraction.MET_PF.redoMetPhiCorrection   = True
  process.PATextraction.MET_PF.redoMetTypeICorrection = False # Automatically true if redoJetCorrection is True


  # Multijet analysis configuration
  process.PATextraction.plugins = cms.PSet(
    multijetExtractorAnalysis = cms.PSet(
      PUJets = cms.PSet(
	removePUJets = cms.bool(True)
	),
    
      firstJet = cms.PSet(
	eta_max = cms.double(1.3)
	),
	
      jet = cms.PSet(
	pt_min = cms.double(25),
	eta_max = cms.double(2.8),
	number_min = cms.double(3)
	),

      recoilJets = cms.PSet(
	pt_min = cms.double(25.),
	eta_max = cms.double(5.0),
	),
	
      vertex = cms.PSet(
	#number_min = cms.double(1),
	tracks_min = cms.double(5)
	),
	
      secondJet = cms.PSet(
	pt_max = cms.double(750),
	a_max = cms.double(0.6)
	),
	
      alpha = cms.PSet(
	alpha_max = cms.double(0.3)
	),
	
      beta = cms.PSet(
	beta_min = cms.double(1)
	),
	
      muons = cms.PSet(
	isolation_max = cms.double(0.20)
	),
	
      electrons = cms.PSet(
	isolation_max = cms.double(0.15)
	),
	
      photons = cms.PSet(
	isolation_max = cms.double(0.15)
	)
	
      ) 

      )



  #########################################
  #
  # Launch the job
  #
  #########################################


  process.p = cms.Path(process.PATextraction)
  process.MessageLogger.cerr.FwkReport.reportEvery = 1000

  return process
