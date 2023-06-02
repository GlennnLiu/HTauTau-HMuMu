
# DATA_TAG = "ULAPV" # Change to PromptReco for Run2016 period H
LEPTON_SETUP = 2016  # current default = 2017 = Moriond2017
APPLYMUCORR = True  # Switch off muon scale corrections
APPLYJEC = True     #
APPLYJER = True     #
RECORRECTMET = True #

PROCESS_CR = True
PROCESS_CRTT = True
PROCESS_CRWJ = True

#DOMETRECOIL = True

#APPLY_QCD_GGF_UNCERT = True
#For MC:
#PD = ""
#MCFILTER = ""
#IsMC = True
#OLDPATTRIGGER = False#True
#For DATA: 
IsMC = False
PD = "SingleMuon"

# Get absolute path
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/HTauTauHMuMu/AnalysisStep/test/"

### ----------------------------------------------------------------------
### Standard sequence
### ----------------------------------------------------------------------

execfile(PyFilePath + "analyzer.py")
#execfile(PyFilePath + "prod/pyFragments/RecoProbabilities.py")

if not IsMC:
	process.source.inputCommands = cms.untracked.vstring("keep *", "drop LHERunInfoProduct_*_*_*", "drop LHEEventProduct_*_*_*") ###FIXME In 9X this removes all collections for MC

### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------

process.source.fileNames = cms.untracked.vstring(
#"/store/mc/RunIISummer20UL16MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v17-v1/2520000/04A698D5-2AF9-B548-9A6D-DB5AFE92F0A6.root"
#"/store/mc/RunIISummer20UL16MiniAODv2/GluGluHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v17-v2/2520000/005D2E78-273E-D241-BF7B-295CEB6B444D.root"
#"/store/mc/RunIISummer20UL16MiniAODv2/WJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v17-v2/2520000/26778B66-EF37-4F45-9261-55914F5FAE44.root"
#'/store/mc/RunIISummer20UL16MiniAODv2/GluGluHToWWTo2L2Nu_M125_TuneCP5_13TeV_powheg2_JHUGenV714_pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v17-v2/260000/F24F2632-D967-8F46-A5F5-7AB27188625C.root'
#'/store/mc/RunIISummer20UL16MiniAODv2/VBFHToWWTo2L2Nu_M125_TuneCP5_13TeV_powheg2_JHUGenV714_pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v17-v2/2530000/3B2B31D9-593B-F845-9412-E7B6E3A23C63.root'
#'/store/mc/RunIISummer20UL16MiniAODv2/GluGluHToMuMu_M125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v17-v2/70000/009EBD07-8910-0F4C-A962-70F3E7E260DB.root'
#"/store/mc/RunIISummer20UL16MiniAODv2/EWKZ2Jets_ZToLL_M-50_TuneCP5_withDipoleRecoil_13TeV-madgraph-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v17-v1/70000/1CB2BD9A-308E-E544-B69C-627458DDA424.root"
#"/store/mc/RunIISummer20UL16MiniAODv2/VBFHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v17-v2/230000/00C065D5-4429-EB4B-9F43-28AE1A146705.root"
#"/store/mc/RunIISummer20UL16MiniAODv2/VBFHToMuMu_M125_TuneCP5_withDipoleRecoil_13TeV-powheg-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v17-v1/70000/19F81AB8-D7F9-9F4D-8EE5-33E37F67E10D.root"
#"/store/mc/RunIISummer20UL16MiniAODv2/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v17-v1/120000/038D9F5E-C340-FD45-B613-DD495181C237.root"
'/store/data/Run2016H/SingleMuon/MINIAOD/UL2016_MiniAODv2-v2/120000/DA9B9166-C273-E54C-8957-7E54DACC6999.root'
### APV ###
# "/store/mc/RunIISummer20UL16MiniAODAPVv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v11-v1/120000/008E4139-6019-CE4C-B83C-A849F56F57B3.root"
)

#process.calibratedPatElectrons.isSynchronization = cms.bool(True) #process.calibratedPatElectrons.isSynchronization = cms.bool(True) # Not needed anymore since new EGamma smearing is event deterministic
#process.calibratedMuons.isSynchronization = cms.bool(True)

process.maxEvents.input = 10000
#process.source.skipEvents = cms.untracked.uint32(5750)

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


### ----------------------------------------------------------------------
### Analyzer for Plots
### ----------------------------------------------------------------------


#process.source.eventsToProcess = cms.untracked.VEventRange("1:8670")

# Debug
process.dumpUserData =  cms.EDAnalyzer("dumpUserData",
     dumpTrigger = cms.untracked.bool(True),
     muonSrcs = cms.PSet(
        muons = cms.InputTag("slimmedMuons"),
#       slimmedMuons = cms.InputTag("slimmedMuons"),
#        muons = cms.InputTag("appendPhotons:muons"),
     ),
     electronSrcs = cms.PSet(
#       slimmedElectron = cms.InputTag("slimmedElectrons"),
        electrons = cms.InputTag("appendPhotons:electrons"),
#        RSE = cms.InputTag("appendPhotons:looseElectrons"),
#        TLE = cms.InputTag("appendPhotons:electronstle"), #These are actually photons, should add a photonSrcs section for them.
     ),
     candidateSrcs = cms.PSet(
        Z     = cms.InputTag("ZCand"),
#        ZRSE     = cms.InputTag("ZCandlooseEle"),
#        ZTLE     = cms.InputTag("ZCandtle"),
#        ZZ  = cms.InputTag("ZZCand"),
#        ZZRSE     = cms.InputTag("ZZCandlooseEle"),
#        ZZTLE     = cms.InputTag("ZZCandtle"),
#        ZLL  = cms.InputTag("ZLLCand"),
#        ZL  = cms.InputTag("ZlCand"),
     ),
     jetSrc = cms.InputTag("cleanJets"),
)

# Create lepton sync file
#process.PlotsZZ.dumpForSync = True;
#process.p = cms.EndPath( process.PlotsZZ)

# Keep all events in the tree, even if no candidate is selected
#process.ZZTree.skipEmptyEvents = False

# replace the paths in analyzer.py
#process.trees = cms.EndPath(process.ZZTree)

#Dump reconstructed variables
#process.appendPhotons.debug = cms.untracked.bool(True)
#process.fsrPhotons.debug = cms.untracked.bool(True)
#process.dump = cms.Path(process.dumpUserData)

#Print MC history
#process.mch = cms.EndPath(process.printTree)


#Monitor memory usage
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1)
#)
