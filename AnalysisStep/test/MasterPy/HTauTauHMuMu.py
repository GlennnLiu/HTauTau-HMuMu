import FWCore.ParameterSet.Config as cms
from HTauTauHMuMu.AnalysisStep.defaults import *
import os, sys

process = cms.Process("HttHmm")

### ----------------------------------------------------------------------
### Flags that need to be set
### ----------------------------------------------------------------------

#Set defaults for variables used in this file (in case they are not defined by a caller script)
declareDefault("IsMC", True, globals())

# Set of effective areas, rho corrections, etc. (can be 2011, 2012, 2015 or 2016)
declareDefault("LEPTON_SETUP", 2018, globals())

# Flag that reflects the actual sqrts of the sample (can be 2011, 2012, 2015 or 2016)
# Can differ from SAMPLE_TYPE for samples that are rescaled to a different sqrts.
declareDefault("SAMPLE_TYPE", LEPTON_SETUP, globals())

declareDefault("YEAR", LEPTON_SETUP, globals())

# Control global tag to be used for 2018 data to distinguish between ReReco (period A, B, C) and PromptReco (period D)
declareDefault("DATA_TAG", "ReReco", globals())

#Optional name of the sample/dataset being analyzed
declareDefault("SAMPLENAME", "", globals())

#Apply muon scale correction
declareDefault("APPLYMUCORR", True, globals())

#Reapply JEC
declareDefault("APPLYJEC", True, globals())

#Apply JER
declareDefault("APPLYJER", True, globals())

#Recorrect MET
declareDefault("RECORRECTMET", True, globals())

#FSR mode 
declareDefault("FSRMODE", "RunII", globals())

#Bunch spacing (can be 25 or 50)
declareDefault("BUNCH_SPACING", 25, globals())

#Best candidate comparator (see interface/Comparators.h)
declareDefault("BESTCANDCOMPARATOR", "byBestZ1bestZ2", globals())

# Set to True to make candidates with the full combinatorial of loose leptons (for debug; much slower)
declareDefault("KEEPLOOSECOMB", False, globals())

# Activate paths for loose electron categories
declareDefault("ADDLOOSEELE", False, globals())

# Activate trigger paths in MC; note that for 2016, only reHLT samples have the correct triggers!!!
declareDefault("APPLYTRIG", True, globals())

# Set to True to re-activate the now-deprecated PATMuonCleanerBySegments
UseMuonCleanerBySegments = False 

# CMSSW version
CMSSW_VERSION = os.environ['CMSSW_VERSION']
CMSSWVERSION = int(CMSSW_VERSION.split("_")[1])

# if SELSETUP=="Legacy" and not BESTCANDCOMPARATOR=="byBestZ1bestZ2":
#     print "WARNING: In ZZ4lAnalysis.py the SELSETUP=\"Legacy\" flag is meant to reproduce the Legacy results, ignoring the setting of the BESTCANDCOMPARATOR: ",BESTCANDCOMPARATOR
#     BESTCANDCOMPARATOR = "byBestZ1bestZ2"

# The isolation cuts for electrons and muons. FIXME: there is an hardcoded instance of these values in src/LeptonIsoHelper.cc !!
ELEISOCUT = 0.15 # [FIXME] Remove isolation cuts from the code completely
MUISOCUT  = 0.15 # [FIXME] Remove isolation cuts from the code completely

### ----------------------------------------------------------------------
### Set the GT
### ----------------------------------------------------------------------
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag

if (SAMPLE_TYPE == 2016):
    if IsMC:
        if (DATA_TAG == "ULAPV"):
            # preVFP samples have different GT
            # Use DATA_TAG (include in MC csv files) to distinguish between pre/post VFP
            process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun2_asymptotic_preVFP_v11', '')
        else:
            process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun2_asymptotic_v17', '')
    else:
        process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v35', '')

elif (SAMPLE_TYPE == 2017):
    if IsMC:
        process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mc2017_realistic_v9', '')
    else:
        process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v35', '')

elif (SAMPLE_TYPE == 2018):
    if IsMC:
        process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v16_L1v1', '')
    else:
        if (DATA_TAG == "PromptReco"):
            #In UL probably not needed anymore. Leaving it here for consistency with ReReco
            #This is not changed wrt to ReReco, probably for 2018 there is not anymore the difference between ReReco and PromptReco
            process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v35', '')
        else:
            process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v35', '')

print '\t',process.GlobalTag.globaltag


### ----------------------------------------------------------------------
### Standard stuff
### ----------------------------------------------------------------------
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


### ----------------------------------------------------------------------
### Source
### ----------------------------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      '/store/cmst3/user/cmgtools/CMG/GluGluToHToZZTo4L_M-130_7TeV-powheg-pythia6/Fall11-PU_S6_START42_V14B-v1/AODSIM/V5/PAT_CMG_V5_2_0/patTuple_1.root'
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


### ----------------------------------------------------------------------
### Trigger bit Requests
### ----------------------------------------------------------------------
import HLTrigger.HLTfilters.hltHighLevel_cfi

process.hltFilterSingleMu  = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilterSingleEle = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilterDiMu  = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilterDiTau = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilterMuTau = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilterEleTau = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilterMuEle = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()

process.hltFilterSingleMu.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilterSingleEle.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilterDiMu.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilterDiTau.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilterMuTau.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilterEleTau.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilterMuEle.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")

process.hltFilterSingleMu.throw = cms.bool(False) #FIXME: beware of this!
process.hltFilterSingleEle.throw = cms.bool(False) #FIXME: beware of this!
process.hltFilterDiMu.throw = cms.bool(False) #FIXME: beware of this!
process.hltFilterDiTau.throw = cms.bool(False) #FIXME: beware of this!
process.hltFilterMuTau.throw = cms.bool(False) #FIXME: beware of this!
process.hltFilterEleTau.throw = cms.bool(False) #FIXME: beware of this!
process.hltFilterMuEle.throw = cms.bool(False) #FIXME: beware of this!

### 2016 triggers - final
if (LEPTON_SETUP == 2016):
    process.hltFilterSingleMu.HLTPaths = ["HLT_IsoMu22_v*","HLT_IsoMu22_eta2p1_v*","HLT_IsoTkMu22_v*","HLT_IsoTkMu22_eta2p1_v*","HLT_IsoMu24_v*","HLT_IsoTkMu24_v*"]
    process.hltFilterSingleEle.HLTPaths = ["HLT_Ele25_eta2p1_WPTight_Gsf_v*","HLT_Ele27_WPTight_Gsf_v*","HLT_Ele27_eta2p1_WPLoose_Gsf_v*", "HLT_Ele32_eta2p1_WPTight_Gsf_v*"]
    process.hltFilterDiMu.HLTPaths = ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*","HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*"]
    process.hltFilterDiTau.HLTPaths = ["HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg_v*","HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg_v*"]
    process.hltFilterMuTau.HLTPaths = ["HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v*","HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v*"]
    process.hltFilterEleTau.HLTPaths = ["HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_SingleL1_v*","HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau20_v*","HLT_Ele24_eta2p1_WPLoose_Gsf_LooseIsoPFTau30_v*"]
    process.hltFilterMuEle.HLTPaths = ["HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v*","HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*"]
    
    process.triggerSingleMu  = cms.Path(process.hltFilterSingleMu)
    process.triggerSingleEle = cms.Path(process.hltFilterSingleEle)
    process.triggerDiMu = cms.Path(process.hltFilterDiMu)
    process.triggerDiTau = cms.Path(process.hltFilterDiTau)
    process.triggerMuTau = cms.Path(process.hltFilterMuTau)
    process.triggerEleTau = cms.Path(process.hltFilterEleTau)
    process.triggerMuEle = cms.Path(process.hltFilterMuEle)

### 2017 triggers - final
elif (LEPTON_SETUP == 2017):
    process.hltFilterSingleMu.HLTPaths = ["HLT_IsoMu24_v*","HLT_IsoMu27_v*"]
    process.hltFilterSingleEle.HLTPaths = ["HLT_Ele32_WPTight_Gsf_v*","HLT_Ele35_WPTight_Gsf_v*"]
    process.hltFilterDiMu.HLTPaths = ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*","HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8_v*"]
    process.hltFilterDiTau.HLTPaths = ["HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v*","HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v*","HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v*"]
    process.hltFilterMuTau.HLTPaths = ["HLT_IsoMu20_eta2p1_LooseChargedIsoPFTau27_eta2p1_CrossL1_v*"]
    process.hltFilterEleTau.HLTPaths = ["HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTau30_eta2p1_CrossL1_v*"]
    process.hltFilterMuEle.HLTPaths = ["HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*"]
   
    process.triggerSingleMu  = cms.Path(process.hltFilterSingleMu)
    process.triggerSingleEle = cms.Path(process.hltFilterSingleEle)
    process.triggerDiMu = cms.Path(process.hltFilterDiMu)
    process.triggerDiTau = cms.Path(process.hltFilterDiTau)
    process.triggerMuTau = cms.Path(process.hltFilterMuTau)
    process.triggerEleTau = cms.Path(process.hltFilterEleTau)
    process.triggerMuEle = cms.Path(process.hltFilterMuEle)

### 2018 triggers - FIXME: to be updated (26/6/18)
elif (LEPTON_SETUP == 2018):
    process.hltFilterSingleMu.HLTPaths = ["HLT_IsoMu24_v*","HLT_IsoMu27_v*"]
    process.hltFilterSingleEle.HLTPaths = ["HLT_Ele32_WPTight_Gsf_v*","HLT_Ele35_WPTight_Gsf_v*"]
    process.hltFilterDiMu.HLTPaths = ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v*"]
    process.hltFilterDiTau.HLTPaths = ["HLT_DoubleMediumChargedIsoPFTau40_Trk1_TightID_eta2p1_Reg_v*","HLT_DoubleTightChargedIsoPFTau40_Trk1_eta2p1_Reg_v*","HLT_DoubleTightChargedIsoPFTau35_Trk1_TightID_eta2p1_Reg_v*","HLT_DoubleTightChargedIsoPFTauHPS40_Trk1_eta2p1_Reg_v*"]
    process.hltFilterMuTau.HLTPaths = ["HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v*","HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_TightID_CrossL1_v*"]
    process.hltFilterEleTau.HLTPaths = ["HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_v*","HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_TightID_CrossL1_v*"]
    process.hltFilterMuEle.HLTPaths = ["HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*","HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*"]
   
    process.triggerSingleMu  = cms.Path(process.hltFilterSingleMu)
    process.triggerSingleEle = cms.Path(process.hltFilterSingleEle)
    process.triggerDiMu = cms.Path(process.hltFilterDiMu)
    process.triggerDiTau = cms.Path(process.hltFilterDiTau)
    process.triggerMuTau = cms.Path(process.hltFilterMuTau)
    process.triggerEleTau = cms.Path(process.hltFilterEleTau)
    process.triggerMuEle = cms.Path(process.hltFilterMuEle)

TRIGGERSET="slimmedPatTrigger"

### ----------------------------------------------------------------------
### MET FILTERS
### ----------------------------------------------------------------------
# process.METFilters  = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
# process.METFilters.TriggerResultsTag  = cms.InputTag("TriggerResults","","RECO")
# if (IsMC):
#   process.METFilters.TriggerResultsTag  = cms.InputTag("TriggerResults","","PAT")

# if (LEPTON_SETUP == 2017):#MET Filters available in miniAOD as described here https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
#   if (IsMC):
#      process.METFilters.HLTPaths = ["Flag_goodVertices","Flag_globalSuperTightHalo2016Filter","Flag_HBHENoiseFilter","Flag_HBHENoiseIsoFilter","Flag_EcalDeadCellTriggerPrimitiveFilter","Flag_BadPFMuonFilter","Flag_BadChargedCandidateFilter"]
#   else:
#      process.METFilters.HLTPaths = ["Flag_goodVertices","Flag_globalSuperTightHalo2016Filter","Flag_HBHENoiseFilter","Flag_HBHENoiseIsoFilter","Flag_EcalDeadCellTriggerPrimitiveFilter","Flag_BadPFMuonFilter","Flag_BadChargedCandidateFilter","Flag_eeBadScFilter"]

# process.triggerMETFilters = cms.Path(process.METFilters)

### ----------------------------------------------------------------------
### MC Filters and tools
### ----------------------------------------------------------------------

process.heavyflavorfilter = cms.EDFilter('HeavyFlavorFilter2',
#                                 src= cms.InputTag("genParticles"), # genParticles available only in PAT
                                 src= cms.InputTag("prunedGenParticles"),
                                 status2 = cms.bool(True),
                                 status3 = cms.bool(False),
                                 hDaughterVeto = cms.bool(False),
                                 zDaughterVeto = cms.bool(True),
                                 ptcut=cms.double(0)
                                 )


process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.drawTree = cms.EDAnalyzer("ParticleTreeDrawer",
                                   src = cms.InputTag("prunedGenParticles"),
                                   printP4 = cms.untracked.bool(False),
                                   printPtEtaPhi = cms.untracked.bool(False),
                                   printVertex = cms.untracked.bool(False),
                                   printStatus = cms.untracked.bool(True),
                                   printIndex = cms.untracked.bool(False) )


process.printTree = cms.EDAnalyzer("ParticleListDrawer",
                                   maxEventsToPrint = cms.untracked.int32(-1),
                                   printVertex = cms.untracked.bool(False),
                                   src = cms.InputTag("prunedGenParticles")
                                   )


process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   calibratedPatElectrons = cms.PSet(
                                                       initialSeed = cms.untracked.uint32(1),
                                                       engineName = cms.untracked.string('TRandom3')
                                                       ),
                                                   )

# FIXME Add total kinematics filter for MC

process.goodPrimaryVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag("offlineSlimmedPrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
  filter = cms.bool(True),
)

### ----------------------------------------------------------------------
### HTXS categorisation
### ----------------------------------------------------------------------
if(IsMC):
    process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
    process.mergedGenParticles = cms.EDProducer("MergedGenParticleProducer",
                                               inputPruned = cms.InputTag("prunedGenParticles"),
                                               inputPacked = cms.InputTag("packedGenParticles"),
                                              )
    process.myGenerator = cms.EDProducer("GenParticles2HepMCConverter",
                                        genParticles = cms.InputTag("mergedGenParticles"),
                                        genEventInfo = cms.InputTag("generator"),
                                        signalParticlePdgIds = cms.vint32(25), ## for the Higgs analysis
                                       )
    process.rivetProducerHTXS = cms.EDProducer('HTXSRivetProducer',
                                              HepMCCollection = cms.InputTag('myGenerator','unsmeared'),
                                              LHERunInfo = cms.InputTag('externalLHEProducer'),
                                              ProductionMode = cms.string('AUTO'),
                                             )
    process.htxs = cms.Path(process.mergedGenParticles*process.myGenerator*process.rivetProducerHTXS)

### ----------------------------------------------------------------------
### ----------------------------------------------------------------------
### Loose lepton selection + cleaning + embedding of user data
### ----------------------------------------------------------------------
### ----------------------------------------------------------------------

DXY_DZ = "abs(dB('PV2D'))<0.045 && abs(dB('PVDZ'))<0.2" #dxy, dz cuts
SIP =  "abs(dB('PV3D')/edB('PV3D')) < 4"
GOODELECTRON = "userFloat('isEleNoIsoID90') && " + SIP
GOODMUON     = "userFloat('mediumID') && " + SIP

APPLYTESCORRECTION = True

#------- MUONS -------

#--- Set correct identifier for muon corrections
#--- UL Rochester from: https://gitlab.cern.ch/akhukhun/roccor/-/tree/Run2.v5
#--- Corresponding TWiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/RochcorMuon
if LEPTON_SETUP == 2016: # Rochester corrections for 2016 data
    if ( DATA_TAG == 'ULAPV' ):
        process.calibratedMuons = cms.EDProducer("RochesterPATMuonCorrector",
                                            src = cms.InputTag("slimmedMuons"),
                                            identifier = cms.string("RoccoR2016aUL"),
                                            isMC = cms.bool(IsMC),
                                            isSynchronization = cms.bool(False),
                                            )
    else:
        process.calibratedMuons = cms.EDProducer("RochesterPATMuonCorrector",
                                            src = cms.InputTag("slimmedMuons"),
                                            identifier = cms.string("RoccoR2016bUL"),
                                            isMC = cms.bool(IsMC),
                                            isSynchronization = cms.bool(False),
                                            )
elif LEPTON_SETUP == 2017:# Rochester corrections for 2017 data
     process.calibratedMuons = cms.EDProducer("RochesterPATMuonCorrector",
                                         src = cms.InputTag("slimmedMuons"),
                                         identifier = cms.string("RoccoR2017UL"),
                                         isMC = cms.bool(IsMC),
                                         isSynchronization = cms.bool(False),
                                         )
elif LEPTON_SETUP == 2018:# Rochester corrections for 2018 data
     process.calibratedMuons = cms.EDProducer("RochesterPATMuonCorrector",
                                         src = cms.InputTag("slimmedMuons"),
                                         identifier = cms.string("RoccoR2018UL"),
                                         isMC = cms.bool(IsMC),
                                         isSynchronization = cms.bool(False),
                                         )
else:
    if APPLYMUCORR:
        print "APPLYMUCORR not configured for LEPTON_SETUP =", LEPTON_SETUP
        sys.exit()



process.bareSoftMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("calibratedMuons"),
    cut = cms.string("pt>10 && abs(eta)<2.4 && (isGlobalMuon || (isTrackerMuon && numberOfMatchedStations>0))")
#    Lowering pT cuts
#    cut = cms.string("(isGlobalMuon || (isTrackerMuon && numberOfMatchedStations>0)) && pt>3 && p>3.5 && abs(eta)<2.4")
)


# MC matching. As the genParticles are no more available in cmg, we re-match with genParticlesPruned.
process.muonMatch = cms.EDProducer("MCMatcher", # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                   src     = cms.InputTag("softMuons"), # RECO objects to match
                                   matched = cms.InputTag("prunedGenParticles"),   # mc-truth particle collection
                                   mcPdgId     = cms.vint32(13), # one or more PDG ID (13 = muon); absolute values (see below)
                                   checkCharge = cms.bool(True), # True = require RECO and MC objects to have the same charge
                                   mcStatus = cms.vint32(1),     # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
                                   maxDeltaR = cms.double(0.5),  # Minimum deltaR for the match
                                   maxDPtRel = cms.double(0.5),  # Minimum deltaPt/Pt for the match
                                   resolveAmbiguities = cms.bool(True),     # Forbid two RECO objects to match to the same GEN object
                                   resolveByMatchQuality = cms.bool(False), # False = just match input in order; True = pick lowest deltaR pair first
                                   )

process.softMuons = cms.EDProducer("MuFiller",
    src = cms.InputTag("bareSoftMuons"),
    sampleType = cms.int32(SAMPLE_TYPE),
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
    cut = cms.string('pt>10 && abs(eta) < 2.5 &&'+ DXY_DZ), #dxy, dz cuts
    # TriggerResultsLabel = cms.InputTag('TriggerResults','','HLT'),
    # TriggerSet = cms.InputTag(TRIGGERSET),
    flags = cms.PSet(
        isSIP = cms.string(SIP),
        isGood = cms.string(GOODMUON),
        isIsoFSRUncorr  = cms.string("userFloat('combRelIsoPF')<" + str(MUISOCUT)),
    )
)


process.muons =  cms.Sequence(process.calibratedMuons + process.bareSoftMuons + process.softMuons)

if not APPLYMUCORR :
    process.muons.replace(process.calibratedMuons, None)
    process.bareSoftMuons.src = cms.InputTag("slimmedMuons")

    
#--- Derecated muon cleaner; keep this option for future reference. 
if UseMuonCleanerBySegments:
    process.cleanedMu = cms.EDProducer("PATMuonCleanerBySegments",
                                       src = cms.InputTag("calibratedMuons"),
                                       preselection = cms.string("track.isNonnull"),
                                       passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
                                       fractionOfSharedSegments = cms.double(0.499))
    process.muons.replace(process.bareSoftMuons,cms.Sequence(process.cleanedMu+process.bareSoftMuons))    
    process.bareSoftMuons.src = "cleanedMu"
    if not APPLYMUCORR:
        process.cleanedMu.src = "slimmedMuons"


#------- ELECTRONS -------

#--- Run2 electron momentum scale and resolution corrections

process.selectedSlimmedElectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("slimmedElectrons"),
    cut = cms.string("pt>15 && abs(eta)<2.5")
)

#--- Photon ID modules seem to be OK also for UL cf: 
#--- https://github.com/cms-egamma/EgammaPostRecoTools/blob/master/python/EgammaPostRecoTools.py#L63

if (LEPTON_SETUP == 2016):
    from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
    if ( DATA_TAG == 'ULAPV'):
        setupEgammaPostRecoSeq(process,
                            #   runEnergyCorrections=True,
                            #   runVID=True,
                            #   eleIDModules=['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Summer16UL_ID_ISO_cff','RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff'],
                            #   phoIDModules=['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V2_cff'],
                              era='2016preVFP-UL')
    else:
        setupEgammaPostRecoSeq(process,
                            #   runEnergyCorrections=True,
                            #   runVID=True,
                            #   eleIDModules=['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Summer16UL_ID_ISO_cff','RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff'],
                            #   phoIDModules=['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V2_cff'],
                              era='2016postVFP-UL')
if (LEPTON_SETUP == 2017):
    from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
    setupEgammaPostRecoSeq(process,
                        #   runEnergyCorrections=True,
                        #   runVID=True,
                        #   eleIDModules=['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Summer17UL_ID_ISO_cff', 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff'],
                        #   phoIDModules=['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V2_cff'],
                          era='2017-UL')

if (LEPTON_SETUP == 2018):
    from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
    setupEgammaPostRecoSeq(process,
                        #   runEnergyCorrections=True,
                        #   runVID=True,
                        #   eleIDModules=['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Summer18UL_ID_ISO_cff','RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff','RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Fall17_iso_V2_cff'],
                        #   phoIDModules=['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Fall17_94X_V2_cff'],
                          era='2018-UL')

process.bareSoftElectrons = cms.EDFilter("PATElectronRefSelector",
   src = cms.InputTag("selectedSlimmedElectrons"),
   cut = cms.string("") #move pt>7 && abs(eta)<2.5 cut to softElectrons so that smear/scale corrections are done before the pT cut
   )

process.softElectrons = cms.EDProducer("EleFiller",
   src    = cms.InputTag("bareSoftElectrons"),
#    TriggerResultsLabel = cms.InputTag('TriggerResults','','HLT'),
#    TriggerSet = cms.InputTag(TRIGGERSET),
   sampleType = cms.int32(SAMPLE_TYPE),
   setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
   cut = cms.string('pt>10 && abs(eta) < 2.5 &&'+ DXY_DZ),
   flags = cms.PSet(
        isSIP = cms.string(SIP),
        isGood = cms.string(GOODELECTRON),
        isIsoFSRUncorr  = cms.string("userFloat('combRelIsoPF')<"+str(ELEISOCUT))
        ),
   )


process.electrons = cms.Sequence(process.egammaPostRecoSeq + process.selectedSlimmedElectrons + process.bareSoftElectrons + process.softElectrons)



#--- TrackLess Electrons
process.bareSoftPhotons = cms.EDFilter("PATPhotonRefSelector",
   src = cms.InputTag("slimmedPhotons"),
   cut = cms.string("pt>10 && abs(eta)<2.5")
   )

process.softPhotons = cms.EDProducer("Philler",
   src    = cms.InputTag("bareSoftPhotons"),
   srcElectron = cms.InputTag("softElectrons"),
   mvaValuesMap = cms.InputTag("photonMVAValueMapProducer:TLEMVAEstimatorRun2Fall15V1Values"),
   mvaValuesMap2 = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring15NonTrig25nsV2Values"),
   sampleType = cms.int32(SAMPLE_TYPE),
   setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
   cut = cms.string("1"), # dxy, dz not applied
   flags = cms.PSet(
        ID = cms.string("userFloat('isBDT')"),
        isSIP = cms.string(SIP),
        isGood = cms.string(GOODELECTRON),
        ),
   )


process.electronMatch = cms.EDProducer("MCMatcher",       # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                       src         = cms.InputTag("bareSoftElectrons"), # RECO objects to match
                                       matched     = cms.InputTag("prunedGenParticles"), # mc-truth particle collection
                                       mcPdgId     = cms.vint32(11),               # one or more PDG ID (11 = electron); absolute values (see below)
                                       checkCharge = cms.bool(True),               # True = require RECO and MC objects to have the same charge
                                       mcStatus    = cms.vint32(1),                # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
                                       maxDeltaR   = cms.double(0.5),              # Minimum deltaR for the match
                                       maxDPtRel   = cms.double(0.5),              # Minimum deltaPt/Pt for the match
                                       resolveAmbiguities    = cms.bool(True),     # Forbid two RECO objects to match to the same GEN object
                                       resolveByMatchQuality = cms.bool(False),    # False = just match input in order; True = pick lowest deltaR pair first
                                       )


### ----------------------------------------------------------------------
### Lepton Cleaning (clean electrons collection from muons)
### ----------------------------------------------------------------------

process.cleanSoftElectrons = cms.EDProducer("PATElectronCleaner",
    # pat electron input source
    src = cms.InputTag("softElectrons"),
    # preselection (any string-based cut for pat::Electron)
    preselection = cms.string(''),
    # overlap checking configurables
    checkOverlaps = cms.PSet(
        muons = cms.PSet(
           src       = cms.InputTag("softMuons"), # Start from loose lepton def
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string("userFloat('isGood')"),
           deltaR              = cms.double(0.05),
           checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
           pairCut             = cms.string(""),
           requireNoOverlaps   = cms.bool(True), # overlaps don't cause the electron to be discared
        )
    ),
    # finalCut (any string-based cut for pat::Electron)
    finalCut = cms.string(''),
)


#------- TAU LEPTONS -------

TAUCUT       = "pt>30 & abs(eta)<2.3"#"tauID('byCombinedIsolationDeltaBetaCorrRaw3Hits') < 1000.0 && pt>18"
SOSOTAU      = "decayMode()!=5 && decayMode()!=6 && tauID('decayModeFindingNewDMs') == 1 && userFloat('dz') < 1"
GOODTAU      = SOSOTAU + " && tauID('byVVVLooseDeepTau2017v2p1VSjet') == 1 && tauID('byVVVLooseDeepTau2017v2p1VSe') == 1 && tauID('byVLooseDeepTau2017v2p1VSmu') == 1"
GOODTAU_MU   = SOSOTAU + " && tauID('byTightDeepTau2017v2p1VSmu') == 1 && tauID('byVVLooseDeepTau2017v2p1VSe') == 1"
GOODTAU_ELE  = SOSOTAU + " && tauID('byVLooseDeepTau2017v2p1VSmu') == 1 && tauID('byTightDeepTau2017v2p1VSe') == 1"
GOODTAU_TAU  = SOSOTAU + " && tauID('byVLooseDeepTau2017v2p1VSmu') == 1 && tauID('byVVLooseDeepTau2017v2p1VSe') == 1"
ISOTAU       = "tauID('byMediumDeepTau2017v2p1VSjet') == 1"

# GOODTAU_MU   = SOSOTAU + " && tauID('byTightDeepTau2017v2p1VSmu') == 1 && tauID('byVLooseDeepTau2017v2p1VSe') == 1 && tauID('byTightDeepTau2017v2p1VSjet') == 1"
# GOODTAU_ELE  = SOSOTAU + " && tauID('byTightDeepTau2017v2p1VSmu') == 1 && tauID('byMediumDeepTau2017v2p1VSe') == 1 && tauID('byMediumDeepTau2017v2p1VSjet') == 1"
# GOODTAU_TAU  = SOSOTAU + " && tauID('byTightDeepTau2017v2p1VSmu') == 1 && tauID('byVLooseDeepTau2017v2p1VSe') == 1 && tauID('byTightDeepTau2017v2p1VSjet') == 1"


process.bareTaus = cms.EDFilter("PATTauRefSelector",
    src = cms.InputTag("slimmedTaus"),
    cut = cms.string(TAUCUT)
    )


##NOT USED FOR NOW, TBD Later
process.cleanTaus = cms.EDProducer("PATTauCleaner",
    src = cms.InputTag("bareTaus"),
    # preselection (any string-based cut on pat::Tau)
    preselection = cms.string(
            'tauID("decayModeFinding") > 0.5 &'
            ' tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits") > 0.5 &'
            ' tauID("againstMuonTight") > 0.5 &'
            ' tauID("againstElectronMedium") > 0.5'
        ),
 
   
   # overlap checking configurables
   checkOverlaps = cms.PSet(
      muons = cms.PSet(
          src       = cms.InputTag("cleanPatMuons"),
          algorithm = cms.string("byDeltaR"),
          preselection        = cms.string(""),
          deltaR              = cms.double(0.3),
          checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
          pairCut             = cms.string(""),
          requireNoOverlaps   = cms.bool(False), # overlaps don't cause the electron to be discared
          ),
      electrons = cms.PSet(
          src       = cms.InputTag("cleanPatElectrons"),
          algorithm = cms.string("byDeltaR"),
          preselection        = cms.string(""),
          deltaR              = cms.double(0.3),
          checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
          pairCut             = cms.string(""),
          requireNoOverlaps   = cms.bool(False), # overlaps don't cause the electron to be discared
          ),
      ),
        # finalCut (any string-based cut on pat::Tau)
        finalCut = cms.string(' '),
)

# TES: https://github.com/cms-tau-pog/TauIDSFs
# NominalTESCorrection=-1#in percent\
APPLYTESCORRECTION = APPLYTESCORRECTION if IsMC else False # always false if data


TESyear = "UL2016_preVFP"

if YEAR=='postVFP':
    TESyear = 'UL2016_postVFP'

if YEAR == 2017:
    TESyear = "UL2017"

if YEAR == 2018:
    TESyear = "UL2018"


process.softTaus = cms.EDProducer("TauFiller",
   src = cms.InputTag("bareTaus"),
   genCollection = cms.InputTag("prunedGenParticles"),
   vtxCollection = cms.InputTag("goodPrimaryVertices"),
   cut = cms.string(TAUCUT),
   discriminator = cms.string("byIsolationMVA3oldDMwoLTraw"),

   ApplyTESCentralCorr = cms.bool(APPLYTESCORRECTION),
   flags = cms.PSet(
        isSIP = cms.string(""),
        isGood = cms.string(GOODTAU),
        passCombRelIsoPFFSRCorr = cms.string(ISOTAU),
        isGood_Mu  = cms.string(GOODTAU_MU),
        isGood_Ele = cms.string(GOODTAU_ELE),
        isGood_Tau = cms.string(GOODTAU_TAU)
        ),

   year = cms.string(TESyear)
   )


process.taus=cms.Sequence(process.bareTaus + process.softTaus)


### ----------------------------------------------------------------------
### L1 Prefiring issue for 2016 and 2017 data
### ----------------------------------------------------------------------

# Recipe taken from https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1PrefiringWeightRecipe

if(IsMC and LEPTON_SETUP == 2016):
    from PhysicsTools.PatUtils.l1PrefiringWeightProducer_cfi import l1PrefiringWeightProducer
    if( DATA_TAG == 'ULAPV' ):
        process.prefiringweight = l1PrefiringWeightProducer.clone(
                                                                DataEraECAL = cms.string("UL2016preVFP"),
                                                                DataEraMuon = cms.string("2016preVFP"),
                                                                UseJetEMPt = cms.bool(False),
                                                                PrefiringRateSystematicUnctyECAL = cms.double(0.2),
                                                                PrefiringRateSystematicUnctyMuon = cms.double(0.2))
    else:
        process.prefiringweight = l1PrefiringWeightProducer.clone(
                                                                DataEraECAL = cms.string("UL2016postVFP"),
                                                                DataEraMuon = cms.string("2016postVFP"),
                                                                UseJetEMPt = cms.bool(False),
                                                                PrefiringRateSystematicUnctyECAL = cms.double(0.2),
                                                                PrefiringRateSystematicUnctyMuon = cms.double(0.2))

if(IsMC and LEPTON_SETUP == 2017):
    from PhysicsTools.PatUtils.l1PrefiringWeightProducer_cfi import l1PrefiringWeightProducer
    process.prefiringweight = l1PrefiringWeightProducer.clone(
                                                            TheJets = cms.InputTag("patJetsReapplyJEC"),
                                                            DataEraECAL = cms.string("UL2017BtoF"),
                                                            DataEraMuon = cms.string("20172018"),
                                                            UseJetEMPt = cms.bool(False),
                                                            PrefiringRateSystematicUnctyECAL = cms.double(0.2),
                                                            PrefiringRateSystematicUnctyMuon = cms.double(0.2))
    
if(IsMC and LEPTON_SETUP == 2018):
    from PhysicsTools.PatUtils.l1PrefiringWeightProducer_cfi import l1PrefiringWeightProducer
    process.prefiringweight = l1PrefiringWeightProducer.clone(
                                                             TheJets = cms.InputTag("patJetsReapplyJEC"),
                                                             DataEraECAL = cms.string("None"),
                                                             DataEraMuon = cms.string("20172018"),
                                                             UseJetEMPt = cms.bool(False),
                                                             PrefiringRateSystematicUnctyECAL = cms.double(0.2),
                                                             PrefiringRateSystematicUnctyMuon = cms.double(0.2))

if(IsMC):
    process.Prefiring = cms.Path(process.prefiringweight)



### ----------------------------------------------------------------------
### Search for FSR candidates
### ----------------------------------------------------------------------

# Create a photon collection; cfg extracted from "UFHHTauTauHMuMuRun2.FSRPhotons.fsrPhotons_cff"
process.fsrPhotons = cms.EDProducer("PhotonFiller",
    electronSrc = cms.InputTag("slimmedElectrons"),
    sampleType = cms.int32(SAMPLE_TYPE),
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
    photonSel = cms.string(FSRMODE)  # "skip", "passThrough", "RunII"
)

import PhysicsTools.PatAlgos.producersLayer1.pfParticleProducer_cfi
process.boostedFsrPhotons = PhysicsTools.PatAlgos.producersLayer1.pfParticleProducer_cfi.patPFParticles.clone(
    pfCandidateSource = 'fsrPhotons'
)

process.appendPhotons = cms.EDProducer("LeptonPhotonMatcher",
    muonSrc = cms.InputTag("softMuons"),
    electronSrc = cms.InputTag("cleanSoftElectrons"),
    photonSrc = cms.InputTag("boostedFsrPhotons"),
    sampleType = cms.int32(SAMPLE_TYPE),
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
    photonSel = cms.string(FSRMODE),  # "skip", "passThrough", "RunII"
    muon_iso_cut = cms.double(MUISOCUT),
    electron_iso_cut = cms.double(ELEISOCUT),
    debug = cms.untracked.bool(False),
    )

if ADDLOOSEELE:
    process.appendPhotons.looseElectronSrc = cms.InputTag("cleanSoftLooseElectrons")
#    process.appendPhotons.tleSrc = cms.InputTag("softPhotons")
#    process.appendPhotons.TLEMinPt = cms.double(25.)

# All leptons, any F/C.
# CAVEAT: merging creates copies of the objects, so that CandViewShallowCloneCombiner is not able to find
# overlaps between merged collections and the original ones.
process.softLeptons = cms.EDProducer("CandViewMerger",
#    src = cms.VInputTag(cms.InputTag("softMuons"), cms.InputTag("cleanSoftElectrons"))
    src = cms.VInputTag(cms.InputTag("appendPhotons:muons"), cms.InputTag("appendPhotons:electrons"), cms.InputTag("softTaus"))
)


        
#----------------- MET, also used for SVFit -------------------------

process.METSequence = cms.Sequence()

PFMetName = "slimmedMETs"
uncorrPFMetTag = cms.InputTag(PFMetName)

    # Shift met due to central corrections of TES, EES, JER and JES
process.ShiftMETcentral = cms.EDProducer ("ShiftMETcentral",
					srcMET = uncorrPFMetTag,
					tauUncorrected = cms.InputTag("bareTaus"),
					tauCorrected = cms.InputTag("softTaus"),
                    jetCollection = cms.InputTag("dressedJets")
					)

srcMETTag = None
srcMETTag = cms.InputTag("ShiftMETcentral")

process.METSignificance = cms.EDProducer ("ExtractMETSignificance", 
					srcMET=uncorrPFMetTag 
					)

    # add variables with MET shifted for TES corrections
process.ShiftMETforTES = cms.EDProducer ("ShiftMETforTES", 
					srcMET  = srcMETTag, 
					tauCollection = cms.InputTag("softTaus")
					)

    # add variables with MET shifted for EES corrections (E->tau ES)
process.ShiftMETforEES = cms.EDProducer ("ShiftMETforEES",
					srcMET  = srcMETTag,
					tauCollection = cms.InputTag("softTaus")
					)

    # add variables with MET shifted for MES corrections (Mu->tau ES)
process.ShiftMETforMES = cms.EDProducer ("ShiftMETforMES",
					srcMET  = srcMETTag,
					tauCollection = cms.InputTag("softTaus")
					)

    # add variables with MET shifted for JES corrections
process.ShiftMETforJES = cms.EDProducer ("ShiftMETforJES", 
					srcMET  = srcMETTag, 
					jetCollection = cms.InputTag("dressedJets")
					)

    # add variables with MET shifted for EES corrections (E->tau ES)
process.ShiftMETforJER = cms.EDProducer ("ShiftMETforJER",
					srcMET  = srcMETTag,
					jetCollection = cms.InputTag("dressedJets")
					)

    # Get a standalone Puppi MET significance collection
#process.PuppiMETSignificance = cms.EDProducer ("ExtractMETSignificance",
#					srcMET=cms.InputTag("slimmedMETsPuppi")
#					)

    # Shift PUPPI met due to central corrections of TES and EES
#process.ShiftPuppiMETcentral = cms.EDProducer ("ShiftMETcentral",
#					srcMET = cms.InputTag("slimmedMETsPuppi"),
#					tauUncorrected = cms.InputTag("bareTaus"),
#					tauCorrected = cms.InputTag("softTaus")
#					)

process.METSequence += process.METSignificance
process.METSequence += process.ShiftMETcentral
# process.METSequence += process.ShiftMETforTES
# process.METSequence += process.ShiftMETforEES
# process.METSequence += process.ShiftMETforMES
process.METSequence += process.ShiftMETforJES
# process.METSequence += process.ShiftMETforJER
#process.METSequence += process.PuppiMETSignificance
#process.METSequence += process.ShiftPuppiMETcentral
metTag=uncorrPFMetTag

process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')

# In 2017 and 2018 some problematic crystals --> pass list of crystals
if YEAR == 2017 or YEAR == 2018:
    baddetEcallist = cms.vuint32(
            [872439604,872422825,872420274,872423218,
            872423215,872416066,872435036,872439336,
            872420273,872436907,872420147,872439731,
            872436657,872420397,872439732,872439339,
            872439603,872422436,872439861,872437051,
            872437052,872420649,872422436,872421950,
            872437185,872422564,872421566,872421695,
            872421955,872421567,872437184,872421951,
            872421694,872437056,872437057,872437313])

# In 2016 no problem --> pass empty list
if YEAR == 2016:
    baddetEcallist = cms.vuint32([])

process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
        "EcalBadCalibFilter",
        EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
        ecalMinEt        = cms.double(50.),
        baddetEcal    = baddetEcallist,
        taggingMode = cms.bool(True),
        debug = cms.bool(False)
    )


### ----------------------------------------------------------------------
### ----------------------------------------------------------------------
### BUILD CANDIDATES
### ----------------------------------------------------------------------
### ----------------------------------------------------------------------



### ----------------------------------------------------------------------
### Dileptons: combine/merge leptons into intermediate (bare) collections;
###            Embed additional user variables into final collections
### ----------------------------------------------------------------------
TWOGOODLEPTONS = "( userFloat('d0.isGood') && userFloat('d1.isGood') && userFloat('isGoodTau') )" # Z made of 2 good leptons (ISO not yet applied)
ISOLEPTON1 = "userFloat('d0.passCombRelIsoPFFSRCorr')"
ISOLEPTON2 = "userFloat('d1.passCombRelIsoPFFSRCorr')"
OS = "daughter(0).pdgId()*daughter(1).pdgId() < 0"
SS = "daughter(0).pdgId()*daughter(1).pdgId() > 0"
TWOISOLEPTONS = "( userFloat('d0.passCombRelIsoPFFSRCorr') && userFloat('d1.passCombRelIsoPFFSRCorr') )"
TWOSFLEPTONS = "(abs(daughter(0).pdgId())==15 || abs(daughter(1).pdgId())==15 || abs(daughter(0).pdgId()*daughter(1).pdgId())==143 || abs(daughter(0).pdgId()*daughter(1).pdgId())==169)"
### NOTE: Isolation cut has been moved to ZZ candidates as we now correct for FSR of all four photons.
### Because if this, isBestZ flags are no longer correct; BESTZ_AMONG is set to "" for safety

### ----------------------------------------------------------------------
### Dileptons (Z->mm, Z->etau, Z->mutau, Z->tautau, Z->emu)
### ----------------------------------------------------------------------

# l+l- (SFOS, e and mu and tau)
# process.bareZCand = cms.EDProducer("PATCandViewShallowCloneCombiner",
#     decay = cms.string('softLeptons@+ softLeptons@-'),
#     cut = cms.string("True"), # see below
#     checkCharge = cms.bool(True)
# )
process.bareZCand = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('softLeptons softLeptons'),
    #cut = cms.string('deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi)>0.02'), # protect against ghosts
    cut = cms.string("True"),
    checkCharge = cms.bool(False)
)

if KEEPLOOSECOMB:
    process.bareZCand.cut = cms.string('mass > 0 && '+TWOSFLEPTONS) # Propagate also combinations of loose leptons (for debugging)
else:
    if FSRMODE == "RunII" : # Just keep combinations of tight leptons (passing ID, SIP and ISO)
        process.bareZCand.cut = cms.string("deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi)>0.02 && mass > 0 && "+TWOSFLEPTONS+" && daughter(0).masterClone.userFloat('isGood') && daughter(1).masterClone.userFloat('isGood') && daughter(0).masterClone.userFloat('passLooseCombRelIsoPFFSRCorr') && daughter(1).masterClone.userFloat('passLooseCombRelIsoPFFSRCorr')")
    else : # Just keep combinations of tight leptons (passing ID and SIP; iso cannot be required at this point if FSRMode is "skip")
        process.bareZCand.cut = cms.string("mass > 0 && "+TWOSFLEPTONS+" && daughter(0).masterClone.userFloat('isGood') && daughter(1).masterClone.userFloat('isGood')")
        

#FLAVOUR     = "((daughter(0).pdgId()*daughter(1).pdgId()==-121 && userFloat('eleHLTMatch')) ||  (daughter(0).pdgId()*daughter(1).pdgId()==-169 && userFloat('muHLTMatch')))"

# BESTZ_AMONG = ( TWOGOODLEPTONS + "&&" + TWOISOLEPTONS + "&& userFloat('DR')>0.5" )
BESTZ_AMONG = ( TWOGOODLEPTONS + "&& userFloat('DR')>0.5" )

TWOGOODISOLEPTONS = ( TWOGOODLEPTONS + "&&" + TWOISOLEPTONS )

process.ZCand = cms.EDProducer("ZCandidateFiller",
    src	       = cms.InputTag("bareZCand"),
                               
    TriggerResultsLabel = cms.InputTag('TriggerResults','','HLT'),
    TriggerSet = cms.InputTag(TRIGGERSET),
                               
    sampleType = cms.int32(SAMPLE_TYPE),
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
    bestZAmong = cms.string(BESTZ_AMONG),
    FSRMode = cms.string(FSRMODE), # "skip", "RunII"
    flags = cms.PSet(
        GoodLeptons = cms.string(TWOGOODLEPTONS),
        IsoLeptons = cms.string(TWOISOLEPTONS),
        GoodIsoLeptons = cms.string(TWOGOODISOLEPTONS),
        OS = cms.string(OS),
        SS = cms.string(SS),
        IsoLepton1 = cms.string(ISOLEPTON1),
        IsoLepton2 = cms.string(ISOLEPTON2)
        # Z1Presel = cms.string(Z1PRESEL),
    )
)


# ll, same flavour/any charge, for control regions only
# process.bareLLCand = cms.EDProducer("CandViewShallowCloneCombiner",
#     decay = cms.string('softLeptons softLeptons'),
#     #cut = cms.string('deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi)>0.02'), # protect against ghosts
#     cut = cms.string("deltaR(daughter(0).eta, daughter(0).phi, daughter(1).eta, daughter(1).phi)>0.02 && " + TWOSFLEPTONS), # protect against ghosts && same flavour
#     checkCharge = cms.bool(False)
# )


# process.LLCand = cms.EDProducer("ZCandidateFiller",
#     src = cms.InputTag("bareLLCand"),
#     sampleType = cms.int32(SAMPLE_TYPE),
#     setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
#     bestZAmong = cms.string(BESTZ_AMONG),
#     FSRMode = cms.string(FSRMODE), # "skip", "RunII"

#     TriggerResults = cms.InputTag('TriggerResults','','HLT'),
#     TriggerSet = cms.InputTag(TRIGGERSET),

#     flags = cms.PSet(
#         GoodLeptons = cms.string(ZLEPTONSEL),
#         Z1Presel = cms.string(Z1PRESEL),
#     )
# )


### ----------------------------------------------------------------------
### TriLeptons (for fake rate)
### ----------------------------------------------------------------------
# Z_PLUS_LEP_MIJ=("sqrt(pow(daughter(0).daughter({0}).energy+daughter(1).energy, 2) - " +
#                 "     pow(daughter(0).daughter({0}).px    +daughter(1).px, 2) -" +
#                 "     pow(daughter(0).daughter({0}).py    +daughter(1).py, 2) -" +
#                 "     pow(daughter(0).daughter({0}).pz    +daughter(1).pz, 2))")

# process.ZlCand = cms.EDProducer("PATCandViewShallowCloneCombiner",
#     decay = cms.string('ZCand softLeptons'),
#     cut = cms.string("deltaR(daughter(0).daughter(0).eta, daughter(0).daughter(0).phi, daughter(1).eta, daughter(1).phi)>0.02 &&" + # Ghost suppression
#                      "deltaR(daughter(0).daughter(1).eta, daughter(0).daughter(1).phi, daughter(1).eta, daughter(1).phi)>0.02 &&" +
#                      ("(daughter(0).daughter(0).charge == daughter(1).charge || %s > 4) && " % ( Z_PLUS_LEP_MIJ.format(0))) +       # mLL>4 for the OS pair (Giovanni's impl)
#                      ("(daughter(0).daughter(1).charge == daughter(1).charge || %s > 4) && " % ( Z_PLUS_LEP_MIJ.format(1))) +
#                      "daughter(0).masterClone.userFloat('isBestZ') &&" + # This includes the Z1 isolation requirement
#                      "daughter(0).masterClone.userFloat('Z1Presel')"
#                      ),
#     checkCharge = cms.bool(False)
# )


### ----------------------------------------------------------------------
### Jets
### ----------------------------------------------------------------------

# DEFAULT 2016 jet pileup ID training: _chsalgos_81x
# 2017-2018: jet pileup ID trainings to be updated
# Updated for UL Training: _chsalgos_106X_UL17, _chsalgos_106X_UL18
from RecoJets.JetProducers.PileupJetID_cfi import _chsalgos_106X_UL16, _chsalgos_106X_UL17, _chsalgos_106X_UL18
from RecoJets.JetProducers.PileupJetIDCutParams_cfi import full_81x_chs_wp
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

process.load("CondCore.CondDB.CondDB_cfi")

# Update jet collection in 2016 MC ntuples to include DeepCSV and DeepFlavour b tag algorithms
# NB: In 2016 Data we do not want to update jet collection to include DeepCSV info (it is already there!)

if (SAMPLE_TYPE == 2016):
    process.load("RecoJets.JetProducers.PileupJetID_cfi")
    process.pileupJetIdUpdated = process.pileupJetId.clone(
        jets=cms.InputTag("slimmedJets"),
        inputIsCorrected=False,
        applyJec=True,
        vertexes=cms.InputTag("offlineSlimmedPrimaryVertices"),
        algos=cms.VPSet(_chsalgos_106X_UL16)
    )  
elif (SAMPLE_TYPE == 2017):
    process.load("RecoJets.JetProducers.PileupJetID_cfi")
    process.pileupJetIdUpdated = process.pileupJetId.clone(
        jets=cms.InputTag("slimmedJets"),
        inputIsCorrected=False,
        applyJec=True,
        vertexes=cms.InputTag("offlineSlimmedPrimaryVertices"),
        algos=cms.VPSet(_chsalgos_106X_UL17)
    )
elif (SAMPLE_TYPE == 2018):
    process.load("RecoJets.JetProducers.PileupJetID_cfi")
    process.pileupJetIdUpdated = process.pileupJetId.clone(
        jets=cms.InputTag("slimmedJets"),
        inputIsCorrected=False,
        applyJec=True,
        vertexes=cms.InputTag("offlineSlimmedPrimaryVertices"),
        algos=cms.VPSet(_chsalgos_106X_UL18)
    )
else:
    process.load("RecoJets.JetProducers.PileupJetID_cfi")
    process.pileupJetIdUpdated = process.pileupJetId.clone(
        jets=cms.InputTag("slimmedJets"),
        inputIsCorrected=False,
        applyJec=True,
        vertexes=cms.InputTag("offlineSlimmedPrimaryVertices"),
        # Make sure correct algo is used
        algos=cms.VPSet(_chsalgos_106X_UL17)
    )

#-- BTagging updated following UL recipes:
#-- https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation
#-- 2016APV DeepCSV, medium WP: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16preVFP
#-- 2016 DeepCSV, medium WP: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL16postVFP
#-- 2018 DeepCSV, medium WP: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL18
#-- 2017 DeepCSV, medium WP: https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation106XUL17

theBTagger=""
theBTaggerThr=0
theBTagSFFile=""
theBTagMCEffFile=""
if (LEPTON_SETUP == 2016):
    if( DATA_TAG == 'ULAPV'):
        theBTagger="pfDeepCSVJetTags:probb"
        theBTaggerThr=0.6001
        theBTagSFFile="HTauTauHMuMu/AnalysisStep/data/BTagging/DeepCSV_106XUL16preVFPSF_v1_hzz.csv"
        theBTagMCEffFile="HTauTauHMuMu/AnalysisStep/data/BTagging/bTagEfficiencies_2016_LegacyPaper.root"
    else:
        theBTagger="pfDeepCSVJetTags:probb"
        theBTaggerThr=0.5847
        theBTagSFFile="HTauTauHMuMu/AnalysisStep/data/BTagging/DeepCSV_106XUL16postVFPSF_v2_hzz.csv"
        theBTagMCEffFile="HTauTauHMuMu/AnalysisStep/data/BTagging/bTagEfficiencies_2016_LegacyPaper.root"
elif (LEPTON_SETUP == 2017):
    theBTagger="pfDeepCSVJetTags:probb"
    theBTaggerThr=0.4506
    theBTagSFFile="HTauTauHMuMu/AnalysisStep/data/BTagging/wp_deepCSV_106XUL17_v3_hzz.csv"
    theBTagMCEffFile="HTauTauHMuMu/AnalysisStep/data/BTagging/bTagEfficiencies_2017_LegacyPaper.root"
elif (LEPTON_SETUP == 2018):
    theBTagger="pfDeepCSVJetTags:probb"
    theBTaggerThr=0.4168
    theBTagSFFile="HTauTauHMuMu/AnalysisStep/data/BTagging/wp_deepCSV_106XUL18_v2_hzz.csv"
    theBTagMCEffFile="HTauTauHMuMu/AnalysisStep/data/BTagging/bTagEfficiencies_2018_LegacyPaper.root"
else:
    sys.exit("HTauTauHMuMu.py: Need to define the btagging for the new setup!")

### q/g likelihood
qgDatabaseVersion = 'cmssw8020_v2'
process.QGPoolDBESSource = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(messageLevel = cms.untracked.int32(1)),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
        cms.PSet(
            record = cms.string('QGLikelihoodRcd'),
            tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_AK4PFchs'),
            label  = cms.untracked.string('QGL_AK4PFchs')
        ),
      ),
      connect = cms.string('sqlite_fip:HTauTauHMuMu/AnalysisStep/data/QGTagging/QGL_'+qgDatabaseVersion+'.db')
)
process.es_prefer_qg = cms.ESPrefer('PoolDBESSource','QGPoolDBESSource')
process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets = cms.InputTag( 'slimmedJets' )
process.QGTagger.jetsLabel = cms.string('QGL_AK4PFchs')

### DRESSED JETS: b tag applied
process.dressedJets = cms.EDProducer("JetFiller",
    src = cms.InputTag("slimmedJets"),
    sampleType = cms.int32(SAMPLE_TYPE),
    setup = cms.int32(LEPTON_SETUP),
    cut = cms.string(""),#("pt>20 && abs(eta)<4.7 && userFloat('JetID') && (userFloat('PUjetID') || pt>50)"),
    isMC = cms.bool(IsMC),
    bTaggerName = cms.string(theBTagger),
    bTaggerThreshold = cms.double(theBTaggerThr),
    applyJEC = cms.bool(APPLYJEC),
    jecType = cms.string("AK4PFchs"),
    applyJER = cms.bool(APPLYJER),
    jerType = cms.string("AK4PFchs"),
    bTagSFFile = cms.string(theBTagSFFile),
    bTagMCEffFile = cms.string(theBTagMCEffFile),
    flags = cms.PSet()
    )


### Load JEC
#Taken from https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC#Recommended_for_MC

if (APPLYJEC and SAMPLE_TYPE == 2016):
    if IsMC:
        if ( DATA_TAG == 'ULAPV' ):
            process.jec = cms.ESSource("PoolDBESSource",
                DBParameters = cms.PSet(
                    messageLevel = cms.untracked.int32(1)
                    ),
                timetype = cms.string('runnumber'),
                toGet = cms.VPSet(
                    cms.PSet(
                        record = cms.string('JetCorrectionsRecord'),
                        tag    = cms.string('JetCorrectorParametersCollection_Summer19UL16APV_V7_MC_AK4PFchs'),
                        label  = cms.untracked.string('AK4PFchs')
                        ),
                    ),
            connect = cms.string('sqlite_fip:HTauTauHMuMu/AnalysisStep/data/JEC/Summer19UL16APV_V7_MC.db'),             )
        else:
            process.jec = cms.ESSource("PoolDBESSource",
                DBParameters = cms.PSet(
                    messageLevel = cms.untracked.int32(1)
                    ),
                timetype = cms.string('runnumber'),
                toGet = cms.VPSet(
                    cms.PSet(
                        record = cms.string('JetCorrectionsRecord'),
                        tag    = cms.string('JetCorrectorParametersCollection_Summer19UL16_V7_MC_AK4PFchs'),
                        label  = cms.untracked.string('AK4PFchs')
                        ),
                    ),
            connect = cms.string('sqlite_fip:HTauTauHMuMu/AnalysisStep/data/JEC/Summer19UL16_V7_MC.db'),             )
    else:
        # For Data no distinction between pre and postVFP periods: https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC#2016_Data
        process.jec = cms.ESSource("PoolDBESSource",
            DBParameters = cms.PSet(
                messageLevel = cms.untracked.int32(1)
                ),
            timetype = cms.string('runnumber'),
            toGet = cms.VPSet(
                cms.PSet(
                    record = cms.string('JetCorrectionsRecord'),
                    tag    = cms.string('JetCorrectorParametersCollection_Summer19UL16_RunBCDEFGH_Combined_V7_DATA_AK4PFchs'),
                    label  = cms.untracked.string('AK4PFchs')
                    ),
                ),
            connect = cms.string('sqlite_fip:HTauTauHMuMu/AnalysisStep/data/JEC/Summer19UL16_RunBCDEFGH_Combined_V7_DATA.db'),
            )

    ## Add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
    process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

    ### REAPPLY JEC
    #--- In principle now in UL we have DeepCSV already in the samples

    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJetCorrFactors
    process.patJetCorrFactorsReapplyJEC = updatedPatJetCorrFactors.clone(
        src = cms.InputTag("slimmedJets"),
        levels = ['L1FastJet','L2Relative','L3Absolute'],
        payload = 'AK4PFchs' )
    if not IsMC:
        process.patJetCorrFactorsReapplyJEC.levels = ['L1FastJet','L2Relative','L3Absolute','L2L3Residual']

    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJets
    process.patJetsReapplyJEC = updatedPatJets.clone(
        jetSource = cms.InputTag("slimmedJets"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
        )

    ### Add pileup id and discriminant to patJetsReapplyJEC
    process.patJetsReapplyJEC.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']
    process.patJetsReapplyJEC.userData.userInts.src += ['pileupJetIdUpdated:fullId']
    
    ### Replace inputs in QGTagger and dressedJets
    process.QGTagger.srcJets = cms.InputTag( 'patJetsReapplyJEC')
    process.dressedJets.src = cms.InputTag('patJetsReapplyJEC')


if (APPLYJEC and SAMPLE_TYPE == 2017):
    if IsMC:
        process.jec = cms.ESSource("PoolDBESSource",
            DBParameters = cms.PSet(
                messageLevel = cms.untracked.int32(1)
                ),
            timetype = cms.string('runnumber'),
            toGet = cms.VPSet(
                cms.PSet(
                    record = cms.string('JetCorrectionsRecord'),
                    tag    = cms.string('JetCorrectorParametersCollection_Summer19UL17_V5_MC_AK4PFchs'),
                    label  = cms.untracked.string('AK4PFchs')
                    ),
                ),
              connect = cms.string('sqlite_fip:HTauTauHMuMu/AnalysisStep/data/JEC/Summer19UL17_V5_MC.db'),
            )
    else:
        process.jec = cms.ESSource("PoolDBESSource",
            DBParameters = cms.PSet(
                messageLevel = cms.untracked.int32(1)
                ),
            timetype = cms.string('runnumber'),
            toGet = cms.VPSet(
                cms.PSet(
                    record = cms.string('JetCorrectionsRecord'),
                    tag    = cms.string('JetCorrectorParametersCollection_Summer19UL17_RunBCDEF_V5_DATA_AK4PFchs'),
                    label  = cms.untracked.string('AK4PFchs')
                    ),
                ),
            connect = cms.string('sqlite_fip:HTauTauHMuMu/AnalysisStep/data/JEC/Summer19UL17_RunBCDEF_V5_DATA.db'),
            )

    ## Add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
    process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

    ### REAPPLY JEC
    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJetCorrFactors
    process.patJetCorrFactorsReapplyJEC = updatedPatJetCorrFactors.clone(
        src = cms.InputTag("slimmedJets"),
        levels = ['L1FastJet','L2Relative','L3Absolute'],
        payload = 'AK4PFchs' )
    if not IsMC:
        process.patJetCorrFactorsReapplyJEC.levels = ['L1FastJet','L2Relative','L3Absolute','L2L3Residual']

    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJets
    process.patJetsReapplyJEC = updatedPatJets.clone(
        jetSource = cms.InputTag("slimmedJets"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
        )

    ### Add pileup id and discriminant to patJetsReapplyJEC
    process.patJetsReapplyJEC.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']
    process.patJetsReapplyJEC.userData.userInts.src += ['pileupJetIdUpdated:fullId']

    ### Replace inputs in QGTagger and dressedJets
    process.QGTagger.srcJets = cms.InputTag( 'patJetsReapplyJEC')
    process.dressedJets.src = cms.InputTag('patJetsReapplyJEC')


if (APPLYJEC and SAMPLE_TYPE == 2018):
    if IsMC:
        process.jec = cms.ESSource("PoolDBESSource",
            DBParameters = cms.PSet(
                messageLevel = cms.untracked.int32(1)
                ),
            timetype = cms.string('runnumber'),
            toGet = cms.VPSet(
                cms.PSet(
                    record = cms.string('JetCorrectionsRecord'),
                    tag    = cms.string('JetCorrectorParametersCollection_Summer19UL18_V5_MC_AK4PFchs'), #for 10_2_X MC
                    label  = cms.untracked.string('AK4PFchs')
                    ),
                ),
              connect = cms.string('sqlite_fip:HTauTauHMuMu/AnalysisStep/data/JEC/Summer19UL18_V5_MC.db'), #for Summer19UL MC
            )
    else:
        process.jec = cms.ESSource("PoolDBESSource",
            DBParameters = cms.PSet(
                messageLevel = cms.untracked.int32(1)
                ),
            timetype = cms.string('runnumber'),
            toGet = cms.VPSet(
                cms.PSet(
                    record = cms.string('JetCorrectionsRecord'),
                    tag    = cms.string('JetCorrectorParametersCollection_Summer19UL18_V5_DATA_AK4PFchs'),
                    label  = cms.untracked.string('AK4PFchs')
                    ),
                ),
            connect = cms.string('sqlite_fip:HTauTauHMuMu/AnalysisStep/data/JEC/Summer19UL18_V5_DATA.db'),
            )

    ## Add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
    process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

    ### REAPPLY JEC
    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJetCorrFactors
    process.patJetCorrFactorsReapplyJEC = updatedPatJetCorrFactors.clone(
        src = cms.InputTag("slimmedJets"),
        levels = ['L1FastJet','L2Relative','L3Absolute'],
        payload = 'AK4PFchs' )
    if not IsMC:
        process.patJetCorrFactorsReapplyJEC.levels = ['L1FastJet','L2Relative','L3Absolute','L2L3Residual']

    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJets
    process.patJetsReapplyJEC = updatedPatJets.clone(
        jetSource = cms.InputTag("slimmedJets"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
        )

    ### Add pileup id and discriminant to patJetsReapplyJEC
    process.patJetsReapplyJEC.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']
    process.patJetsReapplyJEC.userData.userInts.src += ['pileupJetIdUpdated:fullId']

    ### Replace inputs in QGTagger and dressedJets
    process.QGTagger.srcJets = cms.InputTag( 'patJetsReapplyJEC')
    process.dressedJets.src = cms.InputTag('patJetsReapplyJEC')


# JER from https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution#Smearing_procedures
# Updated for UL (all three years)
if (APPLYJER and SAMPLE_TYPE == 2016):
    process.load('Configuration.StandardSequences.Services_cff')
    process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
    from CondCore.DBCommon.CondDBSetup_cfi import *
    if ( DATA_TAG == 'ULAPV'):
        process.jer = cms.ESSource("PoolDBESSource",
                                     CondDBSetup,
                                     connect = cms.string('sqlite_fip:HTauTauHMuMu/AnalysisStep/data/JER/Summer20UL16APV_JRV3_MC.db'),
                                     toGet = cms.VPSet(
                                                       cms.PSet(
                                                                record = cms.string('JetResolutionRcd'),
                                                                tag    = cms.string('JR_Summer20UL16APV_JRV3_MC_PtResolution_AK4PFchs'),
                                                                label  = cms.untracked.string('AK4PFchs_pt')
                                                                ),
                                                       cms.PSet(
                                                                record = cms.string('JetResolutionRcd'),
                                                                tag    = cms.string('JR_Summer20UL16APV_JRV3_MC_PhiResolution_AK4PFchs'),
                                                                label  = cms.untracked.string('AK4PFchs_phi')
                                                                ),
                                                       cms.PSet(
                                                                record = cms.string('JetResolutionScaleFactorRcd'),
                                                                tag    = cms.string('JR_Summer20UL16APV_JRV3_MC_SF_AK4PFchs'),
                                                                label  = cms.untracked.string('AK4PFchs')
                                                                )
                                                       )
                                     )
        process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')
    else:
        process.jer = cms.ESSource("PoolDBESSource",
                                     CondDBSetup,
                                     connect = cms.string('sqlite_fip:HTauTauHMuMu/AnalysisStep/data/JER/Summer20UL16_JRV3_MC.db'),
                                     toGet = cms.VPSet(
                                                       cms.PSet(
                                                                record = cms.string('JetResolutionRcd'),
                                                                tag    = cms.string('JR_Summer20UL16_JRV3_MC_PtResolution_AK4PFchs'),
                                                                label  = cms.untracked.string('AK4PFchs_pt')
                                                                ),
                                                       cms.PSet(
                                                                record = cms.string('JetResolutionRcd'),
                                                                tag    = cms.string('JR_Summer20UL16_JRV3_MC_PhiResolution_AK4PFchs'),
                                                                label  = cms.untracked.string('AK4PFchs_phi')
                                                                ),
                                                       cms.PSet(
                                                                record = cms.string('JetResolutionScaleFactorRcd'),
                                                                tag    = cms.string('JR_Summer20UL16_JRV3_MC_SF_AK4PFchs'),
                                                                label  = cms.untracked.string('AK4PFchs')
                                                                )
                                                       )
                                     )
        process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')

if (APPLYJER and SAMPLE_TYPE == 2017):
    process.load('Configuration.StandardSequences.Services_cff')
    process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
    from CondCore.DBCommon.CondDBSetup_cfi import *
    process.jer = cms.ESSource("PoolDBESSource",
                                 CondDBSetup,
                                 connect = cms.string('sqlite_fip:HTauTauHMuMu/AnalysisStep/data/JER/Summer19UL17_JRV3_MC.db'),
                                 toGet = cms.VPSet(
                                                   cms.PSet(
                                                            record = cms.string('JetResolutionRcd'),
                                                            tag    = cms.string('JR_Summer19UL17_JRV3_MC_PtResolution_AK4PFchs'),
                                                            label  = cms.untracked.string('AK4PFchs_pt')
                                                            ),
                                                   cms.PSet(
                                                            record = cms.string('JetResolutionRcd'),
                                                            tag    = cms.string('JR_Summer19UL17_JRV3_MC_PhiResolution_AK4PFchs'),
                                                            label  = cms.untracked.string('AK4PFchs_phi')
                                                            ),
                                                   cms.PSet(
                                                            record = cms.string('JetResolutionScaleFactorRcd'),
                                                            tag    = cms.string('JR_Summer19UL17_JRV3_MC_SF_AK4PFchs'),
                                                            label  = cms.untracked.string('AK4PFchs')
                                                            )
                                                   )
                                 )
    process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')


if (APPLYJER and SAMPLE_TYPE == 2018):
    process.load('Configuration.StandardSequences.Services_cff')
    process.load("JetMETCorrections.Modules.JetResolutionESProducer_cfi")
    from CondCore.DBCommon.CondDBSetup_cfi import *
    process.jer = cms.ESSource("PoolDBESSource",
                               CondDBSetup,
                               #connect = cms.string('sqlite_fip:HTauTauHMuMu/AnalysisStep/data/JER/Autumn18_V1_MC.db'),
                               connect = cms.string('sqlite_fip:HTauTauHMuMu/AnalysisStep/data/JER/Summer19UL18_JRV2_MC.db'),
                               toGet = cms.VPSet(
                                                 cms.PSet(
                                                          record = cms.string('JetResolutionRcd'),
                                                          tag    = cms.string('JR_Summer19UL18_JRV2_MC_PtResolution_AK4PFchs'),
                                                          label  = cms.untracked.string('AK4PFchs_pt')
                                                          ),
                                                 cms.PSet(
                                                          record = cms.string('JetResolutionRcd'),
                                                          tag    = cms.string('JR_Summer19UL18_JRV2_MC_PhiResolution_AK4PFchs'),
                                                          label  = cms.untracked.string('AK4PFchs_phi')
                                                          ),
                                                 cms.PSet(
                                                          record = cms.string('JetResolutionScaleFactorRcd'),
                                                          tag    = cms.string('JR_Summer19UL18_JRV2_MC_SF_AK4PFchs'),
                                                          label  = cms.untracked.string('AK4PFchs')
                                                          )
                                                 )
                               )
    process.es_prefer_jer = cms.ESPrefer('PoolDBESSource', 'jer')

### Clean jets wrt. good (preFSR-)isolated leptons
process.cleanJets = cms.EDProducer("JetsWithLeptonsRemover",
                                   Jets      = cms.InputTag("dressedJets"),
                                   Muons     = cms.InputTag("appendPhotons:muons"),
                                   Electrons = cms.InputTag("appendPhotons:electrons"),
                                   Diboson   = cms.InputTag(""),
                                   JetPreselection      = cms.string("pt>20 && abs(eta)<4.7 && userFloat('JetID') && (userFloat('PUjetID') || pt>50)"),
                                   MuonPreselection = cms.string("userFloat('isGood') && userFloat('passCombRelIsoPFFSRCorr')"),
                                   ElectronPreselection = cms.string("userFloat('isGood') && userFloat('passCombRelIsoPFFSRCorr')"),
                                   DiBosonPreselection  = cms.string(""),
                                   MatchingType = cms.string("byDeltaR"),
                                   cleanFSRFromLeptons = cms.bool(True),
                                   DebugPlots = cms.untracked.bool(False),
                                   DebugPrintOuts = cms.untracked.bool(False)
                                   )


### ----------------------------------------------------------------------
### Missing ET
### ----------------------------------------------------------------------

#metTag = cms.InputTag("slimmedMETs")

# NB: b tag UPDATE DOES NOT WORK including this part related to MET => Updated info in jets get lost
if (RECORRECTMET and SAMPLE_TYPE == 2016):

#    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
#    runMetCorAndUncFromMiniAOD(process,
#                               isData=(not IsMC),
#    )
#    metTag = cms.InputTag("slimmedMETs","","ZZ")

    # NB: removed to use properly update jet collection to include DeepCSV in 2016 ntuples
    ### somehow MET recorrection gets this lost again...                          
    process.patJetsReapplyJEC.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']
    process.patJetsReapplyJEC.userData.userInts.src += ['pileupJetIdUpdated:fullId']

#[FIXME] Does not work in CMSSW_10_3_1 currently                                                                                                                                      
### Recorrect MET, cf. https://indico.cern.ch/event/759372/contributions/3149378/attachments/1721436/2779341/metreport.pdf slide 10                                        
###                and https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETUncertaintyPrescription#Instructions_for_9_4_X_X_9_for_2 
if (RECORRECTMET and SAMPLE_TYPE == 2017):                                                                                                                   
#                                                                                                                                                                          
#    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD                                                        
#                                                                                                                                                                         
#    runMetCorAndUncFromMiniAOD(process,                                                                                                                              
#                               isData=(not IsMC),                                                                                                                          
#                               fixEE2017 = True,                                                                                                                 
#                               fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139} ,                          
#                               postfix = "ModifiedMET"                                                                                                            
#                               )                                                                                                                                            
#    metTag = cms.InputTag("slimmedMETsModifiedMET","","ZZ")                                                                                                                 
#                                                                                                                                                                      
#    ### somehow MET recorrection gets this lost again...                                                                                                                 
    process.patJetsReapplyJEC.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']
    process.patJetsReapplyJEC.userData.userInts.src +=['pileupJetIdUpdated:fullId']                                                                              

if (RECORRECTMET and SAMPLE_TYPE == 2018):

#    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
#    runMetCorAndUncFromMiniAOD(process,
#                               isData=(not IsMC),
#                               )
#    metTag = cms.InputTag("slimmedMETs","","ZZ")

    ### somehow MET recorrection gets this lost again...                                                                                                             
    process.patJetsReapplyJEC.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']
    process.patJetsReapplyJEC.userData.userInts.src += ['pileupJetIdUpdated:fullId']



### ----------------------------------------------------------------------
### Paths
### ----------------------------------------------------------------------

process.preSkimCounter = cms.EDProducer("EventCountProducer")
process.PVfilter =  cms.Path(process.preSkimCounter+process.goodPrimaryVertices)

if APPLYJEC:
    process.Jets = cms.Path(process.pileupJetIdUpdated + process.patJetCorrFactorsReapplyJEC + process.patJetsReapplyJEC + process.QGTagger + process.dressedJets )
else:
    process.Jets = cms.Path( process.QGTagger + process.dressedJets + process.cleanJets)


# if (RECORRECTMET and SAMPLE_TYPE == 2016):
#     if IsMC:
#         process.MET = cms.Path(process.fullPatMetSequence)
#     else:
#         process.MET = cms.Path(process.fullPatMetSequence)

# #[FIXME] Does not work in CMSSW_10_3_1 currently                                                                                                         
# #if (RECORRECTMET and SAMPLE_TYPE == 2017):                                                                                                                       
# #    if IsMC:                                                                                                                                                   
# #        process.MET = cms.Path(process.fullPatMetSequenceModifiedMET)                                                                                              
# #    else:                                                                                                                                                                  
# #        process.MET = cms.Path(process.fullPatMetSequenceModifiedMET)                                                                                                     

# if (RECORRECTMET and SAMPLE_TYPE == 2017):
#     if IsMC:
#         process.MET = cms.Path(process.fullPatMetSequence)
#     else:
#         process.MET = cms.Path(process.fullPatMetSequence)

# if (RECORRECTMET and SAMPLE_TYPE == 2018):
#     if IsMC:
#         process.MET = cms.Path(process.fullPatMetSequence)
#     else:
#         process.MET = cms.Path(process.fullPatMetSequence)



### ----------------------------------------------------------------------
### Filters
### ----------------------------------------------------------------------
### Create filter for events with one candidate in the SR

# Prepare lepton collections
process.Candidates = cms.Path(
       process.muons             +
       process.electrons         + process.cleanSoftElectrons +
       process.taus              +
       process.fsrPhotons        + process.boostedFsrPhotons +
       process.appendPhotons     +
       process.softLeptons       +
       process.dressedJets       + process.cleanJets         +
       process.METSequence       +
       process.bareZCand         + process.ZCand
    )

# Optional sequence to build control regions. To get it, add
#process.CRPath = cms.Path(process.CRZl) # only trilepton
#OR
#process.CRPath = cms.Path(process.CR)   # trilep+4lep CRs

# process.CRZl = cms.Sequence(
#        process.bareZCand	 + process.ZCand	+
#        #process.bareZCand         + process.SVZCand       + process.ZCand     +
#        process.ZlCand
#    )



### Skim, triggers and MC filters (Only store filter result, no filter is applied)

# 2012 skim, Reimplementation by Giovanni
#process.load("HTauTauHMuMu.AnalysisStep.Skim2012_cfg")
#process.SkimSequence = cms.Sequence(process.HZZSkim2012)
#process.Skim = cms.Path(process.SkimSequence)

#SkimPaths = cms.vstring('Skim')
SkimPaths = cms.vstring('PVfilter') #Do not apply skim, just require a good PV

# process.HF = cms.Path(process.heavyflavorfilter)

# FIXME total kin filter?


if (ADDLOOSEELE) :
    import os
    execfile(os.environ['CMSSW_BASE'] + "/src/HTauTauHMuMu/AnalysisStep/test/MasterPy/LooseEle.py")
#    execfile(os.environ['CMSSW_BASE'] + "/src/HTauTauHMuMu/AnalysisStep/test/MasterPy/TracklessEle.py")
