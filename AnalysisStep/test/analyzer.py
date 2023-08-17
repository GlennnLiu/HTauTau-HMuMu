from HTauTauHMuMu.AnalysisStep.defaults import *
from HTauTauHMuMu.AnalysisStep.miscenums import *

### ----------------------------------------------------------------------
###
### Example analyzer
###
###----------------------------------------------------------------------

# Set defaults for variables used in this file (in case they are not defined by a caller script)
declareDefault("PD", "", globals()) # "" for MC, "DoubleEle", "DoubleMu", or "MuEG" for data
declareDefault("MCFILTER", "", globals())
declareDefault("XSEC", 1, globals())
declareDefault("GENXSEC", 1, globals())
declareDefault("GENBR", 1, globals())
declareDefault("PROCESS_CR", False, globals())
declareDefault("PROCESS_CRTT", False, globals())
declareDefault("PROCESS_CRWJ", False, globals())

declareDefault("DOMETRECOIL", False, globals())
#declareDefault("ADDZTREE", False, globals())

# LHE info
#  VVDECAYMODE\VVMODE  / ZZ==1 / WW==0  / Yukawa==2 / Zgam=3 / gamgam=4 / Z+nj=5
#                     0: 4l    / lnulnu / 2l        / 2l     / gam      / 2l
#                     1: 4q    / 4q     / 2q        / 2q     / -        / 2q
#                     2: 2l2q  / lnu2q  / -         / -      / -        / -
#                     3: 2l2nu / -      / -         / -      / -        / -
#                     4: 2q2nu / -      / -         / -      / -        / -
#                     5: 4nu   / -      / -         / 2nu    / -        / 2nu
#                    -1: [ Any                                                 ]
#                    -2: [ 2l2X         ]
#                    -3: [ 2nu2X        ]
#                    -4: [ 2q2X         ]
declareDefault("VVMODE", 1, globals())
declareDefault("VVDECAYMODE", 0, globals())
declareDefault("ADDLHEKINEMATICS", False, globals())

# K factors
declareDefault("APPLY_K_NNLOQCD_ZZGG", 0, globals()) # 0: Do not; 1: NNLO/LO; 2: NNLO/NLO; 3: NLO/LO
declareDefault("APPLY_K_NNLOQCD_ZZQQB", False, globals())
declareDefault("APPLY_K_NLOEW_ZZQQB", False, globals())

#failed events
declareDefault("SKIP_EMPTY_EVENTS", True, globals())
declareDefault("FAILED_TREE_LEVEL", 1, globals())

#ggF uncertainties for HTXS
declareDefault("APPLY_QCD_GGF_UNCERT", False, globals() )

declareDefault("OLDPATTRIGGER", False, globals() )

if FAILED_TREE_LEVEL and not SKIP_EMPTY_EVENTS:
    raise ValueError(
                     "Inconsistent options: FAILED_TREE_LEVEL={}, SKIP_EMPTY_EVENTS={}\n"
                     "If you want to write a failed tree, set SKIP_EMPTY_EVENTS=True"
                     .format(FAILED_TREE_LEVEL, SKIP_EMPTY_EVENTS)
                    )

# Get absolute path
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/HTauTauHMuMu/AnalysisStep/test/"

### ----------------------------------------------------------------------
### Standard sequence
### ----------------------------------------------------------------------

execfile(PyFilePath + "MasterPy/HTauTauHMuMu.py")


### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------
process.maxEvents.input = -1
#process.options.wantSummary = False


### ----------------------------------------------------------------------
### Output root file (monitoring histograms)
### ----------------------------------------------------------------------
process.TFileService=cms.Service('TFileService',
                                fileName=cms.string('HTauTauHMuMu.root')
                                )

### ----------------------------------------------------------------------
### Analyzers for Plots
### ----------------------------------------------------------------------

# All events together
# process.PlotsZZ    = cms.EDAnalyzer("ZZ4lAnalyzer",
#                                     channel = cms.untracked.string('ZZ'),
#                                     candCollection = cms.untracked.string('ZZCand'),
#                                     isMC = cms.untracked.bool(IsMC),
#                                     sampleType = cms.int32(SAMPLE_TYPE),
#                                     setup = cms.int32(LEPTON_SETUP),
#                                     skimPaths = cms.vstring(SkimPaths),
#                                     PD = cms.string(PD),
#                                     MCFilterPath = cms.string(MCFILTER),
#                                     sampleName = cms.string(SAMPLENAME),
#                                     dumpForSync = cms.untracked.bool(False),
#                                     )

### Control Region Plots

# PlotCRSetup    = cms.EDAnalyzer("ZZ4lAnalyzerCR",
#                                 channel = cms.untracked.string('aChannel'),
#                                 candCollection = cms.untracked.string('aCand'),
#                                 isMC = cms.untracked.bool(IsMC),
#                                 sampleType = cms.int32(SAMPLE_TYPE),
#                                 setup = cms.int32(LEPTON_SETUP),
#                                 skimPaths = cms.vstring(SkimPaths),
#                                 PD = cms.string(PD),
#                                 MCFilterPath = cms.string(MCFILTER),
#                                 )

# process.PlotsCRZLL = PlotCRSetup.clone()
# process.PlotsCRZLL.channel               = "ZLL"
# process.PlotsCRZLL.candCollection        = 'ZLLCand'


#Count events with at least 1 Z
# process.ZFiltered = cms.EDFilter("PATCompositeCandidateRefSelector",
#     src = cms.InputTag("ZCand"),
#     cut = cms.string("userFloat('GoodLeptons')")
# )
# process.sStep4 = cms.EDFilter("CandViewCountFilter",
#                               src = cms.InputTag("ZFiltered"),
#                               minNumber = cms.uint32(1)
#                               )
#process.step4 = cms.Path(process.SkimSequence + process.ZFiltered + process.sStep4 )


#Count events with at least 1 ZZ
# process.ZZFiltered = cms.EDFilter("PATCompositeCandidateRefSelector",
#     src = cms.InputTag("ZZCand"),
#     cut = cms.string("userFloat('GoodLeptons')")
# )
# process.sStep5 = cms.EDFilter("CandViewCountFilter",
#                                 src = cms.InputTag("ZZFiltered"),
#                                 minNumber = cms.uint32(1)
#                             )
#process.step5 = cms.Path(process.SkimSequence + process.ZZFiltered + process.sStep5 )


### ----------------------------------------------------------------------
### Analyzer for Trees
### ----------------------------------------------------------------------

TreeSetup = cms.EDAnalyzer("LLNtupleMaker",
                           channel = cms.untracked.string('aChannel'),
                           CandCollection = cms.untracked.string('ZCand'),
                           cut = cms.string(""),
                           fileName = cms.untracked.string('candTree'),
                           isMC = cms.untracked.bool(IsMC),
                           sampleType = cms.int32(SAMPLE_TYPE),
                           setup = cms.int32(LEPTON_SETUP),
                           skimPaths = cms.vstring(SkimPaths),
                           PD = cms.string(PD),
                           MCFilterPath = cms.string(MCFILTER),
                           
                           dataTag=cms.string(DATA_TAG), #added for recognizing UL16 pre/post VFP
                           
                           #for MET and SV fit
                           metSrc = cms.InputTag("ShiftMETcentral"),#srcMETTag,#metTag,
                           covSrc = cms.InputTag("METSignificance", "METCovariance"),
                           metJESup = cms.InputTag("ShiftMETforJES","up"),
                           metJESdown = cms.InputTag("ShiftMETforJES","down"),
#                            METdxUPTES = cms.InputTag("ShiftMETforTES", "METdxUP"),
#                            METdyUPTES = cms.InputTag("ShiftMETforTES", "METdyUP"),
#                            METdxDOWNTES = cms.InputTag("ShiftMETforTES", "METdxDOWN"),
#                            METdyDOWNTES = cms.InputTag("ShiftMETforTES", "METdyDOWN"),
#                            METdxUPEES = cms.InputTag("ShiftMETforEES", "METdxUPEES"),
#                            METdyUPEES = cms.InputTag("ShiftMETforEES", "METdyUPEES"),
#                            METdxDOWNEES = cms.InputTag("ShiftMETforEES", "METdxDOWNEES"),
#                            METdyDOWNEES = cms.InputTag("ShiftMETforEES", "METdyDOWNEES"),
#                            METdxUPMES = cms.InputTag("ShiftMETforMES", "METdxUPMES"),
#                            METdyUPMES = cms.InputTag("ShiftMETforMES", "METdyUPMES"),
#                            METdxDOWNMES = cms.InputTag("ShiftMETforMES", "METdxDOWNMES"),
#                            METdyDOWNMES = cms.InputTag("ShiftMETforMES", "METdyDOWNMES"),
#                            METdxUPJES = cms.InputTag("ShiftMETforJES", "METdxUPJES"),
#                            METdyUPJES = cms.InputTag("ShiftMETforJES", "METdyUPJES"),
#                            METdxDOWNJES = cms.InputTag("ShiftMETforJES", "METdxDOWNJES"),
#                            METdyDOWNJES = cms.InputTag("ShiftMETforJES", "METdyDOWNJES"),
#                            METdxUPJER = cms.InputTag("ShiftMETforJER", "METdxUPJER"),
#                            METdyUPJER = cms.InputTag("ShiftMETforJER", "METdyUPJER"),
#                            METdxDOWNJER = cms.InputTag("ShiftMETforJER", "METdxDOWNJER"),
#                            METdyDOWNJER = cms.InputTag("ShiftMETforJER", "METdyDOWNJER"),
                           
                           applyTrigger = cms.bool(APPLYTRIG), #Skip events failing required triggers. They are stored with sel<0 if set to false
                           applyTrigEff = cms.bool(False), #Add trigger efficiency as a weight, for samples where the trigger cannot be applied (obsoltete)
                           doMETRecoil = cms.bool(DOMETRECOIL),
                           skipEmptyEvents = cms.bool(SKIP_EMPTY_EVENTS),
                           failedTreeLevel = cms.int32(FAILED_TREE_LEVEL),
                           sampleName = cms.string(SAMPLENAME),
                           GenXSEC = cms.double(GENXSEC),
                           GenBR = cms.double(GENBR),

                           # Reco MEs to pick from the candidate
                           recoProbabilities = cms.vstring(),

                           # LHE info. parameters
                           lheProbabilities = cms.vstring(),
                           xsec = cms.double(XSEC),
                           VVMode = cms.int32(VVMODE),
                           VVDecayMode = cms.int32(VVDECAYMODE),
                           AddLHEKinematics = cms.bool(ADDLHEKINEMATICS),
                           Apply_K_NNLOQCD_ZZGG = cms.int32(APPLY_K_NNLOQCD_ZZGG),
                           Apply_K_NNLOQCD_ZZQQB = cms.bool(APPLY_K_NNLOQCD_ZZQQB),
                           Apply_K_NLOEW_ZZQQB = cms.bool(APPLY_K_NLOEW_ZZQQB),
                           Apply_QCD_GGF_UNCERT = cms.bool(APPLY_QCD_GGF_UNCERT),
                           )

process.SRTree = TreeSetup.clone()
process.SRTree.channel = 'SR'
process.SRTree.cut = cms.string('userFloat("OS") && userFloat("IsoLeptons")')

process.CRTTTree = TreeSetup.clone()
process.CRTTTree.channel = 'CRTT'
process.CRTTTree.cut = cms.string('userFloat("OS") && userFloat("IsoLepton1") && userFloat("NmedB")>=1')

process.CRQCDTree = TreeSetup.clone()
process.CRQCDTree.channel = 'CRQCD'
process.CRQCDTree.cut = cms.string('userFloat("SS") && userFloat("IsoLepton1") && userFloat("NmedB")<1')

process.CRWJTree = TreeSetup.clone()
process.CRWJTree.channel = 'CRWJ'
process.CRWJTree.cut = cms.string('userFloat("OS") && userFloat("IsoLepton1") && userFloat("NmedB")<1 && userFloat("MtLMET")>70')

process.CRQCDvSRTree = TreeSetup.clone()
process.CRQCDvSRTree.channel = 'CRQCDvSR'
process.CRQCDvSRTree.cut = cms.string('userFloat("SS") && !userFloat("IsoLepton1") && userFloat("NmedB")<1')

process.CRQCDvSROS1Tree = TreeSetup.clone()
process.CRQCDvSROS1Tree.channel = 'CRQCDvSROS1'
process.CRQCDvSROS1Tree.cut = cms.string('userFloat("OS") && !userFloat("IsoLepton1") && userFloat("IsoLepton2") && userFloat("NmedB")<1')

process.CRQCDvSROSTree = TreeSetup.clone()
process.CRQCDvSROSTree.channel = 'CRQCDvSROS'
process.CRQCDvSROSTree.cut = cms.string('userFloat("OS") && !userFloat("IsoLepton1") && !userFloat("IsoLepton2") && userFloat("NmedB")<1')

process.CRQCDvSRSSTree = TreeSetup.clone()
process.CRQCDvSRSSTree.channel = 'CRQCDvSRSS'
process.CRQCDvSRSSTree.cut = cms.string('userFloat("SS") && !userFloat("IsoLepton1") && !userFloat("IsoLepton2") && userFloat("NmedB")<1')

process.CRWJvSRTree = TreeSetup.clone()
process.CRWJvSRTree.channel = 'CRWJvSR'
process.CRWJvSRTree.cut = cms.string('userFloat("OS") && userFloat("IsoLepton1") && userFloat("NmedB")<1')

process.CRAPPOSTree = TreeSetup.clone()
process.CRAPPOSTree.channel = 'CRAPPOS'
process.CRAPPOSTree.cut = cms.string('userFloat("OS") && userFloat("IsoLepton1") && !userFloat("IsoLepton2") && userFloat("NmedB")<1')

process.CRAPPSSTree = TreeSetup.clone()
process.CRAPPSSTree.channel = 'CRAPPSS'
process.CRAPPSSTree.cut = cms.string('userFloat("SS") && userFloat("IsoLepton1") && !userFloat("IsoLepton2") && userFloat("NmedB")<1')

# ### ----------------------------------------------------------------------
# ### Z tree
# ### ----------------------------------------------------------------------
# process.ZTree = cms.EDAnalyzer("ZNtupleMaker",
#                                channel = cms.untracked.string('ZZ'),
#                                CandCollection = cms.untracked.string('ZCand'),
#                                fileName = cms.untracked.string('candTree'),
#                                isMC = cms.untracked.bool(IsMC),
#                                sampleType = cms.int32(SAMPLE_TYPE),
#                                setup = cms.int32(LEPTON_SETUP),
#                                skimPaths = cms.vstring(SkimPaths),
#                                PD = cms.string(PD),
#                                MCFilterPath = cms.string(MCFILTER),
#                                metSrc = srcMETTag,#metTag,
#                                skipEmptyEvents = cms.bool(True),
#                                sampleName = cms.string(SAMPLENAME),
#                                xsec = cms.double(XSEC),
#                                dataTag = cms.string(DATA_TAG)
#                                )



# Debug
#Define candidates to be dumped
process.ZZFiltered = cms.EDFilter("PATCompositeCandidateRefSelector",
                                  src = cms.InputTag("ZZCand"),
                                  cut = cms.string("userFloat('isBestCand')")
                                  )
### Select only events with one such candidate
process.ZZSelection= cms.EDFilter("CandViewCountFilter",
                                  src = cms.InputTag("ZZFiltered"),
                                  minNumber = cms.uint32(1)
                                  )

# ### Select CR events
# process.CRFiltered = cms.EDFilter("PATCompositeCandidateRefSelector",
#                                   src = cms.InputTag("ZLLCand"),
#                                   cut = cms.string("((userFloat('isBestCRZLLss')&&userFloat('CRZLLss')))||(userFloat('isBestCRZLLos_2P2F')&&userFloat('CRZLLos_2P2F'))||(userFloat('isBestCRZLLos_3P1F')&&userFloat('CRZLLos_3P1F'))")
#                                   )
# ### Select only events with one such candidate
# process.CRSelection= cms.EDFilter("CandViewCountFilter",
#                                   src = cms.InputTag("CRFiltered"),
#                                   minNumber = cms.uint32(1)
#                                   )


process.dumpUserData =  cms.EDAnalyzer("dumpUserData",
     dumpTrigger = cms.untracked.bool(True),
     muonSrcs =  cms.PSet(
        muons = cms.InputTag("appendPhotons:muons"),
     ),
     electronSrcs = cms.PSet(
        electrons = cms.InputTag("appendPhotons:electrons"),
     ),
     candidateSrcs = cms.PSet(
        Z     = cms.InputTag("ZCand"),
        # ZZ  = cms.InputTag("ZZCand"),
        # ZLL   =cms.InputTag("ZLLCand"),    # Starting point for all CRs
     ),
    jetSrc = cms.InputTag("cleanJets"),
)

process.trees = cms.EndPath(process.CRQCDvSROS1Tree)
# process.trees = cms.EndPath(process.SRTree)
# if PROCESS_CR:
#    process.trees += cms.Sequence( process.CRQCDTree + process.CRWJTree + process.CRQCDvSRTree + process.CRQCDvSROS1Tree + process.CRQCDvSROSTree + process.CRQCDvSRSSTree + process.CRAPPOSTree + process.CRAPPSSTree)
# if PROCESS_CRTT:
#    process.trees += process.CRTTTree
# if PROCESS_CRWJ:
#    process.trees += process.CRWJvSRTree