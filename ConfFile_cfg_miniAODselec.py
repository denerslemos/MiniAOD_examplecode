import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process('Demo',eras.Run2_2018_pp_on_AA)

process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')

# Limit the output messages
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool( True ),
)

process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
                'file:/eos/cms/store/group/phys_heavyions/caber/MiniAODValidation/MiniAODfiles/step2_HIMiniAOD.root',
                ),
                skipBadFiles=cms.untracked.bool(True),
				duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            )

#define number of events to be processed (-1 means all)
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Set the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '112X_dataRun2_PromptLike_HI_v3', '') #centrality table is included here!

#-> This is needed for EDAnalyzer{

# Add PbPb collision event selection 
process.load('HeavyIonsAnalysis.EventAnalysis.skimanalysis_cfi')
process.load('HeavyIonsAnalysis.EventAnalysis.collisionEventSelection_cff')

# Define the event selection sequence
process.eventFilter_HM = cms.Sequence(
    process.phfCoincFilter2Th4 *
    process.primaryVertexFilter *
    process.clusterCompatibilityFilter
)

#}

# Define the output
process.TFileService = cms.Service("TFileService",fileName = cms.string('PbPb_trk_miniAOD_withevselec.root'))

process.demo = cms.EDAnalyzer('DemoAnalyzer',
                        vertexCollection  = cms.InputTag("offlineSlimmedPrimaryVerticesRecovery"),#"offlinePrimaryVertices"),
                        tracks   = cms.InputTag("packedPFCandidates"),
                        chi2src = cms.InputTag("packedPFCandidateTrackChi2"),
			HFfilters = cms.InputTag("hiHFfilters","hiHFfilters","PAT"),
                        CentralitySrc    = cms.InputTag("hiCentrality","","reRECO"),
                        CentralityBinSrc = cms.InputTag("centralityBin","HFtowers")
)

process.p = cms.Path(process.eventFilter_HM*process.demo) # with selection
process.schedule = cms.Schedule(process.p)

                                                                                  
