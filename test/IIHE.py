# In order to run the code for MC/DATA on lxplus
#[cmsRun IIHE.py grid=False DataProcessing="mc2016" DataFormat="mc" dataset="RunIISummer16MiniAODv2" sample="ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1" address="MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000" file="0E4CD60A-0FC3-E611-BCB5-0CC47A7C3420.root"  ]
#root://eoscms//cms/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/36CDAE89-B3BE-E611-B022-0025905B8604.root

# In order to run the code for DATA on lxplus
#[cmsRun IIHE.py grid=False DataProcessing="rerecodata" DataFormat="data" dataset="Run2016D" sample="DoubleEG" address="MINIAOD/03Feb2017-v1/100000" file="CC7F39AC-F5EA-E611-A125-0CC47A13CDB0.root"  ]
#root://eoscms//cms/store/data/Run2015D/DoubleEG/MINIAOD/16Dec2015-v2/00000/F6E918C9-87A6-E511-B3D3-0CC47A4D76B2.root

import sys
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as opts
import copy
import os

options = opts.VarParsing ("analysis")
options.register("sample",
                 "ST_tW_top_5f_NoFullyHadronicDecays_13TeV-powheg_TuneCUETP8M1", 
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 "Sample to analyze")
options.register("address",
                 "MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000",
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 "address of sample in eos like: MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000")
options.register("file",
                 "0E4CD60A-0FC3-E611-BCB5-0CC47A7C3420.root",
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 "file to analyze")
options.register("DataProcessing",
                 "mc2017",
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 "Data processing types. Options are:mc2016,rerecodata,promptdata")
options.register("dataset",
                 "RunIISummer16MiniAODv2",
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 "datasets to analyze: RunIISummer16MiniAODv2, Run2016B-H")
options.register("DataFormat",
                 "mc",
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.string,
                 "Data format. Options are:mc,data")
options.register('grid',
                 True,
                 opts.VarParsing.multiplicity.singleton,
                 opts.VarParsing.varType.bool,
                 'If you run on grid or localy on eos')
options.parseArguments()

###########################################################################################
#os.system('wget https://github.com/rgoldouz/miniAOD_IIHE/archive/CMSSW_8_0_26_patch1.zip')
#os.system("unzip -a " + "CMSSW_8_0_26_patch1.zip" )
#os.system("mv " + "miniAOD_IIHE-CMSSW_8_0_26_patch1/data" + " ..")
#os.system("rm -rf CMSSW_8_0_26_patch1.zip miniAOD_IIHE-CMSSW_8_0_26_patch1")
##########################################################################################
#                                      Global tags                                       #
##########################################################################################
globalTag = "80"
if options.DataProcessing == "mc2017":
  globalTag = "94X_mc2017_realistic_v10"
if options.DataProcessing == "data2017":
  globalTag = "92X_dataRun2_Jun23ReReco_PixelCommissioning"
##########################################################################################
#                                  Start the sequences                                   #
##########################################################################################

process = cms.Process("IIHEAnalysis")

process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.load('Configuration.StandardSequences.Services_cff')

process.GlobalTag.globaltag = globalTag
print "Global Tag is ", process.GlobalTag.globaltag
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

##########################################################################################
#                                         Files                                          #
##########################################################################################
if options.DataFormat == "mc":
  path = "root://eoscms//eos/cms/store/"+ options.DataFormat + "/" + options.dataset  + "/" + options.sample + "/" + options.address + "/" + options.file

if options.DataFormat == "data":
  path = "root://eoscms//eos/cms/store/"+ options.DataFormat + "/" + options.dataset  + "/" + options.sample + "/" + options.address + "/" + options.file

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
#    eventsToProcess = cms.untracked.VEventRange('1:19792:3958249')
)
#process.source.fileNames.append( path )
#process.source.fileNames.append( "file:03Feb2017data.root" )
#process.source.fileNames.append( "file:TW_80_miniAOD.root" )
#process.source.fileNames.append( "file:2017data.root" )
process.source.fileNames.append( "file:MC2017_tW.root" )
###
filename_out = "outfile.root"
if options.DataFormat == "mc" and not options.grid:
#  filename_out = "file:/tmp/output_%s" % (options.sample + "_" + options.file)
  filename_out = "outfile_MC_"+options.file+".root"
if options.DataFormat == "data" and not options.grid:
  filename_out = "outfile_Data.root"

#filename_out = "outfile.root"
process.out = cms.OutputModule("PoolOutputModule", fileName = cms.untracked.string(filename_out) )
process.TFileService = cms.Service("TFileService", fileName = cms.string(filename_out) )


##########################################################################################
#                                   NEW tau id disc.                                     #
##########################################################################################

from RecoTauTag.RecoTau.TauDiscriminatorTools import noPrediscriminants
process.load('RecoTauTag.Configuration.loadRecoTauTagMVAsFromPrepDB_cfi')
from RecoTauTag.RecoTau.PATTauDiscriminationByMVAIsolationRun2_cff import *
from RecoTauTag.RecoTau.PATTauDiscriminationAgainstElectronMVA6_cfi import *

tauIdDiscrMVA_trainings_run2_2017 = {
  'tauIdMVAIsoDBoldDMwLT2017' : "tauIdMVAIsoDBoldDMwLT2017",
  }
tauIdDiscrMVA_WPs_run2_2017 = {
  'tauIdMVAIsoDBoldDMwLT2017' : {
    'Eff95' : "DBoldDMwLTEff95",
    'Eff90' : "DBoldDMwLTEff90",
    'Eff80' : "DBoldDMwLTEff80",
    'Eff70' : "DBoldDMwLTEff70",
    'Eff60' : "DBoldDMwLTEff60",
    'Eff50' : "DBoldDMwLTEff50",
    'Eff40' : "DBoldDMwLTEff40"
    }
  }
tauIdDiscrMVA_2017_version = "v1"

process.rerunDiscriminationAgainstElectronMVA6 = patTauDiscriminationAgainstElectronMVA6.clone(
   PATTauProducer = cms.InputTag('slimmedTaus'),
   Prediscriminants = noPrediscriminants,
   #Prediscriminants = requireLeadTrack,
   loadMVAfromDB = cms.bool(True),
   returnMVA = cms.bool(True),
   method = cms.string("BDTG"),
   mvaName_NoEleMatch_woGwoGSF_BL = cms.string("RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_BL"),
   mvaName_NoEleMatch_wGwoGSF_BL = cms.string("RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_BL"),
   mvaName_woGwGSF_BL = cms.string("RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_BL"),
   mvaName_wGwGSF_BL = cms.string("RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_BL"),
   mvaName_NoEleMatch_woGwoGSF_EC = cms.string("RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_woGwoGSF_EC"),
   mvaName_NoEleMatch_wGwoGSF_EC = cms.string("RecoTauTag_antiElectronMVA6v1_gbr_NoEleMatch_wGwoGSF_EC"),
   mvaName_woGwGSF_EC = cms.string("RecoTauTag_antiElectronMVA6v1_gbr_woGwGSF_EC"),
   mvaName_wGwGSF_EC = cms.string("RecoTauTag_antiElectronMVA6v1_gbr_wGwGSF_EC"),
   minMVANoEleMatchWOgWOgsfBL = cms.double(0.0),
   minMVANoEleMatchWgWOgsfBL  = cms.double(0.0),
   minMVAWOgWgsfBL            = cms.double(0.0),
   minMVAWgWgsfBL             = cms.double(0.0),
   minMVANoEleMatchWOgWOgsfEC = cms.double(0.0),
   minMVANoEleMatchWgWOgsfEC  = cms.double(0.0),
   minMVAWOgWgsfEC            = cms.double(0.0),
   minMVAWgWgsfEC             = cms.double(0.0),
   srcElectrons = cms.InputTag('slimmedElectrons'),
   usePhiAtEcalEntranceExtrapolation = cms.bool(True)
)

def loadMVA_WPs_run2_2017(process):

      for training, gbrForestName in tauIdDiscrMVA_trainings_run2_2017.items():
         process.loadRecoTauTagMVAsFromPrepDB.toGet.append(
            cms.PSet(
               record = cms.string('GBRWrapperRcd'),
               tag = cms.string("RecoTauTag_%s%s" % (gbrForestName, tauIdDiscrMVA_2017_version)),
               label = cms.untracked.string("RecoTauTag_%s%s" % (gbrForestName, tauIdDiscrMVA_2017_version))
            )
         )

         for WP in tauIdDiscrMVA_WPs_run2_2017[training].keys():
            process.loadRecoTauTagMVAsFromPrepDB.toGet.append(
               cms.PSet(
                  record = cms.string('PhysicsTGraphPayloadRcd'),
                  tag = cms.string("RecoTauTag_%s%s_WP%s" % (gbrForestName, tauIdDiscrMVA_2017_version, WP)),
                  label = cms.untracked.string("RecoTauTag_%s%s_WP%s" % (gbrForestName, tauIdDiscrMVA_2017_version, WP))
               )
            )

         process.loadRecoTauTagMVAsFromPrepDB.toGet.append(
            cms.PSet(
               record = cms.string('PhysicsTFormulaPayloadRcd'),
               tag = cms.string("RecoTauTag_%s%s_mvaOutput_normalization" % (gbrForestName, tauIdDiscrMVA_2017_version)),
               label = cms.untracked.string("RecoTauTag_%s%s_mvaOutput_normalization" % (gbrForestName, tauIdDiscrMVA_2017_version))
            )
)

loadMVA_WPs_run2_2017(process)
process.rerunDiscriminationByIsolationMVArun2v1raw = patDiscriminationByIsolationMVArun2v1raw.clone(
   PATTauProducer = cms.InputTag('slimmedTaus'),
   Prediscriminants = noPrediscriminants,
   loadMVAfromDB = cms.bool(True),
   mvaName = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1"),  # name of the training you want to use: training with 2017 MC_v1 for oldDM
   mvaOpt = cms.string("DBoldDMwLTwGJ"), # option you want to use for your training (i.e., which variables are used to compute the BDT score)
   requireDecayMode = cms.bool(True),
   verbosity = cms.int32(0)
)

process.rerunDiscriminationByIsolationMVArun2v1VLoose = patDiscriminationByIsolationMVArun2v1VLoose.clone(
   PATTauProducer = cms.InputTag('slimmedTaus'),
   Prediscriminants = noPrediscriminants,
   toMultiplex = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1raw'),
   key = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1raw:category'),
   loadMVAfromDB = cms.bool(True),
   mvaOutput_normalization = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_mvaOutput_normalization"), # normalization fo the training you want to use
   mapping = cms.VPSet(
      cms.PSet(
         category = cms.uint32(0),
         cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff90"), # this is the name of the working point you want to use
         variable = cms.string("pt"),
      )
   )
)

#loadMVA_WPs_run2_2017(process)
# here we produce all the other working points for the training
process.rerunDiscriminationByIsolationMVArun2v1VVLoose = process.rerunDiscriminationByIsolationMVArun2v1VLoose.clone()
process.rerunDiscriminationByIsolationMVArun2v1VVLoose.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff95")
process.rerunDiscriminationByIsolationMVArun2v1Loose = process.rerunDiscriminationByIsolationMVArun2v1VLoose.clone()
process.rerunDiscriminationByIsolationMVArun2v1Loose.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff80")
process.rerunDiscriminationByIsolationMVArun2v1Medium = process.rerunDiscriminationByIsolationMVArun2v1VLoose.clone()
process.rerunDiscriminationByIsolationMVArun2v1Medium.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff70")
process.rerunDiscriminationByIsolationMVArun2v1Tight = process.rerunDiscriminationByIsolationMVArun2v1VLoose.clone()
process.rerunDiscriminationByIsolationMVArun2v1Tight.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff60")
process.rerunDiscriminationByIsolationMVArun2v1VTight = process.rerunDiscriminationByIsolationMVArun2v1VLoose.clone()
process.rerunDiscriminationByIsolationMVArun2v1VTight.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff50")
process.rerunDiscriminationByIsolationMVArun2v1VVTight = process.rerunDiscriminationByIsolationMVArun2v1VLoose.clone()
process.rerunDiscriminationByIsolationMVArun2v1VVTight.mapping[0].cut = cms.string("RecoTauTag_tauIdMVAIsoDBoldDMwLT2017v1_WPEff40")

# this sequence has to be included in your cms.Path() before your analyzer which accesses the new variables is called.
process.rerunMvaIsolation2SeqRun2 = cms.Sequence(
   process.rerunDiscriminationByIsolationMVArun2v1raw
   *process.rerunDiscriminationByIsolationMVArun2v1VLoose
   *process.rerunDiscriminationByIsolationMVArun2v1VVLoose
   *process.rerunDiscriminationByIsolationMVArun2v1Loose
   *process.rerunDiscriminationByIsolationMVArun2v1Medium
   *process.rerunDiscriminationByIsolationMVArun2v1Tight
   *process.rerunDiscriminationByIsolationMVArun2v1VTight
   *process.rerunDiscriminationByIsolationMVArun2v1VVTight
)

# embed new id's into new tau collection
embedID = cms.EDProducer("PATTauIDEmbedder",
   src = cms.InputTag('slimmedTaus'),
   tauIDSources = cms.PSet(
      againstElectronMVA6RawNew = cms.InputTag('rerunDiscriminationAgainstElectronMVA6'),
      byIsolationMVArun2v1DBoldDMwLTrawNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1raw'),
      byVLooseIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1VLoose'),
      byVVLooseIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1VVLoose'),
      byLooseIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1Loose'),
      byMediumIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1Medium'),
      byTightIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1Tight'),
      byVTightIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1VTight'),
      byVVTightIsolationMVArun2v1DBoldDMwLTNew = cms.InputTag('rerunDiscriminationByIsolationMVArun2v1VVTight'),
      ),
   )

setattr(process, "NewTauIDsEmbedded", embedID)


##########################################################################################
#                                   IIHETree options                                     #
##########################################################################################

#Track isolation correction value for HEEP v7
#process.load("RecoEgamma.ElectronIdentification.heepIdVarValueMapProducer_cfi")
#EGamma VID for various working points
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
my_id_modules = ["RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff",
                 "RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff",
                 "RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff"]
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)



#Electron energy scale and smearing
process.load("Configuration.StandardSequences.Services_cff")
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                    calibratedPatElectrons  = cms.PSet( initialSeed = cms.untracked.uint32(81),
                                                    engineName = cms.untracked.string("TRandom3")),
                                                  )
process.load("EgammaAnalysis.ElectronTools.calibratedPatElectronsRun2_cfi")
process.calibratedPatElectrons.isMC = cms.bool("mc" in options.DataProcessing)
process.calibratedPatElectrons.electrons = cms.InputTag("slimmedElectrons")
process.calibratedPatElectrons.isSynchronization = cms.bool(False)
##########################################################################################
#                            MY analysis input!                              ####
##########################################################################################
process.load("UserCode.IIHETree.IIHETree_cfi")
process.IIHEAnalysis.globalTag = cms.string(globalTag)
process.IIHEAnalysis.isData  = cms.untracked.bool("data" in options.DataProcessing)
process.IIHEAnalysis.isMC    = cms.untracked.bool("mc" in options.DataProcessing)
#****Collections added before the analysis
# VID output
process.IIHEAnalysis.eleTrkPtIsoLabel                            = cms.InputTag("heepIDVarValueMaps"    ,"eleTrkPtIso"       ,"IIHEAnalysis" )
process.IIHEAnalysis.VIDVeto                                     = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto"  )
process.IIHEAnalysis.VIDLoose                                    = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose" )
process.IIHEAnalysis.VIDMedium                                   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium")
process.IIHEAnalysis.VIDTight                                    = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight" )
process.IIHEAnalysis.eleMediumIdMap                             = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90" )
process.IIHEAnalysis.eleTightIdMap                             = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80" )
process.IIHEAnalysis.mvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values")
process.IIHEAnalysis.mvaCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories")

process.IIHEAnalysis.VIDHEEP7                                    = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV70"                   )
# Collections for DATA only.
process.IIHEAnalysis.particleFlowEGammaGSFixedCollection         = cms.InputTag("particleFlowEGammaGSFixed", "dupECALClusters"              )
process.IIHEAnalysis.ecalMultiAndGSGlobalRecHitEBCollection      = cms.InputTag("ecalMultiAndGSGlobalRecHitEB","dupESClusters"        ,"PAT")
process.IIHEAnalysis.METsMuEGCleanCollection                     = cms.InputTag("slimmedMETsMuEGClean"                                      )
process.IIHEAnalysis.discardedMuonCollection                     = cms.InputTag("packedPFCandidatesDiscarded"                               )

#use 80 regression + scale/smearing for electron
process.IIHEAnalysis.calibratedElectronCollection    = cms.InputTag("calibratedPatElectrons","","IIHEAnalysis")

#jet smeared collection
#process.IIHEAnalysis.JetCollection                   = cms.InputTag("basicJetsForMet" ,"","IIHEAnalysis")
#process.IIHEAnalysis.JetCollectionSmeared            = cms.InputTag("patSmearedJets"               ,"","IIHEAnalysis")
#process.IIHEAnalysis.JetCollectionEnUp               = cms.InputTag("shiftedPatJetEnUp"            ,"","IIHEAnalysis")
#process.IIHEAnalysis.JetCollectionEnDown             = cms.InputTag("shiftedPatJetEnDown"          ,"","IIHEAnalysis")
#process.IIHEAnalysis.JetCollectionSmearedJetResUp    = cms.InputTag("shiftedPatSmearedJetResUp"    ,"","IIHEAnalysis")
#process.IIHEAnalysis.JetCollectionSmearedJetResDown  = cms.InputTag("shiftedPatSmearedJetResDown"  ,"","IIHEAnalysis")

#MET collections
process.IIHEAnalysis.patPFMetCollection                        = cms.InputTag("patPFMet"                  , ""                ,"IIHEAnalysis"  )
process.IIHEAnalysis.patPFMetT1Collection                      = cms.InputTag("patPFMetT1"                , ""                ,"IIHEAnalysis"  )
process.IIHEAnalysis.patPFMetT1JetEnDownCollection             = cms.InputTag("patPFMetT1JetEnDown"       , ""                ,"IIHEAnalysis"  )
process.IIHEAnalysis.patPFMetT1JetEnUpCollection               = cms.InputTag("patPFMetT1JetEnUp"         , ""                ,"IIHEAnalysis"  )
process.IIHEAnalysis.patPFMetT1SmearCollection                 = cms.InputTag("patPFMetT1Smear"           , ""                ,"IIHEAnalysis"  )
process.IIHEAnalysis.patPFMetT1SmearJetEnDownCollection        = cms.InputTag("patPFMetT1SmearJetEnDown"  , ""                ,"IIHEAnalysis"  )
process.IIHEAnalysis.patPFMetT1SmearJetEnUpCollection          = cms.InputTag("patPFMetT1SmearJetEnUp"    , ""                ,"IIHEAnalysis"  )
process.IIHEAnalysis.patPFMetT1SmearJetResDownCollection       = cms.InputTag("patPFMetT1SmearJetResDown" , ""                ,"IIHEAnalysis"  )
process.IIHEAnalysis.patPFMetT1SmearJetResUpCollection         = cms.InputTag("patPFMetT1SmearJetResUp"   , ""                ,"IIHEAnalysis"  )
process.IIHEAnalysis.patPFMetT1TxyCollection                   = cms.InputTag("patPFMetT1Txy"             , ""                ,"IIHEAnalysis"  )
process.IIHEAnalysis.patPFMetFinalCollection                   = cms.InputTag("slimmedMETs"               , ""                ,"IIHEAnalysis"  )

#particle level fiducial collections
process.IIHEAnalysis.particleLevelJetsCollection               = cms.InputTag("pseudoTop"       , "jets"     ,"IIHEAnalysis"  )
process.IIHEAnalysis.particleLevelak1DressedLeptonCollection   = cms.InputTag("pseudoTop"       , "leptons"  ,"IIHEAnalysis"  )
process.IIHEAnalysis.particleLevelMETCollection                = cms.InputTag("pseudoTop"       , "mets","IIHEAnalysis"  )



process.IIHEAnalysis.includeLeptonsAcceptModule  = cms.untracked.bool(True)
process.IIHEAnalysis.includeTriggerModule        = cms.untracked.bool(True)
process.IIHEAnalysis.includeEventModule          = cms.untracked.bool(True)
process.IIHEAnalysis.includeVertexModule         = cms.untracked.bool(True)
process.IIHEAnalysis.includeElectronModule       = cms.untracked.bool(True)
process.IIHEAnalysis.includeMuonModule           = cms.untracked.bool(True)
process.IIHEAnalysis.includeMETModule            = cms.untracked.bool(True)
process.IIHEAnalysis.includeJetModule            = cms.untracked.bool(True)
process.IIHEAnalysis.includeTauModule            = cms.untracked.bool(True)
process.IIHEAnalysis.includeL1Module            = cms.untracked.bool(True)
process.IIHEAnalysis.includeMCTruthModule        = cms.untracked.bool("mc" in options.DataProcessing)
process.IIHEAnalysis.includeLHEWeightModule        = cms.untracked.bool("mc" in options.DataProcessing)
#process.IIHEAnalysis.includeDataModule            = cms.untracked.bool("data" in options.DataProcessing)


process.IIHEAnalysis.includeAutoAcceptEventModule                = cms.untracked.bool(True)
##########################################################################################
#                            Woohoo!  We"re ready to start!                              #
##########################################################################################
#process.p1 = cms.Path(process.kt6PFJetsForIsolation+process.IIHEAnalysis)
#process.out = cms.OutputModule(
#    "PoolOutputModule",
#    fileName = cms.untracked.string("EDM.root")
#    )
process.p1 = cms.Path(
    process.rerunMvaIsolation2SeqRun2 *
    process.rerunDiscriminationAgainstElectronMVA6 *
    getattr(process, "NewTauIDsEmbedded") *
    process.calibratedPatElectrons *
    process.egmGsfElectronIDSequence * 
    process.IIHEAnalysis
    )

#process.outpath = cms.EndPath(process.out)

