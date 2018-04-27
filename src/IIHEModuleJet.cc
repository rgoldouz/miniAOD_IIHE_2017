#include "UserCode/IIHETree/interface/IIHEModuleJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include <iostream>
#include <TMath.h>
#include <vector>

using namespace std ;
using namespace reco;
using namespace edm ;
//////////////////////////////////////////////////////////////////////////////////////////
////                                  Main IIHEJetModule
////////////////////////////////////////////////////////////////////////////////////////////

IIHEModuleJet::IIHEModuleJet(const edm::ParameterSet& iConfig, edm::ConsumesCollector && iC):IIHEModule(iConfig){
  pfJetLabel_                  =  iConfig.getParameter<edm::InputTag>("JetCollection");
  pfJetToken_                  =  iC.consumes<View<pat::Jet> > (pfJetLabel_);

  ETThreshold_ = iConfig.getUntrackedParameter<double>("jetPtThreshold") ;
  isMC_ = iConfig.getUntrackedParameter<bool>("isMC") ;

}
IIHEModuleJet::~IIHEModuleJet(){}

// ------------ method called once each job just before starting event loop  ------------
void IIHEModuleJet::beginJob(){
  addBranch("jet_n", kUInt);

  setBranchType(kVectorFloat);
  addBranch("jet_px");
  addBranch("jet_py");
  addBranch("jet_pz");
  addBranch("jet_pt");
  addBranch("jet_eta");
  addBranch("jet_theta");
  addBranch("jet_phi");
  addBranch("jet_energy");
  addBranch("jet_mass");
  addBranch("jet_chargedEmEnergyFraction");
  addBranch("jet_neutralHadronEnergyFraction");
  addBranch("jet_neutralEmEnergyFraction");
  addBranch("jet_chargedHadronEnergyFraction");
  addBranch("jet_muonEnergyFraction");
  setBranchType(kVectorInt);
  addBranch("jet_chargedMultiplicity");
  addBranch("jet_neutralMultiplicity");
  addBranch("jet_partonFlavour");
  addBranch("jet_hadronFlavour");

  setBranchType(kVectorFloat);
  addBranch("jet_CSVv2");
  addBranch("jet_CvsL");
  addBranch("jet_CvsB");
  addBranch("jet_MVA2BJets");
  setBranchType(kVectorInt);
  addBranch("jet_isJetIDLoose");
  addBranch("jet_isJetIDTight");
  addBranch("jet_isJetIDTightLepVeto");

  if (isMC_){

  setBranchType(kVectorFloat);
  addBranch("jet_Smeared_pt");
  addBranch("jet_Smeared_energy");
  addBranch("jet_SmearedJetResUp_pt");
  addBranch("jet_SmearedJetResUp_energy");
  addBranch("jet_SmearedJetResDown_pt");
  addBranch("jet_SmearedJetResDown_energy");
  addBranch("jet_EnUp_pt");
  addBranch("jet_EnUp_energy");
  addBranch("jet_EnDown_pt");
  addBranch("jet_EnDown_energy");

  addBranch("jet_BtagSF_loose");
  addBranch("jet_BtagSFbcUp_loose");
  addBranch("jet_BtagSFbcDown_loose");
  addBranch("jet_BtagSFudsgUp_loose");
  addBranch("jet_BtagSFudsgDown_loose");

  addBranch("jet_BtagSF_medium");
  addBranch("jet_BtagSFbcUp_medium");
  addBranch("jet_BtagSFbcDown_medium");
  addBranch("jet_BtagSFudsgUp_medium");
  addBranch("jet_BtagSFudsgDown_medium");

  addBranch("jet_BtagSF_tight");
  addBranch("jet_BtagSFbcUp_tight");
  addBranch("jet_BtagSFbcDown_tight");
  addBranch("jet_BtagSFudsgUp_tight");
  addBranch("jet_BtagSFudsgDown_tight");
  }

}

// ------------ method called to for each event  ------------
void IIHEModuleJet::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){


  edm::Handle<edm::View<pat::Jet> > pfJetHandle_;
  iEvent.getByToken(pfJetToken_, pfJetHandle_);

  const string ctr = "central";
  const string vup = "up";
  const string vdown = "down";


  for ( unsigned int i = 0; i <pfJetHandle_->size(); ++i) {
    Ptr<pat::Jet> pfjet = pfJetHandle_->ptrAt( i );

    if(pfjet->pt() < ETThreshold_) continue ;

    store("jet_px"    , pfjet->px()) ;
    store("jet_py"    , pfjet->py()) ;
    store("jet_pz"    , pfjet->pz()) ;
    store("jet_pt"    , pfjet->pt()) ;
    store("jet_eta"   , pfjet->eta()) ;
    store("jet_theta" , pfjet->theta()) ;
    store("jet_phi"   , pfjet->phi()) ;
    store("jet_energy", pfjet->energy()) ;
    store("jet_mass"  , pfjet->mass()) ;
    store("jet_neutralHadronEnergyFraction"         ,pfjet->neutralHadronEnergyFraction());
    store("jet_neutralEmEnergyFraction"             ,pfjet->neutralEmEnergyFraction());
    store("jet_chargedHadronEnergyFraction"         ,pfjet->chargedHadronEnergyFraction());
    store("jet_muonEnergyFraction"                  ,pfjet->muonEnergyFraction());
    store("jet_chargedEmEnergyFraction"             ,pfjet->chargedEmEnergyFraction());
    store("jet_chargedMultiplicity"                 ,pfjet->chargedMultiplicity());
    store("jet_neutralMultiplicity"                 ,pfjet->neutralMultiplicity());
    store("jet_partonFlavour"                       ,pfjet->partonFlavour());
    store("jet_hadronFlavour"                       ,pfjet->hadronFlavour());

    store("jet_CSVv2",pfjet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    store("jet_CvsL",pfjet->bDiscriminator("pfCombinedCvsLJetTags"));
    store("jet_CvsB",pfjet->bDiscriminator("pfCombinedCvsBJetTags"));
    store("jet_MVA2BJets", pfjet->bDiscriminator("pfCombinedMVAV2BJetTags"));
    float eta = pfjet->eta();
    float NHF = pfjet->neutralHadronEnergyFraction();
    float NEMF = pfjet->neutralEmEnergyFraction();
    float CHF = pfjet->chargedHadronEnergyFraction();
    float MUF = pfjet->muonEnergyFraction();
    float CEMF = pfjet->chargedEmEnergyFraction();
    int NumConst = pfjet->chargedMultiplicity()+pfjet->neutralMultiplicity();
    int NumNeutralParticles =pfjet->neutralMultiplicity();
    int CHM = pfjet->chargedMultiplicity(); 

    bool isJetIDLoose;
    bool isJetIDTight;
    bool isJetIDTightLepVeto;

    if( fabs(eta) <= 2.7){
      isJetIDLoose = ((NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4));
      isJetIDTight = ((NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(eta)>2.4));
      isJetIDTightLepVeto = ((NHF<0.90 && NEMF<0.90 && NumConst>1 && MUF<0.8) && ((abs(eta)<=2.4 && CHF>0 && CHM>0 && CEMF<0.90) || abs(eta)>2.4));
      }
    else if (abs(eta)>2.7 && abs(eta)<=3.0 ){
      isJetIDLoose = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2 );
      isJetIDTight = (NHF<0.98 && NEMF>0.01 && NumNeutralParticles>2 );
      isJetIDTightLepVeto = false;
      }
    else{
      isJetIDLoose = (NEMF<0.90 && NumNeutralParticles>10 );
      isJetIDTight = (NEMF<0.90 && NumNeutralParticles>10 ) ;
      isJetIDTightLepVeto = false;
    }
    store("jet_isJetIDLoose"                 ,int(isJetIDLoose));
    store("jet_isJetIDTight"                 ,int(isJetIDTight));
    store("jet_isJetIDTightLepVeto"          ,int(isJetIDTightLepVeto));
    if (isMC_){

//JetBTagWeight( edm::View<pat::Jet>&b, size_t ijet, const vector<BTagEntry::OperatingPoinconst string &bc_full_syst, const string &udsg_full_syst,const string &bc_full_syst, const string &udsg_full_syst,const string &bc_fast_syst, const string &udsg_fast_syst,bool do_deep_csv, bool do_by_proc, Runs runs)


   }

  }
}
void IIHEModuleJet::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup){}
void IIHEModuleJet::beginEvent(){}
void IIHEModuleJet::endEvent(){}


// ------------ method called once each job just after ending the event loop  ------------
void IIHEModuleJet::endJob(){
}

DEFINE_FWK_MODULE(IIHEModuleJet);
