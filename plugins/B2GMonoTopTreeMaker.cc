// -*- C++ -*-
//
// Package:    Analysis/B2GMonoTopTreeMaker
// Class:      B2GMonoTopTreeMaker
// 
/**\class B2GMonoTopTreeMaker B2GMonoTopTreeMaker.cc Analysis/B2GMonoTopTreeMaker/plugins/B2GMonoTopTreeMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  James Dolen
// Modified By:      Ryan Mueller
//         Created:  Sat, 30 Apr 2016 17:40:42 GMT
//
//
// system include files
#include <memory>
#include <iostream>    
#include <algorithm>   
#include <bitset>   

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// DataFormats
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

// TFileService
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// Gen particle
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// JEC
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

// JER
#include "JetMETCorrections/Modules/interface/JetResolution.h"

// Electron
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/PatCandidates/interface/VIDCutFlowResult.h"

// Trigger
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

// Vertex
#include "DataFormats/VertexReco/interface/Vertex.h"

// Pileup
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

// LHE weights
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Utilities
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

// root
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <TRandom3.h>

#include "Analysis/B2GMonoTop/interface/B2GMonoTopTreeMaker.h"
#include "Analysis/B2GMonoTop/interface/eventDataStruct.h"
#include "Analysis/B2GMonoTop/interface/eventTTree.h"

#ifdef __CINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif

//RS gluon PDF weights
namespace LHAPDF {
  void initPDFSet(int nset, const std::string& filename, int member=0);
  int numberPDF(int nset);
  void usePDFMember(int nset, int member);
  double xfx(int nset, double x, double Q, int fl);
  double getXmin(int nset, int member);
  double getXmax(int nset, int member);
  double getQ2min(int nset, int member);
  double getQ2max(int nset, int member);
  void extrapolate(bool extrapolate=true);
}

//  888     888                                8888888888                         888    d8b                            
//  888     888                                888                                888    Y8P                            
//  888     888                                888                                888                                   
//  888     888 .d8888b   .d88b.  888d888      8888888 888  888 88888b.   .d8888b 888888 888  .d88b.  88888b.  .d8888b  
//  888     888 88K      d8P  Y8b 888P"        888     888  888 888 "88b d88P"    888    888 d88""88b 888 "88b 88K      
//  888     888 "Y8888b. 88888888 888          888     888  888 888  888 888      888    888 888  888 888  888 "Y8888b. 
//  Y88b. .d88P      X88 Y8b.     888          888     Y88b 888 888  888 Y88b.    Y88b.  888 Y88..88P 888  888      X88 
//   "Y88888P"   88888P'  "Y8888  888          888      "Y88888 888  888  "Y8888P  "Y888 888  "Y88P"  888  888  88888P' 
//                                                                                                                      
//                                                                                                                      
//                                                                                                                      

void setVector0( Float_t vec[4]){
    vec[0] = 0;
    vec[1] = 0;
    vec[2] = 0;
    vec[3] = 0;
    return;
}

void setVectorTL( Float_t vec[4], TLorentzVector  TL4){
    vec[0] = TL4.Pt();
    vec[1] = TL4.Eta();
    vec[2] = TL4.Phi();
    vec[3] = TL4.M();
    return;
}





//   .d8888b.                                                         888      8888888b.                    888            
//  d88P  Y88b                                                        888      888  "Y88b                   888            
//  888    888                                                        888      888    888                   888            
//  888         .d88b.  88888b.  .d8888b        8888b.  88888b.   .d88888      888    888  .d88b.  .d8888b  888888 888d888 
//  888        d88""88b 888 "88b 88K               "88b 888 "88b d88" 888      888    888 d8P  Y8b 88K      888    888P"   
//  888    888 888  888 888  888 "Y8888b.      .d888888 888  888 888  888      888    888 88888888 "Y8888b. 888    888     
//  Y88b  d88P Y88..88P 888  888      X88      888  888 888  888 Y88b 888      888  .d88P Y8b.          X88 Y88b.  888     
//   "Y8888P"   "Y88P"  888  888  88888P'      "Y888888 888  888  "Y88888      8888888P"   "Y8888   88888P'  "Y888 888     
//                                                                                                                         
//                                                                                                                         
//      

B2GMonoTopTreeMaker::B2GMonoTopTreeMaker(const edm::ParameterSet& iConfig):
    ak4jetToken_(consumes<pat::JetCollection>(    edm::InputTag("slimmedJets"))),
    ak8jetToken_(consumes<pat::JetCollection>(    iConfig.getParameter<edm::InputTag>("ak8chsInput"))),  //edm::InputTag("slimmedJetsAK8"))),
    puppijetToken_(consumes<pat::JetCollection>(  iConfig.getParameter<edm::InputTag>("ak8puppiInput"))),
    ca12puppijetToken_(consumes<pat::JetCollection>(  iConfig.getParameter<edm::InputTag>("ca12puppiInput"))),
    ak8CHSSoftDropSubjetsToken_(consumes<pat::JetCollection>(  iConfig.getParameter<edm::InputTag>("ak8chsSubjetsInput"))),
    ak8PuppiSoftDropSubjetsToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("ak8puppiSubjetsInput"))),
    ca12PuppiSoftDropSubjetsToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("ca12puppiSubjetsInput"))),
    ak4genjetToken_(consumes<reco::GenJetCollection>(edm::InputTag("slimmedGenJets"))),
    ak8genjetToken_(consumes<reco::GenJetCollection>(edm::InputTag("slimmedGenJetsAK8"))),
    prunedGenToken_(consumes<edm::View<reco::GenParticle> >(edm::InputTag("prunedGenParticles"))),
    rhoToken_(consumes<double>(edm::InputTag("fixedGridRhoFastjetAll"))),
    vtxToken_(consumes<std::vector<reco::Vertex> >(edm::InputTag("offlineSlimmedPrimaryVertices"))),
    triggerResultsMETFilterToken_(consumes<edm::TriggerResults>(edm::InputTag("TriggerResults", "", "RECO"))),  //"PAT"
    triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerBits"))),//"TriggerResults", "", "HLT2"))),
    triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(edm::InputTag("patTrigger"))),
    triggerPrescales_(consumes<pat::PackedTriggerPrescales>(edm::InputTag("patTrigger"))),  //   selectedPatTrigger))),
    badMuonFilterToken_(consumes<bool>(edm::InputTag("BadPFMuonFilter",""))),
    badChargedCandidateFilterToken_(consumes<bool>(edm::InputTag("BadChargedCandidateFilter",""))),
    muonToken_(consumes<pat::MuonCollection>(edm::InputTag("slimmedMuons"))),
    electronToken_(consumes<edm::View<pat::Electron>>(edm::InputTag("slimmedElectrons"))), //Collection
    metToken_(consumes<pat::METCollection>(edm::InputTag("slimmedMETs","","Ana"))),
    pileupInfoToken_(consumes<std::vector<PileupSummaryInfo>>(edm::InputTag("slimmedAddPileupInfo"))),
    lheSrc_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheSrc"))),
    pdfToken_(consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
    beamSpotToken_(consumes<reco::BeamSpot>(edm::InputTag("offlineBeamSpot"))),
    conversionsToken_(consumes<reco::ConversionCollection>(edm::InputTag("reducedEgamma:reducedConversions"))),
    eleIdFullInfoMapToken_HLTpre_  (consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleIdFullInfoMapToken_HLTpre"))),
    eleIdFullInfoMapToken_Loose_   (consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleIdFullInfoMapToken_Loose"))),
    eleIdFullInfoMapToken_Medium_  (consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleIdFullInfoMapToken_Medium"))),
    eleIdFullInfoMapToken_Tight_   (consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleIdFullInfoMapToken_Tight"))),
    eleIdFullInfoMapToken_HEEP_    (consumes<edm::ValueMap<vid::CutFlowResult> >(iConfig.getParameter<edm::InputTag>("eleIdFullInfoMapToken_HEEP"))),
    verbose_         (iConfig.getParameter<bool>  ("verbose")),
    verboseGen_      (iConfig.getParameter<bool>  ("verboseGen")),
    useToolbox_      (iConfig.getParameter<bool>  ("useToolbox")),
    runGenLoop_      (iConfig.getParameter<bool>  ("runGenLoop")),
    runHadTree_      (iConfig.getParameter<bool>  ("runHadTree")),
    runLeptTree_     (iConfig.getParameter<bool>  ("runLeptTree")),
    isZprime_        (iConfig.getParameter<bool>  ("isZprime")),
    isMonoTop_         (iConfig.getParameter<bool>  ("isMonoTop")),
    isRSG_           (iConfig.getParameter<bool>  ("isRSG")),
    isRun2016F_      (iConfig.getParameter<bool>  ("isRun2016F")),
    jecPayloadsAK4chs_ (iConfig.getParameter<std::vector<std::string> >  ("jecPayloadsAK4chs")),
    jecPayloadsAK8chs_ (iConfig.getParameter<std::vector<std::string> >  ("jecPayloadsAK8chs")),
    jecPayloadsAK4pup_ (iConfig.getParameter<std::vector<std::string> >  ("jecPayloadsAK4pup")),
    jecPayloadsAK8pup_ (iConfig.getParameter<std::vector<std::string> >  ("jecPayloadsAK8pup")),
    jecPayloadsAK4chsSecondary_ (iConfig.getParameter<std::vector<std::string> >  ("jecPayloadsAK4chsSecondary")),
    jecPayloadsAK8chsSecondary_ (iConfig.getParameter<std::vector<std::string> >  ("jecPayloadsAK8chsSecondary")),
    jecPayloadsAK4pupSecondary_ (iConfig.getParameter<std::vector<std::string> >  ("jecPayloadsAK4pupSecondary")),
    jecPayloadsAK8pupSecondary_ (iConfig.getParameter<std::vector<std::string> >  ("jecPayloadsAK8pupSecondary")),
    jertextAK4_ (iConfig.getParameter<std::string>  ("jertextAK4")),
    jertextAK8_ (iConfig.getParameter<std::string>  ("jertextAK8")),
    jerSFtext_  (iConfig.getParameter<std::string>  ("jerSFtext"))
{
  std::cout<<"B2GMonoTopTreeMaker::B2GMonoTopTreeMaker"<<std::endl;

  //RS gluon PDF weights
  LHAPDF::initPDFSet(1, "NNPDF30_lo_as_0130");

  usesResource("TFileService");

  edm::Service<TFileService> fs;

  h_cutflow_had                      =  fs->make<TH1D>("h_cutflow_had"                     ,"",20,0,20);
  h_cutflow_lept                     =  fs->make<TH1D>("h_cutflow_lept"                    ,"",20,0,20);

  h_trigger_efficency_1               = fs->make<TH1D>("h_trigger_efficency_1", "trigger efficency vs MET btag + mu", 100, 0, 1200);
  h_trigger_efficency_2               = fs->make<TH1D>("h_trigger_efficency_2", "trigger efficency vs MET pfmet300", 100, 0, 1200);
  h_trigger_efficency_3               = fs->make<TH1D>("h_trigger_efficency_3", "trigger efficency vs MET 120met + 120 MHT", 100, 0, 1200);
  h_trigger_efficency_4               = fs->make<TH1D>("h_trigger_efficency_4", "trigger efficency vs MET 170 ", 100, 0, 1200);

  h_trigger_accept_1               = fs->make<TH1D>("h_trigger_accept_1", "trigger accept vs MET btag + mu", 100, 0, 1200);
  h_trigger_accept_2               = fs->make<TH1D>("h_trigger_accept_2", "trigger accept vs MET pfmet300", 100, 0, 1200);
  h_trigger_accept_3               = fs->make<TH1D>("h_trigger_accept_3", "trigger accept vs MET 120met + 120 MHT", 100, 0, 1200);
  h_trigger_accept_4               = fs->make<TH1D>("h_trigger_accept_4", "trigger accept vs MET 170", 100, 0, 1200);

  h_trigger_reject_1               = fs->make<TH1D>("h_trigger_reject_1", "trigger reject valuesGen Met", 100, 0, 1200);
  h_trigger_reject_2               = fs->make<TH1D>("h_trigger_reject_2", "trigger reject valuesGen Met", 100, 0, 1200);
  h_trigger_reject_3               = fs->make<TH1D>("h_trigger_reject_3", "trigger reject valuesGen Met", 100, 0, 1200);
  h_trigger_reject_4               = fs->make<TH1D>("h_trigger_reject_4", "trigger reject valuesGen Met", 100, 0, 1200);

  h_trigger_efficency_1_topPt               = fs->make<TH1D>("h_trigger_efficency_1_topPt", "trigger efficency vs Gen Top p_{T} btag + mu", 100, 0, 1200);
  h_trigger_efficency_2_topPt               = fs->make<TH1D>("h_trigger_efficency_2_topPt", "trigger efficency vs Gen Top p_{T} pfmet300", 100, 0, 1200);
  h_trigger_efficency_3_topPt               = fs->make<TH1D>("h_trigger_efficency_3_topPt", "trigger efficency vs Gen Top p_{T} 120met + 120 MHT", 100, 0, 1200);
  h_trigger_efficency_4_topPt               = fs->make<TH1D>("h_trigger_efficency_4_topPt", "trigger efficency vs Gen Top p_{T} 120met + 120 MHT", 100, 0, 1200);
  h_trigger_accept_1_topPt               = fs->make<TH1D>("h_trigger_accept_1_topPt", "trigger accept vs Gen Top p_{T} btag + mu", 100, 0, 1200);
  h_trigger_accept_2_topPt               = fs->make<TH1D>("h_trigger_accept_2_topPt", "trigger accept vs Gen Top p_{T} pfmet300", 100, 0, 1200);
  h_trigger_accept_3_topPt               = fs->make<TH1D>("h_trigger_accept_3_topPt", "trigger accept vs Gen Top p_{T} 120met + 120 MHT", 100, 0, 1200);
  h_trigger_accept_4_topPt               = fs->make<TH1D>("h_trigger_accept_4_topPt", "trigger accept vs Gen Top p_{T} 120met + 120 MHT", 100, 0, 1200);
  h_trigger_reject_1_topPt               = fs->make<TH1D>("h_trigger_reject_1_topPt", "trigger reject valuesGen MTop p_{T}", 100, 0, 1200);
  h_trigger_reject_2_topPt               = fs->make<TH1D>("h_trigger_reject_2_topPt", "trigger reject valuesGen MTop p_{T}", 100, 0, 1200);
  h_trigger_reject_3_topPt               = fs->make<TH1D>("h_trigger_reject_3_topPt", "trigger reject valuesGen MTop p_{T}", 100, 0, 1200);
  h_trigger_reject_4_topPt               = fs->make<TH1D>("h_trigger_reject_4_topPt", "trigger reject valuesGen MTop p_{T}", 100, 0, 1200);
  h_NtrueIntPU                       =  fs->make<TH1D>("h_NtrueIntPU"                      ,"",200,0,200);
  h_NPV                              =  fs->make<TH1D>("h_NPV"                             ,"",200,0,200);
  h_NPVreweighted                    =  fs->make<TH1D>("h_NPVreweighted"                   ,"",200,0,200);
  h_NPVgood                          =  fs->make<TH1D>("h_NPVgood"                         ,"",200,0,200);
  h_NPVgoodreweighted                =  fs->make<TH1D>("h_NPVgoodreweighted"               ,"",200,0,200);

  trigsToRunMu.push_back("HLT_PFMET120_NoiseCleaned_BTagCSV07_v");
  trigsToRunMu.push_back("HLT_PFMET120_BTagCSV_p067_v");
  trigsToRunMu.push_back("HLT_PFMET120_Mu5_v");

  trigsToRunHad.push_back("HLT_PFMET120_NoiseCleaned_BTagCSV07_v");
  trigsToRunHad.push_back("HLT_PFMET300_v");
  trigsToRunHad.push_back("HLT_PFMET400_v");
  trigsToRunHad.push_back("HLT_PFMET500_v");
  trigsToRunHad.push_back("HLT_PFMET600_v");

  trigsToRun.push_back("HLT_PFMET120_NoiseCleaned_BTagCSV07_v");
  trigsToRun.push_back("HLT_PFMET120_BTagCSV_p067_v");
  trigsToRun.push_back("HLT_PFMET120_Mu5_v");
  trigsToRun.push_back("HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV_p13_v");
  trigsToRun.push_back("HLT_PFHT300_v");
  trigsToRun.push_back("HLT_PFHT350_v");
  trigsToRun.push_back("HLT_PFHT400_v");
  trigsToRun.push_back("HLT_PFHT475_v");
  trigsToRun.push_back("HLT_PFHT600_v");
  trigsToRun.push_back("HLT_PFHT650_v");
  trigsToRun.push_back("HiLT_PFHT800_v");
  trigsToRun.push_back("HLT_PFHT900_v");
  trigsToRun.push_back("HLT_PFHT650_WideJetMJJ900"); //HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v6
  trigsToRun.push_back("HLT_PFHT650_WideJetMJJ950"); //HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v6
//
  //// Single jet
  trigsToRun.push_back("HLT_CaloJet500_NoJetID_v");
  trigsToRun.push_back("HLT_PFJet320_v");
  trigsToRun.push_back("HLT_PFJet400_v");
  trigsToRun.push_back("HLT_PFJet450_v");
  trigsToRun.push_back("HLT_PFJet500_v");
  trigsToRun.push_back("HLT_AK8PFJet450_v");
  trigsToRun.push_back("HLT_AK8PFJet500_v");
//
  //// Substructure
  trigsToRun.push_back("HLT_AK8PFJet360_TrimMass30_v");
  trigsToRun.push_back("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_v");
  trigsToRun.push_back("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v");
//
  //// Substructure + b-tag
  trigsToRun.push_back("HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_v");
  trigsToRun.push_back("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v");
//
  //// Muon
  trigsToRun.push_back("HLT_Mu45_eta2p1_v");
  trigsToRun.push_back("HLT_Mu50_v");
  trigsToRun.push_back("HLT_Mu55_v");
  trigsToRun.push_back("HLT_TkMu50_v");
  trigsToRun.push_back("HLT_IsoMu22_eta2p1_v");
  trigsToRun.push_back("HLT_IsoMu24_v");
  trigsToRun.push_back("HLT_IsoMu27_v");
//
  //// Muon + jet
  trigsToRun.push_back("HLT_Mu30_eta2p1_PFJet150_PFJet50_v");
  trigsToRun.push_back("HLT_Mu40_eta2p1_PFJet200_PFJet50_v");
//
  //// Electron
  trigsToRun.push_back("HLT_Ele32_eta2p1_WPTight_Gsf_v");
  trigsToRun.push_back("HLT_Ele35_WPLoose_Gsf_v");
  trigsToRun.push_back("HLT_Ele105_CaloIdVT_GsfTrkIdT_v");
  trigsToRun.push_back("HLT_Ele115_CaloIdVT_GsfTrkIdT_v");
//
  // Electron + jet
  trigsToRun.push_back("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v");
  trigsToRun.push_back("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet140_v");
  trigsToRun.push_back("HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v");


//  888    888               888                          d8b               88888888888                       
//  888    888               888                          Y8P                   888                           
//  888    888               888                                                888                           
//  8888888888  8888b.   .d88888 888d888 .d88b.  88888b.  888  .d8888b          888  888d888 .d88b.   .d88b.  
//  888    888     "88b d88" 888 888P"  d88""88b 888 "88b 888 d88P"             888  888P"  d8P  Y8b d8P  Y8b 
//  888    888 .d888888 888  888 888    888  888 888  888 888 888               888  888    88888888 88888888 
//  888    888 888  888 Y88b 888 888    Y88..88P 888  888 888 Y88b.             888  888    Y8b.     Y8b.     
//  888    888 "Y888888  "Y88888 888     "Y88P"  888  888 888  "Y8888P          888  888     "Y8888   "Y8888  
//                                                                                                            
//                                                                                                            
//          http://patorjk.com/software/taag/#p=display&v=1&c=c%2B%2B&f=Colossal&t=Leptonic%20Tree          




    TTree *example; 
    eventDataStruct * example2;
    //bookEventTree(example, example2);



  TreeHad = new TTree("TreeHad","TreeHad"); 

  TreeHad->Branch("HLTtriggers"        , "vector<std::string>", &HLTtriggers);
  TreeHad->Branch("HLTtriggersPass"        , "vector<bool>", &HLTtriggersPass);
  TreeHad->Branch("HLTtriggersPrescales"        , "vector<int>", &HLTtriggersPrescales);

  TreeHad->Branch("HadTrigPrescales"   , "vector<int>", &HadTrigPrescales);
  TreeHad->Branch("HadTrigPass"        , "vector<bool>", &HadTrigPass);
  TreeHad->Branch("HadTrigAcceptBits"  , &HadTrigAcceptBits);

  TreeHad->Branch("HadTrigPrescalesMu"   , "vector<int>", &HadTrigPrescalesMu);
  TreeHad->Branch("HadTrigPassMu"        , "vector<bool>", &HadTrigPassMu);
  TreeHad->Branch("HadTrigAcceptBitsMu"  , &HadTrigAcceptBitsMu);

  TreeHad->Branch("HadTrigPrescalesHad"   , "vector<int>", &HadTrigPrescalesHad);
  TreeHad->Branch("HadTrigPassHad"        , "vector<bool>", &HadTrigPassHad);
  TreeHad->Branch("HadTrigAcceptBitsHad"  , &HadTrigAcceptBitsHad);
                     


  TreeHad->Branch("Gen_MET_pT",         & Gen_MET_pT          , "Gen_MET_pT/F"        );
  TreeHad->Branch("Gen_MET_phi",        & Gen_MET_phi         , "Gen_MET_phi/F"       );
  TreeHad->Branch("Gen_MET_eta",        & Gen_MET_eta         , "Gen_MET_eta/F"       );

  TreeHad->Branch("Gen_array_t_p4"                   , & Gen_array_t_p4                   ,  "Gen_array_t_p4[4]/F"                      );
  TreeHad->Branch("Gen_array_final_t_p4"             , & Gen_array_final_t_p4             ,  "Gen_array_final_t_p4[4]/F"                );
  TreeHad->Branch("Gen_array_b_p4"                   , & Gen_array_b_p4                   ,  "Gen_array_b_p4[4]/F"                      );
  TreeHad->Branch("Gen_array_W_p4"                   , & Gen_array_W_p4                   ,  "Gen_array_W_p4[4]/F"                      );
  TreeHad->Branch("Gen_array_Wd1_p4"                 , & Gen_array_Wd1_p4                 ,  "Gen_array_Wd1_p4[4]/F"                    );
  TreeHad->Branch("Gen_array_Wd2_p4"                 , & Gen_array_Wd2_p4                 ,  "Gen_array_Wd2_p4[4]/F"                    );
  TreeHad->Branch("Gen_array_hardest_parton_hardScatterOutgoing_p4", & Gen_array_hardest_parton_hardScatterOutgoing_p4 ,  "Gen_array_hardest_parton_hardScatterOutgoing_p4[4]/F");
  TreeHad->Branch("Gen_array_second_hardest_parton_hardScatterOutgoing_p4", & Gen_array_second_hardest_parton_hardScatterOutgoing_p4 ,  "Gen_array_second_hardest_parton_hardScatterOutgoing_p4[4]/F" );
  TreeHad->Branch("tophadronic"                         , & tophadronic                         ,  "tophadronic/O"                      );
  TreeHad->Branch("topleptonic"                         , & topleptonic                         ,  "topleptonic/O"                      );
  TreeHad->Branch("parton1id"                           , & parton1id                           ,  "parton1id/I"                        );
  TreeHad->Branch("parton2id"                           , & parton2id                           ,  "parton2id/I"                        );
  TreeHad->Branch("Wd1_id"                              , & Wd1_id                              ,  "Wd1_id/I"                           );
  TreeHad->Branch("Wd2_id"                              , & Wd2_id                              ,  "Wd2_id/I"                           ); 

  TreeHad->Branch("AK4JetLV_pt","vector<float>",&AK4JetLV_pt);
  TreeHad->Branch("AK4JetLV_eta","vector<float>",&AK4JetLV_eta);
  TreeHad->Branch("AK4JetLV_phi","vector<float>",&AK4JetLV_phi);
  TreeHad->Branch("AK4JetLV_mass","vector<float>",&AK4JetLV_mass);

  TreeHad->Branch("AK4JetLV_corr", "vector<float>", &AK4JetLV_corr);
  TreeHad->Branch("AK4JetLV_corrUp", "vector<float>", &AK4JetLV_corrUp);
  TreeHad->Branch("AK4JetLV_corrDn", "vector<float>", &AK4JetLV_corrDn);
  TreeHad->Branch("AK4JetLV_SF", "vector<float>", &AK4JetLV_SF);
  TreeHad->Branch("AK4JetLV_SFUp", "vector<float>", &AK4JetLV_SFUp);
  TreeHad->Branch("AK4JetLV_SFDn", "vector<float>", &AK4JetLV_SFDn);
  TreeHad->Branch("AK4JetLV_ptsmear", "vector<float>", &AK4JetLV_ptsmear);
  TreeHad->Branch("AK4JetLV_ptsmearUp", "vector<float>", &AK4JetLV_ptsmearUp);
  TreeHad->Branch("AK4JetLV_ptsmearDn", "vector<float>", &AK4JetLV_ptsmearDn);

  TreeHad->Branch("AK8JetLV_corr", "vector<float>", &AK8JetLV_corr);
  TreeHad->Branch("AK8JetLV_corrUp", "vector<float>", &AK8JetLV_corrUp);
  TreeHad->Branch("AK8JetLV_corrDn", "vector<float>", &AK8JetLV_corrDn);
  TreeHad->Branch("AK8JetLV_SF", "vector<float>", &AK8JetLV_SF);
  TreeHad->Branch("AK8JetLV_SFUp", "vector<float>", &AK8JetLV_SFUp);
  TreeHad->Branch("AK8JetLV_SFDn", "vector<float>", &AK8JetLV_SFDn);
  TreeHad->Branch("AK8JetLV_ptsmear", "vector<float>", &AK8JetLV_ptsmear);
  TreeHad->Branch("AK8JetLV_ptsmearUp", "vector<float>", &AK8JetLV_ptsmearUp);
  TreeHad->Branch("AK8JetLV_ptsmearDn", "vector<float>", &AK8JetLV_ptsmearDn);

  TreeHad->Branch("AK8JetLV_pt","vector<float>",&AK8JetLV_pt);
  TreeHad->Branch("AK8JetLV_eta","vector<float>",&AK8JetLV_eta);
  TreeHad->Branch("AK8JetLV_phi","vector<float>",&AK8JetLV_phi);
  TreeHad->Branch("AK8JetLV_mass","vector<float>",&AK8JetLV_mass);

  TreeHad->Branch("AK8SubjetLV_pt","vector<float>",&AK8SubjetLV_pt);
  TreeHad->Branch("AK8SubjetLV_eta","vector<float>",&AK8SubjetLV_eta);
  TreeHad->Branch("AK8SubjetLV_phi","vector<float>",&AK8SubjetLV_phi);
  TreeHad->Branch("AK8SubjetLV_mass","vector<float>",&AK8SubjetLV_mass);

  TreeHad->Branch("AK4JetBtag_p","vector<float>",&AK4JetBtag_p);
  TreeHad->Branch("AK8JetTau1_p","vector<float>",&AK8JetTau1_p);
  TreeHad->Branch("AK8JetTau2_p","vector<float>",&AK8JetTau2_p);
  TreeHad->Branch("AK8JetTau3_p","vector<float>",&AK8JetTau3_p);
  TreeHad->Branch("AK8JetSoftdropMass_p","vector<float>",&AK8JetSoftdropMass_p);



  TreeHad->Branch("JetPtRaw"                             , & JetPtRaw                          ,    "JetPtRaw/F"                               );                                  
  TreeHad->Branch("JetEtaRaw"                            , & JetEtaRaw                         ,    "JetEtaRaw/F"                              );                                   
  TreeHad->Branch("JetPhiRaw"                            , & JetPhiRaw                         ,    "JetPhiRaw/F"                              );                                   
  TreeHad->Branch("JetMassRaw"                           , & JetMassRaw                        ,    "JetMassRaw/F"                             );                                                                                      
  TreeHad->Branch("JetArea"                              , & JetArea                           ,    "JetArea/F"                                );                                 
                                    
  TreeHad->Branch("JetSDmassRaw"                         , & JetSDmassRaw                      ,    "JetSDmassRaw/F"                           );                                               
  TreeHad->Branch("JetSDmassSubjetCorrL23"                     , & JetSDmassSubjetCorrL23                  ,    "JetSDmassSubjetCorrL23/F"                       );                                                                                                     
  TreeHad->Branch("JetSDmassSubjetCorrL123"                    , & JetSDmassSubjetCorrL123                 ,    "JetSDmassSubjetCorrL123/F"                      );                                                      
  TreeHad->Branch("JetSDptRaw"                           , & JetSDptRaw                        ,    "JetSDptRaw/F"                             );                                                                                                 
  TreeHad->Branch("JetSDetaRaw"                          , & JetSDetaRaw                       ,    "JetSDetaRaw/F"                            );                                               
  TreeHad->Branch("JetSDphiRaw"                          , & JetSDphiRaw                       ,    "JetSDphiRaw/F"                            );  

  TreeHad->Branch("JetMassPruned"                        , & JetMassPruned                     ,    "JetMassPruned/F"                          );                                       
  TreeHad->Branch("JetMassTrimmed"                       , & JetMassTrimmed                    ,    "JetMassTrimmed/F"                         );                                       
  TreeHad->Branch("JetTau1"                              , & JetTau1                           ,    "JetTau1/F"                                );                                 
  TreeHad->Branch("JetTau2"                              , & JetTau2                           ,    "JetTau2/F"                                );                                 
  TreeHad->Branch("JetTau3"                              , & JetTau3                           ,    "JetTau3/F"                                );                                 
  TreeHad->Branch("JetTau4"                              , & JetTau4                           ,    "JetTau4/F"                                );                                 
  TreeHad->Branch("JetTau32"                             , & JetTau32                          ,    "JetTau32/F"                               );                                  
  TreeHad->Branch("JetTau21"                             , & JetTau21                          ,    "JetTau21/F"                               );                                  
  TreeHad->Branch("JetSDmaxbdisc"                        , & JetSDmaxbdisc                     ,    "JetSDmaxbdisc/F"                          );                                       
  TreeHad->Branch("JetSDmaxbdiscflavHadron"              , & JetSDmaxbdiscflavHadron           ,    "JetSDmaxbdiscflavHadron/F"                );                                           
  TreeHad->Branch("JetSDmaxbdiscflavParton"              , & JetSDmaxbdiscflavParton           ,    "JetSDmaxbdiscflavParton/F"                );  
                                         
  TreeHad->Branch("JetSDsubjet0pt"                       , & JetSDsubjet0pt                    ,    "JetSDsubjet0pt/F"                         );    
  TreeHad->Branch("JetSDsubjet0mass"                     , & JetSDsubjet0mass                  ,    "JetSDsubjet0mass/F"                       );
  TreeHad->Branch("JetSDsubjet0eta"                      , & JetSDsubjet0eta                   ,    "JetSDsubjet0eta/F"                        );
  TreeHad->Branch("JetSDsubjet0phi"                      , & JetSDsubjet0phi                   ,    "JetSDsubjet0phi/F"                        );
  TreeHad->Branch("JetSDsubjet0area"                     , & JetSDsubjet0area                  ,    "JetSDsubjet0area/F"                       );
  TreeHad->Branch("JetSDsubjet0flavHadron"               , & JetSDsubjet0flavHadron            ,    "JetSDsubjet0flavHadron/F"                 );
  TreeHad->Branch("JetSDsubjet0flavParton"               , & JetSDsubjet0flavParton            ,    "JetSDsubjet0flavParton/F"                 );
  TreeHad->Branch("JetSDsubjet0matchedgenjetpt"          , & JetSDsubjet0matchedgenjetpt       ,    "JetSDsubjet0matchedgenjetpt/F"            ); 
  TreeHad->Branch("JetSDsubjet0tau1"                     , & JetSDsubjet0tau1                  ,    "JetSDsubjet0tau1/F"                       );
  TreeHad->Branch("JetSDsubjet0tau2"                     , & JetSDsubjet0tau2                  ,    "JetSDsubjet0tau2/F"                       );
  TreeHad->Branch("JetSDsubjet0tau3"                     , & JetSDsubjet0tau3                  ,    "JetSDsubjet0tau3/F"                       ); 
  TreeHad->Branch("JetSDsubjet0bdisc"                    , & JetSDsubjet0bdisc                 ,    "JetSDsubjet0bdisc/F"                      );                                          
  TreeHad->Branch("JetSDsubjet1pt"                       , & JetSDsubjet1pt                    ,    "JetSDsubjet1pt/F"                         );    
  TreeHad->Branch("JetSDsubjet1mass"                     , & JetSDsubjet1mass                  ,    "JetSDsubjet1mass/F"                       );
  TreeHad->Branch("JetSDsubjet1eta"                      , & JetSDsubjet1eta                   ,    "JetSDsubjet1eta/F"                        );
  TreeHad->Branch("JetSDsubjet1phi"                      , & JetSDsubjet1phi                   ,    "JetSDsubjet1phi/F"                        );  
  TreeHad->Branch("JetSDsubjet1area"                     , & JetSDsubjet1area                  ,    "JetSDsubjet1area/F"                       );
  TreeHad->Branch("JetSDsubjet1flavHadron"               , & JetSDsubjet1flavHadron            ,    "JetSDsubjet1flavHadron/F"                 );
  TreeHad->Branch("JetSDsubjet1flavParton"               , & JetSDsubjet1flavParton            ,    "JetSDsubjet1flavParton/F"                 );
  TreeHad->Branch("JetSDsubjet1matchedgenjetpt"          , & JetSDsubjet1matchedgenjetpt       ,    "JetSDsubjet1matchedgenjetpt/F"            ); 
  TreeHad->Branch("JetSDsubjet1tau1"                     , & JetSDsubjet1tau1                  ,    "JetSDsubjet1tau1/F"                       );
  TreeHad->Branch("JetSDsubjet1tau2"                     , & JetSDsubjet1tau2                  ,    "JetSDsubjet1tau2/F"                       );
  TreeHad->Branch("JetSDsubjet1tau3"                     , & JetSDsubjet1tau3                  ,    "JetSDsubjet1tau3/F"                       );                                           
  TreeHad->Branch("JetSDsubjet1bdisc"                    , & JetSDsubjet1bdisc                 ,    "JetSDsubjet1bdisc/F"                      );                                     

  TreeHad->Branch("JetPuppiPtRaw"                           , & JetPuppiPtRaw                        ,    "JetPuppiPtRaw/F"                             );                                    
  TreeHad->Branch("JetPuppiEtaRaw"                          , & JetPuppiEtaRaw                       ,    "JetPuppiEtaRaw/F"                            );                                     
  TreeHad->Branch("JetPuppiPhiRaw"                          , & JetPuppiPhiRaw                       ,    "JetPuppiPhiRaw/F"                            );                                     
  TreeHad->Branch("JetPuppiMassRaw"                         , & JetPuppiMassRaw                      ,    "JetPuppiMassRaw/F"                           );                                      
  TreeHad->Branch("JetPuppiArea"                         , & JetPuppiArea                      ,    "JetPuppiArea/F"                           );                                      

  TreeHad->Branch("JetPuppiMassPruned"                    , & JetPuppiMassPruned               ,   "JetPuppiMassPruned/F"                          );
  TreeHad->Branch("JetPuppiMassTrimmed"                   , & JetPuppiMassTrimmed              ,   "JetPuppiMassTrimmed/F"                         );


  
  TreeHad->Branch("JetPuppiSDmassRaw"                         , & JetPuppiSDmassRaw                    ,    "JetPuppiSDmassRaw/F"                          );
  TreeHad->Branch("JetPuppiSDmassSubjetCorr"                     , & JetPuppiSDmassSubjetCorr                ,    "JetPuppiSDmassSubjetCorr/F"                      );
  TreeHad->Branch("JetPuppiSDptRaw"                           , & JetPuppiSDptRaw                      ,    "JetPuppiSDptRaw/F"                            );
  TreeHad->Branch("JetPuppiSDetaRaw"                          , & JetPuppiSDetaRaw                     ,    "JetPuppiSDetaRaw/F"                           );
  TreeHad->Branch("JetPuppiSDphiRaw"                          , & JetPuppiSDphiRaw                     ,    "JetPuppiSDphiRaw/F"                           );
                         

  TreeHad->Branch("JetPuppiTau1"                         , & JetPuppiTau1                      ,    "JetPuppiTau1/F"                           );                                      
  TreeHad->Branch("JetPuppiTau2"                         , & JetPuppiTau2                      ,    "JetPuppiTau2/F"                           );                                      
  TreeHad->Branch("JetPuppiTau3"                         , & JetPuppiTau3                      ,    "JetPuppiTau3/F"                           );                                      
  TreeHad->Branch("JetPuppiTau4"                         , & JetPuppiTau4                      ,    "JetPuppiTau4/F"                           );                                      
  TreeHad->Branch("JetPuppiTau32"                        , & JetPuppiTau32                     ,    "JetPuppiTau32/F"                          );                                       
  TreeHad->Branch("JetPuppiTau21"                        , & JetPuppiTau21                     ,    "JetPuppiTau21/F"                          );                                       

  TreeHad->Branch("JetPuppiSDmaxbdisc"                   , & JetPuppiSDmaxbdisc                ,    "JetPuppiSDmaxbdisc/F"                     );                                            
  TreeHad->Branch("JetPuppiSDmaxbdiscflavHadron"         , & JetPuppiSDmaxbdiscflavHadron      ,    "JetPuppiSDmaxbdiscflavHadron/F"           );                                                
  TreeHad->Branch("JetPuppiSDmaxbdiscflavParton"         , & JetPuppiSDmaxbdiscflavParton      ,    "JetPuppiSDmaxbdiscflavParton/F"           );                                                
  TreeHad->Branch("JetPuppiSDsubjet0pt"                  , & JetPuppiSDsubjet0pt               ,    "JetPuppiSDsubjet0pt/F"                    );    
  TreeHad->Branch("JetPuppiSDsubjet0mass"                , & JetPuppiSDsubjet0mass             ,    "JetPuppiSDsubjet0mass/F"                  );
  TreeHad->Branch("JetPuppiSDsubjet0eta"                 , & JetPuppiSDsubjet0eta              ,    "JetPuppiSDsubjet0eta/F"                   );
  TreeHad->Branch("JetPuppiSDsubjet0phi"                 , & JetPuppiSDsubjet0phi              ,    "JetPuppiSDsubjet0phi/F"                   );
  TreeHad->Branch("JetPuppiSDsubjet0area"                , & JetPuppiSDsubjet0area             ,    "JetPuppiSDsubjet0area/F"                  );
  TreeHad->Branch("JetPuppiSDsubjet0flavHadron"          , & JetPuppiSDsubjet0flavHadron       ,    "JetPuppiSDsubjet0flavHadron/F"            );
  TreeHad->Branch("JetPuppiSDsubjet0flavParton"          , & JetPuppiSDsubjet0flavParton       ,    "JetPuppiSDsubjet0flavParton/F"            );
  TreeHad->Branch("JetPuppiSDsubjet0matchedgenjetpt"     , & JetPuppiSDsubjet0matchedgenjetpt  ,    "JetPuppiSDsubjet0matchedgenjetpt/F"       );
  TreeHad->Branch("JetPuppiSDsubjet0tau1"                , & JetPuppiSDsubjet0tau1             ,    "JetPuppiSDsubjet0tau1/F"                  );
  TreeHad->Branch("JetPuppiSDsubjet0tau2"                , & JetPuppiSDsubjet0tau2             ,    "JetPuppiSDsubjet0tau2/F"                  );
  TreeHad->Branch("JetPuppiSDsubjet0tau3"                , & JetPuppiSDsubjet0tau3             ,    "JetPuppiSDsubjet0tau3/F"                  );
  TreeHad->Branch("JetPuppiSDsubjet0bdisc"               , & JetPuppiSDsubjet0bdisc            ,    "JetPuppiSDsubjet0bdisc/F"                 );                                                
  TreeHad->Branch("JetPuppiSDsubjet1pt"                  , & JetPuppiSDsubjet1pt               ,    "JetPuppiSDsubjet1pt/F"                    );    
  TreeHad->Branch("JetPuppiSDsubjet1mass"                , & JetPuppiSDsubjet1mass             ,    "JetPuppiSDsubjet1mass/F"                  );
  TreeHad->Branch("JetPuppiSDsubjet1eta"                 , & JetPuppiSDsubjet1eta              ,    "JetPuppiSDsubjet1eta/F"                   );
  TreeHad->Branch("JetPuppiSDsubjet1phi"                 , & JetPuppiSDsubjet1phi              ,    "JetPuppiSDsubjet1phi/F"                   );  
  TreeHad->Branch("JetPuppiSDsubjet1area"                , & JetPuppiSDsubjet1area             ,    "JetPuppiSDsubjet1area/F"                  );
  TreeHad->Branch("JetPuppiSDsubjet1flavHadron"          , & JetPuppiSDsubjet1flavHadron       ,    "JetPuppiSDsubjet1flavHadron/F"            );
  TreeHad->Branch("JetPuppiSDsubjet1flavParton"          , & JetPuppiSDsubjet1flavParton       ,    "JetPuppiSDsubjet1flavParton/F"            );
  TreeHad->Branch("JetPuppiSDsubjet1matchedgenjetpt"     , & JetPuppiSDsubjet1matchedgenjetpt  ,    "JetPuppiSDsubjet1matchedgenjetpt/F"       );
  TreeHad->Branch("JetPuppiSDsubjet1tau1"                , & JetPuppiSDsubjet1tau1             ,    "JetPuppiSDsubjet1tau1/F"                  );
  TreeHad->Branch("JetPuppiSDsubjet1tau2"                , & JetPuppiSDsubjet1tau2             ,    "JetPuppiSDsubjet1tau2/F"                  );
  TreeHad->Branch("JetPuppiSDsubjet1tau3"                , & JetPuppiSDsubjet1tau3             ,    "JetPuppiSDsubjet1tau3/F"                  );  
  TreeHad->Branch("JetPuppiSDECF1"                       , & JetPuppiSDECF1                    ,    "JetPuppiSDECF1/F"                         );
  TreeHad->Branch("JetPuppiSDECF2"                       , & JetPuppiSDECF2                    ,    "JetPuppiSDECF2/F"                         );
  TreeHad->Branch("JetPuppiSDECF3"                       , & JetPuppiSDECF3                    ,    "JetPuppiSDECF3/F"                         );
  TreeHad->Branch("JetPuppiSDECF4"                       , & JetPuppiSDECF4                    ,    "JetPuppiSDECF4/F"                         );
  TreeHad->Branch("JetPuppiSDECF5"                       , & JetPuppiSDECF5                    ,    "JetPuppiSDECF5/F"                         );
  TreeHad->Branch("JetPuppiSDC_2"                       , & JetPuppiSDC_2                    ,    "JetPuppiSDC_2/F"                         );
  TreeHad->Branch("JetPuppiSDD_2"                       , & JetPuppiSDD_2                    ,    "JetPuppiSDD_2/F"                         );
  TreeHad->Branch("JetPuppiSDC_3"                       , & JetPuppiSDC_3                    ,    "JetPuppiSDC_3/F"                         );
  TreeHad->Branch("JetPuppiSDD_3"                       , & JetPuppiSDD_3                    ,    "JetPuppiSDD_3/F"                         );
  TreeHad->Branch("JetPuppiSDsubjet1bdisc"               , & JetPuppiSDsubjet1bdisc            ,    "JetPuppiSDsubjet1bdisc/F"                 );                                                                                                                        


  TreeHad->Branch("JetCHF"                               , & JetCHF                            ,    "JetCHF/F"                                 );                                
  TreeHad->Branch("JetNHF"                               , & JetNHF                            ,    "JetNHF/F"                                 );                                
  TreeHad->Branch("JetCM"                                , & JetCM                             ,    "JetCM/F"                                  );                               
  TreeHad->Branch("JetNM"                                , & JetNM                             ,    "JetNM/F"                                  );                               
  TreeHad->Branch("JetNEF"                               , & JetNEF                            ,    "JetNEF/F"                                 );                                
  TreeHad->Branch("JetCEF"                               , & JetCEF                            ,    "JetCEF/F"                                 );                                
  TreeHad->Branch("JetMF"                                , & JetMF                             ,    "JetMF/F"                                  );                               
  TreeHad->Branch("JetMult"                              , & JetMult                           ,    "JetMult/F"                                );
  TreeHad->Branch("JetPuppiCHF"                          , & JetPuppiCHF                       ,    "JetPuppiCHF/F"                            );                                
  TreeHad->Branch("JetPuppiNHF"                          , & JetPuppiNHF                       ,    "JetPuppiNHF/F"                            );                                
  TreeHad->Branch("JetPuppiCM"                           , & JetPuppiCM                        ,    "JetPuppiCM/F"                             );                               
  TreeHad->Branch("JetPuppiNM"                           , & JetPuppiNM                        ,    "JetPuppiNM/F"                             );                               
  TreeHad->Branch("JetPuppiNEF"                          , & JetPuppiNEF                       ,    "JetPuppiNEF/F"                            );                                
  TreeHad->Branch("JetPuppiCEF"                          , & JetPuppiCEF                       ,    "JetPuppiCEF/F"                            );                                
  TreeHad->Branch("JetPuppiMF"                           , & JetPuppiMF                        ,    "JetPuppiMF/F"                             );                               
  TreeHad->Branch("JetPuppiMult"                         , & JetPuppiMult                      ,    "JetPuppiMult/F"                           );                                  
  TreeHad->Branch("JetMassCorrFactor"                    , & JetMassCorrFactor                 ,    "JetMassCorrFactor/F"                      );                                           
  TreeHad->Branch("JetMassCorrFactorUp"                  , & JetMassCorrFactorUp               ,    "JetMassCorrFactorUp/F"                    );                                             
  TreeHad->Branch("JetMassCorrFactorDn"                  , & JetMassCorrFactorDn               ,    "JetMassCorrFactorDn/F"                    );                                             
  TreeHad->Branch("JetCorrFactor"                        , & JetCorrFactor                     ,    "JetCorrFactor/F"                          );                                       
  TreeHad->Branch("JetCorrFactorUp"                      , & JetCorrFactorUp                   ,    "JetCorrFactorUp/F"                        );                                         
  TreeHad->Branch("JetCorrFactorDn"                      , & JetCorrFactorDn                   ,    "JetCorrFactorDn/F"                        );                                         
  TreeHad->Branch("JetPtSmearFactor"                     , & JetPtSmearFactor                  ,    "JetPtSmearFactor/F"                       );                                          
  TreeHad->Branch("JetPtSmearFactorUp"                   , & JetPtSmearFactorUp                ,    "JetPtSmearFactorUp/F"                     );                                            
  TreeHad->Branch("JetPtSmearFactorDn"                   , & JetPtSmearFactorDn                ,    "JetPtSmearFactorDn/F"                     );                                            
  TreeHad->Branch("JetPuppiMassCorrFactor"               , & JetPuppiMassCorrFactor            ,    "JetPuppiMassCorrFactor/F"                 );                                                
  TreeHad->Branch("JetPuppiMassCorrFactorUp"             , & JetPuppiMassCorrFactorUp          ,    "JetPuppiMassCorrFactorUp/F"               );                                                  
  TreeHad->Branch("JetPuppiMassCorrFactorDn"             , & JetPuppiMassCorrFactorDn          ,    "JetPuppiMassCorrFactorDn/F"               );                                                  
  TreeHad->Branch("JetPuppiCorrFactor"                   , & JetPuppiCorrFactor                ,    "JetPuppiCorrFactor/F"                     );                                            
  TreeHad->Branch("JetPuppiCorrFactorUp"                 , & JetPuppiCorrFactorUp              ,    "JetPuppiCorrFactorUp/F"                   );                                              
  TreeHad->Branch("JetPuppiCorrFactorDn"                 , & JetPuppiCorrFactorDn              ,    "JetPuppiCorrFactorDn/F"                   );                                              
  TreeHad->Branch("JetPuppiPtSmearFactor"                , & JetPuppiPtSmearFactor             ,    "JetPuppiPtSmearFactor/F"                  );                                               
  TreeHad->Branch("JetPuppiPtSmearFactorUp"              , & JetPuppiPtSmearFactorUp           ,    "JetPuppiPtSmearFactorUp/F"                );                                                 
  TreeHad->Branch("JetPuppiPtSmearFactorDn"              , & JetPuppiPtSmearFactorDn           ,    "JetPuppiPtSmearFactorDn/F"                );                                                 
  TreeHad->Branch("JetMatchedGenJetPt"                   , & JetMatchedGenJetPt                ,    "JetMatchedGenJetPt/F"                     );                                            
  TreeHad->Branch("JetMatchedGenJetMass"                 , & JetMatchedGenJetMass              ,    "JetMatchedGenJetMass/F"                   ); 
  TreeHad->Branch("JetPuppiMatchedGenJetPt"              , & JetPuppiMatchedGenJetPt           ,    "JetPuppiMatchedGenJetPt/F"                );                                            
  TreeHad->Branch("JetPuppiMatchedGenJetMass"            , & JetPuppiMatchedGenJetMass         ,    "JetPuppiMatchedGenJetMass/F"              ); 
                           
  TreeHad->Branch("JetGenMatched_TopHadronic"            , & JetGenMatched_TopHadronic         ,    "JetGenMatched_TopHadronic/I"              );      
  TreeHad->Branch("JetGenMatched_TopPt"                  , & JetGenMatched_TopPt               ,    "JetGenMatched_TopPt/F"                    );      
  TreeHad->Branch("JetGenMatched_TopEta"                 , & JetGenMatched_TopEta              ,    "JetGenMatched_TopEta/F"                   );      
  TreeHad->Branch("JetGenMatched_TopPhi"                 , & JetGenMatched_TopPhi              ,    "JetGenMatched_TopPhi/F"                   );      
  TreeHad->Branch("JetGenMatched_TopMass"                , & JetGenMatched_TopMass             ,    "JetGenMatched_TopMass/F"                  );      
  TreeHad->Branch("JetGenMatched_bPt"                    , & JetGenMatched_bPt                 ,    "JetGenMatched_bPt/F"                      );      
  TreeHad->Branch("JetGenMatched_WPt"                    , & JetGenMatched_WPt                 ,    "JetGenMatched_WPt/F"                      );      
  TreeHad->Branch("JetGenMatched_Wd1Pt"                  , & JetGenMatched_Wd1Pt               ,    "JetGenMatched_Wd1Pt/F"                    );      
  TreeHad->Branch("JetGenMatched_Wd2Pt"                  , & JetGenMatched_Wd2Pt               ,    "JetGenMatched_Wd2Pt/F"                    );      
  TreeHad->Branch("JetGenMatched_Wd1ID"                  , & JetGenMatched_Wd1ID               ,    "JetGenMatched_Wd1ID/F"                    );      
  TreeHad->Branch("JetGenMatched_Wd2ID"                  , & JetGenMatched_Wd2ID               ,    "JetGenMatched_Wd2ID/F"                    );      
  TreeHad->Branch("JetGenMatched_MaxDeltaRPartonTop"     , & JetGenMatched_MaxDeltaRPartonTop  ,    "JetGenMatched_MaxDeltaRPartonTop/F"       );      
  TreeHad->Branch("JetGenMatched_MaxDeltaRWPartonTop"    , & JetGenMatched_MaxDeltaRWPartonTop ,    "JetGenMatched_MaxDeltaRWPartonTop/F"      );      
  TreeHad->Branch("JetGenMatched_MaxDeltaRWPartonW"      , & JetGenMatched_MaxDeltaRWPartonW   ,    "JetGenMatched_MaxDeltaRWPartonW/F"        );      
  TreeHad->Branch("JetGenMatched_DeltaR_t_b"             , & JetGenMatched_DeltaR_t_b          ,    "JetGenMatched_DeltaR_t_b/F"               );      
  TreeHad->Branch("JetGenMatched_DeltaR_t_W"             , & JetGenMatched_DeltaR_t_W          ,    "JetGenMatched_DeltaR_t_W/F"               );      
  TreeHad->Branch("JetGenMatched_DeltaR_t_Wd1"           , & JetGenMatched_DeltaR_t_Wd1        ,    "JetGenMatched_DeltaR_t_Wd1/F"             );      
  TreeHad->Branch("JetGenMatched_DeltaR_t_Wd2"           , & JetGenMatched_DeltaR_t_Wd2        ,    "JetGenMatched_DeltaR_t_Wd2/F"             );      
  TreeHad->Branch("JetGenMatched_DeltaR_W_b1"            , & JetGenMatched_DeltaR_W_b1         ,    "JetGenMatched_DeltaR_W_b1/F"              );      
  TreeHad->Branch("JetGenMatched_DeltaR_W_Wd1"           , & JetGenMatched_DeltaR_W_Wd1        ,    "JetGenMatched_DeltaR_W_Wd1/F"             );      
  TreeHad->Branch("JetGenMatched_DeltaR_W_Wd2"           , & JetGenMatched_DeltaR_W_Wd2        ,    "JetGenMatched_DeltaR_W_Wd2/F"             );      
  TreeHad->Branch("JetGenMatched_DeltaR_Wd1_Wd2"         , & JetGenMatched_DeltaR_Wd1_Wd2      ,    "JetGenMatched_DeltaR_Wd1_Wd2/F"           );      
  TreeHad->Branch("JetGenMatched_DeltaR_Wd1_b"           , & JetGenMatched_DeltaR_Wd1_b        ,    "JetGenMatched_DeltaR_Wd1_b/F"             );      
  TreeHad->Branch("JetGenMatched_DeltaR_Wd2_b"           , & JetGenMatched_DeltaR_Wd2_b        ,    "JetGenMatched_DeltaR_Wd2_b/F"             );      
  TreeHad->Branch("JetGenMatched_DeltaR_jet_t"           , & JetGenMatched_DeltaR_jet_t        ,    "JetGenMatched_DeltaR_jet_t/F"             );      
  TreeHad->Branch("JetGenMatched_DeltaR_jet_W"           , & JetGenMatched_DeltaR_jet_W        ,    "JetGenMatched_DeltaR_jet_W/F"             );      
  TreeHad->Branch("JetGenMatched_DeltaR_jet_b"           , & JetGenMatched_DeltaR_jet_b        ,    "JetGenMatched_DeltaR_jet_b/F"             );      
  TreeHad->Branch("JetGenMatched_DeltaR_jet_Wd1"         , & JetGenMatched_DeltaR_jet_Wd1      ,    "JetGenMatched_DeltaR_jet_Wd1/F"           );      
  TreeHad->Branch("JetGenMatched_DeltaR_jet_Wd2"         , & JetGenMatched_DeltaR_jet_Wd2      ,    "JetGenMatched_DeltaR_jet_Wd2/F"           );      
  TreeHad->Branch("JetGenMatched_DeltaR_pup0_b"          , & JetGenMatched_DeltaR_pup0_b       ,    "JetGenMatched_DeltaR_pup0_b/F"            );      
  TreeHad->Branch("JetGenMatched_DeltaR_pup0_Wd1"        , & JetGenMatched_DeltaR_pup0_Wd1     ,    "JetGenMatched_DeltaR_pup0_Wd1/F"          );      
  TreeHad->Branch("JetGenMatched_DeltaR_pup0_Wd2"        , & JetGenMatched_DeltaR_pup0_Wd2     ,    "JetGenMatched_DeltaR_pup0_Wd2/F"          );      
  TreeHad->Branch("JetGenMatched_DeltaR_pup1_b"          , & JetGenMatched_DeltaR_pup1_b       ,    "JetGenMatched_DeltaR_pup1_b/F"            );      
  TreeHad->Branch("JetGenMatched_DeltaR_pup1_Wd1"        , & JetGenMatched_DeltaR_pup1_Wd1     ,    "JetGenMatched_DeltaR_pup1_Wd1/F"          );      
  TreeHad->Branch("JetGenMatched_DeltaR_pup1_Wd2"        , & JetGenMatched_DeltaR_pup1_Wd2     ,    "JetGenMatched_DeltaR_pup1_Wd2/F"          );               
  TreeHad->Branch("JetGenMatched_partonPt"               , & JetGenMatched_partonPt            ,    "JetGenMatched_partonPt/F"                 );      
  TreeHad->Branch("JetGenMatched_partonEta"              , & JetGenMatched_partonEta           ,    "JetGenMatched_partonEta/F"                );      
  TreeHad->Branch("JetGenMatched_partonPhi"              , & JetGenMatched_partonPhi           ,    "JetGenMatched_partonPhi/F"                );      
  TreeHad->Branch("JetGenMatched_partonMass"             , & JetGenMatched_partonMass          ,    "JetGenMatched_partonMass/F"               );      
  TreeHad->Branch("JetGenMatched_partonID"               , & JetGenMatched_partonID            ,    "JetGenMatched_partonID/F"                 );      
  TreeHad->Branch("JetGenMatched_DeltaRjetParton"        , & JetGenMatched_DeltaRjetParton     ,    "JetGenMatched_DeltaRjetParton/F"          );      
  std::cout<<"Setup semi-lept jets in tree"<<std::endl;


  TreeHad->Branch("HadMETpx"                        , & HadMETpx                     , "HadMETpx/F"                  );
  TreeHad->Branch("HadMETpy"                        , & HadMETpy                     , "HadMETpy/F"                  );
  TreeHad->Branch("HadMETpt"                        , & HadMETpt                     , "HadMETpt/F"                  );
  TreeHad->Branch("HadMETphi"                       , & HadMETphi                    , "HadMETphi/F"                 );
  TreeHad->Branch("HadMETsumET"                     , & HadMETsumET                  , "HadMETsumET/F"               );
  TreeHad->Branch("HadMETgenMET"                    , & HadMETgenMET                 , "HadMETgenMET/F"              );
  TreeHad->Branch("HadMETuncorPt"                   , & HadMETuncorPt                , "HadMETuncorPt/F"             );

  TreeHad->Branch("HadMETshiftedPtJetEnUp"      , & HadMETshiftedPtJetEnUp   , "HadMETshiftedPtJetEnUp/F"     );
  TreeHad->Branch("HadMETshiftedPtJetEnDn"      , & HadMETshiftedPtJetEnDn   , "HadMETshiftedPtJetEnDn/F"     );
  TreeHad->Branch("HadMETshiftedPtElEnUp"       , & HadMETshiftedPtElEnUp    , "HadMETshiftedPtElEnUp/F"      );
  TreeHad->Branch("HadMETshiftedPtElEnDn"       , & HadMETshiftedPtElEnDn    , "HadMETshiftedPtElEnDn/F"      );
  TreeHad->Branch("HadMETshiftedPtMuEnUp"       , & HadMETshiftedPtMuEnUp    , "HadMETshiftedPtMuEnUp/F"      );
  TreeHad->Branch("HadMETshiftedPtMuEnDn"       , & HadMETshiftedPtMuEnDn    , "HadMETshiftedPtMuEnDn/F"      );
  TreeHad->Branch("HadMETshiftedPtJetResUp"     , & HadMETshiftedPtJetResUp  , "HadMETshiftedPtJetResUp/F"    );
  TreeHad->Branch("HadMETshiftedPtJetResDn"     , & HadMETshiftedPtJetResDn  , "HadMETshiftedPtJetResDn/F"    );
  TreeHad->Branch("HadMETshiftedPtUnclEnUp"     , & HadMETshiftedPtUnclEnUp  , "HadMETshiftedPtUnclEnUp/F"    );
  TreeHad->Branch("HadMETshiftedPtUnclEnDn"     , & HadMETshiftedPtUnclEnDn  , "HadMETshiftedPtUnclEnDn/F"    );

  TreeHad->Branch("HadNvtx"                         , & HadNvtx                      , "HadNvtx/F"                   );
  TreeHad->Branch("HadNvtxGood"                     , & HadNvtxGood                  , "HadNvtxGood/F"               );
  TreeHad->Branch("HadRho"                          , & HadRho                       , "HadRho/F"                    );
  TreeHad->Branch("HadEventWeight"                  , & HadEventWeight               , "HadEventWeight/F"            );
  TreeHad->Branch("HadPUweight"                     , & HadPUweight                  , "HadPUweight/F"            );
  TreeHad->Branch("HadPUweight_MBup"                , & HadPUweight_MBup             , "HadPUweight_MBup/F"            );
  TreeHad->Branch("HadPUweight_MBdn"                , & HadPUweight_MBdn             , "HadPUweight_MBdn/F"            );
       
  

  TreeHad->Branch("HadGenTTmass"                    , & HadGenTTmass                 , "HadGenTTmass/F"              );
  TreeHad->Branch("HadGenCountHadTop"               , & HadGenCountHadTop            , "HadGenCountHadTop/I"         );
  
  TreeHad->Branch("HTlep"                                , & HTlep                             , "HTlep/F"                  );
  TreeHad->Branch("ST"                                   , & ST                                , "ST/F"                     );
  TreeHad->Branch("ST_CorrDn"                            , & ST_CorrDn                         , "ST_CorrDn/F"              );
  TreeHad->Branch("ST_CorrUp"                            , & ST_CorrUp                         , "ST_CorrUp/F"              );
  TreeHad->Branch("ST_PtSmearNom"                        , & ST_PtSmearNom                     , "ST_PtSmearNom/F"          );
  TreeHad->Branch("ST_PtSmearUp"                         , & ST_PtSmearUp                      , "ST_PtSmearUp/F"           );
  TreeHad->Branch("ST_PtSmearDn"                         , & ST_PtSmearDn                      , "ST_PtSmearDn/F"           );
  
  TreeHad->Branch("HadQ2weight_CorrDn"              , & HadQ2weight_CorrDn           , "HadQ2weight_CorrDn/F"        );
  TreeHad->Branch("HadQ2weight_CorrUp"              , & HadQ2weight_CorrUp           , "HadQ2weight_CorrUp/F"        );
  TreeHad->Branch("HadNNPDF3weight_CorrDn"          , & HadNNPDF3weight_CorrDn       , "HadNNPDF3weight_CorrDn/F"    );
  TreeHad->Branch("HadNNPDF3weight_CorrUp"          , & HadNNPDF3weight_CorrUp       , "HadNNPDF3weight_CorrUp/F"    );
  TreeHad->Branch("HadRunNum"                       , & HadRunNum                    , "HadRunNum/I"                 );
  TreeHad->Branch("HadLumiBlock"                    , & HadLumiBlock                 , "HadLumiBlock/I"              );
  TreeHad->Branch("HadEventNum"                     , & HadEventNum                  , "HadEventNum/I"               );
  TreeHad->Branch("HadPassMETFilters"               , & HadPassMETFilters            , "HadPassMETFilters/I"         );


  TreeHad->Branch("AK4_uncorr_pt",          & AK4_uncorr_pt               , "AK4_uncorr_pt/F");      
  TreeHad->Branch("AK4_corr_pt",          & AK4_corr_pt               , "AK4_corr_pt/F");  

  TreeHad->Branch("CA12JetPtRaw",          & CA12JetPtRaw               , "CA12JetPtRaw/F");                
  TreeHad->Branch("CA12JetEtaRaw",          & CA12JetEtaRaw               , "CA12JetEtaRaw/F");                 
  TreeHad->Branch("CA12JetPhiRaw",          & CA12JetPhiRaw               , "CA12JetPhiRaw/F");                 
  TreeHad->Branch("CA12JetMassRaw",          & CA12JetMassRaw               , "CA12JetMassRaw/F");                

  TreeHad->Branch("CA12JetTau1",          & CA12JetTau1               , "CA12JetTau1/F");                   
  TreeHad->Branch("CA12JetTau2",          & CA12JetTau2               , "CA12JetTau2/F");                   
  TreeHad->Branch("CA12JetTau3",          & CA12JetTau3               , "CA12JetTau3/F");                   
  TreeHad->Branch("CA12JetTau4",          & CA12JetTau4               , "CA12JetTau4/F");                   
  TreeHad->Branch("CA12JetTau32",          & CA12JetTau32               , "CA12JetTau32/F");                  
  TreeHad->Branch("CA12JetTau21",          & CA12JetTau21               , "CA12JetTau21/F");                  

  TreeHad->Branch("CA12Jetsubjet0bdisc",          & CA12Jetsubjet0bdisc               , "CA12Jetsubjet0bdisc/F");           
  TreeHad->Branch("CA12Jetsubjet1bdisc",          & CA12Jetsubjet1bdisc               , "CA12Jetsubjet1bdisc/F");           
  TreeHad->Branch("CA12Jetmaxbdisc",          & CA12Jetmaxbdisc               , "CA12Jetmaxbdisc/F");               

  TreeHad->Branch("CA12Jetsubjet0pt",          & CA12Jetsubjet0pt               , "CA12Jetsubjet0pt/F");              
  TreeHad->Branch("CA12Jetsubjet0mass",          & CA12Jetsubjet0mass               , "CA12Jetsubjet0mass/F");            
  TreeHad->Branch("CA12Jetsubjet0eta",          & CA12Jetsubjet0eta               , "CA12Jetsubjet0eta/F");             
  TreeHad->Branch("CA12Jetsubjet0phi",          & CA12Jetsubjet0phi               , "CA12Jetsubjet0phi/F");             
  TreeHad->Branch("CA12Jetsubjet0area",          & CA12Jetsubjet0area               , "CA12Jetsubjet0area/F");            

  TreeHad->Branch("CA12Jetsubjet1pt",          & CA12Jetsubjet1pt               , "CA12Jetsubjet1pt/F");              
  TreeHad->Branch("CA12Jetsubjet1mass",          & CA12Jetsubjet1mass               , "CA12Jetsubjet1mass/F");            
  TreeHad->Branch("CA12Jetsubjet1eta",          & CA12Jetsubjet1eta               , "CA12Jetsubjet1eta/F");             
  TreeHad->Branch("CA12Jetsubjet1phi",          & CA12Jetsubjet1phi               , "CA12Jetsubjet1phi/F");             
  TreeHad->Branch("CA12Jetsubjet1area",          & CA12Jetsubjet1area               , "CA12Jetsubjet1area/F");            



  TreeHad->Branch("AK4ReconstructedJetPt",          & AK4ReconstructedJetPt               , "AK4ReconstructedJetPt/F");         
  TreeHad->Branch("AK4ReconstructedJetEta",          & AK4ReconstructedJetEta               , "AK4ReconstructedJetEta/F");        
  TreeHad->Branch("AK4ReconstructedJetPhi",          & AK4ReconstructedJetPhi               , "AK4ReconstructedJetPhi/F");        
  TreeHad->Branch("AK4ReconstructedJetMass",          & AK4ReconstructedJetMass               , "AK4ReconstructedJetMass/F");       


  TreeHad->Branch("AK4bJetPtRaw",      & AK4bJetPtRaw ,             "AK4bJetPtRaw/F");       
  TreeHad->Branch("AK4bJetEtaRaw",      & AK4bJetEtaRaw ,             "AK4bJetEtaRaw/F");      
  TreeHad->Branch("AK4bJetPhiRaw",      & AK4bJetPhiRaw ,             "AK4bJetPhiRaw/F");      
  TreeHad->Branch("AK4bJetMassRaw",      & AK4bJetMassRaw ,             "AK4bJetMassRaw/F");     

  TreeHad->Branch("AK4bJet_PtSmear",      & AK4bJet_PtSmear ,             "AK4bJet_PtSmear/F");      
  TreeHad->Branch("AK4bJet_PtSmearUp",      & AK4bJet_PtSmearUp ,             "AK4bJet_PtSmearUp/F");    
  TreeHad->Branch("AK4bJet_PtSmearDn",      & AK4bJet_PtSmearDn ,             "AK4bJet_PtSmearDn/F");    
  TreeHad->Branch("AK4bJet_PtUncorr",      & AK4bJet_PtUncorr ,             "AK4bJet_PtUncorr/F");   
  TreeHad->Branch("AK4bJet_Corr",      & AK4bJet_Corr ,             "AK4bJet_Corr/F");       
  TreeHad->Branch("AK4bJet_CorrUp",      & AK4bJet_CorrUp ,             "AK4bJet_CorrUp/F");     
  TreeHad->Branch("AK4bJet_CorrDn",      & AK4bJet_CorrDn ,             "AK4bJet_CorrDn/F");     
  TreeHad->Branch("AK4bJet_bDisc",      & AK4bJet_bDisc ,             "AK4bJet_bDisc/F");   

  TreeHad->Branch("AK4WJetPtRaw",      & AK4WJetPtRaw ,             "AK4WJetPtRaw/F");       
  TreeHad->Branch("AK4WJetEtaRaw",      & AK4WJetEtaRaw ,             "AK4WJetEtaRaw/F");      
  TreeHad->Branch("AK4WJetPhiRaw",      & AK4WJetPhiRaw ,             "AK4WJetPhiRaw/F");      
  TreeHad->Branch("AK4WJetMassRaw",      & AK4WJetMassRaw ,             "AK4WJetMassRaw/F");     

  TreeHad->Branch("AK4WJet_PtSmear",      & AK4WJet_PtSmear ,             "AK4WJet_PtSmear/F");      
  TreeHad->Branch("AK4WJet_PtSmearUp",      & AK4WJet_PtSmearUp ,             "AK4WJet_PtSmearUp/F");    
  TreeHad->Branch("AK4WJet_PtSmearDn",      & AK4WJet_PtSmearDn ,             "AK4WJet_PtSmearDn/F");    
  TreeHad->Branch("AK4WJet_PtUncorr",      & AK4WJet_PtUncorr ,             "AK4WJet_PtUncorr/F");   
  TreeHad->Branch("AK4WJet_Corr",      & AK4WJet_Corr ,             "AK4WJet_Corr/F");       
  TreeHad->Branch("AK4WJet_CorrUp",      & AK4WJet_CorrUp ,             "AK4WJet_CorrUp/F");     
  TreeHad->Branch("AK4WJet_CorrDn",      & AK4WJet_CorrDn ,             "AK4WJet_CorrDn/F");     
  TreeHad->Branch("AK4WJet_bDisc",      & AK4WJet_bDisc ,             "AK4WJet_bDisc/F");      


  TreeHad->Branch("AK4W2JetPtRaw",      & AK4W2JetPtRaw ,             "AK4W2JetPtRaw/F");       
  TreeHad->Branch("AK4W2JetEtaRaw",      & AK4W2JetEtaRaw ,             "AK4W2JetEtaRaw/F");      
  TreeHad->Branch("AK4W2JetPhiRaw",      & AK4W2JetPhiRaw ,             "AK4W2JetPhiRaw/F");      
  TreeHad->Branch("AK4W2JetMassRaw",      & AK4W2JetMassRaw ,             "AK4W2JetMassRaw/F");     

  TreeHad->Branch("AK4W2Jet_PtSmear",      & AK4W2Jet_PtSmear ,             "AK4W2Jet_PtSmear/F");      
  TreeHad->Branch("AK4W2Jet_PtSmearUp",      & AK4W2Jet_PtSmearUp ,             "AK4W2Jet_PtSmearUp/F");    
  TreeHad->Branch("AK4W2Jet_PtSmearDn",      & AK4W2Jet_PtSmearDn ,             "AK4W2Jet_PtSmearDn/F");    
  TreeHad->Branch("AK4W2Jet_PtUncorr",      & AK4W2Jet_PtUncorr ,             "AK4W2Jet_PtUncorr/F");   
  TreeHad->Branch("AK4W2Jet_Corr",      & AK4W2Jet_Corr ,             "AK4W2Jet_Corr/F");       
  TreeHad->Branch("AK4W2Jet_CorrUp",      & AK4W2Jet_CorrUp ,             "AK4W2Jet_CorrUp/F");     
  TreeHad->Branch("AK4W2Jet_CorrDn",      & AK4W2Jet_CorrDn ,             "AK4W2Jet_CorrDn/F");     
  TreeHad->Branch("AK4W2Jet_bDisc",      & AK4W2Jet_bDisc ,             "AK4W2Jet_bDisc/F");      

  TreeHad->Branch("MuPhi", "vector<float>",&MuPhi);
  TreeHad->Branch("MuPt", "vector<float>",&MuPt);
  TreeHad->Branch("MuEta", "vector<float>",&MuEta);
  TreeHad->Branch("MuMass", "vector<float>",&MuMass);
  TreeHad->Branch("MuIso", "vector<float>",&MuIso);
  TreeHad->Branch("MuIsoTrk", "vector<float>",&MuIsoTrk);
  TreeHad->Branch("MuMedium", "vector<int>", &MuMedium);
  TreeHad->Branch("MuTight", "vector<int>", &MuTight);

  TreeHad->Branch("Electron_Phi", "vector<float>",&Electron_Phi);
  TreeHad->Branch("Electron_Pt", "vector<float>",&Electron_Pt);
  TreeHad->Branch("Electron_Eta", "vector<float>",&Electron_Eta);
  TreeHad->Branch("Electron_Mass", "vector<float>",&Electron_Mass);
  TreeHad->Branch("Elecron_absiso", "vector<float>",&Elecron_absiso);
  TreeHad->Branch("Elecron_relIsoWithDBeta", "vector<float>",&Elecron_relIsoWithDBeta);
  TreeHad->Branch("Elecron_absiso_EA", "vector<float>",&Elecron_absiso_EA);
  TreeHad->Branch("Elecron_relIsoWithEA", "vector<float>",&Elecron_relIsoWithEA);



TreeHad->Branch("Electron_iso_passHLTpre", "vector<int>", &Electron_iso_passHLTpre);
TreeHad->Branch("Electron_iso_passLoose", "vector<int>", &Electron_iso_passLoose);
TreeHad->Branch("Electron_iso_passMedium", "vector<int>", &Electron_iso_passMedium);
TreeHad->Branch("Electron_iso_passTight", "vector<int>", &Electron_iso_passTight);
TreeHad->Branch("Electron_iso_passHEEP", "vector<int>", &Electron_iso_passHEEP);
TreeHad->Branch("Electron_noiso_passLoose", "vector<int>", &Electron_noiso_passLoose);
TreeHad->Branch("Electron_noiso_passMedium", "vector<int>", &Electron_noiso_passMedium);
TreeHad->Branch("Electron_noiso_passTight", "vector<int>", &Electron_noiso_passTight);
TreeHad->Branch("Electron_noiso_passHEEP", "vector<int>", &Electron_noiso_passHEEP);

//  888                       888                     d8b               88888888888                       
//  888                       888                     Y8P                   888                           
//  888                       888                                           888                           
//  888      .d88b.  88888b.  888888 .d88b.  88888b.  888  .d8888b          888  888d888 .d88b.   .d88b.  
//  888     d8P  Y8b 888 "88b 888   d88""88b 888 "88b 888 d88P"             888  888P"  d8P  Y8b d8P  Y8b 
//  888     88888888 888  888 888   888  888 888  888 888 888               888  888    88888888 88888888 
//  888     Y8b.     888 d88P Y88b. Y88..88P 888  888 888 Y88b.             888  888    Y8b.     Y8b.     
//  88888888 "Y8888  88888P"   "Y888 "Y88P"  888  888 888  "Y8888P          888  888     "Y8888   "Y8888  
//                   888                                                                                  
//                   888                                                                                  
//                   888          
    
  
  TreeLept = new TTree("TreeLept","TreeLept"); 
       

  TreeLept->Branch("LeptTrigPrescales" , "vector<int>",  &LeptTrigPrescales);
  TreeLept->Branch("LeptTrigPass"      , "vector<bool>", &LeptTrigPass);
  TreeLept->Branch("LeptTrigAcceptBits", &LeptTrigAcceptBits);

  TreeLept->Branch("Gen_array_t_p4"                   , & Gen_array_t_p4                   ,  "Gen_array_t_p4[4]/F"                      );
  TreeLept->Branch("Gen_array_final_t_p4"             , & Gen_array_final_t_p4             ,  "Gen_array_final_t_p4[4]/F"                );
  TreeLept->Branch("Gen_array_b_p4"                   , & Gen_array_b_p4                   ,  "Gen_array_b_p4[4]/F"                      );
  TreeLept->Branch("Gen_array_W_p4"                   , & Gen_array_W_p4                   ,  "Gen_array_W_p4[4]/F"                      );
  TreeLept->Branch("Gen_array_Wd1_p4"                 , & Gen_array_Wd1_p4                 ,  "Gen_array_Wd1_p4[4]/F"                    );
  TreeLept->Branch("Gen_array_Wd2_p4"                 , & Gen_array_Wd2_p4                 ,  "Gen_array_Wd2_p4[4]/F"                    );
  TreeLept->Branch("Gen_array_hardest_parton_hardScatterOutgoing_p4", & Gen_array_hardest_parton_hardScatterOutgoing_p4 ,  "Gen_array_hardest_parton_hardScatterOutgoing_p4[4]/F");
  TreeLept->Branch("Gen_array_second_hardest_parton_hardScatterOutgoing_p4", & Gen_array_second_hardest_parton_hardScatterOutgoing_p4 ,  "Gen_array_second_hardest_parton_hardScatterOutgoing_p4[4]/F" );
  TreeLept->Branch("tophadronic"                         , & tophadronic                         ,  "tophadronic/O"                      );
  TreeLept->Branch("topleptonic"                         , & topleptonic                         ,  "topleptonic/O"                      );
  TreeLept->Branch("parton1id"                           , & parton1id                           ,  "parton1id/I"                        );
  TreeLept->Branch("parton2id"                           , & parton2id                           ,  "parton2id/I"                        );
  TreeLept->Branch("Wd1_id"                              , & Wd1_id                              ,  "Wd1_id/I"                           );
  TreeLept->Branch("Wd2_id"                              , & Wd2_id                              ,  "Wd2_id/I"                           ); 

  std::cout<<"Finished constructor"<<std::endl;        
}        



B2GMonoTopTreeMaker::~B2GMonoTopTreeMaker()
{
}

//  888b     d888                        888                            8888888888                         888    d8b                            
//  8888b   d8888                        888                            888                                888    Y8P                            
//  88888b.d88888                        888                            888                                888                                   
//  888Y88888P888  .d88b.  88888b.d88b.  88888b.   .d88b.  888d888      8888888 888  888 88888b.   .d8888b 888888 888  .d88b.  88888b.  .d8888b  
//  888 Y888P 888 d8P  Y8b 888 "888 "88b 888 "88b d8P  Y8b 888P"        888     888  888 888 "88b d88P"    888    888 d88""88b 888 "88b 88K      
//  888  Y8P  888 88888888 888  888  888 888  888 88888888 888          888     888  888 888  888 888      888    888 888  888 888  888 "Y8888b. 
//  888   "   888 Y8b.     888  888  888 888 d88P Y8b.     888          888     Y88b 888 888  888 Y88b.    Y88b.  888 Y88..88P 888  888      X88 
//  888       888  "Y8888  888  888  888 88888P"   "Y8888  888          888      "Y88888 888  888  "Y8888P  "Y888 888  "Y88P"  888  888  88888P' 
//                                                                                                                                               
//                                                                                                                                               
//                                                                                                                                               

//  8888888888                           888         888                                
//  888                                  888         888                                
//  888                                  888         888                                
//  8888888   888  888  .d88b.  88888b.  888888      888      .d88b.   .d88b.  88888b.  
//  888       888  888 d8P  Y8b 888 "88b 888         888     d88""88b d88""88b 888 "88b 
//  888       Y88  88P 88888888 888  888 888         888     888  888 888  888 888  888 
//  888        Y8bd8P  Y8b.     888  888 Y88b.       888     Y88..88P Y88..88P 888 d88P 
//  8888888888  Y88P    "Y8888  888  888  "Y888      88888888 "Y88P"   "Y88P"  88888P"  
//                                                                             888      
//                                                                             888      
//                                                                             888      

void
B2GMonoTopTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  using namespace pat;
  using namespace LHAPDF;
  // using namespace fastjet;


      if (nEvents%100==0) std::cout << "nEvents: " << nEvents << std::endl;
      nEvents = nEvents +1;
      
      PFMET120_BTagCSV_Mu5_Trigger = false;
      PFMET300_Trigger = false;
      HLT_PFMET120_PFMHT120_Trigger = false;
      HLT_PFMET170_Trigger =  false;

       MuPhi->clear();
       MuPt->clear();
       MuEta->clear();
       MuMass->clear();
       MuIso->clear();
       MuIsoTrk->clear();
       MuTight->clear();
       MuMedium->clear();
      

       Electron_Phi->clear();
       Electron_Pt->clear();
       Electron_Eta->clear();
       Electron_Mass->clear();
       Elecron_absiso->clear();
       Elecron_relIsoWithDBeta->clear();
       Elecron_absiso_EA->clear();
       Elecron_relIsoWithEA->clear();

       Electron_iso_passHLTpre->clear();
       Electron_iso_passLoose->clear();
       Electron_iso_passMedium->clear();
       Electron_iso_passTight->clear();
       Electron_iso_passHEEP->clear();
       Electron_noiso_passLoose->clear();
       Electron_noiso_passMedium->clear();
       Electron_noiso_passTight->clear();
       Electron_noiso_passHEEP->clear();


  AK4JetLV_pt->clear();
  AK4JetLV_eta->clear();
  AK4JetLV_phi->clear();
  AK4JetLV_mass->clear();
  AK4JetLV_corr->clear();
  AK4JetLV_corrUp->clear();
  AK4JetLV_corrDn->clear();
  AK4JetLV_SF->clear();
  AK4JetLV_SFUp->clear();
  AK4JetLV_SFDn->clear();
  AK4JetLV_ptsmear->clear();
  AK4JetLV_ptsmearUp->clear();
  AK4JetLV_ptsmearDn->clear();
  AK8JetLV_pt->clear();
  AK8JetLV_eta->clear();
  AK8JetLV_phi->clear();
  AK8JetLV_mass->clear();
  AK8JetLV_corr->clear();
  AK8JetLV_corrUp->clear();
  AK8JetLV_corrDn->clear();
  AK8JetLV_SF->clear();
  AK8JetLV_SFUp->clear();
  AK8JetLV_SFDn->clear();
  AK8JetLV_ptsmear->clear();
  AK8JetLV_ptsmearUp->clear();
  AK8JetLV_ptsmearDn->clear();
  AK8SubjetLV_pt->clear();
  AK8SubjetLV_eta->clear();
  AK8SubjetLV_phi->clear();
  AK8SubjetLV_mass->clear();
  AK4JetBtag_p->clear();
  AK8JetTau1_p->clear();
  AK8JetTau2_p->clear();
  AK8JetTau3_p->clear();
  AK8JetSoftdropMass_p->clear();


  if (verbose_) {
    cout<<"----------------------------------------------------------------------------------------------------------"<<endl;
    cout<<"Analyze event "<<iEvent.id().event()<<" run "<<iEvent.id().run()<<" lumiblock "<<iEvent.id().luminosityBlock()<<endl;
  }

  h_cutflow_had   ->Fill(0.5);
  h_cutflow_lept ->Fill(0.5);

  //
  //    .d8888b.  8888888888 888b    888     8888888b.                   888    d8b          888                   
  //   d88P  Y88b 888        8888b   888     888   Y88b                  888    Y8P          888                   
  //   888    888 888        88888b  888     888    888                  888                 888                   
  //   888        8888888    888Y88b 888     888   d88P  8888b.  888d888 888888 888  .d8888b 888  .d88b.  .d8888b  
  //   888  88888 888        888 Y88b888     8888888P"      "88b 888P"   888    888 d88P"    888 d8P  Y8b 88K      
  //   888    888 888        888  Y88888     888        .d888888 888     888    888 888      888 88888888 "Y8888b. 
  //   Y88b  d88P 888        888   Y8888     888        888  888 888     Y88b.  888 Y88b.    888 Y8b.          X88 
  //    "Y8888P88 8888888888 888    Y888     888        "Y888888 888      "Y888 888  "Y8888P 888  "Y8888   88888P' 
  //                                                          

  TLorentzVector t_p4;
  TLorentzVector final_t_p4;
  TLorentzVector b_p4;
  TLorentzVector W_p4;

  TLorentzVector Wd1_p4;
  TLorentzVector Wd2_p4;

  TLorentzVector hardest_parton_hardScatterOutgoing_p4;
  TLorentzVector second_hardest_parton_hardScatterOutgoing_p4;


  tophadronic=false;
  topleptonic=false;


  int hardest_parton_hardScatterOutgoing_pt        = 0;
  int second_hardest_parton_hardScatterOutgoing_pt = 0;

  parton1id = 0;
  parton2id = 0;
  Wd1_id = 0 ;
  Wd2_id = 0 ;

  double counttop = 0;
  if (!iEvent.isRealData() and runGenLoop_) {
    Handle<edm::View<reco::GenParticle> > genpart;
    iEvent.getByToken(prunedGenToken_,genpart);  


    // Classify the event based on the number of top quarks
    for(size_t i=0; i<genpart->size();i++){
      if (fabs((*genpart)[i].pdgId())==6 && (*genpart)[i].status()<30 && (*genpart)[i].status()>=20) counttop++;  // Z' events: status 22 = top from Z', status 52 with 2 daughters = the top that decays (after radiating a bunch of times)
    }
    if (verboseGen_) cout<<"counttop "<<counttop<<endl;
   
    // Loop over all pruned gen particles and find the 4-vectors of the top, W, B, W duaghters
    double countW = 0;
    double countb = 0;
    for(size_t i=0; i<genpart->size();i++){
      int id        = (*genpart)[i].pdgId();
      int status    = (*genpart)[i].status();
      int ndau      = (*genpart)[i].numberOfDaughters();
      double px     = (*genpart)[i].px();
      double py     = (*genpart)[i].py();
      double pz     = (*genpart)[i].pz();
      double e      = (*genpart)[i].energy();
      double m      = (*genpart)[i].mass();
      double pt     = (*genpart)[i].pt();
      double eta    = (*genpart)[i].eta();
      double phi    = (*genpart)[i].phi();
      

      // Find the particles from the hard scatter (for QCD samples)
      if (status==23 && counttop==0){
        if (pt>hardest_parton_hardScatterOutgoing_pt){
          second_hardest_parton_hardScatterOutgoing_pt = hardest_parton_hardScatterOutgoing_pt;
          second_hardest_parton_hardScatterOutgoing_p4 = hardest_parton_hardScatterOutgoing_p4;
          hardest_parton_hardScatterOutgoing_pt = pt;
          hardest_parton_hardScatterOutgoing_p4.SetPxPyPzE( px, py, pz, e );
          parton1id = id;
          if (verboseGen_) cout<<"---------- pt>hardest_parton_hardScatterOutgoing_pt - parton1id = "<<parton1id<<endl;
        }
        else if (pt>second_hardest_parton_hardScatterOutgoing_pt){
          second_hardest_parton_hardScatterOutgoing_pt = pt;
          second_hardest_parton_hardScatterOutgoing_p4.SetPxPyPzE( px, py, pz, e ); 
          parton2id = id;
          if (verboseGen_) cout<<"---------- pt>second_hardest_parton_hardScatterOutgoing_pt - parton2id = "<<parton2id<<endl;
        }
      }
      
      // Get tops from hard subprocess (for MonoTop samples)
      if (fabs(id)==6 && status<30 && status>=20) {
        t_p4.SetPxPyPzE( px, py, pz, e ); 
        parton1id = id;
        if (verboseGen_) cout<<"..top (hard)"<< " " << id <<endl;//" with pt "<<pt<<" status "<<status<<" ndau "<< ndau <<" pt "<<pt<<" eta "<<eta<<" phi "<<phi<<" parton1id = "<<parton1id<<endl;
      }


      // Get the tops which decay - record b and W information
      if (ndau==2 && fabs(id)==6){
        final_t_p4.SetPxPyPzE( px, py, pz, e ); 
        if (verboseGen_) cout<<"....two daughters top pt "<<pt<<" status "<<status<<" ndau "<< ndau <<" pt "<<pt<<" eta "<<eta<<" phi "<<phi<<endl;
        for (int daught =0; daught<2; daught++)
        {
          if ( fabs((*genpart)[i].daughter( daught )->pdgId())==5 )  b_p4.SetPxPyPzE( (*genpart)[i].daughter( daught )->px(), (*genpart)[i].daughter( daught )->py(), (*genpart)[i].daughter( daught )->pz(), (*genpart)[i].daughter( daught )->energy() );
          if ( fabs((*genpart)[i].daughter( daught )->pdgId())==24 ) W_p4.SetPxPyPzE( (*genpart)[i].daughter( daught )->px(), (*genpart)[i].daughter( daught )->py(), (*genpart)[i].daughter( daught )->pz(), (*genpart)[i].daughter( daught )->energy() );
          if (verboseGen_) cout<<"......top daughter ID "<< (*genpart)[i].daughter( daught )->pdgId() <<" pt "<< (*genpart)[i].daughter( daught )->pt()  <<endl;
        }
      }

      // Get the Ws which decay - record their daughter information
      if (ndau==2 && fabs(id)==24){
        if (verboseGen_) cout<<"....W+ with 2 daughters  id "<<id<<" status "<<status<<" ndau "<<ndau<<" pt "<<pt<<" eta "<<eta<<" phi "<<phi<<endl;
        if (verboseGen_) cout<<"......dd0 "<<(*genpart)[i].daughter( 0 )->pdgId()<<" ndau "<<(*genpart)[i].daughter( 0 )->numberOfDaughters()<<endl;
        if (verboseGen_) cout<<"......dd1 "<<(*genpart)[i].daughter( 1 )->pdgId()<<" ndau "<<(*genpart)[i].daughter( 1 )->numberOfDaughters()<<endl;
        Wd1_p4.SetPxPyPzE( (*genpart)[i].daughter( 0 )->px(), (*genpart)[i].daughter( 0 )->py(), (*genpart)[i].daughter( 0 )->pz(), (*genpart)[i].daughter( 0 )->energy() );
        Wd2_p4.SetPxPyPzE( (*genpart)[i].daughter( 1 )->px(), (*genpart)[i].daughter( 1 )->py(), (*genpart)[i].daughter( 1 )->pz(), (*genpart)[i].daughter( 1 )->energy() );
        if ( fabs( (*genpart)[i].daughter( 0 )->pdgId() ) < 6 && fabs( (*genpart)[i].daughter( 1 )->pdgId() ) < 6) tophadronic = true;  
        if ( fabs( (*genpart)[i].daughter( 0 )->pdgId() ) <= 18 && fabs( (*genpart)[i].daughter( 0 )->pdgId() ) >= 11) topleptonic = true;  
        Wd1_id = (*genpart)[i].daughter( 0 )->pdgId();
        Wd2_id = (*genpart)[i].daughter( 1 )->pdgId();
      }

    } // end genParticle loop


    if (verboseGen_)
    {
      cout<<"second_hardest_parton_hardScatterOutgoing_pt "<<second_hardest_parton_hardScatterOutgoing_pt      <<endl;                
      cout<<"second_hardest_parton_hardScatterOutgoing_p4pt "<<second_hardest_parton_hardScatterOutgoing_p4.Pt() <<endl;                      
      cout<<"second_hardest_parton_hardScatterOutgoing_eta "<<second_hardest_parton_hardScatterOutgoing_p4.Eta() <<endl;                      
      cout<<"hardest_parton_hardScatterOutgoing_pt        "<<hardest_parton_hardScatterOutgoing_pt             <<endl;  
      cout<<"hardest_parton_hardScatterOutgoing_p4pt        "<<hardest_parton_hardScatterOutgoing_p4.Pt()        <<endl;       
      cout<<"hardest_parton_hardScatterOutgoing_eta        "<<hardest_parton_hardScatterOutgoing_p4.Eta()        <<endl;       
      cout<<"parton1id = "<<parton1id<<endl;
      cout<<"parton2id = "<<parton1id<<endl;

      cout<<"tophadronic "<<tophadronic<<endl;

      cout<<"topleptonic "<<topleptonic<<endl;

      cout<<"Wd1_id "<<Wd1_id<<endl;
      cout<<"Wd2_id "<<Wd2_id<<endl;
    

      cout<<"tophadronic "<<tophadronic<<endl;
     
      cout<<"topleptonic "<<topleptonic<<endl;
    


      cout<<"t_p4   Pt "<<t_p4  .Pt()<<" Eta "<<t_p4  .Eta()<<" Phi "<<t_p4  .Phi()<<" M "<<t_p4  .M()<<endl;
  
      cout<<"b_p4   Pt "<<b_p4  .Pt()<<" Eta "<<b_p4  .Eta()<<" Phi "<<b_p4  .Phi()<<" M "<<b_p4  .M()<<endl;
      
      cout<<"W_p4   Pt "<<W_p4  .Pt()<<" Eta "<<W_p4  .Eta()<<" Phi "<<W_p4  .Phi()<<" M "<<W_p4  .M()<<endl;
     
      cout<<"Wd1_p4 Pt "<<Wd1_p4.Pt()<<" Eta "<<Wd1_p4.Eta()<<" Phi "<<Wd1_p4.Phi()<<" M "<<Wd1_p4.M()<<endl;
      cout<<"Wd2_p4 Pt "<<Wd2_p4.Pt()<<" Eta "<<Wd2_p4.Eta()<<" Phi "<<Wd2_p4.Phi()<<" M "<<Wd2_p4.M()<<endl;
     
      cout<<"counttop "<<counttop<<" countW "<<countW<<" countb "<<countb<<endl;
    }

  }



  // 
  // 888    888 888      88888888888 
  // 888    888 888          888     
  // 888    888 888          888     
  // 8888888888 888          888     
  // 888    888 888          888     
  // 888    888 888          888     
  // 888    888 888          888     
  // 888    888 88888888     888     
  // 
  
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;

  iEvent.getByToken(triggerBits_, triggerBits);
  iEvent.getByToken(triggerObjects_, triggerObjects);
  iEvent.getByToken(triggerPrescales_, triggerPrescales);


  

  const int ntrigs = trigsToRun.size();
  if (verbose_) cout<<"trigsToRun size "<<ntrigs<<endl;

  // do the same thing two different ways ( to test)
  std::bitset<38> hltbit;
  vector<bool> trigAccept;

  HadTrigPrescales   ->clear();
  LeptTrigPrescales ->clear();
  HadTrigPass        ->clear();
  LeptTrigPass      ->clear();

  HadTrigPrescalesMu   ->clear();
  HadTrigPassMu        ->clear();

  HadTrigPrescalesHad   ->clear();
  HadTrigPassHad        ->clear();
  HLTtriggers->clear();
  HLTtriggersPass->clear();
  HLTtriggersPrescales->clear();

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  if (verbose_) std::cout << "\n === TRIGGER PATHS === " << std::endl;
  int counttrigs =0;

  string name = "";
  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
     name = names.triggerName(i);
     HLTtriggers->push_back(name);

    bool pass = false;
    int accept = triggerBits->accept(i) ;
    int prescale = 0;
    if (accept ==1 ) pass = true;
    prescale = triggerPrescales->getPrescaleForIndex(i)  ;

    if(pass && (name.find("HLT_PFMET120_BTagCSV_p067_v") !=std::string::npos or name.find("HLT_PFMET120_Mu5_v") !=std::string::npos )  ) {
      PFMET120_BTagCSV_Mu5_Trigger = true;
      
    } else if(pass && (name.find("HLT_PFMET120_PFMHT120_IDTight_v") !=std::string::npos)  ) {
      HLT_PFMET120_PFMHT120_Trigger = true;
      
    } else if(pass && (name.find("HLT_PFMET300_v") !=std::string::npos)  ) {
      PFMET300_Trigger = true;
      
    }else if(pass && ( name.find("HLT_PFMET170_NoiseCleaned_v") !=std::string::npos or name.find("HLT_PFMET170_HBHECleaned_v") !=std::string::npos or name.find("HLT_PFMET170_JetIdCleaned_v") !=std::string::npos or name.find("HLT_PFMET170_BeamHaloCleaned_v") !=std::string::npos )  ) {
      HLT_PFMET170_Trigger = true;
      
    }

    HLTtriggersPass->push_back(pass);
    HLTtriggersPrescales->push_back(prescale);
  }


  // Loop over the list of triggers to save
  for (unsigned int j=0; j<trigsToRun.size(); j++){ 
    if (verbose_) cout<<"try to find "<<setw(50)<< trigsToRun[j];
    
    bool foundtrig  = false;
    bool pass = false;
    int prescale = 0;
    string name = "";
    // Loop over all triggers in the event
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
      
      name = names.triggerName(i);
      std::size_t found = name.find( trigsToRun[j] );

      // If the trigger from the trigger list is found in the event check if it passed
      if ( found !=std::string::npos ) {
        foundtrig = true;
        int accept = triggerBits->accept(i) ;
        if (accept ==1 ) pass = true;
        prescale = triggerPrescales->getPrescaleForIndex(i)  ;
        break;
      }
    }// end loop over all triggers in event

    if (verbose_ && foundtrig)  cout<<"  -> found. pass = "<<pass << ", prescale = " << prescale<<", name = "<< name  << std::endl; 
    if (verbose_ && !foundtrig) cout<<"  -> did not find "<< trigsToRun[j]<<endl;

      trigAccept.push_back(pass);
      HadTrigPrescales   ->push_back(prescale);
      LeptTrigPrescales ->push_back(prescale);
      HadTrigPass        ->push_back(pass);
      LeptTrigPass      ->push_back(pass);
      if (pass)  hltbit[counttrigs]=1;  
      counttrigs++;
  
  }// end loop over list of triggers to save in tree

  for (unsigned int j=0; j<trigsToRunMu.size(); j++){ 
    if (verbose_) cout<<"try to find "<<setw(50)<< trigsToRunMu[j];
    
    bool foundtrig  = false;
    bool pass = false;
    int prescale = 0;
    string name = "";
    // Loop over all triggers in the event
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
      
      name = names.triggerName(i);
      std::size_t found = name.find( trigsToRunMu[j] );

      // If the trigger from the trigger list is found in the event check if it passed
      if ( found !=std::string::npos ) {
        foundtrig = true;
        int accept = triggerBits->accept(i) ;
        if (accept ==1 ) pass = true;
        prescale = triggerPrescales->getPrescaleForIndex(i)  ;
        break;
      }
    }// end loop over all triggers in event

    if (verbose_ && foundtrig)  cout<<"  -> found. pass = "<<pass << ", prescale = " << prescale<<", name = "<< name  << std::endl; 
    if (verbose_ && !foundtrig) cout<<"  -> did not find "<< trigsToRunMu[j]<<endl;

      HadTrigPrescalesMu   ->push_back(prescale);
      HadTrigPassMu        ->push_back(pass);
      if (pass)  hltbit[counttrigs]=1;  
      counttrigs++; 
  }

  for (unsigned int j=0; j<trigsToRunHad.size(); j++){ 
    if (verbose_) cout<<"try to find "<<setw(50)<< trigsToRunHad[j];
    
    bool foundtrig  = false;
    bool pass = false;
    int prescale = 0;
    string name = "";
    // Loop over all triggers in the event
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
      
      name = names.triggerName(i);
      std::size_t found = name.find( trigsToRunHad[j] );

      // If the trigger from the trigger list is found in the event check if it passed
      if ( found !=std::string::npos ) {
        foundtrig = true;
        int accept = triggerBits->accept(i) ;
        if (accept ==1 ) pass = true;
        prescale = triggerPrescales->getPrescaleForIndex(i)  ;
        break;
      }
    }// end loop over all triggers in event

    if (verbose_ && foundtrig)  cout<<"  -> found. pass = "<<pass << ", prescale = " << prescale<<", name = "<< name  << std::endl; 
    if (verbose_ && !foundtrig) cout<<"  -> did not find "<< trigsToRunHad[j]<<endl;

      HadTrigPrescalesHad   ->push_back(prescale);
      HadTrigPassHad        ->push_back(pass);
      if (pass)  hltbit[counttrigs]=1;  
      counttrigs++; 
  }



  if (verbose_) {
    cout<<"trig accept vector. size = "<<trigAccept.size()<<" contents = "<<endl;
    for (unsigned int i=0; i< trigAccept.size(); i++){
      cout<<trigAccept[trigAccept.size()-1-i];
    }
    cout<<endl;
    cout<<"hlt bit = "<<endl;
    cout<<hltbit.to_string()<<endl;
  }

  HadTrigAcceptBits   = hltbit.to_string();
  LeptTrigAcceptBits = hltbit.to_string();
  

  // 
  // 888b     d888 8888888888 88888888888     888b    888          d8b                       8888888888 d8b 888 888                              
  // 8888b   d8888 888            888         8888b   888          Y8P                       888        Y8P 888 888                              
  // 88888b.d88888 888            888         88888b  888                                    888            888 888                              
  // 888Y88888P888 8888888        888         888Y88b 888  .d88b.  888 .d8888b   .d88b.      8888888    888 888 888888  .d88b.  888d888 .d8888b  
  // 888 Y888P 888 888            888         888 Y88b888 d88""88b 888 88K      d8P  Y8b     888        888 888 888    d8P  Y8b 888P"   88K      
  // 888  Y8P  888 888            888         888  Y88888 888  888 888 "Y8888b. 88888888     888        888 888 888    88888888 888     "Y8888b. 
  // 888   "   888 888            888         888   Y8888 Y88..88P 888      X88 Y8b.         888        888 888 Y88b.  Y8b.     888          X88 
  // 888       888 8888888888     888         888    Y888  "Y88P"  888  88888P'  "Y8888      888        888 888  "Y888  "Y8888  888      88888P' 

  bool passMETfilters=false;

  edm::Handle < edm::TriggerResults > metFilters;
  iEvent.getByToken(triggerResultsMETFilterToken_, metFilters);
  edm::TriggerNames const& filterTriggerNames = iEvent.triggerNames(*metFilters);

  int nMETfilters = metFilters->size();
  if (verbose_) cout<<"nMETfilters "<<nMETfilters<<endl;

  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2

  vector <string> filterFlags;
  filterFlags.push_back("Flag_goodVertices");
  if ( iEvent.isRealData() ) filterFlags.push_back("Flag_globalTightHalo2016Filter"); //Data only
  filterFlags.push_back("Flag_HBHENoiseFilter");
  filterFlags.push_back("Flag_HBHENoiseIsoFilter");
  filterFlags.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
  // filterFlags.push_back("Flag_eeBadScFilter");   // No longer suggested for Moriond 2017

  unsigned int count_matched_accept = 0;
  for (int itrig = 0; itrig != nMETfilters; ++itrig){
   std::string trigName = filterTriggerNames.triggerName(itrig);
    bool accept = metFilters->accept(itrig);
    if (verbose_) cout<<trigName<<"  "<<accept;
    if ( std::find( filterFlags.begin(), filterFlags.end(), trigName ) != filterFlags.end() ) {
        if (verbose_) cout<<"  -> matches filterFlags list ("<<trigName<<")"<<endl;
        if (accept) count_matched_accept++;
    }
    else { if (verbose_) cout<<endl;}
  }
  if (verbose_) cout<<"filterFlags.size() "<<filterFlags.size()<<" count_matched_accept "<<count_matched_accept<<endl;
  if ( count_matched_accept==filterFlags.size() ){
    passMETfilters=true;
  }
  if (verbose_) cout<<"miniAOD Flags pass? "<<passMETfilters<<endl;
 

  // RECO problemes -> apply to both data and MC. Not sure if we need these for the Summer16 samples. What about the re-miniAOD?
  Handle<bool> ifilterbadChCand;
  iEvent.getByToken(badChargedCandidateFilterToken_ , ifilterbadChCand);
  bool  filterbadChCandidate = *ifilterbadChCand;
  if (verbose_) cout <<"filterbadChCandidate "<<filterbadChCandidate<<endl;

  Handle<bool> ifilterbadPFMuon;
  iEvent.getByToken(badMuonFilterToken_, ifilterbadPFMuon);
  bool filterbadPFMuon = *ifilterbadPFMuon;
  if (verbose_) cout <<"filterbadPFMuon "<<filterbadPFMuon<<endl;

  passMETfilters = passMETfilters && filterbadChCandidate && filterbadPFMuon;
  if (verbose_) cout<<"passMETfilters = "<< passMETfilters <<endl;


  //
  //  888888 8888888888  .d8888b.      8888888b.                    888                        888          
  //    "88b 888        d88P  Y88b     888   Y88b                   888                        888          
  //     888 888        888    888     888    888                   888                        888          
  //     888 8888888    888            888   d88P  8888b.  888  888 888  .d88b.   8888b.   .d88888 .d8888b  
  //     888 888        888            8888888P"      "88b 888  888 888 d88""88b     "88b d88" 888 88K      
  //     888 888        888    888     888        .d888888 888  888 888 888  888 .d888888 888  888 "Y8888b. 
  //     88P 888        Y88b  d88P     888        888  888 Y88b 888 888 Y88..88P 888  888 Y88b 888      X88 
  //     888 8888888888  "Y8888P"      888        "Y888888  "Y88888 888  "Y88P"  "Y888888  "Y88888  88888P' 
  //   .d88P                                                    888                                         
  // .d88P"                                                Y8b d88P                                         
  //888P"                                                   "Y88P"                                          
  //

  // Run 2016 F into two JEC payloads. 
  //   IOV EF: [276831,278801] corresponds to Summer16_23Sep2016EFV3_DATA (For Runs E/early F)
  //   IOV FG: [278802,280385] corresponds to Summer16_23Sep2016GV3_DATA (For Runs lateF/G)
  std::vector<std::string>  jecPayloadsAK4chsFinal;
  std::vector<std::string>  jecPayloadsAK8chsFinal;
  std::vector<std::string>  jecPayloadsAK4pupFinal;
  std::vector<std::string>  jecPayloadsAK8pupFinal;

  int runNumber = iEvent.id().run();

  if (isRun2016F_ && runNumber < 278801 ){
    if (verbose_) cout<<"Using Run2016F early JEC"<<endl;
    jecPayloadsAK4chsFinal = jecPayloadsAK4chs_;
    jecPayloadsAK8chsFinal = jecPayloadsAK8chs_;
    jecPayloadsAK4pupFinal = jecPayloadsAK4pup_;
    jecPayloadsAK8pupFinal = jecPayloadsAK8pup_;
  } 
  else if (isRun2016F_ && runNumber > 278801 ){
    if (verbose_) cout<<"Using Run2016F late JEC"<<endl;
    jecPayloadsAK4chsFinal = jecPayloadsAK4chsSecondary_;
    jecPayloadsAK8chsFinal = jecPayloadsAK8chsSecondary_;
    jecPayloadsAK4pupFinal = jecPayloadsAK4pupSecondary_;
    jecPayloadsAK8pupFinal = jecPayloadsAK8pupSecondary_;
  } 
  else{
    if (verbose_) cout<<"Using Primary JEC from cfg"<<endl;
    jecPayloadsAK4chsFinal = jecPayloadsAK4chs_;
    jecPayloadsAK8chsFinal = jecPayloadsAK8chs_;
    jecPayloadsAK4pupFinal = jecPayloadsAK4pup_;
    jecPayloadsAK8pupFinal = jecPayloadsAK8pup_;
  }
  
  if (verbose_){
    cout<<"jecPayloadsAK4chs_.size()              "<<jecPayloadsAK4chs_.size()<<endl;
    cout<<"jecPayloadsAK4chsSecondary_.size()     "<<jecPayloadsAK4chsSecondary_.size()<<endl;
    cout<<"jecPayloadsAK4chsFinal.size()          "<<jecPayloadsAK4chsFinal.size()<<endl;
  }

  // AK4chs JEC 
  std::vector<JetCorrectorParameters> vParAK4chs;
   for ( std::vector<std::string>::const_iterator ipayload = jecPayloadsAK4chsFinal.begin(),
     ipayloadEnd = jecPayloadsAK4chsFinal.end(); ipayload != ipayloadEnd - 1; ++ipayload ) {
     if (verbose_)cout<<"AK4chs JEC txt: "<<*ipayload<<endl;
     JetCorrectorParameters pars(*ipayload);
     vParAK4chs.push_back(pars);
  }
  JetCorrectorAK4chs   = boost::shared_ptr<FactorizedJetCorrector>  ( new FactorizedJetCorrector(vParAK4chs) );
  JetCorrUncertAK4chs  = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jecPayloadsAK4chsFinal.back()));
  
  // AK8chs JEC 
  std::vector<JetCorrectorParameters> vParAK8chs;
   for ( std::vector<std::string>::const_iterator ipayload = jecPayloadsAK8chsFinal.begin(),
     ipayloadEnd = jecPayloadsAK8chsFinal.end(); ipayload != ipayloadEnd - 1; ++ipayload ) {
     JetCorrectorParameters pars(*ipayload);
     vParAK8chs.push_back(pars);
  }
  JetCorrectorAK8chs   = boost::shared_ptr<FactorizedJetCorrector>  ( new FactorizedJetCorrector(vParAK8chs) );
  JetCorrUncertAK8chs  = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jecPayloadsAK8chsFinal.back()));

  // AK4pup JEC 
  std::vector<JetCorrectorParameters> vParAK4pup;
   for ( std::vector<std::string>::const_iterator ipayload = jecPayloadsAK4pupFinal.begin(),
     ipayloadEnd = jecPayloadsAK4pupFinal.end(); ipayload != ipayloadEnd - 1; ++ipayload ) {
     JetCorrectorParameters pars(*ipayload);
     vParAK4pup.push_back(pars);
  }
  JetCorrectorAK4pup   = boost::shared_ptr<FactorizedJetCorrector>  ( new FactorizedJetCorrector(vParAK4pup) );
  JetCorrUncertAK4pup  = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jecPayloadsAK4pupFinal.back()));
  
  // AK8pup JEC 
  std::vector<JetCorrectorParameters> vParAK8pup;
   for ( std::vector<std::string>::const_iterator ipayload = jecPayloadsAK8pupFinal.begin(),
     ipayloadEnd = jecPayloadsAK8pupFinal.end(); ipayload != ipayloadEnd - 1; ++ipayload ) {
     JetCorrectorParameters pars(*ipayload);
     vParAK8pup.push_back(pars);
  }
  JetCorrectorAK8pup   = boost::shared_ptr<FactorizedJetCorrector>  ( new FactorizedJetCorrector(vParAK8pup) );
  JetCorrUncertAK8pup  = boost::shared_ptr<JetCorrectionUncertainty>( new JetCorrectionUncertainty(jecPayloadsAK8pupFinal.back()));
  
  // jet resolution from text files
  JME::JetResolution jet_resolution_AK4CHS;
  JME::JetResolution jet_resolution_AK8CHS;
  jet_resolution_AK4CHS = JME::JetResolution(jertextAK4_);
  jet_resolution_AK8CHS = JME::JetResolution(jertextAK8_);

  // jet resolution scale factor from text files
  JME::JetResolutionScaleFactor jer_scaler;
  jer_scaler = JME::JetResolutionScaleFactor(jerSFtext_);


  //
  //  888     888                  888    d8b                            
  //  888     888                  888    Y8P                            
  //  888     888                  888                                   
  //  Y88b   d88P  .d88b.  888d888 888888 888  .d8888b  .d88b.  .d8888b  
  //   Y88b d88P  d8P  Y8b 888P"   888    888 d88P"    d8P  Y8b 88K      
  //    Y88o88P   88888888 888     888    888 888      88888888 "Y8888b. 
  //     Y888P    Y8b.     888     Y88b.  888 Y88b.    Y8b.          X88 
  //      Y8P      "Y8888  888      "Y888 888  "Y8888P  "Y8888   88888P' 
  //                                                                     
  //                                                                     
                                                                     
  edm::Handle<std::vector<reco::Vertex> > vertices;
  iEvent.getByToken(vtxToken_, vertices);
  int nvtx = vertices->size();
  if (vertices->empty()) return; // skip the event if no PV found
  const reco::Vertex &PV = vertices->front();  // save PV for tight muon ID

  int nvtxgood =0 ;
  for(std::vector<reco::Vertex>::const_iterator vtx = vertices->begin(); vtx != vertices->end(); ++vtx) {
    bool isFake = (vtx->chi2()==0 && vtx->ndof()==0);   //// bool isFake = vtx->isFake();  // for AOD
    if ( !isFake &&  vtx->ndof()>=4. && vtx->position().Rho()<=2.0 && fabs(vtx->position().Z())<=24.0) {
      nvtxgood++;
    }
  }
  if (verbose_) cout<<"nvtx "<<nvtx<<" nvtxgood "<<nvtxgood<<endl;

  //
  //  8888888b.  888     888     888       888          d8b          888      888    
  //  888   Y88b 888     888     888   o   888          Y8P          888      888    
  //  888    888 888     888     888  d8b  888                       888      888    
  //  888   d88P 888     888     888 d888b 888  .d88b.  888  .d88b.  88888b.  888888 
  //  8888888P"  888     888     888d88888b888 d8P  Y8b 888 d88P"88b 888 "88b 888    
  //  888        888     888     88888P Y88888 88888888 888 888  888 888  888 888    
  //  888        Y88b. .d88P     8888P   Y8888 Y8b.     888 Y88b 888 888  888 Y88b.  
  //  888         "Y88888P"      888P     Y888  "Y8888  888  "Y88888 888  888  "Y888 
  //                                                             888                 
  //                                                        Y8b d88P                 
  //                                                         "Y88P"                  

  edm::Handle<std::vector<PileupSummaryInfo> > pileup;
  iEvent.getByToken(pileupInfoToken_, pileup);
  int nPU = 0;
  if(pileup.isValid()) { // protection for data
    for(std::vector<PileupSummaryInfo>::const_iterator iPV = pileup->begin(); iPV != pileup->end(); ++iPV) {
      if(iPV->getBunchCrossing() == 0) {
        nPU = iPV->getTrueNumInteractions();  
        //  numGenPV = iPV->getPU_NumInteractions();
        break;
      }
    }
  }

  double puweight   = hPUweight     ->GetBinContent( hPUweight     ->GetXaxis()->FindBin( nPU ) )  ;
  double puweightUp = hPUweight_MBup->GetBinContent( hPUweight_MBup->GetXaxis()->FindBin( nPU ) )  ;
  double puweightDn = hPUweight_MBdn->GetBinContent( hPUweight_MBdn->GetXaxis()->FindBin( nPU ) )  ;

  if (verbose_) std::cout<<"nPU true  "<<nPU<<"  PU weight: "<<puweight<<std::endl;

  h_NtrueIntPU         ->Fill(nPU);
  h_NPV                ->Fill(nvtx);
  h_NPVreweighted      ->Fill(nvtx,puweight);
  h_NPVgood            ->Fill(nvtxgood);
  h_NPVgoodreweighted  ->Fill(nvtxgood,puweight);

  //  888      888    888 8888888888     888       888          d8b          888      888             
  //  888      888    888 888            888   o   888          Y8P          888      888             
  //  888      888    888 888            888  d8b  888                       888      888             
  //  888      8888888888 8888888        888 d888b 888  .d88b.  888  .d88b.  88888b.  888888 .d8888b  
  //  888      888    888 888            888d88888b888 d8P  Y8b 888 d88P"88b 888 "88b 888    88K      
  //  888      888    888 888            88888P Y88888 88888888 888 888  888 888  888 888    "Y8888b. 
  //  888      888    888 888            8888P   Y8888 Y8b.     888 Y88b 888 888  888 Y88b.       X88 
  //  88888888 888    888 8888888888     888P     Y888  "Y8888  888  "Y88888 888  888  "Y888  88888P' 
  //                                                                     888                          
  //                                                                Y8b d88P                          
  //                                                                 "Y88P"                           

 double Q2wgt_up = -999;
 double Q2wgt_down = -999;

 double NNPDF3wgt_up = -999;
 double NNPDF3wgt_down = -999;


 edm::Handle<GenEventInfoProduct> genEventInfo;
 iEvent.getByToken(pdfToken_, genEventInfo);
 double evWeight = 1.0 ;
 double qScale = 1.0 ;
 double pthat = 1.0 ;
 if (genEventInfo.isValid())
 {
   evWeight = genEventInfo->weight();
   qScale   = genEventInfo->qScale();
   pthat    = (genEventInfo->hasBinningValues() ? (genEventInfo->binningValues())[0] : 0.0);

   if(verbose_) cout<<"GenEventInfo: qScale = "<<qScale<<" pthat = "<<pthat<<" evWeight = "<<evWeight<<" 1/pow(pthat/15,4.5) "<<1/pow(pthat/15,4.5)<<endl;
 }

  // 
  // 8888888b.  888               
  // 888   Y88b 888               
  // 888    888 888               
  // 888   d88P 88888b.   .d88b.  
  // 8888888P"  888 "88b d88""88b 
  // 888 T88b   888  888 888  888 
  // 888  T88b  888  888 Y88..88P 
  // 888   T88b 888  888  "Y88P"  
  //                

  Handle<double> rhoH;
  iEvent.getByToken(rhoToken_, rhoH);
  double rho = *rhoH;
  if (verbose_) cout<<"rho = "<<rho<<endl;

  //
  // 888b     d888 8888888888 88888888888 
  // 8888b   d8888 888            888     
  // 88888b.d88888 888            888     
  // 888Y88888P888 8888888        888     
  // 888 Y888P 888 888            888     
  // 888  Y8P  888 888            888     
  // 888   "   888 888            888     
  // 888       888 8888888888     888     
  //                                      

  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(metToken_, mets);
  const pat::MET &met = mets->front();
  if (verbose_){
    cout<<"MET pt "<<met.pt()<<endl;
    cout<<"MET phi "<<met.phi()<<endl;
    cout<<"MET sumEt "<<met.sumEt()<<endl;
    if (!iEvent.isRealData() )  cout<<"genMET "<< met.genMET()->pt()<<endl;
  }

//  888b     d888                            
//  8888b   d8888                            
//  88888b.d88888                            
//  888Y88888P888 888  888  .d88b.  88888b.  
//  888 Y888P 888 888  888 d88""88b 888 "88b 
//  888  Y8P  888 888  888 888  888 888  888 
//  888   "   888 Y88b 888 Y88..88P 888  888 
//  888       888  "Y88888  "Y88P"  888  888 
//                                           
//                                           
//                                           
  


  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muonToken_, muons);

  TLorentzVector mu0_p4;
  bool mu0_isTight=false;
  bool mu0_isMedium=false;
  bool mu0_isHighPt = false;
  double mu0_iso04=0;
  int count_mu=0;



  std::vector<reco::CandidatePtr> muFootprint;

  for (const pat::Muon &mu : *muons) {

      // use only loose muons 
      if (mu.pt() < 30 || !mu.isLooseMuon() || fabs( mu.eta() ) > 2.1) continue;
      // only look at 3 leading muons
      if (count_mu< 2){
        mu0_p4.SetPtEtaPhiM( mu.pt(), mu.eta(), mu.phi(), mu.mass() );

        MuPhi->push_back(mu.phi());
        MuPt->push_back(mu.pt());
        MuEta->push_back(mu.eta());
        MuMass->push_back(mu.mass());


        // Moriond 2017 short term instructions for medium muon ID
        //https://twiki.cern.ch/twiki/bin/viewauth/CMS/SWGuideMuonIdRun2
        bool goodGlob   = mu.isGlobalMuon() && 
                          mu.globalTrack()->normalizedChi2() < 3 && 
                          mu.combinedQuality().chi2LocalPosition < 12 && 
                          mu.combinedQuality().trkKink < 20; 
        bool isMediumBF = muon::isLooseMuon(mu) && 
                          mu.innerTrack()->validFraction() > 0.49 && 
                          muon::segmentCompatibility(mu) > (goodGlob ? 0.303 : 0.451); 
        bool isMediumStandard = muon::isMediumMuon(mu);
        bool isMedium   = false;
        if      (iEvent.isRealData() && runNumber <= 278808) isMedium = isMediumBF;       // Data B-F
        else if (iEvent.isRealData() && runNumber  > 278808) isMedium = isMediumStandard; // Data G-H
        else isMedium = isMediumStandard;  // MC 
        if ( isMedium ) mu0_isMedium = true;

        // Tight ID
        if ( mu.isTightMuon(PV) ) mu0_isTight = true;

        // HighPt ID
        bool isHighPt = mu.isHighPtMuon(PV);
        if ( isHighPt ) mu0_isHighPt = true;

        // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#Accessing_PF_Isolation_from_reco
        double sumChargedHadronPt = mu.pfIsolationR04().sumChargedHadronPt;
        double sumNeutralHadronPt = mu.pfIsolationR04().sumNeutralHadronEt;
        double sumPhotonPt        = mu.pfIsolationR04().sumPhotonEt;
        double sumPUPt            = mu.pfIsolationR04().sumPUPt;
        double pt                 = mu.pt();
        double iso04 = (sumChargedHadronPt+TMath::Max(0.,sumNeutralHadronPt+sumPhotonPt-0.5*sumPUPt))/pt;
        double isoForTrk03 = mu.isolationR03().sumPt/mu.pt();

        MuIso->push_back(iso04);
        MuIsoTrk->push_back(isoForTrk03);
        MuTight->push_back( (int) mu0_isTight);
        MuMedium->push_back( (int) mu0_isMedium);

        if (verbose_) cout<<"Muon pT "<<mu.pt()<<" iso04 "<<iso04<<" isMedium "<<mu0_isTight<<" isTight "<<mu0_isTight<<" isHighPt "<<mu0_isHighPt<<endl;
      } 
      count_mu++;
  }


//  8888888888 888                   888                             
//  888        888                   888                             
//  888        888                   888                             
//  8888888    888  .d88b.   .d8888b 888888 888d888 .d88b.  88888b.  
//  888        888 d8P  Y8b d88P"    888    888P"  d88""88b 888 "88b 
//  888        888 88888888 888      888    888    888  888 888  888 
//  888        888 Y8b.     Y88b.    Y88b.  888    Y88..88P 888  888 
//  8888888888 888  "Y8888   "Y8888P  "Y888 888     "Y88P"  888  888 
//                                                                   
//                                                                   
//                                                                   


  edm::Handle<edm::View<pat::Electron>> electrons; //Collection
  iEvent.getByToken(electronToken_, electrons);

  // Get the conversions collection
  edm::Handle<reco::ConversionCollection> conversions;
  iEvent.getByToken(conversionsToken_, conversions);


  edm::Handle<edm::ValueMap<vid::CutFlowResult> > cutflow_eleId_HLTpre  ;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > cutflow_eleId_Loose   ;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > cutflow_eleId_Medium  ;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > cutflow_eleId_Tight   ;
  edm::Handle<edm::ValueMap<vid::CutFlowResult> > cutflow_eleId_HEEP    ;
  iEvent.getByToken(eleIdFullInfoMapToken_HLTpre_   ,   cutflow_eleId_HLTpre   );
  iEvent.getByToken(eleIdFullInfoMapToken_Loose_    ,   cutflow_eleId_Loose   );
  iEvent.getByToken(eleIdFullInfoMapToken_Medium_   ,   cutflow_eleId_Medium  );
  iEvent.getByToken(eleIdFullInfoMapToken_Tight_    ,   cutflow_eleId_Tight   );
  iEvent.getByToken(eleIdFullInfoMapToken_HEEP_     ,   cutflow_eleId_HEEP    );



  TLorentzVector el0_p4;
  Float_t el0_absiso           =0;
  Float_t el0_relIsoWithDBeta  =0;
  Float_t el0_absiso_EA        =0;
  Float_t el0_relIsoWithEA     =0;
  int el0_iso_passHLTpre    = 0;
  int el0_iso_passLoose     = 0;
  int el0_iso_passMedium    = 0;
  int el0_iso_passTight     = 0;
  int el0_iso_passHEEP      = 0;
  int el0_noiso_passLoose   = 0;
  int el0_noiso_passMedium  = 0;
  int el0_noiso_passTight   = 0;
  int el0_noiso_passHEEP    = 0;
  int count_el=0;

  for (size_t i = 0; i < electrons->size(); ++i){   
    const auto el = electrons->ptrAt(i);          // easier if we use ptrs for the id
    if (el->pt() < 30 || fabs(el->eta())>2.4 ) continue;

     // electron ID
    vid::CutFlowResult full_cutflow_HLTpre = (*cutflow_eleId_HLTpre )[el];
    vid::CutFlowResult full_cutflow_Loose  = (*cutflow_eleId_Loose  )[el];
    vid::CutFlowResult full_cutflow_Medium = (*cutflow_eleId_Medium )[el];
    vid::CutFlowResult full_cutflow_Tight  = (*cutflow_eleId_Tight  )[el];
    vid::CutFlowResult full_cutflow_HEEP   = (*cutflow_eleId_HEEP   )[el];

    bool iso_passHLTpre   =  full_cutflow_HLTpre.cutFlowPassed();
    bool iso_passLoose    =  full_cutflow_Loose .cutFlowPassed();
    bool iso_passMedium   =  full_cutflow_Medium.cutFlowPassed();
    bool iso_passTight    =  full_cutflow_Tight .cutFlowPassed();
    bool iso_passHEEP     =  full_cutflow_HEEP  .cutFlowPassed();

    // get electron ID without isolation cuts
    vid::CutFlowResult masked_cutflow_Loose  = full_cutflow_Loose     .getCutFlowResultMasking(7); // 7 = GsfEleEffAreaPFIsoCut_0
    vid::CutFlowResult masked_cutflow_Medium = full_cutflow_Medium    .getCutFlowResultMasking(7); // 7 = GsfEleEffAreaPFIsoCut_0
    vid::CutFlowResult masked_cutflow_Tight  = full_cutflow_Tight     .getCutFlowResultMasking(7); // 7 = GsfEleEffAreaPFIsoCut_0

    std::vector<std::string> maskCuts;
    // maskCuts.push_back("GsfEleTrkPtIsoCut_0");  // OLD HEEP v6 only
    maskCuts.push_back("GsfEleValueMapIsoRhoCut_0"); // new in HEEP v7
    maskCuts.push_back("GsfEleEmHadD1IsoRhoCut_0");
    vid::CutFlowResult masked_cutflow_HEEP   = full_cutflow_HEEP      .getCutFlowResultMasking(maskCuts);

    bool noiso_passLoose    =  masked_cutflow_Loose .cutFlowPassed();
    bool noiso_passMedium   =  masked_cutflow_Medium.cutFlowPassed();
    bool noiso_passTight    =  masked_cutflow_Tight .cutFlowPassed();
    bool noiso_passHEEP     =  masked_cutflow_HEEP  .cutFlowPassed();

    if (verbose_){
      cout<<"full_cutflow_Loose"<<endl;
      printCutFlowResult(full_cutflow_Loose);
      cout<<"masked_cutflow_Loose"<<endl;
      printCutFlowResult(masked_cutflow_Loose);
      cout<<"masked_cutflow_HEEP"<<endl;
      printCutFlowResult(masked_cutflow_HEEP);
    }


    bool electronEvent = noiso_passLoose || noiso_passHEEP;
    if (electronEvent){
      if (count_el==0){
        if (verbose_) cout<<"Electron pT "<<el->pt()<<endl;

        //calculate isolation variables    
        GsfElectron::PflowIsolationVariables pfIso = el->pfIsolationVariables();
        float absiso = pfIso.sumChargedHadronPt + max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt );
        float relIsoWithDBeta = absiso/el->pt();
        float eta = el->eta();

        float effArea = 0.;
        if(abs(eta)>0.0 && abs(eta)<=1.0) effArea = 0.1752;
        if(abs(eta)>1.0 && abs(eta)<=1.479) effArea = 0.1862;
        if(abs(eta)>1.479 && abs(eta)<=2.0) effArea = 0.1411;
        if(abs(eta)>2.0 && abs(eta)<=2.2) effArea = 0.1534;
        if(abs(eta)>2.2 && abs(eta)<=2.3) effArea = 0.1903;
        if(abs(eta)>2.3 && abs(eta)<=2.4) effArea = 0.2243;
        if(abs(eta)>2.4 && abs(eta)<=2.5) effArea = 0.2687;

        float absiso_EA = pfIso.sumChargedHadronPt + max(0.0 , pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - rho * effArea );
        float relIsoWithEA = absiso_EA/el->pt();

       // el0_p4.SetPtEtaPhiM( el->pt(), el->eta(), el->phi(), el->mass() );




        el0_iso_passHLTpre    = (int) iso_passHLTpre   ;
        el0_iso_passLoose     = (int) iso_passLoose    ;
        el0_iso_passMedium    = (int) iso_passMedium   ;
        el0_iso_passTight     = (int) iso_passTight    ;
        el0_iso_passHEEP      = (int) iso_passHEEP     ;
        el0_noiso_passLoose   = (int) noiso_passLoose  ;
        el0_noiso_passMedium  = (int) noiso_passMedium ;
        el0_noiso_passTight   = (int) noiso_passTight  ;
        el0_noiso_passHEEP    = (int) noiso_passHEEP   ;


        Electron_Phi->push_back(el->phi());
        Electron_Pt->push_back(el->pt());
        Electron_Eta->push_back(el->eta());
        Electron_Mass->push_back(el->mass());

        Elecron_absiso->push_back(absiso);
        Elecron_relIsoWithDBeta->push_back(relIsoWithDBeta);
        Elecron_absiso_EA->push_back(absiso_EA);
        Elecron_relIsoWithEA->push_back(relIsoWithEA);
        Electron_iso_passHLTpre->push_back(el0_iso_passHLTpre);
        Electron_iso_passLoose->push_back(el0_iso_passLoose);
        Electron_iso_passMedium->push_back(el0_iso_passMedium);
        Electron_iso_passTight->push_back(el0_iso_passTight);
        Electron_iso_passHEEP->push_back(el0_iso_passHEEP);
        Electron_noiso_passLoose->push_back(el0_noiso_passLoose);
        Electron_noiso_passMedium->push_back(el0_noiso_passMedium);
        Electron_noiso_passTight->push_back(el0_noiso_passTight);
        Electron_noiso_passHEEP->push_back(el0_noiso_passHEEP);

      } 
      count_el++;
    }
  }



Elecron_absiso = new std::vector<float>                         ;
Elecron_relIsoWithDBeta = new std::vector<float>                ;
Elecron_absiso_EA = new std::vector<float>                      ;
Elecron_relIsoWithEA = new std::vector<float>                   ;
Electron_iso_passHLTpre = new std::vector<int>                  ;
Electron_iso_passLoose = new std::vector<int>                   ;
Electron_iso_passMedium = new std::vector<int>                  ;
Electron_iso_passTight = new std::vector<int>                   ;
Electron_iso_passHEEP = new std::vector<int>                    ;
Electron_noiso_passLoose = new std::vector<int>                 ;
Electron_noiso_passMedium = new std::vector<int>                ;
Electron_noiso_passTight = new std::vector<int>                 ;
Electron_noiso_passHEEP = new std::vector<int>                  ;


  //      
  //         d8888 888    d8P   .d8888b.       .d8888b.  888    888  .d8888b.         d8b          888             
  //        d88888 888   d8P   d88P  Y88b     d88P  Y88b 888    888 d88P  Y88b        Y8P          888             
  //       d88P888 888  d8P    Y88b. d88P     888    888 888    888 Y88b.                          888             
  //      d88P 888 888d88K      "Y88888"      888        8888888888  "Y888b.         8888  .d88b.  888888 .d8888b  
  //     d88P  888 8888888b    .d8P""Y8b.     888        888    888     "Y88b.       "888 d8P  Y8b 888    88K      
  //    d88P   888 888  Y88b   888    888     888    888 888    888       "888        888 88888888 888    "Y8888b. 
  //   d8888888888 888   Y88b  Y88b  d88P     Y88b  d88P 888    888 Y88b  d88P        888 Y8b.     Y88b.       X88 
  //  d88P     888 888    Y88b  "Y8888P"       "Y8888P"  888    888  "Y8888P"         888  "Y8888   "Y888  88888P' 
  //                                                                                  888                          
  //                                                                               d88P                          
  //                                                                             888P"                                                                                              

  edm::Handle<pat::JetCollection> AK8CHS;
  iEvent.getByToken(ak8jetToken_, AK8CHS);

  edm::Handle<reco::GenJetCollection> AK8GENJET;  
  iEvent.getByToken(ak8genjetToken_, AK8GENJET);

  edm::Handle<pat::JetCollection> AK8CHSsub;
  edm::Handle<pat::JetCollection> AK8PUPPI;
  edm::Handle<pat::JetCollection> AK8PUPPIsub;

  edm::Handle<pat::JetCollection> CA12PUPPI;
  edm::Handle<pat::JetCollection> CA12PUPPIsub;

  if (useToolbox_){
    iEvent.getByToken( ak8CHSSoftDropSubjetsToken_   , AK8CHSsub);
    iEvent.getByToken( puppijetToken_ , AK8PUPPI );
    iEvent.getByToken( ak8PuppiSoftDropSubjetsToken_ , AK8PUPPIsub);

    iEvent.getByToken( ca12puppijetToken_ , CA12PUPPI );
    iEvent.getByToken( ca12PuppiSoftDropSubjetsToken_ , CA12PUPPIsub);
    //iEvent.getByToken( ca8puppijetToken_ , CA8PUPPI );
    //iEvent.getByToken( ca8PuppiSoftDropSubjetsToken_ , CA8PUPPIsub);
  }

  int count_AK8CHS = 0;
  int count_hadAK8CHS = 0;
  int count_fill_hadTree =0;

  TLorentzVector AK8jet_had_P4corr;
  TLorentzVector AK8jet0_P4corr;
  TLorentzVector AK8jet1_P4corr;
  TLorentzVector PUPPIjet0_P4;
  TLorentzVector PUPPIjet1_P4;
  TLorentzVector PUPPIjet0_P4corr;
  TLorentzVector PUPPIjet1_P4corr;

  TLorentzVector GenJetMatched0;
  TLorentzVector GenJetMatched1;
  TLorentzVector GenJetMatchedPuppi0;
  TLorentzVector GenJetMatchedPuppi1;

  TLorentzVector leading_CA12;
  TLorentzVector leading_CA12_subjet;



  double closestAK8_to_Jet_dR=99;
  TLorentzVector closestAK8_to_Jet_P4;

  CA12JetPtRaw = 0.0;
  CA12JetEtaRaw = 0.0;
  CA12JetPhiRaw = 0.0;
  CA12JetMassRaw = 0.0;
  CA12Jetsubjet0bdisc = 0.0;
  CA12Jetsubjet1bdisc = 0.0;
  CA12Jetsubjet0pt    = 0.0;
  CA12Jetsubjet0mass  = 0.0;
  CA12Jetsubjet0eta   = 0.0;
  CA12Jetsubjet0phi   = 0.0;
  CA12Jetsubjet0area = 0.0;
  CA12Jetsubjet1pt    = 0.0;
  CA12Jetsubjet1mass  = 0.0;
  CA12Jetsubjet1eta   = 0.0;
  CA12Jetsubjet1phi   = 0.0;
  CA12Jetsubjet1area  = 0.0;
  int count_subjets = 0;


  if (verbose_) cout<<"\nAK8 jet loop"<<endl;
  if (verbose_) cout<<"\nAK8 jet size " <<AK8PUPPI->size() <<  " " << AK8CHS->size() <<endl;




if (useToolbox_){
  if (verbose_)  cout<<"   Puppi jet loop (toolbox)"<<endl; 
  int count_puppi_jet = 0;
  for (const pat::Jet &ipup : *AK8PUPPI) {
    if (verbose_)  cout<<"    puppi jet "<<count_puppi_jet<<" uncorr "<<ipup.correctedP4(0).pt()<<" corr "<< ipup.pt() << endl;
    //reco::Candidate::LorentzVector corrJet = ipup.correctedP4(2);

    reco::Candidate::LorentzVector uncorrJet = ipup.correctedP4(0);

    AK8JetLV_pt->push_back(uncorrJet.pt());
    AK8JetLV_eta->push_back(uncorrJet.eta());
    AK8JetLV_phi->push_back(uncorrJet.phi());
    AK8JetLV_mass->push_back(uncorrJet.mass());
    AK8JetTau1_p->push_back(ipup.userFloat("NjettinessAK8Puppi:tau1"));
    AK8JetTau2_p->push_back(ipup.userFloat("NjettinessAK8Puppi:tau2"));
    AK8JetTau3_p->push_back(ipup.userFloat("NjettinessAK8Puppi:tau3"));
    AK8JetSoftdropMass_p->push_back(ipup.userFloat("ak8PFJetsPuppiSoftDropMass"));

    TLorentzVector AK8PUPPI_P4uncorr;
    AK8PUPPI_P4uncorr.SetPtEtaPhiM(uncorrJet.pt(),uncorrJet.eta(),uncorrJet.phi(),uncorrJet.mass());



    //
    // AK8PUPPI JEC L23 correction
    //------------------------------------

    JetCorrectorAK8pup->setJetEta( uncorrJet.eta() );
    JetCorrectorAK8pup->setJetPt ( uncorrJet.pt() );
    JetCorrectorAK8pup->setJetE  ( uncorrJet.energy() );
    JetCorrectorAK8pup->setJetA  ( ipup.jetArea() );
    JetCorrectorAK8pup->setRho   ( rho );
    JetCorrectorAK8pup->setNPV   ( nvtx );
    vector<float> factorsAK8pup = JetCorrectorAK8pup->getSubCorrections();
    float corr_factorAK8pup_L1      = 1.0;
    float corr_factorAK8pup_L12     = 1.0;
    float corr_factorAK8pup_L123    = 1.0;
    float corr_factorAK8pup_L123res = 1.0;
    if (factorsAK8pup.size() > 0) corr_factorAK8pup_L1       = factorsAK8pup[0];
    if (factorsAK8pup.size() > 1) corr_factorAK8pup_L12      = factorsAK8pup[1];
    if (factorsAK8pup.size() > 2) corr_factorAK8pup_L123     = factorsAK8pup[2];
    if (factorsAK8pup.size() > 3) corr_factorAK8pup_L123res  = factorsAK8pup[3];
    double corr_factorAK8pup_L2 = corr_factorAK8pup_L12/corr_factorAK8pup_L1;
    double corr_factorAK8pup_L3 = corr_factorAK8pup_L123/corr_factorAK8pup_L12;
    double corr_factorAK8pup_res = corr_factorAK8pup_L123res/corr_factorAK8pup_L123;
    //double corr_factor_L23 = corr_factor_L2*corr_factor_L3;
    double corr_factorAK8pup_L23res = corr_factorAK8pup_L2*corr_factorAK8pup_L3*corr_factorAK8pup_res;

    TLorentzVector AK8PUPPI_P4corr;
    AK8PUPPI_P4corr = corr_factorAK8pup_L23res *  AK8PUPPI_P4uncorr;

    double puppi_pt = AK8PUPPI_P4corr.Pt();
    double puppi_eta = AK8PUPPI_P4corr.Eta();
    double puppi_phi = AK8PUPPI_P4corr.Phi();
    double puppi_mass = AK8PUPPI_P4corr.M();

    if(count_puppi_jet==0) PUPPIjet0_P4corr = AK8PUPPI_P4corr;
    if(count_puppi_jet==1) PUPPIjet1_P4corr = AK8PUPPI_P4corr;
    count_puppi_jet++;



    //------------------------------------
    // AK8PUPPI JEC uncertainty
    //------------------------------------
    double corrDn_pup_L23  = 1.0;
    double corrUp_pup_L23  = 1.0;

    if (puppi_pt>10){   //make sure puppi jet was found and filled
        JetCorrUncertAK8pup->setJetPhi(  puppi_phi  );
        JetCorrUncertAK8pup->setJetEta(  puppi_eta  );
        JetCorrUncertAK8pup->setJetPt(   puppi_pt * corr_factorAK8pup_L23res   );
        corrDn_pup_L23   = corr_factorAK8pup_L23res - JetCorrUncertAK8pup->getUncertainty(0);
        JetCorrUncertAK8pup->setJetPhi(   puppi_phi  );
        JetCorrUncertAK8pup->setJetEta(   puppi_eta  );
        JetCorrUncertAK8pup->setJetPt(    puppi_pt * corr_factorAK8pup_L23res   );
        corrUp_pup_L23   = corr_factorAK8pup_L23res + JetCorrUncertAK8pup->getUncertainty(1);;
    }
    if (verbose_){
      cout<<"    -> puppi uncorr pt "<<AK8PUPPI_P4uncorr.Perp()<<" corr pt "<<AK8PUPPI_P4corr.Perp() <<endl;
      cout<<"    -> corr L2L3res " <<corr_factorAK8pup_L23res<<" corr_factorAK8pup_L1 "<<corr_factorAK8pup_L1<<" corrDn_pup_L23"<<corrDn_pup_L23<<" corrUp_pup_L23 "<<corrUp_pup_L23<<endl;
    }

    AK8JetLV_corr->push_back(corr_factorAK8pup_L23res);
    AK8JetLV_corrUp->push_back(corrUp_pup_L23);
    AK8JetLV_corrDn->push_back(corrDn_pup_L23);


   TLorentzVector GenJetMatched;


   if (!iEvent.isRealData()) {
     if (verbose_) cout<<"   Get JER SF"<<endl;

     // get genjet
     const reco::GenJet* genJet = ipup.genJet();
//     bool foundgenjet = false;
     if (genJet) {
       //foundgenjet=true;
       //genpt = genJet->pt();
       GenJetMatched.SetPtEtaPhiM( genJet->pt(), genJet->eta(), genJet->phi(), genJet->mass() );
       if (verbose_) cout<<"     -> found ak8 genJet pt "<<genJet->pt()<<" mass "<<genJet->mass()<<endl;
     }
    }


    //------------------------------------
    // AK8PUPPI JER SF
    //------------------------------------
    TLorentzVector GenJetMatched_Pup;

    double pup_ptsmear   = 1;
    double pup_ptsmearUp = 1;
    double pup_ptsmearDn = 1;
    if (!iEvent.isRealData()) {
      double jer_sf    = jer_scaler.getScaleFactor({{JME::Binning::JetEta, puppi_eta  }});
      double jer_sf_up = jer_scaler.getScaleFactor({{JME::Binning::JetEta, puppi_eta  }}, Variation::UP);
      double jer_sf_dn = jer_scaler.getScaleFactor({{JME::Binning::JetEta, puppi_eta  }}, Variation::DOWN);
      if (verbose_) std::cout << "   PUPPI JER Scale factors (Nominal / Up / Down) : " << jer_sf << " / " << jer_sf_up << " / " << jer_sf_dn << std::endl;
      double recopt    = AK8PUPPI_P4corr.Perp();
      double genpt     = GenJetMatched.Perp();
      double deltapt   = (recopt-genpt)*(jer_sf-1.0);
      double deltaptUp = (recopt-genpt)*(jer_sf_up-1.0);
      double deltaptDn = (recopt-genpt)*(jer_sf_dn-1.0);
      pup_ptsmear   = std::max((double)0.0, (recopt+deltapt)/recopt     );
      pup_ptsmearUp = std::max((double)0.0, (recopt+deltaptUp)/recopt   );
      pup_ptsmearDn = std::max((double)0.0, (recopt+deltaptDn)/recopt   );

      AK8JetLV_SF->push_back(jer_sf);
      AK8JetLV_SFUp->push_back(jer_sf_up);
      AK8JetLV_SFDn->push_back(jer_sf_dn);

      AK8JetLV_ptsmear->push_back(pup_ptsmear);
      AK8JetLV_ptsmearUp->push_back(pup_ptsmearUp);
      AK8JetLV_ptsmearDn->push_back(pup_ptsmearDn);

      if (verbose_){
        cout<<"    -> AK8PUPPI JER recopt  "<<recopt <<endl;  
        cout<<"    -> AK8PUPPI JER genpt   "<<genpt  <<endl;  
        cout<<"    -> AK8PUPPI JER deltapt "<<deltapt<<endl;  
        cout<<"    -> AK8PUPPI JER ptsmear "<<pup_ptsmear<<endl;
        cout<<"    -> AK8PUPPI JER pup_ptsmearUp "<<pup_ptsmearUp<<endl;
        cout<<"    -> AK8PUPPI JER pup_ptsmearDn "<<pup_ptsmearDn<<endl;
      }
    }

  }
}


  // 
  //        d8888 888    d8P      d8888       .d8888b.  888    888  .d8888b.         d8b          888             
  //       d88888 888   d8P      d8P888      d88P  Y88b 888    888 d88P  Y88b        Y8P          888             
  //      d88P888 888  d8P      d8P 888      888    888 888    888 Y88b.                          888             
  //     d88P 888 888d88K      d8P  888      888        8888888888  "Y888b.         8888  .d88b.  888888 .d8888b  
  //    d88P  888 8888888b    d88   888      888        888    888     "Y88b.       "888 d8P  Y8b 888    88K      
  //   d88P   888 888  Y88b   8888888888     888    888 888    888       "888        888 88888888 888    "Y8888b. 
  //  d8888888888 888   Y88b        888      Y88b  d88P 888    888 Y88b  d88P        888 Y8b.     Y88b.       X88 
  // d88P     888 888    Y88b       888       "Y8888P"  888    888  "Y8888P"         888  "Y8888   "Y888  88888P' 
  //                                                                                 888                          
  //                                                                                d88P                          

TLorentzVector AK4_p4[5];

edm::Handle<pat::JetCollection> AK4MINI;
iEvent.getByToken(ak4jetToken_, AK4MINI);

edm::Handle<reco::GenJetCollection> AK4GENJET;  
iEvent.getByToken(ak4genjetToken_, AK4GENJET);

TLorentzVector reconstructed_top;
TLorentzVector AK4_b;
TLorentzVector AK4_W;
TLorentzVector AK4_W2;

AK4_uncorr_pt = 0.0;
AK4_corr_pt = 0.0;

AK4ReconstructedJetPt = 0.0;            
AK4ReconstructedJetEta = 0.0;      
AK4ReconstructedJetPhi = 0.0;      
AK4ReconstructedJetMass = 0.0;      

AK4bJetPtRaw = 0.0;              
AK4bJetEtaRaw = 0.0;       
AK4bJetPhiRaw = 0.0;       
AK4bJetMassRaw = 0.0;      

AK4bJet_PtSmear = 0.0;                  
AK4bJet_PtSmearUp = 0.0;                
AK4bJet_PtSmearDn = 0.0;                
AK4bJet_PtUncorr = 0.0;    
AK4bJet_Corr = 0.0;        
AK4bJet_CorrUp = 0.0;      
AK4bJet_CorrDn = 0.0;     
AK4bJet_bDisc = 0.0;     

AK4WJetPtRaw = 0.0;              
AK4WJetEtaRaw = 0.0;       
AK4WJetPhiRaw = 0.0;       
AK4WJetMassRaw = 0.0;      

AK4WJet_PtSmear = 0.0;                  
AK4WJet_PtSmearUp = 0.0;                
AK4WJet_PtSmearDn = 0.0;                
AK4WJet_PtUncorr = 0.0;    
AK4WJet_Corr = 0.0;        
AK4WJet_CorrUp = 0.0;      
AK4WJet_CorrDn = 0.0;     
AK4WJet_bDisc = 0.0;     

AK4W2JetPtRaw = 0.0;              
AK4W2JetEtaRaw = 0.0;       
AK4W2JetPhiRaw = 0.0;       
AK4W2JetMassRaw = 0.0;      

AK4W2Jet_PtSmear = 0.0;                  
AK4W2Jet_PtSmearUp = 0.0;                
AK4W2Jet_PtSmearDn = 0.0;                
AK4W2Jet_PtUncorr = 0.0;    
AK4W2Jet_Corr = 0.0;        
AK4W2Jet_CorrUp = 0.0;      
AK4W2Jet_CorrDn = 0.0;     
AK4W2Jet_bDisc = 0.0;     

  if( abs(AK8jet_had_P4corr.M() -171) < abs(leading_CA12.M() -171) ){
    reconstructed_top = AK8jet_had_P4corr;
  } else {
    reconstructed_top = leading_CA12;
  }
int count_AK4CHS = 0;
for (const pat::Jet &ijet : *AK4MINI) { 
  if (ijet.pt()<15 || fabs(ijet.eta())>3.0) continue; 

  reco::Candidate::LorentzVector corrJet = ijet.correctedP4(2);

    if( AK4_p4[count_AK4CHS].DeltaR(reconstructed_top) < 2*174/leading_CA12.Pt()){

      if (ijet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") >AK4bJet_bDisc &&  ijet.pt() > 30){
         AK4bJetPtRaw = corrJet.pt();       
         AK4bJetEtaRaw = corrJet.eta();      
         AK4bJetPhiRaw = corrJet.phi();      
         AK4bJetMassRaw = corrJet.mass();     
         AK4bJet_bDisc = ijet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
         AK4_b.SetPtEtaPhiM(corrJet.pt(),corrJet.eta(),corrJet.phi(),corrJet.mass());
       } else {
         break;
       }


    }
    count_AK4CHS++;
  }

count_AK4CHS = 0;

double CSVv2L = 0.605;

if (verbose_) cout<<"AK4 jet loop"<<endl;
for (const pat::Jet &ijet : *AK4MINI) { 
  if (ijet.pt()<30 || fabs(ijet.eta())>3.0) continue; 

    reco::Candidate::LorentzVector uncorrJet = ijet.correctedP4(0);
    //reco::Candidate::LorentzVector corrJet = ijet.correctedP4(3);

    AK4JetLV_pt->push_back(uncorrJet.pt());
    AK4JetLV_eta->push_back(uncorrJet.eta());
    AK4JetLV_phi->push_back(uncorrJet.phi());
    AK4JetLV_mass->push_back(uncorrJet.mass());

    double CSVv2 = ijet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
    AK4JetBtag_p->push_back(CSVv2);

    //------------------------------------
    // AK4CHS JEC correction 
    //------------------------------------
    JetCorrectorAK4chs->setJetEta( uncorrJet.eta() );
    JetCorrectorAK4chs->setJetPt ( uncorrJet.pt() );
    JetCorrectorAK4chs->setJetE  ( uncorrJet.phi() );
    JetCorrectorAK4chs->setJetA  ( ijet.jetArea() );
    JetCorrectorAK4chs->setRho   ( rho );
    JetCorrectorAK4chs->setNPV   ( nvtx );
    double corr = JetCorrectorAK4chs->getCorrection();
    reco::Candidate::LorentzVector corrJet = corr * uncorrJet;

    AK4JetLV_corr->push_back(corr);

    if (CSVv2 > CSVv2L && corrJet.pt() > 30 && fabs(corrJet.eta())<2.5){
          bJet_count++;
    }

    if ( verbose_ ) cout << "   -> after JEC pt,eta,phi,m = " << corrJet.pt() << ", " << corrJet.eta() << ", " << corrJet.phi() << ", " << corrJet.mass() << endl;
    
    if (corrJet.pt()<15 ) continue;  

    //------------------------------------
    // AK4CHS JEC uncertainty
    //------------------------------------
    double corrDn = 1.0;
    JetCorrUncertAK4chs->setJetPhi(  corrJet.phi()  );
    JetCorrUncertAK4chs->setJetEta(  corrJet.eta()  );
    JetCorrUncertAK4chs->setJetPt(   corrJet.pt()   );
    corrDn = corr - JetCorrUncertAK4chs->getUncertainty(0);
    double corrUp = 1.0;
    JetCorrUncertAK4chs->setJetPhi(  corrJet.phi()  );
    JetCorrUncertAK4chs->setJetEta(  corrJet.eta()  );
    JetCorrUncertAK4chs->setJetPt(   corrJet.pt()   );
    corrUp = corr + JetCorrUncertAK4chs->getUncertainty(1);

    AK4JetLV_corrUp->push_back(corrUp);
    AK4JetLV_corrDn->push_back(corrDn);


    if (verbose_) cout<<"   -> corr "<<corr<<" corrDn "<<corrDn<<" corrUp "<< corrUp<<endl;

    //------------------------------------
    // AK4 JER SF
    //------------------------------------
   
    double ptsmear   = 1;
    double ptsmearUp = 1;
    double ptsmearDn = 1;

    if (!iEvent.isRealData()) {
      if (verbose_) cout<<"   Get JER SF"<<endl;

      // get genjet
      double genpt = 0;
      TLorentzVector GenJetMatched;
      const reco::GenJet* genJet = ijet.genJet();
      bool foundgenjet = false;
      if (genJet) {
        foundgenjet=true;
        genpt = genJet->pt();
        GenJetMatched.SetPtEtaPhiM( genJet->pt(), genJet->eta(), genJet->phi(), genJet->mass() );
        reco::Candidate::LorentzVector corrJet = ijet.correctedP4(2);

       if(verbose_) cout << "cor vs un corr " << corrJet.pt() - GenJetMatched.Pt() << " " << uncorrJet.pt() - GenJetMatched.Pt() << endl;

        AK4_uncorr_pt = uncorrJet.pt() - GenJetMatched.Pt();
        AK4_corr_pt = corrJet.pt() - GenJetMatched.Pt();

        if (verbose_) cout<<"      -> Found ak4 genJet pt "<<genJet->pt()<<" mass "<<genJet->mass()<<endl;
      }
      else{ if(verbose_)cout<<"      -> Did not find genJet"<<endl;}
    
      // Set parameters needed for jet resolution and scale factors
      JME::JetParameters jer_parameters;
      jer_parameters.setJetPt ( corrJet.pt()  );
      jer_parameters.setJetEta( corrJet.eta() );
      jer_parameters.setRho   ( rho           );

      // Get resolution
      double res = jet_resolution_AK4CHS.getResolution(jer_parameters); 

      // Get scale factors
      double jer_sf    = jer_scaler.getScaleFactor(jer_parameters                   );
      double jer_sf_up = jer_scaler.getScaleFactor(jer_parameters , Variation::UP   );
      double jer_sf_dn = jer_scaler.getScaleFactor(jer_parameters , Variation::DOWN );
      if (verbose_) std::cout << "      -> JER Scale factors (Nominal / Up / Down) : " << jer_sf << " / " << jer_sf_up << " / " << jer_sf_dn <<"    & Resolution :"<<res<< std::endl;

      AK4JetLV_SF->push_back(jer_sf);
      AK4JetLV_SFUp->push_back(jer_sf_up);
      AK4JetLV_SFDn->push_back(jer_sf_dn);
     
      // Get Smearings  
      // --- If well matched, smear based on GenJet, If not well matched,  gaussian smear based on resolution
      TLorentzVector AK4JetP4;
      AK4JetP4.SetPtEtaPhiM( corrJet.pt(), corrJet.eta(), corrJet.phi(), corrJet.mass() );
      double DeltaR_gen_reco  = AK4JetP4.DeltaR( GenJetMatched );
      double DeltaPt_gen_reco = AK4JetP4.Pt() - GenJetMatched.Pt()  ;
      double jet_distance_param = 0.4; 
      if (verbose_) cout<<"      -> gen pt "<<GenJetMatched.Pt()<<" reco pt "<<AK4JetP4.Pt()<<"  delta "<<DeltaPt_gen_reco<<endl;

      if (genJet && (DeltaR_gen_reco<jet_distance_param/2.0) && (std::abs(DeltaPt_gen_reco)<(3*res*AK4JetP4.Pt())) ) {
        if (verbose_) cout<<"      -> Well matched (recojet,genjet)"<<endl;
        double recopt    = corrJet.pt();
        // double genpt     = GenJetMatched.Perp();
        double deltapt   = (recopt-genpt)*(jer_sf-1.0);
        double deltaptUp = (recopt-genpt)*(jer_sf_up-1.0);
        double deltaptDn = (recopt-genpt)*(jer_sf_dn-1.0);

        ptsmear   = std::max((double)0.0, (recopt+deltapt)/recopt     );
        ptsmearUp = std::max((double)0.0, (recopt+deltaptUp)/recopt   );
        ptsmearDn = std::max((double)0.0, (recopt+deltaptDn)/recopt   );
      }
      else{
        if (verbose_){
          cout<<"      -> Not well matched. DeltaR_gen_reco "<<DeltaR_gen_reco<<" DeltaPt_gen_reco "<<DeltaPt_gen_reco<<" 3*res*AK4JetP4.Pt()) "<<3*res*AK4JetP4.Pt();
          if (!foundgenjet) cout<<". Did not find genjet"<<endl;
          else cout<<endl;
        }
        double sigma   = std::sqrt(jer_sf * jer_sf - 1)       * res ;  
        double sigmaUp = std::sqrt(jer_sf_up * jer_sf_up - 1) * res ;
        double sigmaDn = std::sqrt(jer_sf_dn * jer_sf_dn - 1) * res ;

        TRandom3 rand1(0);
        ptsmear   = std::max( (double)0.0, 1 + rand1.Gaus(0, sigma  ) );
        ptsmearUp = std::max( (double)0.0, 1 + rand1.Gaus(0, sigmaUp) );
        ptsmearDn = std::max( (double)0.0, 1 + rand1.Gaus(0, sigmaDn) );
      }
    }

    AK4JetLV_ptsmear->push_back(ptsmear);
    AK4JetLV_ptsmearUp->push_back(ptsmearUp);
    AK4JetLV_ptsmearDn->push_back(ptsmearDn);

    if (verbose_) cout<<"   -> ptsmear "<<ptsmear<<" ptsmearUp "<<ptsmearUp<<" ptsmearDn "<< ptsmearDn<<endl;
 
  if(count_AK4CHS < 5){
    reco::Candidate::LorentzVector corrJet = ijet.correctedP4(2);
    AK4_p4[count_AK4CHS].SetPtEtaPhiM(corrJet.pt(), corrJet.eta(), corrJet.phi(), corrJet.mass() );


   if(AK4_p4[count_AK4CHS].DeltaR(reconstructed_top) < 2*174/leading_CA12.Pt()){


     if(corrJet.pt() > 30  && AK4WJetPtRaw == 0.0 && AK4bJetPtRaw!=corrJet.pt()){
        AK4WJetPtRaw = corrJet.pt();       
        AK4WJetEtaRaw = corrJet.eta();      
        AK4WJetPhiRaw = corrJet.phi();      
        AK4WJetMassRaw = corrJet.mass();     
        AK4WJet_bDisc = ijet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        AK4_W.SetPtEtaPhiM(corrJet.pt(),corrJet.eta(),corrJet.phi(),corrJet.mass());
     } else if(corrJet.pt() > 30 &&  AK4W2JetPtRaw == 0.0 && AK4bJetPtRaw!=corrJet.pt()){
        AK4W2JetPtRaw = corrJet.pt();       
        AK4W2JetEtaRaw = corrJet.eta();      
        AK4W2JetPhiRaw = corrJet.phi();      
        AK4W2JetMassRaw = corrJet.mass();     
        AK4W2Jet_bDisc = ijet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        AK4_W2.SetPtEtaPhiM(corrJet.pt(),corrJet.eta(),corrJet.phi(),corrJet.mass());

     }


   }

  } 

  count_AK4CHS++;
}
if (count_AK4CHS > 0){

  reconstructed_top = (AK4_b+AK4_W+AK4_W2);

  AK4ReconstructedJetPt = reconstructed_top.Pt();
  AK4ReconstructedJetEta = reconstructed_top.Eta();
  AK4ReconstructedJetPhi = reconstructed_top.Phi();
  AK4ReconstructedJetMass = reconstructed_top.M();
 
}

  setVector0(Gen_array_t_p4);
  setVector0(Gen_array_t_p4);
  setVector0(Gen_array_final_t_p4);
  setVector0(Gen_array_b_p4);
  setVector0(Gen_array_W_p4);
  setVector0(Gen_array_Wd1_p4);
  setVector0(Gen_array_Wd2_p4);
  setVector0(Gen_array_hardest_parton_hardScatterOutgoing_p4);
  setVector0(Gen_array_second_hardest_parton_hardScatterOutgoing_p4);

  if(runHadTree_){
    if(runGenLoop_){
      setVectorTL(Gen_array_t_p4, t_p4);
      setVectorTL(Gen_array_t_p4, t_p4);
      setVectorTL(Gen_array_final_t_p4, final_t_p4);
      setVectorTL(Gen_array_b_p4, b_p4);
      setVectorTL(Gen_array_W_p4, W_p4);
      setVectorTL(Gen_array_Wd1_p4, Wd1_p4);
      setVectorTL(Gen_array_Wd2_p4, Wd2_p4);
      setVectorTL(Gen_array_hardest_parton_hardScatterOutgoing_p4, hardest_parton_hardScatterOutgoing_p4);
      setVectorTL(Gen_array_second_hardest_parton_hardScatterOutgoing_p4, second_hardest_parton_hardScatterOutgoing_p4);

    }

    h_cutflow_had   ->Fill(1.5);

    if ( !iEvent.isRealData() ){

      Gen_MET_pT = met.genMET()->pt();

      if(PFMET120_BTagCSV_Mu5_Trigger){
        h_trigger_efficency_1->Fill(met.pt());
        h_trigger_accept_1->Fill(met.pt());
        h_trigger_efficency_1_topPt->Fill(t_p4.Pt());
        h_trigger_accept_1_topPt->Fill(t_p4.Pt());
      } else {
        h_trigger_reject_1->Fill(met.pt());
        h_trigger_reject_1_topPt->Fill(t_p4.Pt());
      }
      if(PFMET300_Trigger){
        h_trigger_efficency_2->Fill(met.pt());
        h_trigger_accept_2->Fill(met.pt());
        h_trigger_efficency_2_topPt->Fill(t_p4.Pt());
        h_trigger_accept_2_topPt->Fill(t_p4.Pt());
      } else {
        h_trigger_reject_2->Fill(met.pt());
        h_trigger_reject_2_topPt->Fill(t_p4.Pt());
      }
      if(HLT_PFMET120_PFMHT120_Trigger){
        h_trigger_efficency_3->Fill(met.pt());
        h_trigger_accept_3->Fill(met.pt());
        h_trigger_efficency_3_topPt->Fill(t_p4.Pt());
        h_trigger_accept_3_topPt->Fill(t_p4.Pt());
      } else {
        h_trigger_reject_3->Fill(met.pt());
        h_trigger_reject_3_topPt->Fill(t_p4.Pt());
      } 
      if(HLT_PFMET170_Trigger){
        h_trigger_efficency_4->Fill(met.pt());
        h_trigger_accept_4->Fill(met.pt());
        h_trigger_efficency_4_topPt->Fill(t_p4.Pt());
        h_trigger_accept_4_topPt->Fill(t_p4.Pt());
      } else {
        h_trigger_reject_4->Fill(met.pt());
        h_trigger_reject_4_topPt->Fill(t_p4.Pt());
      }

    }

    if ( bJet_count > 0){ 
      h_cutflow_had   ->Fill(2.5);

      if ( met.pt() > 75){ 
        h_cutflow_had   ->Fill(3.5);

        if (verbose_) cout<<"Write Had Tree"<<endl;


        if (!iEvent.isRealData() )  {
          Gen_MET_pT = met.genMET()->pt();
          Gen_MET_phi = met.genMET()->phi();
          Gen_MET_eta = met.genMET()->eta();
        }

        HadMETpx                = met.px();                   
        HadMETpy                = met.py();                   
        HadMETpt                = met.pt();                   
        HadMETphi               = met.phi();                   
        HadMETsumET             = met.sumEt();   
    
        if ( !iEvent.isRealData() )  HadMETgenMET            = met.genMET()->pt();                   
        HadMETuncorPt           = met.uncorPt();                    
           
        HadMETshiftedPtJetEnUp  = met.shiftedPt(pat::MET::JetEnUp            ) ;                    
        HadMETshiftedPtJetEnDn  = met.shiftedPt(pat::MET::JetEnDown          ) ;                    
        HadMETshiftedPtElEnUp   = met.shiftedPt(pat::MET::ElectronEnUp       ) ;                    
        HadMETshiftedPtElEnDn   = met.shiftedPt(pat::MET::ElectronEnDown     ) ;                    
        HadMETshiftedPtMuEnUp   = met.shiftedPt(pat::MET::MuonEnUp           ) ;                    
        HadMETshiftedPtMuEnDn   = met.shiftedPt(pat::MET::MuonEnDown         ) ;                    
        HadMETshiftedPtJetResUp = met.shiftedPt(pat::MET::JetResUp           ) ;                    
        HadMETshiftedPtJetResDn = met.shiftedPt(pat::MET::JetResDown         ) ;                    
        HadMETshiftedPtUnclEnUp = met.shiftedPt(pat::MET::UnclusteredEnUp    ) ;                    
        HadMETshiftedPtUnclEnDn = met.shiftedPt(pat::MET::UnclusteredEnDown  ) ;                    
        HadNvtx                 = nvtx;     
        HadNvtxGood             = nvtxgood;     
        HadNPUtrue              = nPU;     
        HadRho                  = rho ;               
        if ( !iEvent.isRealData() ){
          HadEventWeight          = evWeight ;              
          HadPUweight             = puweight  ; 
          HadPUweight_MBup        = puweightUp ;
          HadPUweight_MBdn        = puweightDn  ;
         
          HadQ2weight_CorrDn      = Q2wgt_down ;              
          HadQ2weight_CorrUp      = Q2wgt_up ;              
          HadNNPDF3weight_CorrDn  = NNPDF3wgt_down ;              
          HadNNPDF3weight_CorrUp  = NNPDF3wgt_up ;   
          
        }  
        else{ 
          HadEventWeight          = 1;    
          HadPUweight             = 1;
          HadPUweight_MBup        = 1;
          HadPUweight_MBdn        = 1;
          HadGenTTmass            = 0;
          HadGenCountHadTop       = 0;            
          HadQ2weight_CorrDn      = 1;       
          HadQ2weight_CorrUp      = 1;     
          HadNNPDF3weight_CorrDn  = 1;           
          HadNNPDF3weight_CorrUp  = 1;
        }
       

              
        HadRunNum               = iEvent.id().run() ;              
        HadLumiBlock            = iEvent.id().luminosityBlock() ;              
        HadEventNum             = iEvent.id().event() ; 
        HadPassMETFilters       = (int) passMETfilters;              

        TreeHad -> Fill();
      }// end met selection
    }// end bjet selection

 
    
  }

  jecPayloadsAK4chsFinal  .clear();
  jecPayloadsAK8chsFinal  .clear();
  jecPayloadsAK4pupFinal  .clear();
  jecPayloadsAK8pupFinal  .clear();

} //end Event loop

// ------------ method called once each job just before starting event loop  ------------
void 
B2GMonoTopTreeMaker::beginJob()
{
  fPUweight = new TFile("PUweight_FinalJSON2016_PileupJSONFeb2017_Xsec69200nominal_MCRunIISummer16MiniAODv2_PUMoriond17.root") ;
  hPUweight      = (TH1D*) fPUweight->Get("PUweight_true");
  hPUweight_MBup = (TH1D*) fPUweight->Get("PUweight_true_MBup");
  hPUweight_MBdn = (TH1D*) fPUweight->Get("PUweight_true_MBdn");
       
  std::cout<<"Test PU reweight file weights: "<<std::endl;
  std::cout<<std::setw(5)<<" NPU "<<std::setw(10)<<" weight "<<std::endl;
  for (int i=1; i<=50; i++ ){
      std::cout<<std::setw(5)<<hPUweight->GetBinLowEdge( i )<<"   "<<std::setw(10)<<hPUweight->GetBinContent( i ) <<std::endl;
  } 
}

void computeEfficencyHist(TH1D* accept, TH1D* reject){
  if (accept->GetNbinsX() != reject->GetNbinsX()) {
    std::cout << "not equal bins" << std::endl;
    return;
  }
  double n_accept =0;
  double n_reject =0;
  for( int i=-1; i < accept->GetNbinsX()+1; i++){
    n_accept = accept->GetBinContent(i);
    n_reject = reject->GetBinContent(i);

    if (n_accept == 0){
      accept->SetBinContent(i, 0);
    } else {
      n_accept = n_accept/(n_accept+n_reject);
      accept->SetBinContent(i, n_accept);
    }
  }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
B2GMonoTopTreeMaker::endJob() 
{
  computeEfficencyHist(h_trigger_efficency_1, h_trigger_reject_1);
  computeEfficencyHist(h_trigger_efficency_2, h_trigger_reject_2);
  computeEfficencyHist(h_trigger_efficency_3, h_trigger_reject_3);
  computeEfficencyHist(h_trigger_efficency_4, h_trigger_reject_4);
  computeEfficencyHist(h_trigger_efficency_1_topPt, h_trigger_reject_1_topPt);
  computeEfficencyHist(h_trigger_efficency_2_topPt, h_trigger_reject_2_topPt);
  computeEfficencyHist(h_trigger_efficency_3_topPt, h_trigger_reject_3_topPt);
  computeEfficencyHist(h_trigger_efficency_4_topPt, h_trigger_reject_4_topPt);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
B2GMonoTopTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void 
B2GMonoTopTreeMaker::printCutFlowResult(vid::CutFlowResult &cutflow){

  printf("    CutFlow name= %s    decision is %d\n", 
   cutflow.cutFlowName().c_str(),
   (int) cutflow.cutFlowPassed());
  int ncuts = cutflow.cutFlowSize();
  printf(" Index                               cut name              isMasked    value-cut-upon     pass?\n");
  for(int icut = 0; icut<ncuts; icut++){
    printf("  %2d      %50s    %d        %f          %d\n", icut,
     cutflow.getNameAtIndex(icut).c_str(),
     (int)cutflow.isCutMasked(icut),
     cutflow.getValueCutUpon(icut),
     (int)cutflow.getCutResultByIndex(icut));
  }
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(B2GMonoTopTreeMaker);

//  LocalWords:  NNPDF


