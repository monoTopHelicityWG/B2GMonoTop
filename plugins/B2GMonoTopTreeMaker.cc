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

// // CMS Top Tagger
// #include "DataFormats/BTauReco/interface/CATopJetTagInfo.h"
// #include "RecoJets/JetAlgorithms/interface/CATopJetHelper.h"

// Fastjet
// #include <fastjet/JetDefinition.hh>
// #include <fastjet/PseudoJet.hh>
// #include "fastjet/tools/Filter.hh"
// #include <fastjet/ClusterSequence.hh>
// #include <fastjet/ClusterSequenceArea.hh>
// #include "fastjet/contrib/SoftDrop.hh"
// #include "fastjet/contrib/EnergyCorrelator.hh"
// #include <fastjet/Selector.hh>
// #include "fastjet/tools/Pruner.hh"
// #include "Nsubjettiness.hh"
// #include "Njettiness.hh"

// root
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include <TRandom3.h>

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

//
// class declaration
//

class B2GMonoTopTreeMaker : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit B2GMonoTopTreeMaker(const edm::ParameterSet&);
      ~B2GMonoTopTreeMaker();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;


      void printCutFlowResult(vid::CutFlowResult &cutflow);

      // ----------member data ---------------------------
      edm::EDGetTokenT<pat::JetCollection>                     ak4jetToken_                     ;
      edm::EDGetTokenT<pat::JetCollection>                     ak8jetToken_                     ;
      edm::EDGetTokenT<pat::JetCollection>                     puppijetToken_                   ;
      edm::EDGetTokenT<pat::JetCollection>                     ak8CHSSoftDropSubjetsToken_      ;
      edm::EDGetTokenT<pat::JetCollection>                     ak8PuppiSoftDropSubjetsToken_    ;
      edm::EDGetTokenT<reco::GenJetCollection>                 ak4genjetToken_                  ;
      edm::EDGetTokenT<reco::GenJetCollection>                 ak8genjetToken_                  ;               
      edm::EDGetTokenT<edm::View<reco::GenParticle>>           prunedGenToken_                  ;               
      edm::EDGetTokenT<double>                                 rhoToken_                        ;
      edm::EDGetTokenT<std::vector<reco::Vertex>>              vtxToken_                        ;
      edm::EDGetTokenT<edm::TriggerResults>                    triggerResultsMETFilterToken_    ;
      edm::EDGetTokenT<edm::TriggerResults>                    triggerBits_                     ;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_                  ;
      edm::EDGetTokenT<pat::PackedTriggerPrescales>            triggerPrescales_                ;
      edm::EDGetTokenT<bool>                                   badMuonFilterToken_              ;
      edm::EDGetTokenT<bool>                                   badChargedCandidateFilterToken_  ;
      edm::EDGetTokenT<pat::MuonCollection>                    muonToken_                       ;
      edm::EDGetTokenT<edm::View<pat::Electron>>               electronToken_                   ;  
      edm::EDGetTokenT<pat::METCollection>                     metToken_                        ;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo>>         pileupInfoToken_                 ;
      edm::EDGetTokenT<LHEEventProduct>                        lheSrc_                          ;
      edm::EDGetTokenT<GenEventInfoProduct>                    pdfToken_                        ;
      edm::EDGetTokenT<reco::BeamSpot>                         beamSpotToken_                   ;                     
      edm::EDGetTokenT<reco::ConversionCollection>             conversionsToken_                ;

      edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult>>      eleIdFullInfoMapToken_HLTpre_    ;
      edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult>>      eleIdFullInfoMapToken_Loose_     ;
      edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult>>      eleIdFullInfoMapToken_Medium_    ;
      edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult>>      eleIdFullInfoMapToken_Tight_     ;
      edm::EDGetTokenT<edm::ValueMap<vid::CutFlowResult>>      eleIdFullInfoMapToken_HEEP_      ;
  
      // Grab from configuration file
    
      bool verbose_         ;
      bool verboseGen_      ;
      
      bool useToolbox_      ;

      bool runGenLoop_      ;
      bool runHadTree_   ;
      bool runLeptTree_ ;

      bool isZprime_        ;
      bool isMonoTop_         ;
      bool isRSG_           ;
      bool isRun2016F_      ;

      std::vector<std::string>  jecPayloadsAK4chs_;
      std::vector<std::string>  jecPayloadsAK8chs_;
      std::vector<std::string>  jecPayloadsAK4pup_;
      std::vector<std::string>  jecPayloadsAK8pup_;

      std::vector<std::string>  jecPayloadsAK4chsSecondary_;
      std::vector<std::string>  jecPayloadsAK8chsSecondary_;
      std::vector<std::string>  jecPayloadsAK4pupSecondary_;
      std::vector<std::string>  jecPayloadsAK8pupSecondary_;

      std::string jertextAK4_;   // jet resolution AK4 jets
      std::string jertextAK8_;   // jet resolution AK8 jets
      std::string jerSFtext_ ;   // jer SF

      // JEC
      boost::shared_ptr<FactorizedJetCorrector>   JetCorrectorAK4chs;
      boost::shared_ptr<FactorizedJetCorrector>   JetCorrectorAK8chs;
      boost::shared_ptr<FactorizedJetCorrector>   JetCorrectorAK4pup;
      boost::shared_ptr<FactorizedJetCorrector>   JetCorrectorAK8pup;
      boost::shared_ptr<JetCorrectionUncertainty> JetCorrUncertAK4chs;
      boost::shared_ptr<JetCorrectionUncertainty> JetCorrUncertAK8chs;
      boost::shared_ptr<JetCorrectionUncertainty> JetCorrUncertAK4pup;
      boost::shared_ptr<JetCorrectionUncertainty> JetCorrUncertAK8pup;
      
      // PU weight
      TFile* fPUweight;
      TH1D* hPUweight;
      TH1D* hPUweight_MBup;
      TH1D* hPUweight_MBdn;

      // A few basic histograms
      TH1D * h_cutflow_had ;
      TH1D * h_cutflow_lept;
      // TH1D * h_ak8puppi_softDropMass          ;
      // TH1D * h_ak8chs_softDropMass            ;
      // TH1D * h_ak8chs_softDropMass_reweighted ;
      // TH1D * h_ak8chs_pt                      ;
      // TH1D * h_ak8chs_pt_reweighted           ;
      TH1D * h_NtrueIntPU        ;
      TH1D * h_NPV               ;               
      TH1D * h_NPVreweighted     ;     
      TH1D * h_NPVgood           ;               
      TH1D * h_NPVgoodreweighted ;     

      // -- Triggers to be saved in tree
      std::vector<std::string> trigsToRun;


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
                                                                 
                                                                 
    
      TTree *TreeHad; 

      std::vector<int>  *HadTrigPrescales = new std::vector<int>  ;
      std::vector<bool> *HadTrigPass      = new std::vector<bool> ;
      std::string HadTrigAcceptBits;

      Float_t Jet0PtRaw                                 ;
      Float_t Jet0EtaRaw                                ;
      Float_t Jet0PhiRaw                                ;
      Float_t Jet0MassRaw                               ;

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

      TTree *TreeLept;
      // std::vector<std::string> *LeptTrigNames     = new std::vector<std::string>;
      std::vector<int> *LeptTrigPrescales = new std::vector<int>;
      std::vector<bool> *LeptTrigPass    = new std::vector<bool>;     
      std::string LeptTrigAcceptBits;

      Float_t JetPtRaw                               ;      
      Float_t JetEtaRaw                              ;
      Float_t JetPhiRaw                              ;
      Float_t JetMassRaw                             ;     


};



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
    ak8CHSSoftDropSubjetsToken_(consumes<pat::JetCollection>(  iConfig.getParameter<edm::InputTag>("ak8chsSubjetsInput"))),
    ak8PuppiSoftDropSubjetsToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("ak8puppiSubjetsInput"))),
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
  // h_ak8puppi_softDropMass            =  fs->make<TH1D>("h_ak8puppi_softDropMass"           ,"",200,0,400);
  // h_ak8chs_softDropMass              =  fs->make<TH1D>("h_ak8chs_softDropMass"             ,"",200,0,400);
  // h_ak8chs_softDropMass_reweighted   =  fs->make<TH1D>("h_ak8chs_softDropMass_reweighted"  ,"",200,0,400);
  // h_ak8chs_pt                        =  fs->make<TH1D>("h_ak8chs_pt"                       ,"",200,0,4000);
  // h_ak8chs_pt_reweighted             =  fs->make<TH1D>("h_ak8chs_pt_reweighted"            ,"",200,0,4000);
  h_NtrueIntPU                       =  fs->make<TH1D>("h_NtrueIntPU"                      ,"",200,0,200);
  h_NPV                              =  fs->make<TH1D>("h_NPV"                             ,"",200,0,200);
  h_NPVreweighted                    =  fs->make<TH1D>("h_NPVreweighted"                   ,"",200,0,200);
  h_NPVgood                          =  fs->make<TH1D>("h_NPVgood"                         ,"",200,0,200);
  h_NPVgoodreweighted                =  fs->make<TH1D>("h_NPVgoodreweighted"               ,"",200,0,200);


  trigsToRun.push_back("HLT_PFHT300_v");
  trigsToRun.push_back("HLT_PFHT350_v");
  trigsToRun.push_back("HLT_PFHT400_v");
  trigsToRun.push_back("HLT_PFHT475_v");
  trigsToRun.push_back("HLT_PFHT600_v");
  trigsToRun.push_back("HLT_PFHT650_v");
  trigsToRun.push_back("HLT_PFHT800_v");
  trigsToRun.push_back("HLT_PFHT900_v");
  trigsToRun.push_back("HLT_PFHT650_WideJetMJJ900"); //HLT_PFHT650_WideJetMJJ900DEtaJJ1p5_v6
  trigsToRun.push_back("HLT_PFHT650_WideJetMJJ950"); //HLT_PFHT650_WideJetMJJ950DEtaJJ1p5_v6

  // Single jet
  trigsToRun.push_back("HLT_CaloJet500_NoJetID_v");
  trigsToRun.push_back("HLT_PFJet320_v");
  trigsToRun.push_back("HLT_PFJet400_v");
  trigsToRun.push_back("HLT_PFJet450_v");
  trigsToRun.push_back("HLT_PFJet500_v");
  trigsToRun.push_back("HLT_AK8PFJet450_v");
  trigsToRun.push_back("HLT_AK8PFJet500_v");

  // Substructure
  trigsToRun.push_back("HLT_AK8PFJet360_TrimMass30_v");
  trigsToRun.push_back("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_v");
  trigsToRun.push_back("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_v");

  // Substructure + b-tag
  trigsToRun.push_back("HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV_p20_v");
  trigsToRun.push_back("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV_p20_v");

  // Muon
  trigsToRun.push_back("HLT_Mu45_eta2p1_v");
  trigsToRun.push_back("HLT_Mu50_v");
  trigsToRun.push_back("HLT_Mu55_v");
  trigsToRun.push_back("HLT_TkMu50_v");
  trigsToRun.push_back("HLT_IsoMu22_eta2p1_v");
  trigsToRun.push_back("HLT_IsoMu24_v");
  trigsToRun.push_back("HLT_IsoMu27_v");

  // Muon + jet
  trigsToRun.push_back("HLT_Mu30_eta2p1_PFJet150_PFJet50_v");
  trigsToRun.push_back("HLT_Mu40_eta2p1_PFJet200_PFJet50_v");

  // Electron
  trigsToRun.push_back("HLT_Ele32_eta2p1_WPTight_Gsf_v");
  trigsToRun.push_back("HLT_Ele35_WPLoose_Gsf_v");
  trigsToRun.push_back("HLT_Ele105_CaloIdVT_GsfTrkIdT_v");
  trigsToRun.push_back("HLT_Ele115_CaloIdVT_GsfTrkIdT_v");

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



  TreeHad = new TTree("TreeHad","TreeHad"); 


  TreeHad->Branch("HadTrigPrescales"   , "vector<int>", &HadTrigPrescales);
  TreeHad->Branch("HadTrigPass"        , "vector<bool>", &HadTrigPass);
  TreeHad->Branch("HadTrigAcceptBits"  , &HadTrigAcceptBits);
                         
  TreeHad->Branch("Jet0PtRaw"                             , & Jet0PtRaw                          ,    "Jet0PtRaw/F"                               );                                  
  TreeHad->Branch("Jet0EtaRaw"                            , & Jet0EtaRaw                         ,    "Jet0EtaRaw/F"                              );                                   
  TreeHad->Branch("Jet0PhiRaw"                            , & Jet0PhiRaw                         ,    "Jet0PhiRaw/F"                              );                                   
  TreeHad->Branch("Jet0MassRaw"                           , & Jet0MassRaw                        ,    "Jet0MassRaw/F"                             );    


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
 

  TreeLept->Branch("JetPtRaw"                             , & JetPtRaw                          ,    "JetPtRaw/F"                               );                                  
  TreeLept->Branch("JetEtaRaw"                            , & JetEtaRaw                         ,    "JetEtaRaw/F"                              );                                   
  TreeLept->Branch("JetPhiRaw"                            , & JetPhiRaw                         ,    "JetPhiRaw/F"                              );                                   
  TreeLept->Branch("JetMassRaw"                           , & JetMassRaw                        ,    "JetMassRaw/F"                             );   

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

  bool tophadronic=false;
  bool topleptonic=false;


  double hardest_parton_hardScatterOutgoing_pt        = 0;
  double second_hardest_parton_hardScatterOutgoing_pt = 0;

  int parton1id = 0;
  int parton2id = 0;
  int Wd1_id = 0 ;
  int Wd2_id = 0 ;

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

  // HadTrigNames       ->clear();
  // LeptTrigNames     ->clear();
  HadTrigPrescales   ->clear();
  LeptTrigPrescales ->clear();
  HadTrigPass        ->clear();
  LeptTrigPass      ->clear();

  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  if (verbose_) std::cout << "\n === TRIGGER PATHS === " << std::endl;
  int counttrigs =0;

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
      // HadTrigNames       ->push_back(name);
      // LeptTrigNames     ->push_back(name);
      HadTrigPrescales   ->push_back(prescale);
      LeptTrigPrescales ->push_back(prescale);
      HadTrigPass        ->push_back(pass);
      LeptTrigPass      ->push_back(pass);
      if (pass)  hltbit[counttrigs]=1;  
      counttrigs++;
  
  }// end loop over list of triggers to save in tree

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

  if (isZprime_ || isMonoTop_){
    edm::Handle<LHEEventProduct> EvtHandle;
    iEvent.getByToken(lheSrc_, EvtHandle);

    if (EvtHandle.isValid()){
      double centralWgt = EvtHandle->weights()[0].wgt;

      //Q^2 uncertainties
      if (verbose_) cout << "Calculating Q^2 uncertainties." << endl;
      double maxQ2wgt_frac = 1;
      double minQ2wgt_frac = 1;
      
      for (unsigned int iLHE = 0; iLHE < 9; ++iLHE){
        if (iLHE != 5 && iLHE != 7){
          double Q2wgt = EvtHandle->weights()[iLHE].wgt;
          double Q2wgt_frac = Q2wgt/(centralWgt);
          if (verbose_) cout << "Q^2 Weight: " << Q2wgt << "  Fractional Q^2 Weight: " << Q2wgt_frac << endl;
          maxQ2wgt_frac = max(maxQ2wgt_frac, Q2wgt_frac);
          minQ2wgt_frac = min(minQ2wgt_frac, Q2wgt_frac);
        }
      }
      
      Q2wgt_up = maxQ2wgt_frac;
      Q2wgt_down = minQ2wgt_frac;

      //NNPDF3 uncertainties
      if (verbose_) cout << "Calculating NNPDF3 uncertainties." << endl;
      double NNPDF3wgtAvg = 0.0;
      double NNPDF3wgtRMS = 0.0;
      double NNPDF3wgt = 0.0;
      double NNPDF3wgt_frac = 0.0;

      //MonoTop
      unsigned int PDFstart = 9;
      unsigned int PDFend = 109;

      //Zprime
      if (isZprime_){
        PDFstart = 10;
        PDFend = 110;
      }

      //Making sure central PDF isn't zero                                                                                              
      if (centralWgt == 0){
        NNPDF3wgt_up = 0.0;
        NNPDF3wgt_down = 0.0;
        if (verbose_) cout << "Unphysical: central PDF weight is zero!" << endl;
      }
      else{
        if (verbose_){cout << "PDF Fractional weights: "; 
        for (unsigned int i_lhePDF = PDFstart; i_lhePDF < PDFend; ++i_lhePDF){
          NNPDF3wgt = EvtHandle->weights()[i_lhePDF].wgt;
          NNPDF3wgt_frac = NNPDF3wgt/(centralWgt);
          NNPDF3wgtAvg += NNPDF3wgt_frac;
          
            // cout << "-----" << endl;
            // cout << i_lhePDF - PDFstart
            if (verbose_)cout<< NNPDF3wgt_frac<<", ";
            // cout << "-----" << endl;
            // cout << "" << endl;
          }
        }
        if (verbose_) cout<<endl;

        NNPDF3wgtAvg = NNPDF3wgtAvg/(PDFend - PDFstart);
            
        for (unsigned int i_lhePDF = PDFstart; i_lhePDF < PDFend; ++i_lhePDF){
          NNPDF3wgt = EvtHandle->weights()[i_lhePDF].wgt;
          NNPDF3wgt_frac = NNPDF3wgt/(centralWgt);
          NNPDF3wgtRMS += (NNPDF3wgt_frac - NNPDF3wgtAvg)*(NNPDF3wgt_frac - NNPDF3wgtAvg);
        }

        NNPDF3wgtRMS = sqrt(NNPDF3wgtRMS/(PDFend - PDFstart - 1));
        NNPDF3wgt_up = 1.0 + NNPDF3wgtRMS;
        NNPDF3wgt_down = 1.0 - NNPDF3wgtRMS;

        if (verbose_) cout <<"NNPDF3wgtAvg = "<< NNPDF3wgtAvg<<" NNPDF3wgtRMS = "<< NNPDF3wgtRMS<<endl;

      }
    }
  }

  else if (isRSG_){
    edm::Handle<GenEventInfoProduct> pdfstuff;
    iEvent.getByToken(pdfToken_, pdfstuff);

    LHAPDF::usePDFMember(1,0);

    float q = pdfstuff->pdf()->scalePDF;

    int id1 = pdfstuff->pdf()->id.first;
    double x1 = pdfstuff->pdf()->x.first;
    //double pdf1 = pdfstuff->pdf()->xPDF.first;                                                                                           

    int id2 = pdfstuff->pdf()->id.second;
    double x2 = pdfstuff->pdf()->x.second;
    //double pdf2 = pdfstuff->pdf()->xPDF.second;                                                                                          

    double xpdf1 = LHAPDF::xfx(1, x1, q, id1);
    double xpdf2 = LHAPDF::xfx(1, x2, q, id2);
    double w0 = xpdf1 * xpdf2;
    double sumsq = 0.0;
    for(int i=1; i <=100; ++i){
      LHAPDF::usePDFMember(1,i);
      double xpdf1_new = LHAPDF::xfx(1, x1, q, id1);
      double xpdf2_new = LHAPDF::xfx(1, x2, q, id2);
      double weight = xpdf1_new * xpdf2_new / w0;
      sumsq += ( weight - w0 ) * (weight - w0);
    }

    double rmsWt = sqrt( (1./99.)*sumsq );

    if ( rmsWt > 1.0){
      NNPDF3wgt_up = rmsWt;
      NNPDF3wgt_down = 2 - rmsWt;
    }

    if (rmsWt < 1.0){
      NNPDF3wgt_down = rmsWt;
      NNPDF3wgt_up = 2 - rmsWt;
    }
    
  }

  edm::Handle<GenEventInfoProduct> genEventInfo;
  iEvent.getByToken(pdfToken_, genEventInfo);
  double evWeight = 1.0 ;
  double qScale = 1.0 ;
  double pthat = 1.0 ;
  if (genEventInfo.isValid())
  {
    evWeight = genEventInfo->weight();
    qScale   = genEventInfo->qScale();
    // const std::vector<double>& binningValues = genEventInfo->binningValues(); // in case of Pythia6, this will be pypars/pari(17)
    pthat    = (genEventInfo->hasBinningValues() ? (genEventInfo->binningValues())[0] : 0.0);
    // std::vector<double>& evtWeights = genEventInfo->weights();
    // if (evWeight < 0 ) cout<<"NEGATIVE"<<endl;

    // I'll do this at the tree reading stage
    // if (evWeight < 0 ){
    //   if (verbose_) cout<<"evWeight < 0. evWeight *= -1.0"<<endl;;  
    //   evWeight *= -1.0
    // }
    if(verbose_) cout<<"GenEventInfo: qScale = "<<qScale<<" pthat = "<<pthat<<" evWeight = "<<evWeight<<" 1/pow(pthat/15,4.5) "<<1/pow(pthat/15,4.5)<<endl;
  }


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
  if (useToolbox_){
    iEvent.getByToken( ak8CHSSoftDropSubjetsToken_   , AK8CHSsub);
    iEvent.getByToken( puppijetToken_ , AK8PUPPI );
    iEvent.getByToken( ak8PuppiSoftDropSubjetsToken_ , AK8PUPPIsub);
  }

  TLorentzVector AK8_p4;
  TLorentzVector AK8_p4_2;

  double count_AK8CHS = 0;
  for (const pat::Jet &ijet : *AK8CHS) {
      // if (count_AK8CHS>1) break;
      if (count_AK8CHS==0 && ijet.pt()<30) break;
      if (verbose_) cout<<"\nJet "<<count_AK8CHS<<" with pT "<<ijet.pt()<<" sdMass "<<ijet.userFloat("ak8PFJetsCHSSoftDropMass")<<endl;
  
      reco::Candidate::LorentzVector corrJet = ijet.correctedP4(0);

      if(count_AK8CHS == 0) AK8_p4.SetPtEtaPhiM(corrJet.pt(), corrJet.eta(), corrJet.phi(), corrJet.mass());
      if(count_AK8CHS == 1) AK8_p4_2.SetPtEtaPhiM(corrJet.pt(), corrJet.eta(), corrJet.phi(), corrJet.mass());

      if(verboseGen_  && count_AK8CHS==0) {

        cout << "Gen top pt " << t_p4.Pt() << " " << t_p4.Eta() << " " << t_p4.Phi() << " " << t_p4.M() << endl;
        cout << "leading Ak8 top pt " << AK8_p4.Pt() << " " << AK8_p4.Eta() << " " << AK8_p4.Phi() << " " << AK8_p4.M() << endl;
        
        cout << "deltaR gen, leading ak8 " << AK8_p4.DeltaR(t_p4) << endl;

      }

      if(count_AK8CHS == 1) cout << "deltaR gen, leading ak8 " << AK8_p4.DeltaR(t_p4) << " " << (AK8_p4+AK8_p4_2).M() <<   endl;
      count_AK8CHS++;

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

  int count_AK4CHS = 0;
  if (verbose_) cout<<"AK4 jet loop"<<endl;
  for (const pat::Jet &ijet : *AK4MINI) { 
    if (ijet.pt()<15 || fabs(ijet.eta())>3.0) continue; 

    reco::Candidate::LorentzVector corrJet = ijet.correctedP4(0);
    if(count_AK8CHS < 5){
      AK4_p4[count_AK4CHS].SetPtEtaPhiM(corrJet.pt(), corrJet.eta(), corrJet.phi(), corrJet.mass() );

      cout << "deltaR Ak4, leading ak8 " << AK8_p4.DeltaR(AK4_p4[count_AK4CHS]) << endl;

      cout << "deltaR Gen t, leading ak8 " << t_p4.DeltaR(AK4_p4[count_AK4CHS]) << endl;
      cout << "deltaR Gen b, leading ak8 " << b_p4.DeltaR(AK4_p4[count_AK4CHS]) << endl;
      cout << "deltaR Gen W, leading ak8 " << W_p4.DeltaR(AK4_p4[count_AK4CHS]) << endl;

      cout << "mass, ak4+ak8 " <<  AK4_p4[count_AK4CHS].M() << " " <<   (AK4_p4[count_AK4CHS]+ AK8_p4).M() << endl;

    } 

    count_AK4CHS++;
  }






  if(1==2){
    std::cout << NNPDF3wgt_down<< NNPDF3wgt_up<<Q2wgt_up<<Q2wgt_down<<std::endl;
  }

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

// ------------ method called once each job just after ending the event loop  ------------
void 
B2GMonoTopTreeMaker::endJob() 
{
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


