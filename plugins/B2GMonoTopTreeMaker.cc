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
      edm::EDGetTokenT<pat::JetCollection>                     ca12puppijetToken_                   ;
      edm::EDGetTokenT<pat::JetCollection>                     ak8CHSSoftDropSubjetsToken_      ;
      edm::EDGetTokenT<pat::JetCollection>                     ak8PuppiSoftDropSubjetsToken_    ;
      edm::EDGetTokenT<pat::JetCollection>                     ca12PuppiSoftDropSubjetsToken_    ;
      //edm::EDGetTokenT<pat::JetCollection>                     ca8puppijetToken_                 ;
      //edm::EDGetTokenT<pat::JetCollection>                     ca8PuppiSoftDropSubjetsToken_    ;

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


      double uncorrected_err = 0;
      double uncorrected_se= 0;
      double corrected_err = 0;
      double corrected_se= 0;

      Float_t Gen_array_t_p4[4];
      Float_t Gen_array_final_t_p4[4];
      Float_t Gen_array_b_p4[4];
      Float_t Gen_array_W_p4[4];
    
      Float_t Gen_array_Wd1_p4[4];
      Float_t Gen_array_Wd2_p4[4];
    
      Float_t Gen_array_hardest_parton_hardScatterOutgoing_p4[4];
      Float_t Gen_array_second_hardest_parton_hardScatterOutgoing_p4[4];
    
      bool tophadronic=false;
      bool topleptonic=false;
    
      int parton1id = 0;
      int parton2id = 0;
      int Wd1_id = 0 ;
      int Wd2_id = 0 ;

      Float_t CA12JetPtRaw                               ;      
      Float_t CA12JetEtaRaw                              ;
      Float_t CA12JetPhiRaw                              ;
      Float_t CA12JetMassRaw                             ;

      Float_t CA12JetTau1                           ;
      Float_t CA12JetTau2                           ;
      Float_t CA12JetTau3                           ;
      Float_t CA12JetTau4                           ;
      Float_t CA12JetTau32                          ;
      Float_t CA12JetTau21                          ;

      Float_t CA12Jetsubjet0bdisc                 ;
      Float_t CA12Jetsubjet1bdisc                 ;
      Float_t CA12Jetmaxbdisc                     ;

      Float_t CA12Jetsubjet0pt                    ;
      Float_t CA12Jetsubjet0mass                  ;
      Float_t CA12Jetsubjet0eta                   ;
      Float_t CA12Jetsubjet0phi                   ;
      Float_t CA12Jetsubjet0area                  ;

      Float_t CA12Jetsubjet1pt                    ;
      Float_t CA12Jetsubjet1mass                  ;
      Float_t CA12Jetsubjet1eta                   ;
      Float_t CA12Jetsubjet1phi                   ;
      Float_t CA12Jetsubjet1area                  ;



      Float_t AK4ReconstructedJetPt                             ;      
      Float_t AK4ReconstructedJetEta                             ;
      Float_t AK4ReconstructedJetPhi                             ;
      Float_t AK4ReconstructedJetMass                             ;

      Float_t AK4bJetPtRaw                               ;      
      Float_t AK4bJetEtaRaw                              ;
      Float_t AK4bJetPhiRaw                              ;
      Float_t AK4bJetMassRaw                             ;

      Float_t AK4bJet_PtSmear                           ;              
      Float_t AK4bJet_PtSmearUp                         ;              
      Float_t AK4bJet_PtSmearDn                         ;              
      Float_t AK4bJet_PtUncorr                          ; 
      Float_t AK4bJet_Corr                              ; 
      Float_t AK4bJet_CorrUp                            ; 
      Float_t AK4bJet_CorrDn                            ;
      Float_t AK4bJet_bDisc                            ;

      Float_t AK4WJetPtRaw                               ;      
      Float_t AK4WJetEtaRaw                              ;
      Float_t AK4WJetPhiRaw                              ;
      Float_t AK4WJetMassRaw                             ;

      Float_t AK4WJet_PtSmear                           ;              
      Float_t AK4WJet_PtSmearUp                         ;              
      Float_t AK4WJet_PtSmearDn                         ;              
      Float_t AK4WJet_PtUncorr                          ; 
      Float_t AK4WJet_Corr                              ; 
      Float_t AK4WJet_CorrUp                            ; 
      Float_t AK4WJet_CorrDn                            ;
      Float_t AK4WJet_bDisc                            ;


      Float_t AK4W2JetPtRaw                               ;      
      Float_t AK4W2JetEtaRaw                              ;
      Float_t AK4W2JetPhiRaw                              ;
      Float_t AK4W2JetMassRaw                             ;

      Float_t AK4W2Jet_PtSmear                           ;              
      Float_t AK4W2Jet_PtSmearUp                         ;              
      Float_t AK4W2Jet_PtSmearDn                         ;              
      Float_t AK4W2Jet_PtUncorr                          ; 
      Float_t AK4W2Jet_Corr                              ; 
      Float_t AK4W2Jet_CorrUp                            ; 
      Float_t AK4W2Jet_CorrDn                            ;
      Float_t AK4W2Jet_bDisc                            ;




      Float_t JetPtRaw                               ;      
      Float_t JetEtaRaw                              ;
      Float_t JetPhiRaw                              ;
      Float_t JetMassRaw                             ;
      // Float_t JetP                                   ;
      // Float_t JetPt                                  ;
      // Float_t JetEta                                 ;
      // Float_t JetPhi                                 ;
      // Float_t JetRap                                 ;
      // Float_t JetEnergy                              ;
      // Float_t JetMass                                ;
      Float_t JetArea                                ;
      // Float_t JetSDmass                              ;
      Float_t JetSDmassRaw                           ;
      Float_t JetSDmassSubjetCorrL23                       ;
      // Float_t JetSDmassCorrL23Up                     ;
      // Float_t JetSDmassCorrL23Dn                     ;
      Float_t JetSDmassSubjetCorrL123                      ;
      // Float_t JetSDmassCorrL123Up                    ;
      // Float_t JetSDmassCorrL123Dn                    ;
      // Float_t JetSDmassCorrL23Smear                  ;
      // Float_t JetSDmassCorrL23SmearUp                ;
      // Float_t JetSDmassCorrL23SmearDn                ;
      Float_t JetSDptRaw                             ;
      // Float_t JetSDptCorrL23                         ;
      // Float_t JetSDptCorrL23Up                       ;
      // Float_t JetSDptCorrL23Dn                       ;
      // Float_t JetSDptCorrL123                        ;
      // Float_t JetSDptCorrL123Up                      ;
      // Float_t JetSDptCorrL123Dn                      ;
      // Float_t JetSDptCorrL23Smear                    ;
      // Float_t JetSDptCorrL23SmearUp                  ;
      // Float_t JetSDptCorrL23SmearDn                  ;
      Float_t JetSDetaRaw                            ;
      Float_t JetSDphiRaw                            ;
      Float_t JetMassPruned                          ;
      Float_t JetMassTrimmed                         ;
      Float_t JetTau1                                ;
      Float_t JetTau2                                ;
      Float_t JetTau3                                ;
      Float_t JetTau4                                ;
      Float_t JetTau32                               ;
      Float_t JetTau21                               ;
      Float_t JetSDsubjet0bdisc                      ;
      Float_t JetSDsubjet1bdisc                      ;
      Float_t JetSDmaxbdisc                          ;
      Float_t JetSDmaxbdiscflavHadron                ;
      Float_t JetSDmaxbdiscflavParton                ;
      Float_t JetSDsubjet0pt                         ;
      Float_t JetSDsubjet0mass                       ;
      Float_t JetSDsubjet0eta                        ;
      Float_t JetSDsubjet0phi                        ;
      Float_t JetSDsubjet0area                       ;
      Float_t JetSDsubjet0flavHadron                 ;
      Float_t JetSDsubjet0flavParton                 ;
      Float_t JetSDsubjet0matchedgenjetpt                 ;
      Float_t JetSDsubjet0tau1                       ;
      Float_t JetSDsubjet0tau2                       ;
      Float_t JetSDsubjet0tau3                       ;
      Float_t JetSDsubjet1pt                         ;
      Float_t JetSDsubjet1mass                       ;
      Float_t JetSDsubjet1eta                        ;
      Float_t JetSDsubjet1phi                        ;
      Float_t JetSDsubjet1area                       ;
      Float_t JetSDsubjet1flavHadron                 ;
      Float_t JetSDsubjet1flavParton                 ;
      Float_t JetSDsubjet1matchedgenjetpt                 ;
      Float_t JetSDsubjet1tau1                       ;
      Float_t JetSDsubjet1tau2                       ;
      Float_t JetSDsubjet1tau3                       ;
      // Float_t JetPuppiP                              ;
      Float_t JetPuppiPtRaw                             ;
      Float_t JetPuppiEtaRaw                            ;
      Float_t JetPuppiPhiRaw                            ;
      Float_t JetPuppiMassRaw                           ;
      Float_t JetPuppiArea                           ;
      Float_t JetPuppiSDmassRaw                         ;
      Float_t JetPuppiSDmassSubjetCorr                     ;
      // Float_t JetPuppiSDmassSubjetCorrUp                   ;
      // Float_t JetPuppiSDmassSubjetCorrDn                   ;
      // Float_t JetPuppiSDmassSubjetCorrL23Smear             ;
      // Float_t JetPuppiSDmassSubjetCorrL23SmearUp           ;
      // Float_t JetPuppiSDmassSubjetCorrL23SmearDn           ;
      Float_t JetPuppiSDptRaw                           ;
      // Float_t JetPuppiSDptSubjetCorr                       ;
      // Float_t JetPuppiSDptSubjetCorrUp                     ;
      // Float_t JetPuppiSDptSubjetCorrDn                     ;
      // Float_t JetPuppiSDptSubjetCorrL23Smear               ;
      // Float_t JetPuppiSDptSubjetCorrL23SmearUp             ;
      // Float_t JetPuppiSDptSubjetCorrL23SmearDn             ;
      Float_t JetPuppiSDetaRaw                          ;
      Float_t JetPuppiSDphiRaw                          ;
      Float_t JetPuppiTau1                           ;
      Float_t JetPuppiTau2                           ;
      Float_t JetPuppiTau3                           ;
      Float_t JetPuppiTau4                           ;
      Float_t JetPuppiTau32                          ;
      Float_t JetPuppiTau21                          ;


      Float_t JetPuppiSDsubjet0bdisc                 ;
      Float_t JetPuppiSDsubjet1bdisc                 ;
      Float_t JetPuppiSDmaxbdisc                     ;
      Float_t JetPuppiSDmaxbdiscflavHadron           ;
      Float_t JetPuppiSDmaxbdiscflavParton           ;
      Float_t JetPuppiSDsubjet0pt                    ;
      Float_t JetPuppiSDsubjet0mass                  ;
      Float_t JetPuppiSDsubjet0eta                   ;
      Float_t JetPuppiSDsubjet0phi                   ;
      Float_t JetPuppiSDsubjet0area                  ;
      Float_t JetPuppiSDsubjet0flavHadron            ;
      Float_t JetPuppiSDsubjet0flavParton            ;
      Float_t JetPuppiSDsubjet0matchedgenjetpt            ;
      Float_t JetPuppiSDsubjet0tau1                  ;
      Float_t JetPuppiSDsubjet0tau2                  ;
      Float_t JetPuppiSDsubjet0tau3                  ;
      Float_t JetPuppiSDsubjet1pt                    ;
      Float_t JetPuppiSDsubjet1mass                  ;
      Float_t JetPuppiSDsubjet1eta                   ;
      Float_t JetPuppiSDsubjet1phi                   ;
      Float_t JetPuppiSDsubjet1area                  ;
      Float_t JetPuppiSDsubjet1flavHadron            ;
      Float_t JetPuppiSDsubjet1flavParton            ;
      Float_t JetPuppiSDsubjet1matchedgenjetpt            ;
      Float_t JetPuppiSDsubjet1tau1                  ;
      Float_t JetPuppiSDsubjet1tau2                  ;
      Float_t JetPuppiSDsubjet1tau3                  ;
      Float_t JetPuppiSDECF1                         ;
      Float_t JetPuppiSDECF2                         ;
      Float_t JetPuppiSDECF3                         ;
      Float_t JetPuppiSDECF4                         ;
      Float_t JetPuppiSDECF5                         ;
      Float_t JetPuppiSDC_2                          ;
      Float_t JetPuppiSDD_2                          ;
      Float_t JetPuppiSDC_3                          ;
      Float_t JetPuppiSDD_3                          ;


      // Float_t JetPuppiSDmassUserFloat   ;  
      Float_t JetPuppiMassPruned        ;  
      Float_t JetPuppiMassTrimmed       ;  


      Int_t JetNsubjetsSD ;
      Int_t JetNsubjetsSDPuppi ;

      Float_t JetCHF                                 ;
      Float_t JetNHF                                 ;
      Float_t JetCM                                  ;
      Float_t JetNM                                  ;
      Float_t JetNEF                                 ;
      Float_t JetCEF                                 ;
      Float_t JetMF                                  ;
      Float_t JetMult                                ;
      Float_t JetPuppiCHF                            ;
      Float_t JetPuppiNHF                            ;
      Float_t JetPuppiCM                             ;
      Float_t JetPuppiNM                             ;
      Float_t JetPuppiNEF                            ;
      Float_t JetPuppiCEF                            ;
      Float_t JetPuppiMF                             ;
      Float_t JetPuppiMult                           ;
      Float_t JetMassCorrFactor                      ;
      Float_t JetMassCorrFactorUp                    ;
      Float_t JetMassCorrFactorDn                    ;
      Float_t JetCorrFactor                          ;
      Float_t JetCorrFactorUp                        ;
      Float_t JetCorrFactorDn                        ;
      Float_t JetPtSmearFactor                       ;
      Float_t JetPtSmearFactorUp                     ;
      Float_t JetPtSmearFactorDn                     ;
      Float_t JetPuppiMassCorrFactor                 ;
      Float_t JetPuppiMassCorrFactorUp               ;
      Float_t JetPuppiMassCorrFactorDn               ;
      Float_t JetPuppiCorrFactor                     ;
      Float_t JetPuppiCorrFactorUp                   ;
      Float_t JetPuppiCorrFactorDn                   ;
      Float_t JetPuppiPtSmearFactor                  ;
      Float_t JetPuppiPtSmearFactorUp                ;
      Float_t JetPuppiPtSmearFactorDn                ;
      // Float_t JetEtaScaleFactor                      ;
      // Float_t JetPhiScaleFactor                      ;
      // // Float_t JetMatchedGenJetDR                     ;
      //Float_t JetMatchedGenJetPt                     ;
      //Float_t JetMatchedGenJetMass                   ;
      //Float_t JetPuppiMatchedGenJetPt                     ;
      //Float_t JetPuppiMatchedGenJetMass                   ;
      Float_t JetMatchedGenJetPt                     ;
      Float_t JetMatchedGenJetMass                   ;
      Float_t JetPuppiMatchedGenJetPt                     ;
      Float_t JetPuppiMatchedGenJetMass                   ;

      Int_t   JetGenMatched_TopHadronic              ;
      Float_t JetGenMatched_TopPt                    ;
      Float_t JetGenMatched_TopEta                   ;
      Float_t JetGenMatched_TopPhi                   ;
      Float_t JetGenMatched_TopMass                  ;
      Float_t JetGenMatched_bPt                      ;
      Float_t JetGenMatched_WPt                      ;
      Float_t JetGenMatched_Wd1Pt                    ;
      Float_t JetGenMatched_Wd2Pt                    ;
      Float_t JetGenMatched_Wd1ID                    ;
      Float_t JetGenMatched_Wd2ID                    ;
      Float_t JetGenMatched_MaxDeltaRPartonTop       ;
      Float_t JetGenMatched_MaxDeltaRWPartonTop      ;
      Float_t JetGenMatched_MaxDeltaRWPartonW        ;
      Float_t JetGenMatched_DeltaR_t_b               ;
      Float_t JetGenMatched_DeltaR_t_W               ;
      Float_t JetGenMatched_DeltaR_t_Wd1             ;
      Float_t JetGenMatched_DeltaR_t_Wd2             ;
      Float_t JetGenMatched_DeltaR_W_b1              ;
      Float_t JetGenMatched_DeltaR_W_Wd1             ;
      Float_t JetGenMatched_DeltaR_W_Wd2             ;
      Float_t JetGenMatched_DeltaR_Wd1_Wd2           ;
      Float_t JetGenMatched_DeltaR_Wd1_b             ;
      Float_t JetGenMatched_DeltaR_Wd2_b             ;
      Float_t JetGenMatched_DeltaR_jet_t             ;
      Float_t JetGenMatched_DeltaR_jet_W             ;
      Float_t JetGenMatched_DeltaR_jet_b             ;
      Float_t JetGenMatched_DeltaR_jet_Wd1           ;
      Float_t JetGenMatched_DeltaR_jet_Wd2           ;
      Float_t JetGenMatched_DeltaR_pup0_b            ;
      Float_t JetGenMatched_DeltaR_pup0_Wd1          ;
      Float_t JetGenMatched_DeltaR_pup0_Wd2          ;
      Float_t JetGenMatched_DeltaR_pup1_b            ;
      Float_t JetGenMatched_DeltaR_pup1_Wd1          ;
      Float_t JetGenMatched_DeltaR_pup1_Wd2          ;
      Float_t JetGenMatched_partonPt                 ;
      Float_t JetGenMatched_partonEta                ;
      Float_t JetGenMatched_partonPhi                ;
      Float_t JetGenMatched_partonMass               ;
      Float_t JetGenMatched_partonID                 ;
      Float_t JetGenMatched_DeltaRjetParton          ;

      Float_t HadMETpx                          ;
      Float_t HadMETpy                          ;
      Float_t HadMETpt                          ;
      Float_t HadMETphi                         ;
      Float_t HadMETsumET                       ;
      Float_t HadMETgenMET                      ;
      Float_t HadMETuncorPt                     ;

      Float_t HadMETshiftedPtJetEnUp            ;
      Float_t HadMETshiftedPtJetEnDn            ;
      Float_t HadMETshiftedPtElEnUp             ;
      Float_t HadMETshiftedPtElEnDn             ;
      Float_t HadMETshiftedPtMuEnUp             ;
      Float_t HadMETshiftedPtMuEnDn             ;
      Float_t HadMETshiftedPtJetResUp           ;
      Float_t HadMETshiftedPtJetResDn           ;
      Float_t HadMETshiftedPtUnclEnUp           ;
      Float_t HadMETshiftedPtUnclEnDn           ;


      Float_t HadNvtx                           ;
      Float_t HadNvtxGood                       ;
      Float_t HadNPUtrue                        ;
      Float_t HadRho                            ;
      Float_t HadEventWeight                    ;
      Float_t HadPUweight                       ;
      Float_t HadPUweight_MBup                  ;
      Float_t HadPUweight_MBdn                  ;


      Float_t HadGenTTmass                      ;
      Int_t   HadGenCountHadTop                 ;

      Float_t HTlep                                  ;
      Float_t ST                                     ;                
      Float_t ST_CorrDn                              ;                
      Float_t ST_CorrUp                              ;                
      Float_t ST_PtSmearNom                          ;                
      Float_t ST_PtSmearUp                           ;                
      Float_t ST_PtSmearDn                           ;   

      Float_t HadQ2weight_CorrDn                ;
      Float_t HadQ2weight_CorrUp                ;
      Float_t HadNNPDF3weight_CorrDn            ;
      Float_t HadNNPDF3weight_CorrUp            ;
      Int_t   HadRunNum                         ;
      Int_t   HadLumiBlock                      ;
      Int_t   HadEventNum                       ;
      Int_t   HadPassMETFilters                 ;  



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

      //Float_t JetPtRaw                               ;      
      //Float_t JetEtaRaw                              ;
      //Float_t JetPhiRaw                              ;
      //Float_t JetMassRaw                             ;     


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
    ca12puppijetToken_(consumes<pat::JetCollection>(  iConfig.getParameter<edm::InputTag>("ca12puppiInput"))),
    ak8CHSSoftDropSubjetsToken_(consumes<pat::JetCollection>(  iConfig.getParameter<edm::InputTag>("ak8chsSubjetsInput"))),
    ak8PuppiSoftDropSubjetsToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("ak8puppiSubjetsInput"))),
    ca12PuppiSoftDropSubjetsToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("ca12puppiSubjetsInput"))),
    //ca8puppijetToken_(consumes<pat::JetCollection>(  iConfig.getParameter<edm::InputTag>("ca8puppiInput"))),
    //ca8PuppiSoftDropSubjetsToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("ca8puppiSubjetsInput"))),
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



  TreeHad->Branch("JetPtRaw"                             , & JetPtRaw                          ,    "JetPtRaw/F"                               );                                  
  TreeHad->Branch("JetEtaRaw"                            , & JetEtaRaw                         ,    "JetEtaRaw/F"                              );                                   
  TreeHad->Branch("JetPhiRaw"                            , & JetPhiRaw                         ,    "JetPhiRaw/F"                              );                                   
  TreeHad->Branch("JetMassRaw"                           , & JetMassRaw                        ,    "JetMassRaw/F"                             );                                                       
  // TreeHad->Branch("JetP"                                 , & JetP                              ,    "JetP/F"                                   );                              
  // TreeHad->Branch("JetPt"                                , & JetPt                             ,    "JetPt/F"                                  );                               
  // TreeHad->Branch("JetEta"                               , & JetEta                            ,    "JetEta/F"                                 );                                
  // TreeHad->Branch("JetPhi"                               , & JetPhi                            ,    "JetPhi/F"                                 );                                
  // TreeHad->Branch("JetRap"                               , & JetRap                            ,    "JetRap/F"                                 );                                
  // TreeHad->Branch("JetEnergy"                            , & JetEnergy                         ,    "JetEnergy/F"                              );                                   
  // TreeHad->Branch("JetMass"                              , & JetMass                           ,    "JetMass/F"                                );                                 
  TreeHad->Branch("JetArea"                              , & JetArea                           ,    "JetArea/F"                                );                                 
  
  // TreeHad->Branch("JetSDmass"                            , & JetSDmass                         ,    "JetSDmass/F"                              );                                         
  TreeHad->Branch("JetSDmassRaw"                         , & JetSDmassRaw                      ,    "JetSDmassRaw/F"                           );                                               
  TreeHad->Branch("JetSDmassSubjetCorrL23"                     , & JetSDmassSubjetCorrL23                  ,    "JetSDmassSubjetCorrL23/F"                       );                                                    
  // TreeHad->Branch("JetSDmassSubjetCorrL23Up"                   , & JetSDmassSubjetCorrL23Up                ,    "JetSDmassSubjetCorrL23Up/F"                     );                                                      
  // TreeHad->Branch("JetSDmassSubjetCorrL23Dn"                   , & JetSDmassSubjetCorrL23Dn                ,    "JetSDmassSubjetCorrL23Dn/F"                     );                                                      
  TreeHad->Branch("JetSDmassSubjetCorrL123"                    , & JetSDmassSubjetCorrL123                 ,    "JetSDmassSubjetCorrL123/F"                      );                                                      
  // TreeHad->Branch("JetSDmassCorrL123Up"                  , & JetSDmassCorrL123Up               ,    "JetSDmassCorrL123Up/F"                    );                                                        
  // TreeHad->Branch("JetSDmassCorrL123Dn"                  , & JetSDmassCorrL123Dn               ,    "JetSDmassCorrL123Dn/F"                    );                                                        
  // TreeHad->Branch("JetSDmassCorrL23Smear"                , & JetSDmassCorrL23Smear             ,    "JetSDmassCorrL23Smear/F"                     );                                                     
  // TreeHad->Branch("JetSDmassCorrL23SmearUp"              , & JetSDmassCorrL23SmearUp           ,    "JetSDmassCorrL23SmearUp/F"                   );                                                       
  // TreeHad->Branch("JetSDmassCorrL23SmearDn"              , & JetSDmassCorrL23SmearDn           ,    "JetSDmassCorrL23SmearDn/F"                   );   
  TreeHad->Branch("JetSDptRaw"                           , & JetSDptRaw                        ,    "JetSDptRaw/F"                             );                                               
  // TreeHad->Branch("JetSDptCorrL23"                       , & JetSDptCorrL23                    ,    "JetSDptCorrL23/F"                         );                                                    
  // TreeHad->Branch("JetSDptCorrL23Up"                     , & JetSDptCorrL23Up                  ,    "JetSDptCorrL23Up/F"                       );                                                      
  // TreeHad->Branch("JetSDptCorrL23Dn"                     , & JetSDptCorrL23Dn                  ,    "JetSDptCorrL23Dn/F"                       );                                                      
  // TreeHad->Branch("JetSDptCorrL123"                      , & JetSDptCorrL123                   ,    "JetSDptCorrL123/F"                        );                                                      
  // TreeHad->Branch("JetSDptCorrL123Up"                    , & JetSDptCorrL123Up                 ,    "JetSDptCorrL123Up/F"                      );                                                        
  // TreeHad->Branch("JetSDptCorrL123Dn"                    , & JetSDptCorrL123Dn                 ,    "JetSDptCorrL123Dn/F"                      );                                                        
  // TreeHad->Branch("JetSDptCorrL23Smear"                  , & JetSDptCorrL23Smear               ,    "JetSDptCorrL23Smear/F"                       );                                                     
  // TreeHad->Branch("JetSDptCorrL23SmearUp"                , & JetSDptCorrL23SmearUp             ,    "JetSDptCorrL23SmearUp/F"                     );                                                       
  // TreeHad->Branch("JetSDptCorrL23SmearDn"                , & JetSDptCorrL23SmearDn             ,    "JetSDptCorrL23SmearDn/F"                     );                                                     
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

  // TreeHad->Branch("JetPuppiP"                            , & JetPuppiP                         ,    "JetPuppiP/F"                              );                                    
  TreeHad->Branch("JetPuppiPtRaw"                           , & JetPuppiPtRaw                        ,    "JetPuppiPtRaw/F"                             );                                    
  TreeHad->Branch("JetPuppiEtaRaw"                          , & JetPuppiEtaRaw                       ,    "JetPuppiEtaRaw/F"                            );                                     
  TreeHad->Branch("JetPuppiPhiRaw"                          , & JetPuppiPhiRaw                       ,    "JetPuppiPhiRaw/F"                            );                                     
  TreeHad->Branch("JetPuppiMassRaw"                         , & JetPuppiMassRaw                      ,    "JetPuppiMassRaw/F"                           );                                      
  TreeHad->Branch("JetPuppiArea"                         , & JetPuppiArea                      ,    "JetPuppiArea/F"                           );                                      

  // TreeHad->Branch("JetPuppiSDmassUserFloat"               , & JetPuppiSDmassUserFloat          ,   "JetPuppiSDmassUserFloat/F"                     );
  TreeHad->Branch("JetPuppiMassPruned"                    , & JetPuppiMassPruned               ,   "JetPuppiMassPruned/F"                          );
  TreeHad->Branch("JetPuppiMassTrimmed"                   , & JetPuppiMassTrimmed              ,   "JetPuppiMassTrimmed/F"                         );


  
  TreeHad->Branch("JetPuppiSDmassRaw"                         , & JetPuppiSDmassRaw                    ,    "JetPuppiSDmassRaw/F"                          );
  TreeHad->Branch("JetPuppiSDmassSubjetCorr"                     , & JetPuppiSDmassSubjetCorr                ,    "JetPuppiSDmassSubjetCorr/F"                      );
  // TreeHad->Branch("JetPuppiSDmassSubjetCorrUp"                   , & JetPuppiSDmassSubjetCorrUp              ,    "JetPuppiSDmassSubjetCorrUp/F"                    );
  // TreeHad->Branch("JetPuppiSDmassSubjetCorrDn"                   , & JetPuppiSDmassSubjetCorrDn              ,    "JetPuppiSDmassSubjetCorrDn/F"                    );
  // TreeHad->Branch("JetPuppiSDmassSubjetCorrL23Smear"             , & JetPuppiSDmassSubjetCorrL23Smear        ,    "JetPuppiSDmassSubjetCorrL23Smear/F"              );
  // TreeHad->Branch("JetPuppiSDmassSubjetCorrL23SmearUp"           , & JetPuppiSDmassSubjetCorrL23SmearUp      ,    "JetPuppiSDmassSubjetCorrL23SmearUp/F"            );
  // TreeHad->Branch("JetPuppiSDmassSubjetCorrL23SmearDn"           , & JetPuppiSDmassSubjetCorrL23SmearDn      ,    "JetPuppiSDmassSubjetCorrL23SmearDn/F"            );
  TreeHad->Branch("JetPuppiSDptRaw"                           , & JetPuppiSDptRaw                      ,    "JetPuppiSDptRaw/F"                            );
  // TreeHad->Branch("JetPuppiSDptSubjetCorr"                       , & JetPuppiSDptSubjetCorr                  ,    "JetPuppiSDptSubjetCorr/F"                        );
  // TreeHad->Branch("JetPuppiSDptSubjetCorrUp"                     , & JetPuppiSDptSubjetCorrUp                ,    "JetPuppiSDptSubjetCorrUp/F"                      );
  // TreeHad->Branch("JetPuppiSDptSubjetCorrDn"                     , & JetPuppiSDptSubjetCorrDn                ,    "JetPuppiSDptSubjetCorrDn/F"                      );
  // TreeHad->Branch("JetPuppiSDptSubjetCorrL23Smear"               , & JetPuppiSDptSubjetCorrL23Smear          ,    "JetPuppiSDptSubjetCorrL23Smear/F"                );
  // TreeHad->Branch("JetPuppiSDptSubjetCorrL23SmearUp"             , & JetPuppiSDptSubjetCorrL23SmearUp        ,    "JetPuppiSDptSubjetCorrL23SmearUp/F"              );
  // TreeHad->Branch("JetPuppiSDptSubjetCorrL23SmearDn"             , & JetPuppiSDptSubjetCorrL23SmearDn        ,    "JetPuppiSDptSubjetCorrL23SmearDn/F"              );
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
  // TreeHad->Branch("JetEtaScaleFactor"                    , & JetEtaScaleFactor                 ,    "JetEtaScaleFactor/F"                      );                                           
  // TreeHad->Branch("JetPhiScaleFactor"                    , & JetPhiScaleFactor                 ,    "JetPhiScaleFactor/F"                      );                                           
  // TreeHad->Branch("JetMatchedGenJetDR"                   , & JetMatchedGenJetDR                ,    "JetMatchedGenJetDR/F"                     );  
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
 

  //TreeLept->Branch("JetPtRaw"                             , & JetPtRaw                          ,    "JetPtRaw/F"                               );                                  
  //TreeLept->Branch("JetEtaRaw"                            , & JetEtaRaw                         ,    "JetEtaRaw/F"                              );                                   
  //TreeLept->Branch("JetPhiRaw"                            , & JetPhiRaw                         ,    "JetPhiRaw/F"                              );                                   
  //TreeLept->Branch("JetMassRaw"                           , & JetMassRaw                        ,    "JetMassRaw/F"                             );   

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
  //edm::Handle<pat::JetCollection> CA8PUPPI;
  //edm::Handle<pat::JetCollection> CA8PUPPIsub;
  if (useToolbox_){
    iEvent.getByToken( ak8CHSSoftDropSubjetsToken_   , AK8CHSsub);
    iEvent.getByToken( puppijetToken_ , AK8PUPPI );
    iEvent.getByToken( ak8PuppiSoftDropSubjetsToken_ , AK8PUPPIsub);

    iEvent.getByToken( ca12puppijetToken_ , CA12PUPPI );
    iEvent.getByToken( ca12PuppiSoftDropSubjetsToken_ , CA12PUPPIsub);
    //iEvent.getByToken( ca8puppijetToken_ , CA8PUPPI );
    //iEvent.getByToken( ca8PuppiSoftDropSubjetsToken_ , CA8PUPPIsub);
  }


  //cout<<"\nJet CA8PUPPI "<<endl;
  //for (const pat::Jet &ijet : *CA8PUPPI) {
      //cout<<"\nJet CA8PUPPI "<<" with pT "<<ijet.pt()<< " mass " << ijet.mass() << endl;
  //}

  int count_AK8CHS = 0;
  int count_hadAK8CHS = 0;
  // int count_AK8CHS_good = 0;
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
  //double closestAK8_to_Jet_bdisc=-10;
  


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

  int count_CA12PUPPI = 0;
  for (const pat::Jet &ijet : *CA12PUPPI) {
    cout<<"\nJet CA12 "<<count_CA12PUPPI<<" with pT "<<ijet.pt()<<" sdMass "<<ijet.userFloat("ca12PFJetsPuppiSoftDropMass") <<endl;
    if (count_CA12PUPPI ==0){

       leading_CA12.SetPtEtaPhiM(ijet.pt(),ijet.eta(),ijet.phi(),ijet.userFloat("ca12PFJetsPuppiSoftDropMass"));
       CA12JetPtRaw = ijet.pt();
       CA12JetEtaRaw = ijet.eta();
       CA12JetPhiRaw = ijet.phi();
       CA12JetMassRaw = ijet.userFloat("ca12PFJetsPuppiSoftDropMass");

       
    for (const pat::Jet &isubjet : *CA12PUPPIsub) {
      leading_CA12_subjet.SetPtEtaPhiM(isubjet.pt(),isubjet.eta(),isubjet.phi(),isubjet.mass());
      if(count_subjets==0){
          CA12Jetsubjet0bdisc = isubjet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
          CA12Jetsubjet0pt    = isubjet.pt();
          CA12Jetsubjet0mass  = isubjet.mass();
          CA12Jetsubjet0eta   = isubjet.eta();
          CA12Jetsubjet0phi   = isubjet.phi();
          CA12Jetsubjet0area = isubjet.jetArea();

      }
      if(count_subjets==1){
          CA12Jetsubjet1bdisc = isubjet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
          CA12Jetsubjet1pt    = isubjet.pt();
          CA12Jetsubjet1mass  = isubjet.mass();
          CA12Jetsubjet1eta   = isubjet.eta();
          CA12Jetsubjet1phi   = isubjet.phi();
          CA12Jetsubjet1area = isubjet.jetArea();   
          break;   
      }
    count_subjets++;
    }

    //   CA12Jetsubjet0bdisc
    //   CA12Jetsubjet1bdisc
//
    //   CA12Jetsubjet0pt   
    //   CA12Jetsubjet0mass 
    //   CA12Jetsubjet0eta  
    //   CA12Jetsubjet0phi  
    //   CA12Jetsubjet0area
//
    //   CA12Jetsubjet1pt   
    //   CA12Jetsubjet1mass 
    //   CA12Jetsubjet1eta  
    //   CA12Jetsubjet1phi  
    //   CA12Jetsubjet1area 
//
    }



        count_CA12PUPPI++;
  }




      Float_t CA12JetPtRaw                               ;      
      Float_t CA12JetEtaRaw                              ;
      Float_t CA12JetPhiRaw                              ;
      Float_t CA12JetMassRaw                             ;

      Float_t CA12JetTau1                           ;
      Float_t CA12JetTau2                           ;
      Float_t CA12JetTau3                           ;
      Float_t CA12JetTau4                           ;
      Float_t CA12JetTau32                          ;
      Float_t CA12JetTau21                          ;

      Float_t CA12Jetsubjet0bdisc                 ;
      Float_t CA12Jetsubjet1bdisc                 ;
      Float_t CA12Jetmaxbdisc                     ;

      Float_t CA12Jetsubjet0pt                    ;
      Float_t CA12Jetsubjet0mass                  ;
      Float_t CA12Jetsubjet0eta                   ;
      Float_t CA12Jetsubjet0phi                   ;
      Float_t CA12Jetsubjet0area                  ;

      Float_t CA12Jetsubjet1pt                    ;
      Float_t CA12Jetsubjet1mass                  ;
      Float_t CA12Jetsubjet1eta                   ;
      Float_t CA12Jetsubjet1phi                   ;
      Float_t CA12Jetsubjet1area                  ;


  for (const pat::Jet &ijet : *AK8CHS) {
    // if (count_AK8CHS>1) break;
    if (count_AK8CHS==0 && ijet.pt()<30) break;
    if (verbose_) cout<<"\nJet "<<count_AK8CHS<<" with pT "<<ijet.pt()<<" sdMass "<<ijet.userFloat("ak8PFJetsCHSSoftDropMass")<<endl;

    //------------------------------------
    // Noise jet ID
    //------------------------------------    
    double NHF       = ijet.neutralHadronEnergyFraction();
    double NEMF      = ijet.neutralEmEnergyFraction();
    double CHF       = ijet.chargedHadronEnergyFraction();
    // double MUF       = ijet.muonEnergyFraction();
    double CEMF      = ijet.chargedEmEnergyFraction();
    double NumConst  = ijet.chargedMultiplicity()+ijet.neutralMultiplicity();
    double NM        = ijet.neutralMultiplicity();
    double CM        = ijet.chargedMultiplicity(); 
    double eta       = ijet.eta();

    bool goodJet_looseJetID =  
         ( fabs(eta) <= 2.4 && NHF < 0.99 && NEMF < 0.99 && NumConst >1 && CHF > 0.0  && CM > 0 && CEMF < 0.99   ) 
      || ( fabs(eta) <= 2.7 && fabs(eta) > 2.4 && NHF < 0.99 && NEMF < 0.99 && NumConst >1 ) 
      || ( fabs(eta) <= 3.0 && fabs(eta) > 2.7 && NHF < 0.98 && NEMF > 0.01 && NM > 2 ) 
      || ( fabs(eta)  > 3.0 && NEMF < 0.9 && NM > 10 );

    if (verbose_ && goodJet_looseJetID) cout<<"   -> goodJet "<<endl;

    if (!goodJet_looseJetID) {
      if(verbose_) cout<<"   -> bad AK8 jet. skip.  ( pt "<<ijet.pt()<<" eta "<<ijet.eta()<<" NumConst "<<NumConst<<" )"<<endl;
      continue;
    }
    // count_AK8CHS_good ++;

    //------------------------------------
    // AK8CHS JEC correction 
    //------------------------------------
    reco::Candidate::LorentzVector uncorrJet = ijet.correctedP4(0);
    JetCorrectorAK8chs->setJetEta( uncorrJet.eta() );
    JetCorrectorAK8chs->setJetPt ( uncorrJet.pt() );
    JetCorrectorAK8chs->setJetE  ( uncorrJet.energy() );
    JetCorrectorAK8chs->setJetA  ( ijet.jetArea() );
    JetCorrectorAK8chs->setRho   ( rho );
    JetCorrectorAK8chs->setNPV   ( nvtx );
    double corr = JetCorrectorAK8chs->getCorrection();

    reco::Candidate::LorentzVector corrJet = corr * uncorrJet;
    if (verbose_) cout<<"   -> uncorrected AK8 jet pt "<<uncorrJet.pt()<<" corrected jet pt "<<corrJet.pt()<<endl;
    
    TLorentzVector jet_p4;
    jet_p4.SetPtEtaPhiM( corrJet.pt(), corrJet.eta(), corrJet.phi(), corrJet.mass() );


    if ( jet_p4.Pt() < 30 ) continue;

    if(count_AK8CHS==0) AK8jet0_P4corr = jet_p4;
    if(count_AK8CHS==1) AK8jet1_P4corr = jet_p4;

    if(count_AK8CHS>1) {
      double DeltaR_Jet0  = AK8jet0_P4corr.DeltaR( jet_p4 );
      double DeltaR_Jet1  = AK8jet1_P4corr.DeltaR( jet_p4 );
      if (DeltaR_Jet0 < closestAK8_to_Jet_dR) {
        closestAK8_to_Jet_dR = DeltaR_Jet0;
        closestAK8_to_Jet_P4 = jet_p4 ;
        //closestAK8_to_Jet_bdisc = ijet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"); ;
      }
    }
    // Only need 3rd jet and beyond to find jet 0 and jet1 nearest neighbors (above)
    if(count_AK8CHS>4) continue; 

    //------------------------------------
    // AK8CHS JEC L23 correction
    //------------------------------------
    JetCorrectorAK8chs->setJetEta( uncorrJet.eta() );
    JetCorrectorAK8chs->setJetPt ( uncorrJet.pt() );
    JetCorrectorAK8chs->setJetE  ( uncorrJet.energy() );
    JetCorrectorAK8chs->setJetA  ( ijet.jetArea() );
    JetCorrectorAK8chs->setRho   ( rho );
    JetCorrectorAK8chs->setNPV   ( nvtx );
    // getSubCorrections member function returns the vector of the subcorrections UP to the given level. For example in the example above, factors[0] is the L1 correction and factors[3] is the L1+L2+L3+Residual correction. 
    vector<float> factors = JetCorrectorAK8chs->getSubCorrections();
    float corr_factor_L1      = 1.0;
    float corr_factor_L12     = 1.0;
    float corr_factor_L123    = 1.0;
    float corr_factor_L123res = 1.0;
    if (factors.size() > 0) corr_factor_L1       = factors[0];
    if (factors.size() > 1) corr_factor_L12      = factors[1];
    if (factors.size() > 2) corr_factor_L123     = factors[2];
    if (factors.size() > 3) corr_factor_L123res  = factors[3];
    double corr_factor_L2 = corr_factor_L12/corr_factor_L1;
    double corr_factor_L3 = corr_factor_L123/corr_factor_L12;
    double corr_factor_res = corr_factor_L123res/corr_factor_L123;
    //double corr_factor_L23 = corr_factor_L2*corr_factor_L3;
    double corr_factor_L23res = corr_factor_L2*corr_factor_L3*corr_factor_res;

    //------------------------------------
    // AK8CHS JEC uncertainty
    //------------------------------------
    double corrDn_L23  = 1.0;
    double corrDn_L123 = 1.0;
    JetCorrUncertAK8chs->setJetPhi(  corrJet.phi()  );
    JetCorrUncertAK8chs->setJetEta(  corrJet.eta()  );
    JetCorrUncertAK8chs->setJetPt(   corrJet.pt()   );
    double corrDn_temp1 = JetCorrUncertAK8chs->getUncertainty(0);
    corrDn_L23   = corr_factor_L23res - corrDn_temp1;
    corrDn_L123 = corr - corrDn_temp1;
    double corrUp_L23  = 1.0;
    double corrUp_L123 = 1.0;
    JetCorrUncertAK8chs->setJetPhi(  corrJet.phi()  );
    JetCorrUncertAK8chs->setJetEta(  corrJet.eta()  );
    JetCorrUncertAK8chs->setJetPt(   corrJet.pt()   );
    double corrUp_temp1 = JetCorrUncertAK8chs->getUncertainty(1);
    corrUp_L23   = corr_factor_L23res + corrUp_temp1;
    corrUp_L123 = corr + corrUp_temp1;

    if (verbose_) cout<<"   -> corr "<<corr<<" corr_factor_L23res "<<corr_factor_L23res<<" corrDn_L123 "<<corrDn_L123<<" corrUp_L123 "<<corrUp_L123<<endl;

    //------------------------------------
    // AK8 JER SF
    //------------------------------------
  
    TLorentzVector GenJetMatched;
    double ptsmear   = 1;
    double ptsmearUp = 1;
    double ptsmearDn = 1;

    if (!iEvent.isRealData()) {
      if (verbose_) cout<<"   Get JER SF"<<endl;

      // get genjet
      double genpt = 0;
      const reco::GenJet* genJet = ijet.genJet();
      bool foundgenjet = false;
      if (genJet) {
        foundgenjet=true;
        genpt = genJet->pt();
        GenJetMatched.SetPtEtaPhiM( genJet->pt(), genJet->eta(), genJet->phi(), genJet->mass() );
        if (verbose_) cout<<"     -> found ak8 genJet pt "<<genJet->pt()<<" mass "<<genJet->mass()<<endl;
      }
      else{ if(verbose_)cout<<"     -> did not find ak8 genJet"<<endl;}
    
      // Set parameters needed for jet resolution and scale factors
      JME::JetParameters jer_parameters;
      jer_parameters.setJetPt ( corrJet.pt()  );
      jer_parameters.setJetEta( corrJet.eta() );
      jer_parameters.setRho   ( rho           );

      // Get resolution
      double res = jet_resolution_AK8CHS.getResolution(jer_parameters); 

      // Get scale factors
      double jer_sf    = jer_scaler.getScaleFactor(jer_parameters                   );
      double jer_sf_up = jer_scaler.getScaleFactor(jer_parameters , Variation::UP   );
      double jer_sf_dn = jer_scaler.getScaleFactor(jer_parameters , Variation::DOWN );
      if (verbose_) std::cout << "     -> JER Scale factors (Nominal / Up / Down) : " << jer_sf << " / " << jer_sf_up << " / " << jer_sf_dn <<"    & Resolution :"<<res<< std::endl;
     
      // Get Smearings  
      // --- If well matched, smear based on GenJet, If not well matched,  gaussian smear based on resolution
      TLorentzVector AK8JetP4;
      AK8JetP4.SetPtEtaPhiM( corrJet.pt(), corrJet.eta(), corrJet.phi(), corrJet.mass() );
      double DeltaR_gen_reco  = AK8JetP4.DeltaR( GenJetMatched );
      double DeltaPt_gen_reco = AK8JetP4.Pt() - GenJetMatched.Pt()  ;
      double jet_distance_param = 0.4; 
      if (verbose_) cout<<"     -> gen pt "<<GenJetMatched.Pt()<<" reco pt "<<AK8JetP4.Pt()<<"  delta "<<DeltaPt_gen_reco<<endl;

      if (genJet && (DeltaR_gen_reco<jet_distance_param/2.0) && (std::abs(DeltaPt_gen_reco)<(3*res*AK8JetP4.Pt())) ) {
        if (verbose_) cout<<"     -> Well matched"<<endl;
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
          cout<<"     -> Not well matched. DeltaR_gen_reco "<<DeltaR_gen_reco<<" DeltaPt_gen_reco "<<DeltaPt_gen_reco<<" 3*res*AK4JetP4.Pt()) "<<3*res*AK8JetP4.Pt();
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

    if (verbose_) cout<<"   -> ptsmear "<<ptsmear<<" ptsmearUp "<< ptsmearDn<<" ptsmearDn "<< ptsmearUp<<endl;

    //------------------------------------
    // AK8CHS variables 
    //------------------------------------
    // double pt           = corrJet.pt();
    // double mass         = corrJet.mass();
    // double eta          = corrJet.eta();
    // double phi          = corrJet.phi();
    // double rapidity     = ijet.rapidity();
    // double ndau         = ijet.numberOfDaughters();

    double tau1         = 99;
    double tau2         = 99;
    double tau3         = 99;
    double tau4         = 99;
    double prunedMass   = ijet.userFloat("ak8PFJetsCHSPrunedMass");
    double softDropMass = ijet.userFloat("ak8PFJetsCHSSoftDropMass");
    double trimmedMass  = -1;


    if (useToolbox_){
      tau1         = ijet.userFloat("NjettinessAK8CHS:tau1");
      tau2         = ijet.userFloat("NjettinessAK8CHS:tau2");
      tau3         = ijet.userFloat("NjettinessAK8CHS:tau3");
      tau4         = ijet.userFloat("NjettinessAK8CHS:tau4");
      trimmedMass  = ijet.userFloat("ak8PFJetsCHSTrimmedMass"); 
    }
    else{
      tau1         = ijet.userFloat("NjettinessAK8:tau1");
      tau2         = ijet.userFloat("NjettinessAK8:tau2");
      tau3         = ijet.userFloat("NjettinessAK8:tau3");
    }
    double tau21        = 99;
    double tau32        = 99;

    if (tau1!=0) tau21 = tau2/tau1;
    if (tau2!=0) tau32 = tau3/tau2;


    // //-----------
    // // get jet constituents
    // //-----------
    // std::vector<fastjet::PseudoJet> FJparticles;
    // for (unsigned i = 0; i < ijet.numberOfDaughters() ; i++){
    //   const reco::PFCandidate* this_constituent = dynamic_cast<const reco::PFCandidate*>(ijet.daughter(i));
    //   FJparticles.push_back( fastjet::PseudoJet( this_constituent->px(),
    //            this_constituent->py(),
    //          this_constituent->pz(),
    //          this_constituent->energy() ) );
    // }

    // //------------------------------------
    // // Recluster jet
    // //------------------------------------
    // double R = 0.8;
    // double maxrap = 5.0;
    // unsigned int n_repeat = 1; // default is 1
    // double ghost_area = 0.01; // this is the default
    // fastjet::GhostedAreaSpec area_spec(maxrap, n_repeat, ghost_area);       
    // //fastjet::AreaDefinition area_def(fastjet::active_area, area_spec);
    // fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, area_spec);
    // fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm, R);
    // fastjet::ClusterSequenceArea cs(FJparticles, jet_def, area_def);
    // vector<fastjet::PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());


    // double z_cut = 0.10;
    // double beta  = 0.0;
    // contrib::SoftDrop sd(beta, z_cut);
    // cout << "SoftDrop groomer is: " << sd.description() << endl;

    // for (unsigned jjet = 0; jjet < jets.size(); jjet++) {
    //   // Run SoftDrop and examine the output
    //   PseudoJet sd_jet = sd(jets[jjet]);
    //   cout <<"sd mass "<<sd_jet.m() << endl;
    //   // cout << "original    jet: " << jets[jjet] << endl;
    //   // cout << "SoftDropped jet: " << sd_jet << endl;

    //   assert(sd_jet != 0); //because soft drop is a groomer (not a tagger), it should always return a soft-dropped jet
       
    //   cout << "  delta_R between subjets: " << sd_jet.structure_of<contrib::SoftDrop>().delta_R() << endl;
    //   cout << "  symmetry measure(z):     " << sd_jet.structure_of<contrib::SoftDrop>().symmetry() << endl;
    //   cout << "  mass drop(mu):           " << sd_jet.structure_of<contrib::SoftDrop>().mu() << endl;
    // }

    //------------------------------------
    // AK8PUPPI variables 
    //------------------------------------

    double puppi_pt             = -99;     
    double puppi_mass           = -99;     
    double puppi_eta            = -99;     
    double puppi_phi            = -99;     
    double puppi_area           = -99;     
    double puppi_tau1           = -99;     
    double puppi_tau2           = -99;     
    double puppi_tau3           = -99;     
    double puppi_tau4           = -99;     
    double puppi_prunedMass     = -1;     
    double puppi_trimmedMass    = -1;     
    // double puppi_softDropMass   = -1;     

    TLorentzVector AK8PUPPI_P4uncorr;
    TLorentzVector GenJetMatchedPuppi; 

    double puppi_CHF    = -99;
    double puppi_NHF    = -99;
    double puppi_CM     = -99;
    double puppi_NM     = -99;
    double puppi_NEF    = -99;
    double puppi_CEF    = -99;
    double puppi_MF     = -99;
    double puppi_Mult   = -99;
    double puppi_ECF1  = 0;
    double puppi_ECF2  = 0;
    double puppi_ECF3  = 0;
    double puppi_ECF4  = 0;
    double puppi_ECF5  = 0;

 
    // If you're getting jet variables from miniAOD
    if (!useToolbox_){
      puppi_pt           = ijet.userFloat("ak8PFJetsPuppiValueMap:pt");
      puppi_mass         = ijet.userFloat("ak8PFJetsPuppiValueMap:mass");
      puppi_eta          = ijet.userFloat("ak8PFJetsPuppiValueMap:eta");
      puppi_phi          = ijet.userFloat("ak8PFJetsPuppiValueMap:phi");
      puppi_tau1         = ijet.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1");
      puppi_tau2         = ijet.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2");
      puppi_tau3         = ijet.userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau3");

    }

    // If you've clustered AK8PUPPI using the toolbox, match the Puppi jets to AK8CHS and save the puppi variables

    double minDR_pup_chs = 99;
    if (useToolbox_){
      if (verbose_)  cout<<"   Puppi jet loop (toolbox)"<<endl; 
      int count_puppi_jet = 0;
      for (const pat::Jet &ipup : *AK8PUPPI) {  
        if (verbose_)  cout<<"    puppi jet "<<count_puppi_jet<<" uncorr "<<ipup.correctedP4(0).pt()<<" corr "<<ipup.pt()<<endl;
        double deltaRpup = deltaR(ijet.eta(), ijet.phi(), ipup.eta(), ipup.phi() );
        if (deltaRpup< minDR_pup_chs){
          minDR_pup_chs = deltaRpup;
          if (verbose_)  cout<<"      -> clostest puppi jet so far: deltaRpup "<<deltaRpup<<endl;
          if (deltaRpup<1.0){
            if (verbose_)  cout<<"      -> passes dR"<<endl;
            puppi_pt           = ipup.correctedP4(0).pt();
            puppi_mass         = ipup.correctedP4(0).mass();
            puppi_eta          = ipup.correctedP4(0).eta();
            puppi_phi          = ipup.correctedP4(0).phi();

            puppi_area         = ipup.jetArea();
            puppi_prunedMass   = ipup.userFloat("ak8PFJetsPuppiPrunedMass");
            puppi_trimmedMass  = ipup.userFloat("ak8PFJetsPuppiTrimmedMass");
            // puppi_softDropMass = ipup.userFloat("ak8PFJetsPuppiSoftDropMass");
            puppi_tau1         = ipup.userFloat("NjettinessAK8Puppi:tau1");
            puppi_tau2         = ipup.userFloat("NjettinessAK8Puppi:tau2");
            puppi_tau3         = ipup.userFloat("NjettinessAK8Puppi:tau3");
            puppi_tau4         = ipup.userFloat("NjettinessAK8Puppi:tau4");
            puppi_ECF1         = ipup.userFloat("ak8PFJetsPuppiECF:ecf1");
            puppi_ECF2         = ipup.userFloat("ak8PFJetsPuppiECF:ecf2");
            puppi_ECF3         = ipup.userFloat("ak8PFJetsPuppiECF:ecf3");
            puppi_ECF4         = ipup.userFloat("ak8PFJetsPuppiECF:ecf4");
            puppi_ECF5         = ipup.userFloat("ak8PFJetsPuppiECF:ecf5");

            puppi_CHF          = ipup.chargedHadronEnergy() / ipup.correctedP4(0).E()  ;  
            puppi_NHF          = ipup.neutralHadronEnergy() / ipup.correctedP4(0).E()  ;  
            puppi_CM           = ipup.chargedMultiplicity()  ;                  
            puppi_NM           = ipup.neutralMultiplicity()  ;                  
            puppi_NEF          = ipup.neutralEmEnergy() / ipup.correctedP4(0).E()  ;      
            puppi_CEF          = ipup.chargedEmEnergy() / ipup.correctedP4(0).E()  ;      
            puppi_MF           = ipup.muonEnergy() / ipup.correctedP4(0).E()  ;           
            puppi_Mult         = ipup.numberOfDaughters() ;   

            if (!iEvent.isRealData()){
              const reco::GenJet* genJet = ipup.genJet();
              if (genJet) {
                GenJetMatchedPuppi.SetPtEtaPhiM( genJet->pt(), genJet->eta(), genJet->phi(), genJet->mass() );
                if (verbose_) cout<<"      -> ak8puppi genJet pt "<<genJet->pt()<<" mass "<<genJet->mass()<<endl;
              }
            }
          }
        }
        count_puppi_jet++;
      }
    }
    AK8PUPPI_P4uncorr.SetPtEtaPhiM(puppi_pt, puppi_eta, puppi_phi, puppi_mass );
    if ( count_AK8CHS==0 ) GenJetMatchedPuppi0 = GenJetMatchedPuppi;
    if ( count_AK8CHS==1 ) GenJetMatchedPuppi1 = GenJetMatchedPuppi;

    if (minDR_pup_chs>1.0 && verbose_) cout<<"   Did not find matching PUPPI jet. Setting PUPPI variables to -99"<<endl;
    if (minDR_pup_chs<1.0 && verbose_) cout<<"   Found matching PUPPI jet with pt "<<puppi_pt<<" and mass " <<puppi_mass<<endl;

    double puppi_tau21        = 99;
    double puppi_tau32        = 99;
    if (puppi_tau1!=0) puppi_tau21 = puppi_tau2/puppi_tau1;
    if (puppi_tau2!=0) puppi_tau32 = puppi_tau3/puppi_tau2;

    //------------------------------------
    // AK8PUPPI JEC L23 correction
    //------------------------------------

    JetCorrectorAK8pup->setJetEta( puppi_eta );
    JetCorrectorAK8pup->setJetPt ( puppi_pt );
    JetCorrectorAK8pup->setJetE  ( AK8PUPPI_P4uncorr.E() );
    JetCorrectorAK8pup->setJetA  ( puppi_area );
    JetCorrectorAK8pup->setRho   ( rho );
    JetCorrectorAK8pup->setNPV   ( nvtx );
    // getSubCorrections member function returns the vector of the subcorrections UP to the given level. For example in the example above, factors[0] is the L1 correction and factors[3] is the L1+L2+L3+Residual correction. 
    vector<float> factorsAK8pup = JetCorrectorAK8pup->getSubCorrections();
    float corr_factorAK8pup_L1      = 1.0;
    float corr_factorAK8pup_L12     = 1.0;
    float corr_factorAK8pup_L123    = 1.0;
    float corr_factorAK8pup_L123res = 1.0;
    if (factors.size() > 0) corr_factorAK8pup_L1       = factorsAK8pup[0];
    if (factors.size() > 1) corr_factorAK8pup_L12      = factorsAK8pup[1];
    if (factors.size() > 2) corr_factorAK8pup_L123     = factorsAK8pup[2];
    if (factors.size() > 3) corr_factorAK8pup_L123res  = factorsAK8pup[3];
    double corr_factorAK8pup_L2 = corr_factorAK8pup_L12/corr_factorAK8pup_L1;
    double corr_factorAK8pup_L3 = corr_factorAK8pup_L123/corr_factorAK8pup_L12;
    double corr_factorAK8pup_res = corr_factorAK8pup_L123res/corr_factorAK8pup_L123;
    //double corr_factor_L23 = corr_factor_L2*corr_factor_L3;
    double corr_factorAK8pup_L23res = corr_factorAK8pup_L2*corr_factorAK8pup_L3*corr_factorAK8pup_res;

    TLorentzVector AK8PUPPI_P4corr;
    AK8PUPPI_P4corr = corr_factorAK8pup_L23res *  AK8PUPPI_P4uncorr;

    if(count_AK8CHS==0) PUPPIjet0_P4corr = AK8PUPPI_P4corr;
    if(count_AK8CHS==1) PUPPIjet1_P4corr = AK8PUPPI_P4corr;


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

      if (verbose_){
        cout<<"    -> AK8PUPPI JER recopt  "<<recopt <<endl;  
        cout<<"    -> AK8PUPPI JER genpt   "<<genpt  <<endl;  
        cout<<"    -> AK8PUPPI JER deltapt "<<deltapt<<endl;  
        cout<<"    -> AK8PUPPI JER ptsmear "<<pup_ptsmear<<endl;
        cout<<"    -> AK8PUPPI JER pup_ptsmearUp "<<pup_ptsmearUp<<endl;
        cout<<"    -> AK8PUPPI JER pup_ptsmearDn "<<pup_ptsmearDn<<endl;
      }
    }

    //------------------------------------
    // SoftDrop subjets
    //------------------------------------
    TLorentzVector sub0_P4_uncorr           ;
    TLorentzVector sub0_P4_L23res           ;
    // TLorentzVector sub0_P4_L23resCorrUp     ;
    // TLorentzVector sub0_P4_L23resCorrDn     ;
    // TLorentzVector sub0_P4_L23resPtSmear    ;
    // TLorentzVector sub0_P4_L23resPtSmearUp  ;
    // TLorentzVector sub0_P4_L23resPtSmearDn  ;
    TLorentzVector sub0_P4_L123res          ;
    // TLorentzVector sub0_P4_L123resCorrUp    ;
    // TLorentzVector sub0_P4_L123resCorrDn    ;

    TLorentzVector sub1_P4_uncorr           ;
    TLorentzVector sub1_P4_L23res           ;
    // TLorentzVector sub1_P4_L23resCorrUp     ;
    // TLorentzVector sub1_P4_L23resCorrDn     ;
    // TLorentzVector sub1_P4_L23resPtSmear    ;
    // TLorentzVector sub1_P4_L23resPtSmearUp  ;
    // TLorentzVector sub1_P4_L23resPtSmearDn  ;
    TLorentzVector sub1_P4_L123res          ;
    // TLorentzVector sub1_P4_L123resCorrUp    ;
    // TLorentzVector sub1_P4_L123resCorrDn    ;

    double sub0_area  = 0;
    double sub0_tau1  = 0;
    double sub0_tau2  = 0;
    double sub0_tau3  = 0;
    double sub0_flav_hadron  = 0;
    double sub0_flav_parton  = 0;
    double sub0_bdisc = 0;
    double sub0_genpt = 0;
    double sub1_area  = 0;
    double sub1_tau1  = 0;
    double sub1_tau2  = 0;
    double sub1_tau3  = 0;
    double sub1_flav_hadron  = 0;
    double sub1_flav_parton  = 0;
    double sub1_bdisc = 0;
    double sub1_genpt = 0;
    double mostMassiveSDsubjetMass = 0;
    int count_SD =0;

    int nsubjets_chs = 0;
    int nsubjets_pup = 0;

    if (!useToolbox_){
      auto const & sdSubjets = ijet.subjets("SoftDrop");
      for ( auto const & it : sdSubjets ) {
        double subjetPt       = it->correctedP4(0).pt();
        double subjetEta      = it->correctedP4(0).eta();
        double subjetPhi      = it->correctedP4(0).phi();
        double subjetMass     = it->correctedP4(0).mass();
        double subjetBdisc    = it->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"); 
        double deltaRsubjetJet = deltaR(ijet.eta(), ijet.phi(), subjetEta, subjetPhi);
      
        //------------------------------------
        // subjet JEC 
        //------------------------------------
        reco::Candidate::LorentzVector uncorrSubjet = it->correctedP4(0);
        JetCorrectorAK4chs -> setJetEta( uncorrSubjet.eta()    );
        JetCorrectorAK4chs -> setJetPt ( uncorrSubjet.pt()     );
        JetCorrectorAK4chs -> setJetE  ( uncorrSubjet.energy() );
        JetCorrectorAK4chs -> setJetA  ( it->jetArea()         );
        JetCorrectorAK4chs -> setRho   ( rho                   );
        JetCorrectorAK4chs -> setNPV   ( nvtx                  );
        double subjet_corr_factor_L123res_full = JetCorrectorAK4chs->getCorrection();
        reco::Candidate::LorentzVector corrSubjetL123res = subjet_corr_factor_L123res_full * uncorrSubjet;

        //------------------------------------
        // subjet L23 JEC 
        //------------------------------------
        JetCorrectorAK4chs->setJetEta( uncorrSubjet.eta()    );
        JetCorrectorAK4chs->setJetPt ( uncorrSubjet.pt()     );
        JetCorrectorAK4chs->setJetE  ( uncorrSubjet.energy() );
        JetCorrectorAK4chs->setJetA  ( it->jetArea()         );
        JetCorrectorAK4chs->setRho   ( rho                   );
        JetCorrectorAK4chs->setNPV   ( nvtx                  );
        // getSubCorrections member function returns the vector of the subcorrections UP to the given level. For example in the example above, factors[0] is the L1 correction and factors[3] is the L1+L2+L3+Residual correction. 
        vector<float> subjet_factors = JetCorrectorAK4chs->getSubCorrections();
        float subjet_corr_factor_L1      = 1.0;
        float subjet_corr_factor_L12     = 1.0;
        float subjet_corr_factor_L123    = 1.0;
        float subjet_corr_factor_L123res = 1.0;
        if (factors.size() > 0) subjet_corr_factor_L1      = subjet_factors[0];
        if (factors.size() > 1) subjet_corr_factor_L12     = subjet_factors[1];
        if (factors.size() > 2) subjet_corr_factor_L123    = subjet_factors[2];
        if (factors.size() > 3) subjet_corr_factor_L123res = subjet_factors[3];
        double subjet_corr_factor_L2     = subjet_corr_factor_L12     / subjet_corr_factor_L1     ;
        double subjet_corr_factor_L3     = subjet_corr_factor_L123    / subjet_corr_factor_L12    ;
        double subjet_corr_factor_res    = subjet_corr_factor_L123res / subjet_corr_factor_L123   ;
        double subjet_corr_factor_L23    = subjet_corr_factor_L2 * subjet_corr_factor_L3     ;
        double subjet_corr_factor_L23res = subjet_corr_factor_L2 * subjet_corr_factor_L3 * subjet_corr_factor_res    ;
        if (verbose_) cout<<"subjet corr: L1 "<<subjet_corr_factor_L1<<" L23 "<<subjet_corr_factor_L23<<" L23res "<<subjet_corr_factor_L23res<<" L123res"<<subjet_corr_factor_L123res<<endl;
        reco::Candidate::LorentzVector corrSubjetL23res   = subjet_corr_factor_L23res * uncorrSubjet;
        
        // //------------------------------------
        // // subjet JEC uncertainty
        // //------------------------------------
        // double subjet_corrDn_L23 =   1.0;
        // double subjet_corrDn_L123 = 1.0;
        // JetCorrUncertAK4chs->setJetPhi(  corrSubjetL123res.phi()  );
        // JetCorrUncertAK4chs->setJetEta(  corrSubjetL123res.eta()  );
        // JetCorrUncertAK4chs->setJetPt(   corrSubjetL123res.pt()   );
        // double corrDn_temp2 = JetCorrUncertAK4chs->getUncertainty(0);
        // subjet_corrDn_L23   = subjet_corr_factor_L23res - corrDn_temp2;
        // subjet_corrDn_L123  = subjet_corr_factor_L123res_full - corrDn_temp2;

        // double subjet_corrUp_L23   = 1.0;
        // double subjet_corrUp_L123 = 1.0;
        // JetCorrUncertAK4chs->setJetPhi(  corrSubjetL123res.phi()  );
        // JetCorrUncertAK4chs->setJetEta(  corrSubjetL123res.eta()  );
        // JetCorrUncertAK4chs->setJetPt(   corrSubjetL123res.pt()   );
        // double corrUp_temp2 = JetCorrUncertAK4chs->getUncertainty(1);
        // subjet_corrUp_L23   = subjet_corr_factor_L23res + corrUp_temp2;
        // subjet_corrUp_L123  = subjet_corr_factor_L123res_full + corrUp_temp2;

        // reco::Candidate::LorentzVector corrSubjetL123resCorrDn  = subjet_corrDn_L123  * uncorrSubjet;
        // reco::Candidate::LorentzVector corrSubjetL123resCorrUp  = subjet_corrUp_L123  * uncorrSubjet;
        // reco::Candidate::LorentzVector corrSubjetL23resCorrDn   = subjet_corrDn_L23   * uncorrSubjet;
        // reco::Candidate::LorentzVector corrSubjetL23resCorrUp   = subjet_corrUp_L23   * uncorrSubjet;
     

        //------------------------------------
        // subjet values for Tree
        //------------------------------------
        if (count_SD==0){
          sub0_P4_uncorr            .SetPtEtaPhiM( subjetPt, subjetEta, subjetPhi, subjetMass);
          sub0_P4_L123res           .SetPtEtaPhiM( corrSubjetL123res.pt()   , corrSubjetL123res.eta()   , corrSubjetL123res.phi()   , corrSubjetL123res.mass()    );
          sub0_P4_L23res            .SetPtEtaPhiM( corrSubjetL23res.pt()    , corrSubjetL23res.eta()    , corrSubjetL23res.phi()    , corrSubjetL23res.mass()     );
          // sub0_P4_L123resCorrUp    .SetPtEtaPhiM( corrSubjetL123resCorrUp.pt() , corrSubjetL123resCorrUp.eta() , corrSubjetL123resCorrUp.phi() , corrSubjetL123resCorrUp.mass()  );
          // sub0_P4_L23resCorrUp     .SetPtEtaPhiM( corrSubjetL23resCorrUp.pt()  , corrSubjetL23resCorrUp.eta()  , corrSubjetL23resCorrUp.phi()  , corrSubjetL23resCorrUp.mass()   );
          // sub0_P4_L123resCorrDn    .SetPtEtaPhiM( corrSubjetL123resCorrDn.pt() , corrSubjetL123resCorrDn.eta() , corrSubjetL123resCorrDn.phi() , corrSubjetL123resCorrUp.mass()  );
          // sub0_P4_L23resCorrDn     .SetPtEtaPhiM( corrSubjetL23resCorrDn.pt()  , corrSubjetL23resCorrDn.eta()  , corrSubjetL23resCorrDn.phi()  , corrSubjetL23res.mass()     );
          sub0_area   = it->jetArea() ;
          sub0_flav_parton   = it->partonFlavour();
          sub0_flav_hadron   = it->hadronFlavour();
          sub0_bdisc  = subjetBdisc;
        }
        if (count_SD==1) {
          sub1_P4_uncorr          .SetPtEtaPhiM( subjetPt, subjetEta, subjetPhi, subjetMass);
          sub1_P4_L123res         .SetPtEtaPhiM( corrSubjetL123res.pt()   , corrSubjetL123res.eta()   , corrSubjetL123res.phi()   , corrSubjetL123res.mass()    );
          sub1_P4_L23res          .SetPtEtaPhiM( corrSubjetL23res.pt()    , corrSubjetL23res.eta()    , corrSubjetL23res.phi()    , corrSubjetL23res.mass()     );
          // sub1_P4_L123resCorrUp  .SetPtEtaPhiM( corrSubjetL123resCorrUp.pt() , corrSubjetL123resCorrUp.eta() , corrSubjetL123resCorrUp.phi() , corrSubjetL123resCorrUp.mass()  );
          // sub1_P4_L23resCorrUp   .SetPtEtaPhiM( corrSubjetL23resCorrUp.pt()  , corrSubjetL23resCorrUp.eta()  , corrSubjetL23resCorrUp.phi()  , corrSubjetL23resCorrUp.mass()   );
          // sub1_P4_L123resCorrDn  .SetPtEtaPhiM( corrSubjetL123resCorrDn.pt() , corrSubjetL123resCorrDn.eta() , corrSubjetL123resCorrDn.phi() , corrSubjetL123resCorrUp.mass()  );
          // sub1_P4_L23resCorrDn   .SetPtEtaPhiM( corrSubjetL23resCorrDn.pt()  , corrSubjetL23resCorrDn.eta()  , corrSubjetL23resCorrDn.phi()  , corrSubjetL23res.mass()     );
          sub1_area   = it->jetArea() ;
          sub1_flav_parton   = it->partonFlavour();
          sub1_flav_hadron   = it->hadronFlavour();
          sub1_bdisc  = subjetBdisc;
        }
        if (subjetMass > mostMassiveSDsubjetMass) mostMassiveSDsubjetMass = subjetMass;

        if (verbose_) cout<<" SD Subjet pt "<<subjetPt<<" Eta "<<subjetEta<<" deltaRsubjetJet "<<deltaRsubjetJet<<" Mass "<<subjetMass<<" Bdisc "<<subjetBdisc<<endl;
        count_SD++;
      }
    }
  
    if (useToolbox_){
      if (verbose_) cout<<"   Toolbox AK8 jets. Find chs softdrop subjets "<<endl;
   
      int count_all_subjets =0;
      int count_matched_subjets =0;
      double closest_DR = 99;
      double closest_i = -1;
      double second_closest_DR = 99;
      double second_closest_i  = -1;

      // Loop once to find the subjets closest to the AK8 jet
      for (const pat::Jet &isub : *AK8CHSsub) {  
  
        double subjetPt       = isub.correctedP4(0).pt();
        double subjetEta      = isub.correctedP4(0).eta();
        double subjetPhi      = isub.correctedP4(0).phi();
        double subjetMass     = isub.correctedP4(0).mass();
        // double subjetBdisc    = isub.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"); 

        double deltaRsubjetJet = deltaR(ijet.eta(), ijet.phi(), subjetEta, subjetPhi);

        if (verbose_) cout<<"     Subjet "<<count_all_subjets<<"   "<<subjetMass<<" "<<subjetPt<<" "<<deltaRsubjetJet<<endl;

        if (deltaRsubjetJet<closest_DR){
          second_closest_DR = closest_DR;
          closest_DR        = deltaRsubjetJet;
          second_closest_i  = closest_i;
          closest_i         = count_all_subjets;
        }
        else if (deltaRsubjetJet<second_closest_DR){
          second_closest_DR = deltaRsubjetJet ;
          second_closest_i  = count_all_subjets;
        }
        count_all_subjets++;
      }
      // Loop a second time. If one of the two closest subjets matches the dR requirement save its infromation. Subjet 0 = hardest.
      count_all_subjets =0;
      for (const pat::Jet &isub : *AK8CHSsub) {  
        
  
        double subjetPt       = isub.correctedP4(0).pt();
        double subjetEta      = isub.correctedP4(0).eta();
        double subjetPhi      = isub.correctedP4(0).phi();
        double subjetMass     = isub.correctedP4(0).mass();
        double subjetBdisc    = isub.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"); 

        double deltaRsubjetJet = deltaR(ijet.eta(), ijet.phi(), subjetEta, subjetPhi);

        if (verbose_) cout<<"     Subjet "<<count_all_subjets<<"   "<<subjetMass<<" "<<subjetPt<<" "<<deltaRsubjetJet<<endl;

        if ( count_all_subjets==closest_i || count_all_subjets==second_closest_i ){
          if (verbose_) cout<<"      -> one of two closest "<<endl;
          if (deltaRsubjetJet<0.8){
            nsubjets_chs++;
            if (verbose_) cout<<"      -> dR matched subjet with mass "<< subjetMass<<endl;

            count_matched_subjets++;
            //------------------------------------
            // subjet JEC 
            //------------------------------------
            reco::Candidate::LorentzVector uncorrSubjet = isub.correctedP4(0);
            JetCorrectorAK4chs -> setJetEta( uncorrSubjet.eta()    );
            JetCorrectorAK4chs -> setJetPt ( uncorrSubjet.pt()     );
            JetCorrectorAK4chs -> setJetE  ( uncorrSubjet.energy() );
            JetCorrectorAK4chs -> setJetA  ( isub.jetArea()        );
            JetCorrectorAK4chs -> setRho   ( rho                   );
            JetCorrectorAK4chs -> setNPV   ( nvtx                  );
            double subjet_corr_factor_L123res_full = JetCorrectorAK4chs->getCorrection();
            reco::Candidate::LorentzVector corrSubjetL123res = subjet_corr_factor_L123res_full * uncorrSubjet;

            //------------------------------------
            // subjet L23 JEC 
            //------------------------------------
            JetCorrectorAK4chs->setJetEta( uncorrSubjet.eta()    );
            JetCorrectorAK4chs->setJetPt ( uncorrSubjet.pt()     );
            JetCorrectorAK4chs->setJetE  ( uncorrSubjet.energy() );
            JetCorrectorAK4chs->setJetA  ( isub.jetArea()         );
            JetCorrectorAK4chs->setRho   ( rho                   );
            JetCorrectorAK4chs->setNPV   ( nvtx                  );
            // getSubCorrections member function returns the vector of the subcorrections UP to the given level. For example in the example above, factors[0] is the L1 correction and factors[3] is the L1+L2+L3+Residual correction. 
            vector<float> subjet_factors = JetCorrectorAK4chs->getSubCorrections();
            float subjet_corr_factor_L1      = 1.0;
            float subjet_corr_factor_L12     = 1.0;
            float subjet_corr_factor_L123    = 1.0;
            float subjet_corr_factor_L123res = 1.0;
            if (factors.size() > 0) subjet_corr_factor_L1      = subjet_factors[0];
            if (factors.size() > 1) subjet_corr_factor_L12     = subjet_factors[1];
            if (factors.size() > 2) subjet_corr_factor_L123    = subjet_factors[2];
            if (factors.size() > 3) subjet_corr_factor_L123res = subjet_factors[3];
            double subjet_corr_factor_L2     = subjet_corr_factor_L12     / subjet_corr_factor_L1     ;
            double subjet_corr_factor_L3     = subjet_corr_factor_L123    / subjet_corr_factor_L12    ;
            double subjet_corr_factor_res    = subjet_corr_factor_L123res / subjet_corr_factor_L123   ;
            double subjet_corr_factor_L23    = subjet_corr_factor_L2 * subjet_corr_factor_L3     ;
            double subjet_corr_factor_L23res = subjet_corr_factor_L2 * subjet_corr_factor_L3 * subjet_corr_factor_res    ;
            if (verbose_) cout<<"        -> subjet corr: L1 "<<subjet_corr_factor_L1<<" L23 "<<subjet_corr_factor_L23<<" L23res "<<subjet_corr_factor_L23res<<" L123res "<<subjet_corr_factor_L123res<<endl;
            reco::Candidate::LorentzVector corrSubjetL23res   = subjet_corr_factor_L23res * uncorrSubjet;
            
            // //------------------------------------
            // // subjet JEC uncertainty
            // //------------------------------------
            // double subjet_corrDn_L23 =   1.0;
            // double subjet_corrDn_L123 = 1.0;
            // JetCorrUncertAK4chs->setJetPhi(  corrSubjetL123res.phi()  );
            // JetCorrUncertAK4chs->setJetEta(  corrSubjetL123res.eta()  );
            // JetCorrUncertAK4chs->setJetPt(   corrSubjetL123res.pt()   );
            // double corrDn_temp2 = JetCorrUncertAK4chs->getUncertainty(0);
            // subjet_corrDn_L23   = subjet_corr_factor_L23res - corrDn_temp2;
            // subjet_corrDn_L123  = subjet_corr_factor_L123res_full - corrDn_temp2;

            // double subjet_corrUp_L23   = 1.0;
            // double subjet_corrUp_L123 = 1.0;
            // JetCorrUncertAK4chs->setJetPhi(  corrSubjetL123res.phi()  );
            // JetCorrUncertAK4chs->setJetEta(  corrSubjetL123res.eta()  );
            // JetCorrUncertAK4chs->setJetPt(   corrSubjetL123res.pt()   );
            // double corrUp_temp2 = JetCorrUncertAK4chs->getUncertainty(1);
            // subjet_corrUp_L23   = subjet_corr_factor_L23res + corrUp_temp2;
            // subjet_corrUp_L123  = subjet_corr_factor_L123res_full + corrUp_temp2;

            // reco::Candidate::LorentzVector corrSubjetL123resCorrDn  = subjet_corrDn_L123  * uncorrSubjet;
            // reco::Candidate::LorentzVector corrSubjetL123resCorrUp  = subjet_corrUp_L123  * uncorrSubjet;
            // reco::Candidate::LorentzVector corrSubjetL23resCorrDn   = subjet_corrDn_L23   * uncorrSubjet;
            // reco::Candidate::LorentzVector corrSubjetL23resCorrUp   = subjet_corrUp_L23   * uncorrSubjet;
           
            // //------------------------------------
            // // subjet JER SF
            // //------------------------------------

            // Doesn't work. Matched gensubjets don't match well

            // TLorentzVector GenSubJet;
            // double ptsmear   = 1;
            // double ptsmearUp = 1;
            // double ptsmearDn = 1;
            // if (!iEvent.isRealData()){
            //   const reco::GenJet* genJet = isub.genJet();
            //   if (genJet) {
            //     GenSubJet.SetPtEtaPhiM( genJet->pt(), genJet->eta(), genJet->phi(), genJet->mass() );
            //     if (verbose_) cout<<"  SD subjet genJet pt "<<genJet->pt()<<" mass "<<genJet->mass()<<" reco pt "<<subjetPt<<" reco mass "<<subjetMass<<endl;
            //   }
            //   double jer_sf    = jer_scaler.getScaleFactor({{JME::Binning::JetEta, corrSubjetL23res.eta()}});
            //   double jer_sf_up = jer_scaler.getScaleFactor({{JME::Binning::JetEta, corrSubjetL23res.eta()}}, Variation::UP);
            //   double jer_sf_dn = jer_scaler.getScaleFactor({{JME::Binning::JetEta, corrSubjetL23res.eta()}}, Variation::DOWN);
            //   if (verbose_) std::cout << " SD subjet JER Scale factors (Nominal / Up / Down) : " << jer_sf << " / " << jer_sf_up << " / " << jer_sf_dn << std::endl;
            //   double recopt    = corrSubjetL23res.pt();
            //   double genpt     = GenJetMatched.Perp();
            //   double deltapt   = (recopt-genpt)*(jer_sf-1.0);
            //   double deltaptUp = (recopt-genpt)*(jer_sf_up-1.0);
            //   double deltaptDn = (recopt-genpt)*(jer_sf_dn-1.0);
            //   ptsmear   = std::max((double)0.0, (recopt+deltapt)/recopt     );
            //   ptsmearUp = std::max((double)0.0, (recopt+deltaptUp)/recopt   );
            //   ptsmearDn = std::max((double)0.0, (recopt+deltaptDn)/recopt   );
            //   if (verbose_) std::cout<<" SD subjet ptsmear "<<ptsmear<<" ptsmearUp "<<ptsmearUp<<" ptsmearDn "<<ptsmearDn<<endl;
            // }
            // reco::Candidate::LorentzVector corrSubjetL23resPtSmear   = ptsmear * corrSubjetL23res ;
            // reco::Candidate::LorentzVector corrSubjetL23resPtSmearUp = ptsmearUp * corrSubjetL23res ;
            // reco::Candidate::LorentzVector corrSubjetL23resPtSmearDn = ptsmearDn * corrSubjetL23res ;

            //------------------------------------
            // subjet values for Tree
            //------------------------------------

            double gensubjetpt = 0;
            if (!iEvent.isRealData()){
              const reco::GenJet* genJet = isub.genJet();
              if (genJet) gensubjetpt = genJet->pt();
            }

            if (count_SD==0){
              sub0_P4_uncorr            .SetPtEtaPhiM( subjetPt, subjetEta, subjetPhi, subjetMass);
              sub0_P4_L23res            .SetPtEtaPhiM( corrSubjetL23res          .pt() , corrSubjetL23res          .eta()  , corrSubjetL23res          .phi()  , corrSubjetL23res          .mass()  );
              // sub0_P4_L23resCorrUp      .SetPtEtaPhiM( corrSubjetL23resCorrUp    .pt() , corrSubjetL23resCorrUp    .eta()  , corrSubjetL23resCorrUp    .phi()  , corrSubjetL23resCorrUp    .mass()  );
              // sub0_P4_L23resCorrDn      .SetPtEtaPhiM( corrSubjetL23resCorrDn    .pt() , corrSubjetL23resCorrDn    .eta()  , corrSubjetL23resCorrDn    .phi()  , corrSubjetL23resCorrDn    .mass()  );
              // sub0_P4_L23resPtSmear     .SetPtEtaPhiM( corrSubjetL23resPtSmear   .pt() , corrSubjetL23resPtSmear   .eta()  , corrSubjetL23resPtSmear   .phi()  , corrSubjetL23resPtSmear   .mass()  );
              // sub0_P4_L23resPtSmearUp   .SetPtEtaPhiM( corrSubjetL23resPtSmearUp .pt() , corrSubjetL23resPtSmearUp .eta()  , corrSubjetL23resPtSmearUp .phi()  , corrSubjetL23resPtSmearUp .mass()  );
              // sub0_P4_L23resPtSmearDn   .SetPtEtaPhiM( corrSubjetL23resPtSmearDn .pt() , corrSubjetL23resPtSmearDn .eta()  , corrSubjetL23resPtSmearDn .phi()  , corrSubjetL23resPtSmearDn .mass()  );
              sub0_P4_L123res           .SetPtEtaPhiM( corrSubjetL123res         .pt() , corrSubjetL123res         .eta()  , corrSubjetL123res         .phi()  , corrSubjetL123res         .mass()  );
              // sub0_P4_L123resCorrUp     .SetPtEtaPhiM( corrSubjetL123resCorrUp   .pt() , corrSubjetL123resCorrUp   .eta()  , corrSubjetL123resCorrUp   .phi()  , corrSubjetL123resCorrUp   .mass()  );
              // sub0_P4_L123resCorrDn     .SetPtEtaPhiM( corrSubjetL123resCorrDn   .pt() , corrSubjetL123resCorrDn   .eta()  , corrSubjetL123resCorrDn   .phi()  , corrSubjetL123resCorrDn   .mass()  );
   
              sub0_area          = isub.jetArea() ;
              sub0_flav_parton   = isub.partonFlavour();
              sub0_flav_hadron   = isub.hadronFlavour();
              sub0_bdisc         = subjetBdisc;
              // available from toolbox only (80X)
              sub0_tau1          = isub.userFloat("NsubjettinessAK8PFCHSSoftDropSubjets:tau1");
              sub0_tau2          = isub.userFloat("NsubjettinessAK8PFCHSSoftDropSubjets:tau2");
              sub0_tau3          = isub.userFloat("NsubjettinessAK8PFCHSSoftDropSubjets:tau3");
              sub0_genpt         = gensubjetpt;
            }
            if (count_SD==1) {
              sub1_P4_uncorr            .SetPtEtaPhiM( subjetPt, subjetEta, subjetPhi, subjetMass);
              sub1_P4_L23res            .SetPtEtaPhiM( corrSubjetL23res          .pt() , corrSubjetL23res          .eta()  , corrSubjetL23res          .phi()  , corrSubjetL23res          .mass()  );
              // sub1_P4_L23resCorrUp      .SetPtEtaPhiM( corrSubjetL23resCorrUp    .pt() , corrSubjetL23resCorrUp    .eta()  , corrSubjetL23resCorrUp    .phi()  , corrSubjetL23resCorrUp    .mass()  );
              // sub1_P4_L23resCorrDn      .SetPtEtaPhiM( corrSubjetL23resCorrDn    .pt() , corrSubjetL23resCorrDn    .eta()  , corrSubjetL23resCorrDn    .phi()  , corrSubjetL23resCorrDn    .mass()  );
              // sub1_P4_L23resPtSmear     .SetPtEtaPhiM( corrSubjetL23resPtSmear   .pt() , corrSubjetL23resPtSmear   .eta()  , corrSubjetL23resPtSmear   .phi()  , corrSubjetL23resPtSmear   .mass()  );
              // sub1_P4_L23resPtSmearUp   .SetPtEtaPhiM( corrSubjetL23resPtSmearUp .pt() , corrSubjetL23resPtSmearUp .eta()  , corrSubjetL23resPtSmearUp .phi()  , corrSubjetL23resPtSmearUp .mass()  );
              // sub1_P4_L23resPtSmearDn   .SetPtEtaPhiM( corrSubjetL23resPtSmearDn .pt() , corrSubjetL23resPtSmearDn .eta()  , corrSubjetL23resPtSmearDn .phi()  , corrSubjetL23resPtSmearDn .mass()  );
              sub1_P4_L123res           .SetPtEtaPhiM( corrSubjetL123res         .pt() , corrSubjetL123res         .eta()  , corrSubjetL123res         .phi()  , corrSubjetL123res         .mass()  );
              // sub1_P4_L123resCorrUp     .SetPtEtaPhiM( corrSubjetL123resCorrUp   .pt() , corrSubjetL123resCorrUp   .eta()  , corrSubjetL123resCorrUp   .phi()  , corrSubjetL123resCorrUp   .mass()  );
              // sub1_P4_L123resCorrDn     .SetPtEtaPhiM( corrSubjetL123resCorrDn   .pt() , corrSubjetL123resCorrDn   .eta()  , corrSubjetL123resCorrDn   .phi()  , corrSubjetL123resCorrDn   .mass()  );
   
              sub1_area          = isub.jetArea() ;
              sub1_flav_parton   = isub.partonFlavour();
              sub1_flav_hadron   = isub.hadronFlavour();
              sub1_bdisc         = subjetBdisc;
              // available from toolbox only (80X)
              sub1_tau1          = isub.userFloat("NsubjettinessAK8PFCHSSoftDropSubjets:tau1");
              sub1_tau2          = isub.userFloat("NsubjettinessAK8PFCHSSoftDropSubjets:tau2");
              sub1_tau3          = isub.userFloat("NsubjettinessAK8PFCHSSoftDropSubjets:tau3");
              sub1_genpt         = gensubjetpt;
            }
            if (subjetMass > mostMassiveSDsubjetMass) mostMassiveSDsubjetMass = subjetMass;

            if (verbose_) cout<<"        -> SD Subjet pt "<<subjetPt<<" Eta "<<subjetEta<<" deltaRsubjetJet "<<deltaRsubjetJet<<" Mass "<<subjetMass<<" corrMass "<<corrSubjetL23res.mass() <<" Bdisc "<<subjetBdisc<<endl;
            if (verbose_) cout<<"        ->    sub0_tau1 "<<sub0_tau1<<" sub0_tau2 "<<sub0_tau2<<" sub0_tau3 "<<sub0_tau3<<endl;
            count_SD++;

          }
        }
        count_all_subjets++;
      }
      if (verbose_) cout<<"     count_matched_subjets "<<count_matched_subjets<<endl;
      if (count_matched_subjets!=2 && verbose_ ) cout<<"     CHECKME"<<endl;
      if (verbose_) cout<<"     closest_DR "<<closest_DR<<endl;;
      if (verbose_) cout<<"     second_closest_DR "<<second_closest_DR<<endl;;
      if (verbose_) cout<<"     closest_i "<<closest_i<<endl;
      if (verbose_) cout<<"     second_closest_i "<<second_closest_i<<endl;

      
    }

    TLorentzVector sumSDsubjets_P4_uncorr           ;
    TLorentzVector sumSDsubjets_P4_L23res           ;
    // TLorentzVector sumSDsubjets_P4_L23resCorrUp     ;
    // TLorentzVector sumSDsubjets_P4_L23resCorrDn     ;
    // TLorentzVector sumSDsubjets_P4_L23resPtSmear    ;
    // TLorentzVector sumSDsubjets_P4_L23resPtSmearUp  ;
    // TLorentzVector sumSDsubjets_P4_L23resPtSmearDn  ;
    TLorentzVector sumSDsubjets_P4_L123res          ;
    // TLorentzVector sumSDsubjets_P4_L123resCorrDn    ;
    // TLorentzVector sumSDsubjets_P4_L123resCorrUp    ;

    if (count_SD>1){ 
      sumSDsubjets_P4_uncorr             = sub0_P4_uncorr              + sub1_P4_uncorr            ; 
      sumSDsubjets_P4_L23res             = sub0_P4_L23res              + sub1_P4_L23res            ; 
      // sumSDsubjets_P4_L23resCorrUp       = sub0_P4_L23resCorrUp        + sub1_P4_L23resCorrUp      ; 
      // sumSDsubjets_P4_L23resCorrDn       = sub0_P4_L23resCorrDn        + sub1_P4_L23resCorrDn      ; 
      // sumSDsubjets_P4_L23resPtSmear      = sub0_P4_L23resPtSmear       + sub1_P4_L23resPtSmear     ;
      // sumSDsubjets_P4_L23resPtSmearUp    = sub0_P4_L23resPtSmearUp     + sub1_P4_L23resPtSmearUp   ;
      // sumSDsubjets_P4_L23resPtSmearDn    = sub0_P4_L23resPtSmearDn     + sub1_P4_L23resPtSmearDn   ;
      sumSDsubjets_P4_L123res            = sub0_P4_L123res             + sub1_P4_L123res           ; 
      // sumSDsubjets_P4_L123resCorrUp      = sub0_P4_L123resCorrUp       + sub1_P4_L123resCorrUp     ; 
      // sumSDsubjets_P4_L123resCorrDn      = sub0_P4_L123resCorrDn       + sub1_P4_L123resCorrDn     ; 
    }  

    double maxbdisc = 0 ;
    double maxbdiscflav_hadron = 0 ;
    double maxbdiscflav_parton = 0 ;
    if (sub0_bdisc>=sub1_bdisc){
      maxbdisc = sub0_bdisc;
      maxbdiscflav_hadron = sub0_flav_hadron;
      maxbdiscflav_parton = sub0_flav_parton;
    } 
    else if (sub1_bdisc>sub0_bdisc){
      maxbdisc = sub1_bdisc;
      maxbdiscflav_hadron = sub1_flav_hadron;
      maxbdiscflav_parton = sub1_flav_parton;
    }  

    //------------------------------------
    // PUPPI SoftDrop subjets
    //------------------------------------ 
    TLorentzVector pup0_P4_uncorr           ;
    TLorentzVector pup0_P4_L23res           ;
    // TLorentzVector pup0_P4_L23resCorrUp     ;
    // TLorentzVector pup0_P4_L23resCorrDn     ;
    // TLorentzVector pup0_P4_L23resPtSmear    ;
    // TLorentzVector pup0_P4_L23resPtSmearUp  ;
    // TLorentzVector pup0_P4_L23resPtSmearDn  ;

    TLorentzVector pup1_P4_uncorr           ;
    TLorentzVector pup1_P4_L23res           ;
    // TLorentzVector pup1_P4_L23resCorrUp     ;
    // TLorentzVector pup1_P4_L23resCorrDn     ;
    // TLorentzVector pup1_P4_L23resPtSmear    ;
    // TLorentzVector pup1_P4_L23resPtSmearUp  ;
    // TLorentzVector pup1_P4_L23resPtSmearDn  ;


    double pup0_area  = 0;
    double pup0_tau1  = 0;
    double pup0_tau2  = 0;
    double pup0_tau3  = 0;
    double pup0_flav_hadron  = 0;
    double pup0_flav_parton  = 0;
    double pup0_bdisc = 0;
    double pup0_genpt = 0;
    double pup1_area  = 0;
    double pup1_tau1  = 0;
    double pup1_tau2  = 0;
    double pup1_tau3  = 0;
    double pup1_flav_hadron  = 0;
    double pup1_flav_parton  = 0;
    double pup1_bdisc = 0;
    double pup1_genpt = 0;
    double mostMassiveSDPUPPIsubjetMass = 0;
    int count_pup=0;

    if (!useToolbox_){
      auto const & sdSubjetsPuppi = ijet.subjets("SoftDropPuppi");
      for ( auto const & it : sdSubjetsPuppi ) {
        double subjetPt       = it->correctedP4(0).pt();
        double subjetEta      = it->correctedP4(0).eta();
        double subjetPhi      = it->correctedP4(0).phi();
        double subjetMass     = it->correctedP4(0).mass();
        double subjetBdisc    = it->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"); 
        double deltaRsubjetJet = deltaR(ijet.eta(), ijet.phi(), subjetEta, subjetPhi);
        if (verbose_) cout<<" SD Subjet pt "<<subjetPt<<" Eta "<<subjetEta<<" deltaRsubjetJet "<<deltaRsubjetJet<<" Mass "<<subjetMass<<" Bdisc "<<subjetBdisc<<endl; 
        
        //------------------------------------
        // PUPPI subjet JEC 
        //------------------------------------
        reco::Candidate::LorentzVector uncorrSubjet = it->correctedP4(0);
        JetCorrectorAK4pup -> setJetEta( uncorrSubjet.eta()    );
        JetCorrectorAK4pup -> setJetPt ( uncorrSubjet.pt()     );
        JetCorrectorAK4pup -> setJetE  ( uncorrSubjet.energy() );
        JetCorrectorAK4pup -> setJetA  ( it->jetArea()         );
        JetCorrectorAK4pup -> setRho   ( rho                   );
        JetCorrectorAK4pup -> setNPV   ( nvtx                  );
        double subjet_corr_factor_L23res_full = JetCorrectorAK4pup->getCorrection();
        reco::Candidate::LorentzVector corrSubjetL23res = subjet_corr_factor_L23res_full * uncorrSubjet;

        // //------------------------------------
        // // PUPPI subjet JEC uncertainty
        // //------------------------------------
        // double subjet_corrDn_L23 =   1.0;
        // JetCorrUncertAK4pup->setJetPhi(  corrSubjetL23res.phi()  );
        // JetCorrUncertAK4pup->setJetEta(  corrSubjetL23res.eta()  );
        // JetCorrUncertAK4pup->setJetPt(   corrSubjetL23res.pt()   );
        // subjet_corrDn_L23   = subjet_corr_factor_L23res_full - JetCorrUncertAK4pup->getUncertainty(0);
        // double subjet_corrUp_L23 =   1.0;
        // JetCorrUncertAK4pup->setJetPhi(  corrSubjetL23res.phi()  );
        // JetCorrUncertAK4pup->setJetEta(  corrSubjetL23res.eta()  );
        // JetCorrUncertAK4pup->setJetPt(   corrSubjetL23res.pt()   );
        // subjet_corrUp_L23   = subjet_corr_factor_L23res_full + JetCorrUncertAK4pup->getUncertainty(1);

        // reco::Candidate::LorentzVector corrSubjetL23resCorrDn   = subjet_corrDn_L23   * uncorrSubjet;
        // reco::Candidate::LorentzVector corrSubjetL23resCorrUp   = subjet_corrUp_L23   * uncorrSubjet;
     
        //------------------------------------
        // subjet values for Tree
        //------------------------------------

        if (count_pup==0){
          pup0_P4_uncorr           .SetPtEtaPhiM( subjetPt, subjetEta, subjetPhi, subjetMass);
          pup0_P4_L23res           .SetPtEtaPhiM( corrSubjetL23res.pt()    , corrSubjetL23res.eta()    , corrSubjetL23res.phi()    , corrSubjetL23res.mass()     );
          // pup0_P4_L23resCorrUp    .SetPtEtaPhiM( corrSubjetL23resCorrUp.pt()  , corrSubjetL23resCorrUp.eta()  , corrSubjetL23resCorrUp.phi()  , corrSubjetL23resCorrUp.mass()   );
          // pup0_P4_L23resCorrDn    .SetPtEtaPhiM( corrSubjetL23resCorrDn.pt()  , corrSubjetL23resCorrDn.eta()  , corrSubjetL23resCorrDn.phi()  , corrSubjetL23res.mass()     );
          pup0_area   = it->jetArea() ;
          if (useToolbox_){
            pup0_tau1   = it->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1");
            pup0_tau2   = it->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2");
            pup0_tau3   = it->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau3");
          }
          pup0_flav_parton   = it->partonFlavour();
          pup0_flav_hadron   = it->hadronFlavour();
          pup0_bdisc  = subjetBdisc;
        }
        if (count_pup==1) {
          pup1_P4_uncorr            .SetPtEtaPhiM( subjetPt, subjetEta, subjetPhi, subjetMass);
          pup1_P4_L23res           .SetPtEtaPhiM( corrSubjetL23res.pt()    , corrSubjetL23res.eta()    , corrSubjetL23res.phi()    , corrSubjetL23res.mass()     );
          // pup1_P4_L23resCorrUp    .SetPtEtaPhiM( corrSubjetL23resCorrUp.pt()  , corrSubjetL23resCorrUp.eta()  , corrSubjetL23resCorrUp.phi()  , corrSubjetL23resCorrUp.mass()   );
          // pup1_P4_L23resCorrDn    .SetPtEtaPhiM( corrSubjetL23resCorrDn.pt()  , corrSubjetL23resCorrDn.eta()  , corrSubjetL23resCorrDn.phi()  , corrSubjetL23res.mass()     );
          pup1_area   = it->jetArea() ;
          if (useToolbox_){
            pup1_tau1   = it->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1");
            pup1_tau2   = it->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2");
            pup1_tau3   = it->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau3");
          }
          pup1_flav_parton   = it->partonFlavour();
          pup1_flav_hadron   = it->hadronFlavour();
          pup1_bdisc  = subjetBdisc;
        }

        if (subjetMass > mostMassiveSDPUPPIsubjetMass) mostMassiveSDPUPPIsubjetMass = subjetMass;
        count_pup++;
      }
    }
    if (useToolbox_){
      if (verbose_) cout<<"   Toolbox AK8 jets. Find puppi softdrop subjets "<<endl;

 
      int count_all_subjets =0;
      int count_matched_subjets =0;
      double closest_DR = 99;
      double closest_i = -1;
      double second_closest_DR = 99;
      double second_closest_i  = -1;

      // Loop once to find the subjets closest to the AK8 jet
      for (const pat::Jet &isub : *AK8PUPPIsub) {  
  
        double subjetPt       = isub.correctedP4(0).pt();
        double subjetEta      = isub.correctedP4(0).eta();
        double subjetPhi      = isub.correctedP4(0).phi();
        double subjetMass     = isub.correctedP4(0).mass();
        // double subjetBdisc    = isub.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"); 

        double deltaRsubjetJet = deltaR(puppi_eta, puppi_phi, subjetEta, subjetPhi);

        if (verbose_) cout<<"     PupSubjet "<<count_all_subjets<<"   "<<subjetMass<<" "<<subjetPt<<" "<<deltaRsubjetJet<<endl;

        if (deltaRsubjetJet<closest_DR){
          second_closest_DR = closest_DR;
          closest_DR        = deltaRsubjetJet;
          second_closest_i  = closest_i;
          closest_i         = count_all_subjets;
        }
        else if (deltaRsubjetJet<second_closest_DR){
          second_closest_DR = deltaRsubjetJet ;
          second_closest_i  = count_all_subjets;
        }
        count_all_subjets++;
      }

      // Loop a second time. If one of the two closest subjets matches the dR requirement save its infromation. Subjet 0 = hardest.
      count_all_subjets =0;
      for (const pat::Jet &isub : *AK8PUPPIsub) {  
        
  
        double subjetPt       = isub.correctedP4(0).pt();
        double subjetEta      = isub.correctedP4(0).eta();
        double subjetPhi      = isub.correctedP4(0).phi();
        double subjetMass     = isub.correctedP4(0).mass();
        double subjetBdisc    = isub.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"); 

        double deltaRsubjetJet = deltaR(puppi_eta, puppi_phi, subjetEta, subjetPhi);

        if (verbose_) cout<<"     PupSubjet "<<count_all_subjets<<"   "<<subjetMass<<" "<<subjetPt<<" "<<deltaRsubjetJet<<endl;

        if ( count_all_subjets==closest_i || count_all_subjets==second_closest_i ){
          if (verbose_) cout<<"      -> one of two closest "<<endl;
          if (deltaRsubjetJet<0.8){
            nsubjets_pup++;
            if (verbose_) cout<<"      -> dR matched subjet with mass "<< subjetMass<<endl;

            count_matched_subjets++;

            //------------------------------------
            // PUPPI subjet JEC 
            //------------------------------------
            reco::Candidate::LorentzVector uncorrSubjet = isub.correctedP4(0);
            JetCorrectorAK4pup -> setJetEta( uncorrSubjet.eta()    );
            JetCorrectorAK4pup -> setJetPt ( uncorrSubjet.pt()     );
            JetCorrectorAK4pup -> setJetE  ( uncorrSubjet.energy() );
            JetCorrectorAK4pup -> setJetA  ( isub.jetArea()         );
            JetCorrectorAK4pup -> setRho   ( rho                   );
            JetCorrectorAK4pup -> setNPV   ( nvtx                  );
            double subjet_corr_factor_L23res_full = JetCorrectorAK4pup->getCorrection();
            reco::Candidate::LorentzVector corrSubjetL23res = subjet_corr_factor_L23res_full * uncorrSubjet;

            // //------------------------------------
            // // PUPPI subjet JEC uncertainty
            // //------------------------------------
            // double subjet_corrDn_L23 =   1.0;
            // JetCorrUncertAK4pup->setJetPhi(  corrSubjetL23res.phi()  );
            // JetCorrUncertAK4pup->setJetEta(  corrSubjetL23res.eta()  );
            // JetCorrUncertAK4pup->setJetPt(   corrSubjetL23res.pt()   );
            // subjet_corrDn_L23   = subjet_corr_factor_L23res_full - JetCorrUncertAK4pup->getUncertainty(0);
            // double subjet_corrUp_L23 =   1.0;
            // JetCorrUncertAK4pup->setJetPhi(  corrSubjetL23res.phi()  );
            // JetCorrUncertAK4pup->setJetEta(  corrSubjetL23res.eta()  );
            // JetCorrUncertAK4pup->setJetPt(   corrSubjetL23res.pt()   );
            // subjet_corrUp_L23   = subjet_corr_factor_L23res_full + JetCorrUncertAK4pup->getUncertainty(1);

            // reco::Candidate::LorentzVector corrSubjetL23resCorrDn   = subjet_corrDn_L23   * uncorrSubjet;
            // reco::Candidate::LorentzVector corrSubjetL23resCorrUp   = subjet_corrUp_L23   * uncorrSubjet;
         
            // //------------------------------------
            // // subjet JER SF
            // //------------------------------------

            // Doesn't seem to work. Gen subjet from isub.genJet()  is not a good match (especially true for 2nd subjet). If i have the time try just smearing the jets with no gen info.

            // TLorentzVector GenSubJet;
            // double ptsmear   = 1;
            // double ptsmearUp = 1;
            // double ptsmearDn = 1;
            // if (!iEvent.isRealData()){
            //   const reco::GenJet* genJet = isub.genJet();
            //   if (genJet) {
            //     GenSubJet.SetPtEtaPhiM( genJet->pt(), genJet->eta(), genJet->phi(), genJet->mass() );
            //     if (verbose_) cout<<"  SD subjet genJet pt "<<genJet->pt()<<" mass "<<genJet->mass()<<" reco pt "<<subjetPt<<" reco mass "<<subjetMass<<endl;
            //   }
            //   double jer_sf    = jer_scaler.getScaleFactor({{JME::Binning::JetEta, corrSubjetL23res.eta()}});
            //   double jer_sf_up = jer_scaler.getScaleFactor({{JME::Binning::JetEta, corrSubjetL23res.eta()}}, Variation::UP);
            //   double jer_sf_dn = jer_scaler.getScaleFactor({{JME::Binning::JetEta, corrSubjetL23res.eta()}}, Variation::DOWN);
            //   if (verbose_) std::cout << " SD subjet JER Scale factors (Nominal / Up / Down) : " << jer_sf << " / " << jer_sf_up << " / " << jer_sf_dn << std::endl;
            //   double recopt    = corrSubjetL23res.pt();
            //   double genpt     = GenJetMatched.Perp();
            //   double deltapt   = (recopt-genpt)*(jer_sf-1.0);
            //   double deltaptUp = (recopt-genpt)*(jer_sf_up-1.0);
            //   double deltaptDn = (recopt-genpt)*(jer_sf_dn-1.0);

            //   cout<<"recopt     "<<recopt    <<endl; 
            //   cout<<"genpt      "<<genpt     <<endl; 
            //   cout<<"deltapt    "<<deltapt   <<endl; 

            //   cout<<"recopt + deltapt     "<<recopt +deltapt   <<endl; 
            //   cout<<"(recopt+deltapt)/recopt     "<<(recopt+deltapt)/recopt   <<endl; 

            //   ptsmear   = std::max((double)0.0, (recopt+deltapt)/recopt     );
            //   ptsmearUp = std::max((double)0.0, (recopt+deltaptUp)/recopt   );
            //   ptsmearDn = std::max((double)0.0, (recopt+deltaptDn)/recopt   );
            //   if (verbose_); std::cout<<" SD subjet ptsmear "<<ptsmear<<" ptsmearUp "<<ptsmearUp<<" ptsmearDn "<<ptsmearDn<<endl;
            // }
            // reco::Candidate::LorentzVector corrSubjetL23resPtSmear   = ptsmear   * corrSubjetL23res ;
            // reco::Candidate::LorentzVector corrSubjetL23resPtSmearUp = ptsmearUp * corrSubjetL23res ;
            // reco::Candidate::LorentzVector corrSubjetL23resPtSmearDn = ptsmearDn * corrSubjetL23res ;

            //------------------------------------
            // subjet values for Tree
            //------------------------------------

            double gensubjetpt = 0;
            if (!iEvent.isRealData()){
              const reco::GenJet* genJet = isub.genJet();
              if (genJet) gensubjetpt = genJet->pt();
            }

            if (count_pup==0){
              pup0_P4_uncorr            .SetPtEtaPhiM( subjetPt, subjetEta, subjetPhi, subjetMass);
              pup0_P4_L23res            .SetPtEtaPhiM( corrSubjetL23res          .pt() , corrSubjetL23res          .eta() , corrSubjetL23res          .phi() , corrSubjetL23res          .mass() );
              // pup0_P4_L23resCorrUp      .SetPtEtaPhiM( corrSubjetL23resCorrUp    .pt() , corrSubjetL23resCorrUp    .eta() , corrSubjetL23resCorrUp    .phi() , corrSubjetL23resCorrUp    .mass() );
              // pup0_P4_L23resCorrDn      .SetPtEtaPhiM( corrSubjetL23resCorrDn    .pt() , corrSubjetL23resCorrDn    .eta() , corrSubjetL23resCorrDn    .phi() , corrSubjetL23resCorrDn    .mass() );
              // pup0_P4_L23resPtSmear     .SetPtEtaPhiM( corrSubjetL23resPtSmear   .pt() , corrSubjetL23resPtSmear   .eta() , corrSubjetL23resPtSmear   .phi() , corrSubjetL23resPtSmear   .mass() );
              // pup0_P4_L23resPtSmearUp   .SetPtEtaPhiM( corrSubjetL23resPtSmearUp .pt() , corrSubjetL23resPtSmearUp .eta() , corrSubjetL23resPtSmearUp .phi() , corrSubjetL23resPtSmearUp .mass() );
              // pup0_P4_L23resPtSmearDn   .SetPtEtaPhiM( corrSubjetL23resPtSmearDn .pt() , corrSubjetL23resPtSmearDn .eta() , corrSubjetL23resPtSmearDn .phi() , corrSubjetL23resPtSmearDn .mass() );

              pup0_tau1   = isub.userFloat("NsubjettinessAK8PFPuppiSoftDropSubjets:tau1");
              pup0_tau2   = isub.userFloat("NsubjettinessAK8PFPuppiSoftDropSubjets:tau2");
              pup0_tau3   = isub.userFloat("NsubjettinessAK8PFPuppiSoftDropSubjets:tau3");
            
              pup0_flav_parton   = isub.partonFlavour();
              pup0_flav_hadron   = isub.hadronFlavour();
              pup0_area          = isub.jetArea() ;
              pup0_bdisc         = subjetBdisc;
              pup0_genpt         = gensubjetpt;
            }
            if (count_pup==1) {
              pup1_P4_uncorr            .SetPtEtaPhiM( subjetPt, subjetEta, subjetPhi, subjetMass);
              pup1_P4_L23res            .SetPtEtaPhiM( corrSubjetL23res          .pt() , corrSubjetL23res          .eta() , corrSubjetL23res          .phi() , corrSubjetL23res          .mass() );
              // pup1_P4_L23resCorrUp      .SetPtEtaPhiM( corrSubjetL23resCorrUp    .pt() , corrSubjetL23resCorrUp    .eta() , corrSubjetL23resCorrUp    .phi() , corrSubjetL23resCorrUp    .mass() );
              // pup1_P4_L23resCorrDn      .SetPtEtaPhiM( corrSubjetL23resCorrDn    .pt() , corrSubjetL23resCorrDn    .eta() , corrSubjetL23resCorrDn    .phi() , corrSubjetL23resCorrDn    .mass() );
              // pup1_P4_L23resPtSmear     .SetPtEtaPhiM( corrSubjetL23resPtSmear   .pt() , corrSubjetL23resPtSmear   .eta() , corrSubjetL23resPtSmear   .phi() , corrSubjetL23resPtSmear   .mass() );
              // pup1_P4_L23resPtSmearUp   .SetPtEtaPhiM( corrSubjetL23resPtSmearUp .pt() , corrSubjetL23resPtSmearUp .eta() , corrSubjetL23resPtSmearUp .phi() , corrSubjetL23resPtSmearUp .mass() );
              // pup1_P4_L23resPtSmearDn   .SetPtEtaPhiM( corrSubjetL23resPtSmearDn .pt() , corrSubjetL23resPtSmearDn .eta() , corrSubjetL23resPtSmearDn .phi() , corrSubjetL23resPtSmearDn .mass() );

              pup1_tau1   = isub.userFloat("NsubjettinessAK8PFPuppiSoftDropSubjets:tau1");
              pup1_tau2   = isub.userFloat("NsubjettinessAK8PFPuppiSoftDropSubjets:tau2");
              pup1_tau3   = isub.userFloat("NsubjettinessAK8PFPuppiSoftDropSubjets:tau3");
            
              pup1_flav_parton   = isub.partonFlavour();
              pup1_flav_hadron   = isub.hadronFlavour();
              pup1_area          = isub.jetArea() ;
              pup1_bdisc         = subjetBdisc;
              pup1_genpt         = gensubjetpt;
            }

            if (subjetMass > mostMassiveSDPUPPIsubjetMass) mostMassiveSDPUPPIsubjetMass = subjetMass;
            count_pup++;
            if (verbose_) cout<<"       -> SD Subjet pt "<<subjetPt<<" Eta "<<subjetEta<<" deltaRsubjetJet "<<deltaRsubjetJet<<" Mass "<<subjetMass<<" Bdisc "<<subjetBdisc<<endl; 

          }
        }
        count_all_subjets++;
      }
      if (verbose_) cout<<"     Puppi count_matched_subjets "<<count_matched_subjets<<endl;
      if (count_matched_subjets!=2 && verbose_ ) cout<<"     CHECKME"<<endl;
      if (verbose_) cout<<"     closest_DR "<<closest_DR<<endl;;
      if (verbose_) cout<<"     second_closest_DR "<<second_closest_DR<<endl;;
      if (verbose_) cout<<"     closest_i "<<closest_i<<endl;
      if (verbose_) cout<<"     second_closest_i "<<second_closest_i<<endl;
    }


    TLorentzVector sumPUPsubjets_P4_uncorr          ;
    TLorentzVector sumPUPsubjets_P4_L23res          ;
    // TLorentzVector sumPUPsubjets_P4_L23resCorrUp    ;
    // TLorentzVector sumPUPsubjets_P4_L23resCorrDn    ;
    // TLorentzVector sumPUPsubjets_P4_L23resPtSmear   ;
    // TLorentzVector sumPUPsubjets_P4_L23resPtSmearUp ;
    // TLorentzVector sumPUPsubjets_P4_L23resPtSmearDn ;
    if (count_SD>1){ 
      sumPUPsubjets_P4_uncorr            = pup0_P4_uncorr            + pup1_P4_uncorr            ; 
      sumPUPsubjets_P4_L23res            = pup0_P4_L23res            + pup1_P4_L23res            ; 
      // sumPUPsubjets_P4_L23resCorrUp      = pup0_P4_L23resCorrUp      + pup1_P4_L23resCorrUp      ; 
      // sumPUPsubjets_P4_L23resCorrDn      = pup0_P4_L23resCorrDn      + pup1_P4_L23resCorrDn      ; 
      // sumPUPsubjets_P4_L23resPtSmear     = pup0_P4_L23resPtSmear     + pup1_P4_L23resPtSmear     ; 
      // sumPUPsubjets_P4_L23resPtSmearUp   = pup0_P4_L23resPtSmearUp   + pup1_P4_L23resPtSmearUp   ; 
      // sumPUPsubjets_P4_L23resPtSmearDn   = pup0_P4_L23resPtSmearDn   + pup1_P4_L23resPtSmearDn   ; 
    } 


    double pup_maxbdisc = 0 ;
    double pup_maxbdiscflav_hadron = 0 ;
    double pup_maxbdiscflav_parton = 0 ;
    if (pup0_bdisc>=pup1_bdisc){
      pup_maxbdisc = pup0_bdisc;
      pup_maxbdiscflav_hadron = pup0_flav_hadron;
      pup_maxbdiscflav_parton = pup0_flav_parton;
    } 
    else if (pup1_bdisc>pup0_bdisc){
      pup_maxbdisc = pup1_bdisc;
      pup_maxbdiscflav_hadron = pup1_flav_hadron;
      pup_maxbdiscflav_parton = pup1_flav_parton;
    }  

    //------------------------------------
    // Gen particle info
    //------------------------------------ 
    double deltaR_jet_t1 = 99;
    double deltaR_jet_t2 = 99;
    double deltaR_jet_p1 = 99;
    double deltaR_jet_p2 = 99;

    bool jet_matched_t1 = false;
    //bool jet_matched_t2 = false;
    bool jet_matched_p1 = false;
    //bool jet_matched_p2 = false;

    if (!iEvent.isRealData() and runGenLoop_) {
      deltaR_jet_t1 = jet_p4.DeltaR(t_p4  );
      deltaR_jet_p1 = jet_p4.DeltaR(hardest_parton_hardScatterOutgoing_p4        );
      deltaR_jet_p2 = jet_p4.DeltaR(second_hardest_parton_hardScatterOutgoing_p4 );
      if (deltaR_jet_t1<deltaR_jet_t2) jet_matched_t1 = true;
      //if (deltaR_jet_t2<deltaR_jet_t1) jet_matched_t2 = true;
      if (deltaR_jet_p1<deltaR_jet_p2) jet_matched_p1 = true;
     // if (deltaR_jet_p2<deltaR_jet_p1) jet_matched_p2 = true;
    } 

    //------------------------------------
    // Tree variables
    //------------------------------------ 

  
    if (runHadTree_){
 
      // second jet oposite lepton

      // AK8 jet should be in opposite hemisphere from lepton. If leading jet matches then use it. If it doensn't then check the second leading jet.
      if ( count_AK8CHS==0 && count_fill_hadTree==0 ){
        if (verbose_) cout<<"Put jet variables in had tree  -> count_AK8CHS "<<count_AK8CHS<<" count_fill_hadTree"<<count_fill_hadTree<<endl;
        count_hadAK8CHS++;

        count_fill_hadTree++;
        AK8jet_had_P4corr.SetPtEtaPhiM( corrJet.pt(), corrJet.eta(), corrJet.phi(), corrJet.mass() );

        // basic kinematic and ID variables
        JetPtRaw                              = uncorrJet.pt()      ;                 
        JetEtaRaw                             = uncorrJet.eta()     ;                  
        JetPhiRaw                             = uncorrJet.phi()     ;   
        JetMassRaw                            = uncorrJet.mass()    ;   
        //std::cout << "not yet passed jet 1 " << JetPtRaw << " " << JetEtaRaw << " " << JetPhiRaw << " " << JetMassRaw << std::endl;                                        
        // JetP                                  = corrJet.P()         ;        
        // JetPt                                 = corrJet.pt()        ;                  
        // JetEta                                = corrJet.eta()       ;                  
        // JetPhi                                = corrJet.phi()       ;                  
        // JetRap                                = corrJet.Rapidity()  ;                  
        // JetEnergy                             = corrJet.energy()    ;                  
        // JetMass                               = corrJet.mass()      ;                    
        JetArea                               = ijet.jetArea()      ;                  
       
        JetCHF                                = ijet.chargedHadronEnergy() / uncorrJet.E()  ;                        
        JetNHF                                = ijet.neutralHadronEnergy() / uncorrJet.E()  ;                         
        JetCM                                 = ijet.chargedMultiplicity()  ;                         
        JetNM                                 = ijet.neutralMultiplicity()  ;                          
        JetNEF                                = ijet.neutralEmEnergy() / uncorrJet.E()  ;                            
        JetCEF                                = ijet.chargedEmEnergy() / uncorrJet.E()  ;                          
        JetMF                                 = ijet.muonEnergy() / uncorrJet.E()  ;                         
        JetMult                               = ijet.numberOfDaughters() ;   

        // soft drop mass calculated from soft drop subjets                          
        JetSDmassRaw                          = sumSDsubjets_P4_uncorr   .M()    ;  
        JetSDetaRaw                           = sumSDsubjets_P4_uncorr   .Eta()  ;                    
        JetSDphiRaw                           = sumSDsubjets_P4_uncorr   .Phi()  ;  
        JetSDptRaw                            = sumSDsubjets_P4_uncorr   .Perp() ;  
    
        // experiment with JEC applied separately to each subjet
        JetSDmassSubjetCorrL23                      = sumSDsubjets_P4_L23res          .M()    ;   
        // JetSDmassSubjetCorrL23Up                    = sumSDsubjets_P4_L23resSubjetCorrUp    .M()    ;   
        // JetSDmassSubjetCorrL23Dn                    = sumSDsubjets_P4_L23resSubjetCorrDn    .M()    ; 
        JetSDmassSubjetCorrL123                     = sumSDsubjets_P4_L123res         .M()    ;  
        // JetSDmassCorrL123Up                   = sumSDsubjets_P4_L123resCorrUp   .M()    ;   
        // JetSDmassCorrL123Dn                   = sumSDsubjets_P4_L123resCorrDn   .M()    ;  
        // JetSDmassCorrL23Smear                 = sumSDsubjets_P4_L23resPtSmear   .M()    ;   // This doesn't work. Subjet genjets are not a good match.
        // JetSDmassCorrL23SmearUp               = sumSDsubjets_P4_L23resPtSmearUp .M()    ;
        // JetSDmassCorrL23SmearDn               = sumSDsubjets_P4_L23resPtSmearDn .M()    ;
        // JetSDptCorrL23                        = sumSDsubjets_P4_L23res          .Perp() ;  
        // JetSDptCorrL23Up                      = sumSDsubjets_P4_L23resCorrUp    .Perp() ;  
        // JetSDptCorrL23Dn                      = sumSDsubjets_P4_L23resCorrDn    .Perp() ;  
        // JetSDptCorrL123                       = sumSDsubjets_P4_L123res         .Perp() ;  
        // JetSDptCorrL123Up                     = sumSDsubjets_P4_L123resCorrUp   .Perp() ;  
        // JetSDptCorrL123Dn                     = sumSDsubjets_P4_L123resCorrDn   .Perp() ;  
        // JetSDptCorrL23Smear                   = sumSDsubjets_P4_L23resPtSmear   .Perp() ;
        // JetSDptCorrL23SmearUp                 = sumSDsubjets_P4_L23resPtSmearUp .Perp() ;
        // JetSDptCorrL23SmearDn                 = sumSDsubjets_P4_L23resPtSmearDn .Perp() ;
            
        // user floats from the toolbox
        // JetSDmass                             = softDropMass  ;   // soft Drop mass from miniAOD                 
        JetMassPruned                         = prunedMass    ;     
        JetMassTrimmed                        = trimmedMass   ;     
        JetTau1                               = tau1          ;  
        JetTau2                               = tau2          ;  
        JetTau3                               = tau3          ;  
        JetTau4                               = tau4          ;  
        JetTau32                              = tau32         ;  
        JetTau21                              = tau21         ;  

        JetNsubjetsSD                         = nsubjets_chs   ;
        JetNsubjetsSDPuppi                    = nsubjets_pup   ;

        // Softdrop subjet variables
        JetSDsubjet0bdisc                     = sub0_bdisc            ;  
        JetSDsubjet1bdisc                     = sub1_bdisc            ;   
        JetSDmaxbdisc                         = maxbdisc              ;
        JetSDmaxbdiscflavHadron               = maxbdiscflav_hadron   ;  
        JetSDmaxbdiscflavParton               = maxbdiscflav_parton   ;  
        JetSDsubjet0pt                        = sub0_P4_uncorr.Pt()   ;               
        JetSDsubjet0mass                      = sub0_P4_uncorr.M()    ;  
        JetSDsubjet0eta                       = sub0_P4_uncorr.Eta()  ;  
        JetSDsubjet0phi                       = sub0_P4_uncorr.Phi()  ;  
        JetSDsubjet0area                      = sub0_area             ;  
        JetSDsubjet0flavHadron                = sub0_flav_hadron      ;  
        JetSDsubjet0flavParton                = sub0_flav_parton      ;  
        JetSDsubjet0matchedgenjetpt           = sub0_genpt            ;  
        JetSDsubjet0tau1                      = sub0_tau1             ;  
        JetSDsubjet0tau2                      = sub0_tau2             ;  
        JetSDsubjet0tau3                      = sub0_tau3             ;  
        JetSDsubjet1pt                        = sub1_P4_uncorr.Pt()   ;                    
        JetSDsubjet1mass                      = sub1_P4_uncorr.M()    ; 
        JetSDsubjet1eta                       = sub1_P4_uncorr.Eta()  ;  
        JetSDsubjet1phi                       = sub1_P4_uncorr.Phi()  ;                     
        JetSDsubjet1area                      = sub1_area             ;                    
        JetSDsubjet1flavHadron                = sub1_flav_hadron      ;     
        JetSDsubjet1flavParton                = sub1_flav_parton      ;
        JetSDsubjet1matchedgenjetpt           = sub1_genpt            ;       
        JetSDsubjet1tau1                      = sub1_tau1             ;  
        JetSDsubjet1tau2                      = sub1_tau2             ;  
        JetSDsubjet1tau3                      = sub1_tau3             ; 

        // Angle between puppi jet and chs jet
        // JetDeltaRPuppi                        = minDR_pup_chs;       

        // Puppi jet kinematics (uncorrected) and ID variables
        // JetPuppiP                             = AK8PUPPI_P4uncorr.P()    ;                  
        JetPuppiPtRaw                            = puppi_pt   ;                  
        JetPuppiEtaRaw                           = puppi_eta  ;                   
        JetPuppiPhiRaw                           = puppi_phi  ;                  
        JetPuppiMassRaw                          = puppi_mass ;                  
        JetPuppiArea                          = puppi_area ;                  

        JetPuppiCHF                           = puppi_CHF   ; 
        JetPuppiNHF                           = puppi_NHF   ; 
        JetPuppiCM                            = puppi_CM    ; 
        JetPuppiNM                            = puppi_NM    ; 
        JetPuppiNEF                           = puppi_NEF   ; 
        JetPuppiCEF                           = puppi_CEF   ; 
        JetPuppiMF                            = puppi_MF    ; 
        JetPuppiMult                          = puppi_Mult  ; 

        // Puppi softdrop mass from puppi subjets ( JEC applied separately to each subjet )
        JetPuppiSDmassRaw                        = sumPUPsubjets_P4_uncorr           .M()   ;
        JetPuppiSDmassSubjetCorr              = sumPUPsubjets_P4_L23res           .M()   ;
        // JetPuppiSDmassSubjetCorrUp            = sumPUPsubjets_P4_L23resCorrUp     .M()   ;
        // JetPuppiSDmassSubjetCorrDn            = sumPUPsubjets_P4_L23resCorrDn     .M()   ;
        // JetPuppiSDmassSubjetCorrL23Smear            = sumPUPsubjets_P4_L23resPtSmear    .M()   ;
        // JetPuppiSDmassSubjetCorrL23SmearUp          = sumPUPsubjets_P4_L23resPtSmearUp  .M()   ;
        // JetPuppiSDmassSubjetCorrL23SmearDn          = sumPUPsubjets_P4_L23resPtSmearDn  .M()   ;
        JetPuppiSDptRaw                          = sumPUPsubjets_P4_uncorr           .Perp();
        // JetPuppiSDptSubjetCorr                = sumPUPsubjets_P4_L23res           .Perp();
        // JetPuppiSDptSubjetCorrUp              = sumPUPsubjets_P4_L23resCorrUp     .Perp();
        // JetPuppiSDptSubjetCorrDn              = sumPUPsubjets_P4_L23resCorrDn     .Perp();
        // JetPuppiSDptSubjetCorrL23Smear              = sumPUPsubjets_P4_L23resPtSmear    .Perp();
        // JetPuppiSDptSubjetCorrL23SmearUp            = sumPUPsubjets_P4_L23resPtSmearUp  .Perp();
        // JetPuppiSDptSubjetCorrL23SmearDn            = sumPUPsubjets_P4_L23resPtSmearDn  .Perp();
        JetPuppiSDetaRaw                         = sumPUPsubjets_P4_uncorr           .Eta() ;
        JetPuppiSDphiRaw                         = sumPUPsubjets_P4_uncorr           .Phi() ;

        // PUPPI user floats from the toolbox
        // JetPuppiSDmassUserFloat               = puppi_softDropMass ;
        JetPuppiMassPruned                    = puppi_prunedMass   ;
        JetPuppiMassTrimmed                   = puppi_trimmedMass  ; 
        JetPuppiTau1                          = puppi_tau1         ;                  
        JetPuppiTau2                          = puppi_tau2         ;                  
        JetPuppiTau3                          = puppi_tau3         ;                  
        JetPuppiTau4                          = puppi_tau4         ;                  
        JetPuppiTau32                         = puppi_tau32        ;                  
        JetPuppiTau21                         = puppi_tau21        ;   

        // PUPPI subjet variables               
        JetPuppiSDsubjet0bdisc                = pup0_bdisc                ;
        JetPuppiSDsubjet1bdisc                = pup1_bdisc                ;
        JetPuppiSDmaxbdisc                    = pup_maxbdisc              ;
        JetPuppiSDmaxbdiscflavHadron          = pup_maxbdiscflav_hadron   ;
        JetPuppiSDmaxbdiscflavParton          = pup_maxbdiscflav_parton   ;
        JetPuppiSDsubjet0pt                   = pup0_P4_uncorr.Pt()       ;
        JetPuppiSDsubjet0mass                 = pup0_P4_uncorr.M()        ;
        JetPuppiSDsubjet0eta                  = pup0_P4_uncorr.Eta()      ;
        JetPuppiSDsubjet0phi                  = pup0_P4_uncorr.Phi()      ;
        JetPuppiSDsubjet0area                 = pup0_area                 ;
        JetPuppiSDsubjet0flavHadron           = pup0_flav_hadron          ; 
        JetPuppiSDsubjet0flavParton           = pup0_flav_parton          ;
        JetPuppiSDsubjet0matchedgenjetpt      = pup0_genpt                ;       
        JetPuppiSDsubjet0tau1                 = pup0_tau1                 ;  
        JetPuppiSDsubjet0tau2                 = pup0_tau2                 ;  
        JetPuppiSDsubjet0tau3                 = pup0_tau3                 ; 
        JetPuppiSDsubjet1pt                   = pup1_P4_uncorr.Pt()       ;                 
        JetPuppiSDsubjet1mass                 = pup1_P4_uncorr.M()        ; 
        JetPuppiSDsubjet1eta                  = pup1_P4_uncorr.Eta()      ;
        JetPuppiSDsubjet1phi                  = pup1_P4_uncorr.Phi()      ;             
        JetPuppiSDsubjet1area                 = pup1_area                 ;              
        JetPuppiSDsubjet1flavHadron           = pup1_flav_hadron          ;   
        JetPuppiSDsubjet1flavParton           = pup1_flav_parton          ;   
        JetPuppiSDsubjet1matchedgenjetpt      = pup1_genpt                ;       
        JetPuppiSDsubjet1tau1                 = pup1_tau1                 ;  
        JetPuppiSDsubjet1tau2                 = pup1_tau2                 ;  
        JetPuppiSDsubjet1tau3                 = pup1_tau3                 ;
        JetPuppiSDECF1                        = puppi_ECF1                 ; 
        JetPuppiSDECF2                        = puppi_ECF2                 ; 
        JetPuppiSDECF3                        = puppi_ECF3                 ; 
        JetPuppiSDECF4                        = puppi_ECF4                 ;
        JetPuppiSDECF5                        = puppi_ECF5                 ;  
        JetPuppiSDC_2                         = puppi_ECF3/(puppi_ECF2*puppi_ECF2);
        JetPuppiSDD_2                         = puppi_ECF3/(puppi_ECF2*puppi_ECF2*puppi_ECF2);
        JetPuppiSDC_3                         = puppi_ECF4/(puppi_ECF3*puppi_ECF3);
        JetPuppiSDD_3                         = puppi_ECF4/(puppi_ECF3*puppi_ECF3*puppi_ECF3);



        // AK8CHS JEC scale nom/up/down      
        JetCorrFactor                         = corr ;        
        JetCorrFactorUp                       = corrUp_L123 ;
        JetCorrFactorDn                       = corrDn_L123;
        // AK8CHS L2L3 JEC scale nom/up/down for groomed mass correction
        JetMassCorrFactor                     = corr_factor_L23res ;        
        JetMassCorrFactorUp                   = corrUp_L23 ;
        JetMassCorrFactorDn                   = corrDn_L23 ;
        // AK8CHS JER
        JetPtSmearFactor                      = ptsmear  ;
        JetPtSmearFactorUp                    = ptsmearUp;
        JetPtSmearFactorDn                    = ptsmearDn;         
        
        // AK8PUPPI JEC scale nom/up/down  (use for both full jet and groomed mass corrections)     
        JetPuppiCorrFactor                    = corr_factorAK8pup_L23res;          
        JetPuppiCorrFactorUp                  = corrUp_pup_L23;          
        JetPuppiCorrFactorDn                  = corrDn_pup_L23;    
        
        // AK8PUPPI JER
        JetPuppiPtSmearFactor                 = pup_ptsmear;          
        JetPuppiPtSmearFactorUp               = pup_ptsmearUp;          
        JetPuppiPtSmearFactorDn               = pup_ptsmearDn;  

        // AK8CHS JAR   
        // JetEtaScaleFactor                     = 1;          
        // JetPhiScaleFactor                     = 1;      

        // AK8CHS GenJet
        JetMatchedGenJetPt                    = GenJetMatched.Perp();       
        JetMatchedGenJetMass                  = GenJetMatched.M();   

        // AK8PUPPI GenJet
        JetPuppiMatchedGenJetPt               = GenJetMatchedPuppi.Perp();       
        JetPuppiMatchedGenJetMass             = GenJetMatchedPuppi.M();   
        // JetMatchedGenJetDR                 = GenJetMatched_dRmin;             


        if (!iEvent.isRealData() and runGenLoop_) {
          if (counttop==1 && jet_matched_t1){
            JetGenMatched_TopHadronic         = (int) tophadronic             ;
            JetGenMatched_TopPt               = t_p4.Perp()                   ;
            JetGenMatched_TopEta              = t_p4.Eta()                    ;
            JetGenMatched_TopPhi              = t_p4.Phi()                    ;
            JetGenMatched_TopMass             = t_p4.M()                      ;
            JetGenMatched_bPt                 = b_p4.Perp()                   ;
            JetGenMatched_WPt                 = W_p4.Perp()                   ;
            JetGenMatched_Wd1Pt               = Wd1_p4.Perp()                 ;
            JetGenMatched_Wd2Pt               = Wd2_p4.Perp()                 ;
            JetGenMatched_Wd1ID               = Wd1_id                        ;
            JetGenMatched_Wd2ID               = Wd2_id                        ;
         //   JetGenMatched_MaxDeltaRPartonTop  = max_deltaR_parton_t1           ;
         //   JetGenMatched_MaxDeltaRWPartonTop = max_deltaR_Wparton_t1          ;
         //   JetGenMatched_MaxDeltaRWPartonW   = max_deltaR_Wparton_W1          ;
         //   JetGenMatched_DeltaR_t_b          = deltaR_t1_b1                   ;
         //   JetGenMatched_DeltaR_t_W          = deltaR_t1_W1                   ;
         //   JetGenMatched_DeltaR_t_Wd1        = deltaR_t1_W1d1                 ;
         //   JetGenMatched_DeltaR_t_Wd2        = deltaR_t1_W1d2                 ;
         //   JetGenMatched_DeltaR_W_b1         = deltaR_W1_b1                   ;
         //   JetGenMatched_DeltaR_W_Wd1        = deltaR_W1_W1d1                 ;
         //   JetGenMatched_DeltaR_W_Wd2        = deltaR_W1_W1d2                 ;
         //   JetGenMatched_DeltaR_Wd1_Wd2      = deltaR_W1d1_W1d2               ;
         //   JetGenMatched_DeltaR_Wd1_b        = deltaR_W1d1_b1                 ;
         //   JetGenMatched_DeltaR_Wd2_b        = deltaR_W1d2_b1                 ;
         //   JetGenMatched_DeltaR_jet_t        = deltaR_jet_t1                  ;
            JetGenMatched_DeltaR_jet_W        = jet_p4.DeltaR(W_p4  )         ;
            JetGenMatched_DeltaR_jet_b        = jet_p4.DeltaR(b_p4  )         ;
            JetGenMatched_DeltaR_jet_Wd1      = jet_p4.DeltaR(Wd1_p4)         ;
            JetGenMatched_DeltaR_jet_Wd2      = jet_p4.DeltaR(Wd2_p4)         ;
            JetGenMatched_DeltaR_pup0_b       = pup0_P4_L23res.DeltaR(b_p4)   ;
            JetGenMatched_DeltaR_pup0_Wd1     = pup0_P4_L23res.DeltaR(Wd1_p4) ;
            JetGenMatched_DeltaR_pup0_Wd2     = pup0_P4_L23res.DeltaR(Wd2_p4) ;
            JetGenMatched_DeltaR_pup1_b       = pup1_P4_L23res.DeltaR(b_p4)   ;
            JetGenMatched_DeltaR_pup1_Wd1     = pup1_P4_L23res.DeltaR(Wd1_p4) ;
            JetGenMatched_DeltaR_pup1_Wd2     = pup1_P4_L23res.DeltaR(Wd2_p4) ;

            if (verbose_){
              cout<<"  JetGenMatched_DeltaR_pup0_b   "<<JetGenMatched_DeltaR_pup0_b  <<endl;
              cout<<"  JetGenMatched_DeltaR_pup0_Wd1 "<<JetGenMatched_DeltaR_pup0_Wd1<<endl;
              cout<<"  JetGenMatched_DeltaR_pup0_Wd2 "<<JetGenMatched_DeltaR_pup0_Wd2<<endl;
            }
          }   
          if (counttop==0 && jet_matched_p1){
            JetGenMatched_partonPt               = hardest_parton_hardScatterOutgoing_p4.Perp()           ;
            JetGenMatched_partonEta              = hardest_parton_hardScatterOutgoing_p4.Eta()            ;
            JetGenMatched_partonPhi              = hardest_parton_hardScatterOutgoing_p4.Phi()            ;
            JetGenMatched_partonMass             = hardest_parton_hardScatterOutgoing_p4.M()              ;
            JetGenMatched_partonID               = parton1id                                              ;
            JetGenMatched_DeltaRjetParton        = deltaR_jet_p1                                          ;
          }

        }
      } // end if this jet is opposite the lepton



    }// end if event has 1 lepton
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

TLorentzVector reconstructed_top;
TLorentzVector AK4_b;
TLorentzVector AK4_W;
TLorentzVector AK4_W2;





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

    reco::Candidate::LorentzVector corrJet = ijet.correctedP4(0);
    AK4_p4[count_AK4CHS].SetPtEtaPhiM(corrJet.pt(), corrJet.eta(), corrJet.phi(), corrJet.mass() );

    if( AK4_p4[count_AK4CHS].DeltaR(reconstructed_top) < 2*174/leading_CA12.Pt()){

      if (ijet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") >AK4bJet_bDisc &&  ijet.pt() > 30){
         AK4bJetPtRaw = ijet.pt();       
         AK4bJetEtaRaw = ijet.eta();      
         AK4bJetPhiRaw = ijet.phi();      
         AK4bJetMassRaw = ijet.mass();     
         AK4bJet_bDisc = ijet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
         AK4_b.SetPtEtaPhiM(ijet.pt(),ijet.eta(),ijet.phi(),ijet.mass());
       } else {
         break;
       }


    }
    count_AK4CHS++;
  }


count_AK4CHS = 0;
if (verbose_) cout<<"AK4 jet loop"<<endl;
for (const pat::Jet &ijet : *AK4MINI) { 
  if (ijet.pt()<15 || fabs(ijet.eta())>3.0) continue; 



    reco::Candidate::LorentzVector uncorrJet = ijet.correctedP4(0);
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
     
      // Get Smearings  
      // --- If well matched, smear based on GenJet, If not well matched,  gaussian smear based on resolution
      TLorentzVector AK4JetP4;
      AK4JetP4.SetPtEtaPhiM( corrJet.pt(), corrJet.eta(), corrJet.phi(), corrJet.mass() );
      double DeltaR_gen_reco  = AK4JetP4.DeltaR( GenJetMatched );
      double DeltaPt_gen_reco = AK4JetP4.Pt() - GenJetMatched.Pt()  ;
      double jet_distance_param = 0.4; 
      if (verbose_) cout<<"      -> gen pt "<<GenJetMatched.Pt()<<" reco pt "<<AK4JetP4.Pt()<<"  delta "<<DeltaPt_gen_reco<<endl;


      //uncorrected_err = uncorrected_err + ijet.pt() - genpt;
      //uncorrected_se = uncorrected_se + abs(ijet.pt() - genpt);
      //corrected_err = corrected_err + corrJet.pt() - genpt;
      //corrected_se = corrected_se + abs(corrJet.pt() - genpt);
      //cout << "corr vs uncorr " << corrJet.pt() - genpt << " " << uncorrJet.pt() - genpt << endl;
      //cout << "corr vs uncorr err " << corrected_err << " " << uncorrected_err << endl;
      //cout << "corr vs uncorr se " << corrected_se << " " << uncorrected_se << endl;


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

    if (verbose_) cout<<"   -> ptsmear "<<ptsmear<<" ptsmearUp "<<ptsmearUp<<" ptsmearDn "<< ptsmearDn<<endl;
 
  if(count_AK4CHS < 5){
    AK4_p4[count_AK4CHS].SetPtEtaPhiM(corrJet.pt(), corrJet.eta(), corrJet.phi(), corrJet.mass() );


   if(AK4_p4[count_AK4CHS].DeltaR(reconstructed_top) < 2*174/leading_CA12.Pt()){


     if(ijet.pt() > 30  && AK4WJetPtRaw == 0.0 && AK4bJetPtRaw!=ijet.pt()){
        AK4WJetPtRaw = ijet.pt();       
        AK4WJetEtaRaw = ijet.eta();      
        AK4WJetPhiRaw = ijet.phi();      
        AK4WJetMassRaw = ijet.mass();     
        AK4WJet_bDisc = ijet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        AK4_W.SetPtEtaPhiM(ijet.pt(),ijet.eta(),ijet.phi(),ijet.mass());
     } else if(ijet.pt() > 30 &&  AK4W2JetPtRaw == 0.0 && AK4bJetPtRaw!=ijet.pt()){
        AK4W2JetPtRaw = ijet.pt();       
        AK4W2JetEtaRaw = ijet.eta();      
        AK4W2JetPhiRaw = ijet.phi();      
        AK4W2JetMassRaw = ijet.mass();     
        AK4W2Jet_bDisc = ijet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        AK4_W2.SetPtEtaPhiM(ijet.pt(),ijet.eta(),ijet.phi(),ijet.mass());

     }



   }





    //cout << "gen, ak8, ca12 pt " <<t_p4.Pt() << " " << AK8jet_had_P4corr.Pt() << " " << leading_CA12.Pt()  << endl;
    //cout << "gen, ak8, ca12 m " <<t_p4.M() << " " << AK8jet_had_P4corr.M() << " " <<leading_CA12.M()  << endl;
    //cout << "gen b and W deltaR , 2*mT/ca12 pt " <<b_p4.DeltaR(W_p4) << " " << 2*174/leading_CA12.Pt() << endl;
//
    //cout << "deltaR Gen t, leading ca12 " << t_p4.DeltaR(leading_CA12) << endl;
    //cout << "deltaR Gen t, leading ak8 " << t_p4.DeltaR(AK8jet_had_P4corr) << endl;
    //cout << "deltaR Gen t, leading ak8 " << t_p4.DeltaR(reconstructed_top) << endl;
//
    //cout << "deltaR Ak4, leading CA12 " << leading_CA12.DeltaR(AK4_p4[count_AK4CHS]) << endl;
//
    //cout << "deltaR Ak4, leading ak8 " << AK8jet_had_P4corr.DeltaR(AK4_p4[count_AK4CHS]) << endl;
//
    //cout << "deltaR Gen t, ak4 " << t_p4.DeltaR(AK4_p4[count_AK4CHS]) << endl;
    //cout << "deltaR Gen b, ak4 " << b_p4.DeltaR(AK4_p4[count_AK4CHS]) << endl;
    //cout << "deltaR Gen W, ak4 " << W_p4.DeltaR(AK4_p4[count_AK4CHS]) << endl;
//
    //cout <<  " bjet " << ijet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags")  << endl;
//


  } 

  count_AK4CHS++;
}
if (count_AK4CHS > 0){
  //cout << "b, w, w2 pt " <<AK4bJetPtRaw << " " << AK4WJetPtRaw << " " << AK4W2JetPtRaw  << endl;
  //cout << "b, w, w2 btag  " <<AK4bJet_bDisc << " " << AK4WJet_bDisc << " " << AK4W2Jet_bDisc  << endl;
  //cout << "b, gen DR E(ak4)/E(gen) " << b_p4.DeltaR(AK4_b) << " " << AK4_b.E()/b_p4.E() << endl;
  //cout << "W, gen DR E(ak4)/E(gen) " << W_p4.DeltaR(AK4_W) << " " << AK4_W.E()/W_p4.E() << endl;
  //cout << "W2, gen DR E(ak4)/E(gen) " << W_p4.DeltaR(AK4_W+AK4_W2) << " " << (AK4_W+AK4_W2).E()/W_p4.E() << endl;
  //cout << "b, w, w2 mass  " << (AK4_b+AK4_W+AK4_W2).M()<< " " << (AK4_b+AK4_W+AK4_W2).Pt() << endl;
  //cout << "b, w mass  " << (AK4_b+AK4_W).M() << " " << (AK4_b+AK4_W).Pt()<< endl;
  //cout << "ca12 mass  pt " << leading_CA12.M() << " " << leading_CA12.Pt() << endl;
  //cout << "ak8 mass  pt " << AK8jet_had_P4corr.M() << " " << AK8jet_had_P4corr.Pt() << endl;
  //cout << "gen mass  pt " << t_p4.M() << " " << t_p4.Pt() << endl;

  reconstructed_top = (AK4_b+AK4_W+AK4_W2);

  AK4ReconstructedJetPt = reconstructed_top.Pt();
  AK4ReconstructedJetEta = reconstructed_top.Eta();
  AK4ReconstructedJetPhi = reconstructed_top.Phi();
  AK4ReconstructedJetMass = reconstructed_top.M();
 
}


  //   Float_t AK4ReconstructedJetPt                             ;      
  //   Float_t AK4ReconstructedJetEta                             ;
  //   Float_t AK4ReconstructedJetPhi                             ;
  //   Float_t AK4ReconstructedJetMass                             ;

  //   Float_t AK4bJetPtRaw                               ;      
  //   Float_t AK4bJetEtaRaw                              ;
  //   Float_t AK4bJetPhiRaw                              ;
  //   Float_t AK4bJetMassRaw                             ;

  //   Float_t AK4bJet_PtSmear                           ;              
  //   Float_t AK4bJet_PtSmearUp                         ;              
  //   Float_t AK4bJet_PtSmearDn                         ;              
  //   Float_t AK4bJet_PtUncorr                          ; 
  //   Float_t AK4bJet_Corr                              ; 
  //   Float_t AK4bJet_CorrUp                            ; 
  //   Float_t AK4bJet_CorrDn                            ;

  //   Float_t AK4WJet_PtSmear                           ;              
  //   Float_t AK4WJet_PtSmearUp                         ;              
  //   Float_t AK4WJet_PtSmearDn                         ;              
  //   Float_t AK4WJet_PtUncorr                          ; 
  //   Float_t AK4WJet_Corr                              ; 
  //   Float_t AK4WJet_CorrUp                            ; 
  //   Float_t AK4WJet_CorrDn                            ;


  //   Float_t AK4W2Jet_PtSmear                           ;              
  //   Float_t AK4W2Jet_PtSmearUp                         ;              
  //   Float_t AK4W2Jet_PtSmearDn                         ;              
  //   Float_t AK4W2Jet_PtUncorr                          ; 
  //   Float_t AK4W2Jet_Corr                              ; 
  //   Float_t AK4W2Jet_CorrUp                            ; 
  //   Float_t AK4W2Jet_CorrDn                            ;





  if(1==2){
    std::cout << NNPDF3wgt_down<< NNPDF3wgt_up<<Q2wgt_up<<Q2wgt_down<<std::endl;
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

    if ( 1==1 /*AK8jet_Had_P4corr.Perp()>200*/){
      h_cutflow_had   ->Fill(2.5);

      if ( 1==1 /*fabs( AK8jet_Had_P4corr.Rapidity() ) <2.4 */ ){
        h_cutflow_had   ->Fill(3.5);

        if (verbose_) cout<<"Write Semi-Lept Tree"<<endl;

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

        // if (verbose_){
        //   cout<<" met.pt() "<<   met.pt() <<endl;
        //   cout<<" met.shiftedPt(pat::MET::JetEnUp)             "<<met.shiftedPt(pat::MET::JetEnUp)            <<endl;       
        //   cout<<" met.shiftedPt(pat::MET::JetEnDown)           "<<met.shiftedPt(pat::MET::JetEnDown)          <<endl;               
        //   cout<<" met.shiftedPt(pat::MET::ElectronEnUp)        "<<met.shiftedPt(pat::MET::ElectronEnUp)       <<endl;            
        //   cout<<" met.shiftedPt(pat::MET::ElectronEnDown)      "<<met.shiftedPt(pat::MET::ElectronEnDown)     <<endl;                            
        //   cout<<" met.shiftedPt(pat::MET::MuonEnUp)            "<<met.shiftedPt(pat::MET::MuonEnUp)           <<endl;        
        //   cout<<" met.shiftedPt(pat::MET::MuonEnDown)          "<<met.shiftedPt(pat::MET::MuonEnDown)         <<endl;                      
        //   cout<<" met.shiftedPt(pat::MET::JetResUp)            "<<met.shiftedPt(pat::MET::JetResUp)           <<endl;        
        //   cout<<" met.shiftedPt(pat::MET::JetResDown)          "<<met.shiftedPt(pat::MET::JetResDown)         <<endl;               
        //   cout<<" met.shiftedPt(pat::MET::UnclusteredEnUp)     "<<met.shiftedPt(pat::MET::UnclusteredEnUp)    <<endl;               
        //   cout<<" met.shiftedPt(pat::MET::UnclusteredEnDown)   "<<met.shiftedPt(pat::MET::UnclusteredEnDown)  <<endl;                 
          
        //   cout<<" met.phi() "<<   met.phi() <<endl;
        //   cout<<" met.shiftedPhi(pat::MET::UnclusteredEnUp)    "<<met.shiftedPhi(pat::MET::UnclusteredEnUp)   <<endl;                
        //   cout<<" met.shiftedPhi(pat::MET::UnclusteredEnDown)  "<<met.shiftedPhi(pat::MET::UnclusteredEnDown) <<endl;                  
        //   cout<<" met.shiftedPhi(pat::MET::JetEnUp)            "<<met.shiftedPhi(pat::MET::JetEnUp)           <<endl;        
        //   cout<<" met.shiftedPhi(pat::MET::JetEnDown)          "<<met.shiftedPhi(pat::MET::JetEnDown)         <<endl; 
        //   cout<<" met.shiftedPhi(pat::MET::ElectronEnUp)       "<<met.shiftedPhi(pat::MET::ElectronEnUp)      <<endl;             
        //   cout<<" met.shiftedPhi(pat::MET::ElectronEnDown)     "<<met.shiftedPhi(pat::MET::ElectronEnDown)    <<endl; 
        //   cout<<" met.shiftedPhi(pat::MET::MuonEnUp)           "<<met.shiftedPhi(pat::MET::MuonEnUp)          <<endl;         
        //   cout<<" met.shiftedPhi(pat::MET::MuonEnDown)         "<<met.shiftedPhi(pat::MET::MuonEnDown)        <<endl;
        //   cout<<" met.shiftedPhi(pat::MET::JetResUp)           "<<met.shiftedPhi(pat::MET::JetResUp)          <<endl;         
        //   cout<<" met.shiftedPhi(pat::MET::JetResDown)         "<<met.shiftedPhi(pat::MET::JetResDown)        <<endl;    
        // }

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

       

        //if      (count_mu==1 && count_el==0) LeptonIsMu  = 1  ; 
        //else if (count_mu==0 && count_el==1) LeptonIsMu  = 0  ; 
        //else                                 LeptonIsMu  = -1  ;
//
        //PtRel  = AK4_dRMinLep_p4.Perp( lep0_p4.Vect() );
        //MuIso  = mu0_iso04;
//
        //Elecron_absiso            = el0_absiso           ;  
        //Elecron_relIsoWithDBeta   = el0_relIsoWithDBeta  ;  
        //Elecron_absiso_EA         = el0_absiso_EA        ;  
        //Elecron_relIsoWithEA      = el0_relIsoWithEA     ;  
//
        //Electron_iso_passHLTpre   = el0_iso_passHLTpre  ;
        //Electron_iso_passLoose    = el0_iso_passLoose   ;
        //Electron_iso_passMedium   = el0_iso_passMedium  ;
        //Electron_iso_passTight    = el0_iso_passTight   ;
        //Electron_iso_passHEEP     = el0_iso_passHEEP    ;
        //Electron_noiso_passLoose  = el0_noiso_passLoose ;
        //Electron_noiso_passMedium = el0_noiso_passMedium;
        //Electron_noiso_passTight  = el0_noiso_passTight ;
        //Electron_noiso_passHEEP   = el0_noiso_passHEEP  ;
//
        //MuMedium = (int) mu0_isMedium   ;
        //MuTight  = (int) mu0_isTight    ;
        //MuHighPt = (int) mu0_isHighPt   ;

        TreeHad -> Fill();
      }// end rapidity selection
    }// end pt selection

 
    
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


