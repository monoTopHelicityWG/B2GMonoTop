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

#include "Analysis/B2GMonoTop/interface/eventDataStruct.h"
//#include "Analysis/B2GMonoTop/interface/eventTTree.h"

#ifdef __CINT__
#pragma link C++ class std::vector<TLorentzVector>+;
#endif

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
      bool runTTree_   ;

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
      TH1D * h_trigger_efficency_1;
      TH1D * h_trigger_efficency_2;
      TH1D * h_trigger_efficency_3;
      TH1D * h_trigger_efficency_4;
      TH1D * h_trigger_accept_1;
      TH1D * h_trigger_accept_2;
      TH1D * h_trigger_accept_3;
      TH1D * h_trigger_accept_4;
      TH1D * h_trigger_reject_1;
      TH1D * h_trigger_reject_2;
      TH1D * h_trigger_reject_3;
      TH1D * h_trigger_reject_4;
      TH1D * h_trigger_efficency_1_topPt;
      TH1D * h_trigger_efficency_2_topPt;
      TH1D * h_trigger_efficency_3_topPt;
      TH1D * h_trigger_efficency_4_topPt;
      TH1D * h_trigger_accept_1_topPt;
      TH1D * h_trigger_accept_2_topPt;
      TH1D * h_trigger_accept_3_topPt;
      TH1D * h_trigger_accept_4_topPt;
      TH1D * h_trigger_reject_1_topPt;
      TH1D * h_trigger_reject_2_topPt;
      TH1D * h_trigger_reject_3_topPt;
      TH1D * h_trigger_reject_4_topPt;

      TH1D * h_NtrueIntPU        ;
      TH1D * h_NPV               ;               
      TH1D * h_NPVreweighted     ;     
      TH1D * h_NPVgood           ;               
      TH1D * h_NPVgoodreweighted ;     

      // -- Triggers to be saved in tree
      std::vector<std::string> trigsToRunHad;
      std::vector<std::string> trigsToRunMu;
      std::vector<std::string> trigsToRun;


      Int_t nEvents = 0;
      Int_t totalEvents = 0;

      bool PFMET120_BTagCSV_Mu5_Trigger;
      bool PFMET300_Trigger;
      bool HLT_PFMET120_PFMHT120_Trigger;
      bool HLT_PFMET170_Trigger;

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
                                                                 
                 

      TTree *EventTTree; 
      eventDataStruct * event_data;


};