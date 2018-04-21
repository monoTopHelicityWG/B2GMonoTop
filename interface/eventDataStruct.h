#ifndef eventDataStruct_h
#define eventDataStruct_h
#define BADVAL -999.0

#include <vector>
#include <string>
#include "Rtypes.h"

typedef struct eventDataStruct {

      //do not reset
      int nEvents = 0;
      int nEventSaved = 0;

      std::vector<int>  *HadTrigPrescalesMu = new std::vector<int>  ;
      std::vector<bool> *HadTrigPassMu      = new std::vector<bool> ;
      std::string HadTrigAcceptBitsMu;

      std::vector<int>  *HadTrigPrescalesHad = new std::vector<int>  ;
      std::vector<bool> *HadTrigPassHad      = new std::vector<bool> ;
      std::string HadTrigAcceptBitsHad;

      std::vector<int>  *HadTrigPrescales = new std::vector<int>  ;
      std::vector<bool> *HadTrigPass      = new std::vector<bool> ;
      std::string HadTrigAcceptBits;

      std::vector<std::string>  *HLTtriggers = new std::vector<std::string>  ;
      std::vector<bool> *HLTtriggersPass      = new std::vector<bool> ;
      std::vector<int>  *HLTtriggersPrescales = new std::vector<int>  ;


      bool mu_preselection = false;
      bool had_preselection = false;

      bool PFMET120_BTagCSV_Mu5_Trigger = false;
      bool PFMET300_Trigger = false;
      bool HLT_PFMET120_PFMHT120_Trigger = false;
      bool HLT_PFMET170_Trigger = false;

      double uncorrected_err = 0;
      double uncorrected_se= 0;
      double corrected_err = 0;
      double corrected_se= 0;

      Float_t Gen_MET_pT;
      Float_t Gen_MET_phi;
      Float_t Gen_MET_eta;
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

      int bJet_count = 0;

      std::vector<float>* AK4JetLV_pt = new std::vector<float>;
      std::vector<float>* AK4JetLV_eta = new std::vector<float>;
      std::vector<float>* AK4JetLV_phi = new std::vector<float>;
      std::vector<float>* AK4JetLV_mass = new std::vector<float>;

      std::vector<float>* AK4JetLV_corr = new std::vector<float>;
      std::vector<float>* AK4JetLV_corrUp = new std::vector<float>;
      std::vector<float>* AK4JetLV_corrDn = new std::vector<float>;
      std::vector<float>* AK4JetLV_SF = new std::vector<float>;
      std::vector<float>* AK4JetLV_SFUp = new std::vector<float>;
      std::vector<float>* AK4JetLV_SFDn = new std::vector<float>;
      std::vector<float>* AK4JetLV_ptsmear = new std::vector<float>;
      std::vector<float>* AK4JetLV_ptsmearUp = new std::vector<float>;
      std::vector<float>* AK4JetLV_ptsmearDn = new std::vector<float>;


      std::vector<float>* AK8JetLV_corr = new std::vector<float>;
      std::vector<float>* AK8JetLV_corrUp = new std::vector<float>;
      std::vector<float>* AK8JetLV_corrDn = new std::vector<float>;
      std::vector<float>* AK8JetLV_SF = new std::vector<float>;
      std::vector<float>* AK8JetLV_SFUp = new std::vector<float>;
      std::vector<float>* AK8JetLV_SFDn = new std::vector<float>;
      std::vector<float>* AK8JetLV_ptsmear = new std::vector<float>;
      std::vector<float>* AK8JetLV_ptsmearUp = new std::vector<float>;
      std::vector<float>* AK8JetLV_ptsmearDn = new std::vector<float>;


      std::vector<float>* AK8JetLV_pt = new std::vector<float>;
      std::vector<float>* AK8JetLV_eta = new std::vector<float>;
      std::vector<float>* AK8JetLV_phi = new std::vector<float>;
      std::vector<float>* AK8JetLV_mass = new std::vector<float>;

      std::vector<float>* AK8SubjetLV_pt = new std::vector<float>;
      std::vector<float>* AK8SubjetLV_eta = new std::vector<float>;
      std::vector<float>* AK8SubjetLV_phi = new std::vector<float>;
      std::vector<float>* AK8SubjetLV_mass = new std::vector<float>;

      std::vector<float>* AK4JetBtag_p = new std::vector<float>;
      std::vector<float>* AK8JetTau1_p = new std::vector<float>;
      std::vector<float>* AK8JetTau2_p = new std::vector<float>;
      std::vector<float>* AK8JetTau3_p = new std::vector<float>;
      std::vector<float>* AK8JetSoftdropMass_p = new std::vector<float>;

      Float_t AK4_corr_pt                               ;  
      Float_t AK4_uncorr_pt                               ;  

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
      Float_t JetArea                                ;
      Float_t JetSDmassRaw                           ;
      Float_t JetSDmassSubjetCorrL23                       ;
      Float_t JetSDmassSubjetCorrL123                      ;
      Float_t JetSDptRaw                             ;
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
      Float_t JetPuppiSDptRaw                           ;
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

      std::vector<float>* MuPhi = new std::vector<float>                              ;
      std::vector<float>* MuPt = new std::vector<float>                               ;
      std::vector<float>* MuEta = new std::vector<float>                              ;
      std::vector<float>* MuMass = new std::vector<float>                             ;
      std::vector<float>* MuIso = new std::vector<float>                                  ;
      std::vector<float>* MuIsoTrk = new std::vector<float>                                  ;
      std::vector<int>* MuTight = new std::vector<int>                                ;
      std::vector<int>* MuMedium = new std::vector<int>                               ;

      std::vector<float>* Electron_Phi = new std::vector<float>                              ;
      std::vector<float>* Electron_Pt = new std::vector<float>                               ;
      std::vector<float>* Electron_Eta = new std::vector<float>                              ;
      std::vector<float>* Electron_Mass = new std::vector<float>                             ;
      std::vector<float>* Elecron_absiso = new std::vector<float>                         ;
      std::vector<float>* Elecron_relIsoWithDBeta = new std::vector<float>                ;
      std::vector<float>* Elecron_absiso_EA = new std::vector<float>                      ;
      std::vector<float>* Elecron_relIsoWithEA = new std::vector<float>                   ;

      std::vector<int>* Electron_iso_passHLTpre = new std::vector<int>                  ;
      std::vector<int>* Electron_iso_passLoose = new std::vector<int>                   ;
      std::vector<int>* Electron_iso_passMedium = new std::vector<int>                  ;
      std::vector<int>* Electron_iso_passTight = new std::vector<int>                   ;
      std::vector<int>* Electron_iso_passHEEP = new std::vector<int>                    ;
      std::vector<int>* Electron_noiso_passLoose = new std::vector<int>                 ;
      std::vector<int>* Electron_noiso_passMedium = new std::vector<int>                ;
      std::vector<int>* Electron_noiso_passTight = new std::vector<int>                 ;
      std::vector<int>* Electron_noiso_passHEEP = new std::vector<int>                  ;

      void resetArray(Float_t (&arr)[4]){
            for(int i =0; i < 4; i++){
                  arr[i] = -999.9;
            }
      }

      void resetStruct(){

            //vectors
            HadTrigPrescalesMu->clear();
            HadTrigPassMu->clear();
            HadTrigPrescalesHad->clear();
            HadTrigPassHad->clear();
            HadTrigPrescales->clear();
            HadTrigPass->clear();
            HLTtriggers->clear();
            HLTtriggersPass->clear();
            HLTtriggersPrescales->clear();
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
            AK8JetLV_corr->clear();
            AK8JetLV_corrUp->clear();
            AK8JetLV_corrDn->clear();
            AK8JetLV_SF->clear();
            AK8JetLV_SFUp->clear();
            AK8JetLV_SFDn->clear();
            AK8JetLV_ptsmear->clear();
            AK8JetLV_ptsmearUp->clear();
            AK8JetLV_ptsmearDn->clear();
            AK8JetLV_pt->clear();
            AK8JetLV_eta->clear();
            AK8JetLV_phi->clear();
            AK8JetLV_mass->clear();
            AK8SubjetLV_pt->clear();
            AK8SubjetLV_eta->clear();
            AK8SubjetLV_phi->clear();
            AK8SubjetLV_mass->clear();
            AK4JetBtag_p->clear();
            AK8JetTau1_p->clear();
            AK8JetTau2_p->clear();
            AK8JetTau3_p->clear();
            AK8JetSoftdropMass_p->clear();
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
            Electron_iso_passLoose->clear();
            Electron_iso_passMedium->clear();
            Electron_iso_passTight->clear();
            Electron_iso_passHEEP->clear();
            Electron_noiso_passLoose->clear();
            Electron_noiso_passMedium->clear();
            Electron_noiso_passTight->clear();
            Electron_noiso_passHEEP->clear();
            
            //strings
            HadTrigAcceptBitsMu = "";
            HadTrigAcceptBitsHad = "";
            HadTrigAcceptBits = "";
            
            //bool
            
            mu_preselection = false;
            had_preselection = false;
            PFMET120_BTagCSV_Mu5_Trigger = false;
            PFMET300_Trigger = false;
            HLT_PFMET120_PFMHT120_Trigger = false;
            HLT_PFMET170_Trigger = false;
            tophadronic=false;
            topleptonic=false;
            
            //double
            
            uncorrected_err = 0;
            uncorrected_se= 0;
            corrected_err = 0;
            corrected_se= 0;
            
            
            ///Float_t
            
            Gen_MET_pT = -999.9;
            Gen_MET_phi = -999.9;
            Gen_MET_eta = -999.9;


            resetArray(Gen_array_t_p4);
            resetArray(Gen_array_final_t_p4);
            resetArray(Gen_array_b_p4);
            resetArray(Gen_array_W_p4);
            resetArray(Gen_array_Wd1_p4);
            resetArray(Gen_array_Wd2_p4);
            resetArray(Gen_array_hardest_parton_hardScatterOutgoing_p4);
            resetArray(Gen_array_second_hardest_parton_hardScatterOutgoing_p4);
            AK4_corr_pt = -999.9;
            AK4_uncorr_pt = -999.9;
            CA12JetPtRaw = -999.9;
            CA12JetEtaRaw = -999.9;
            CA12JetPhiRaw = -999.9;
            CA12JetMassRaw = -999.9;
            CA12JetTau1 = -999.9;
            CA12JetTau2 = -999.9;
            CA12JetTau3 = -999.9;
            CA12JetTau4 = -999.9;
            CA12JetTau32 = -999.9;
            CA12JetTau21 = -999.9;
            CA12Jetsubjet0bdisc = -999.9;
            CA12Jetsubjet1bdisc = -999.9;
            CA12Jetmaxbdisc = -999.9;
            CA12Jetsubjet0pt = -999.9;
            CA12Jetsubjet0mass = -999.9;
            CA12Jetsubjet0eta = -999.9;
            CA12Jetsubjet0phi = -999.9;
            CA12Jetsubjet0area = -999.9;
            CA12Jetsubjet1pt = -999.9;
            CA12Jetsubjet1mass = -999.9;
            CA12Jetsubjet1eta = -999.9;
            CA12Jetsubjet1phi = -999.9;
            CA12Jetsubjet1area = -999.9;
            AK4ReconstructedJetPt = -999.9;
            AK4ReconstructedJetEta = -999.9;
            AK4ReconstructedJetPhi = -999.9;
            AK4ReconstructedJetMass = -999.9;
            AK4bJetPtRaw = -999.9;
            AK4bJetEtaRaw = -999.9;
            AK4bJetPhiRaw = -999.9;
            AK4bJetMassRaw = -999.9;
            AK4bJet_PtSmear = -999.9;
            AK4bJet_PtSmearUp = -999.9;
            AK4bJet_PtSmearDn = -999.9;
            AK4bJet_PtUncorr = -999.9;
            AK4bJet_Corr = -999.9;
            AK4bJet_CorrUp = -999.9;
            AK4bJet_CorrDn = -999.9;
            AK4bJet_bDisc = -999.9;
            AK4WJetPtRaw = -999.9;
            AK4WJetEtaRaw = -999.9;
            AK4WJetPhiRaw = -999.9;
            AK4WJetMassRaw = -999.9;
            AK4WJet_PtSmear = -999.9;
            AK4WJet_PtSmearUp = -999.9;
            AK4WJet_PtSmearDn = -999.9;
            AK4WJet_PtUncorr = -999.9;
            AK4WJet_Corr = -999.9;
            AK4WJet_CorrUp = -999.9;
            AK4WJet_CorrDn = -999.9;
            AK4WJet_bDisc = -999.9;
            AK4W2JetPtRaw = -999.9;
            AK4W2JetEtaRaw = -999.9;
            AK4W2JetPhiRaw = -999.9;
            AK4W2JetMassRaw = -999.9;
            AK4W2Jet_PtSmear = -999.9;
            AK4W2Jet_PtSmearUp = -999.9;
            AK4W2Jet_PtSmearDn = -999.9;
            AK4W2Jet_PtUncorr = -999.9;
            AK4W2Jet_Corr = -999.9;
            AK4W2Jet_CorrUp = -999.9;
            AK4W2Jet_CorrDn = -999.9;
            AK4W2Jet_bDisc = -999.9;
            JetPtRaw = -999.9;
            JetEtaRaw = -999.9;
            JetPhiRaw = -999.9;
            JetMassRaw = -999.9;
            JetArea = -999.9;
            JetSDmassRaw = -999.9;
            JetSDmassSubjetCorrL23 = -999.9;
            JetSDmassSubjetCorrL123 = -999.9;
            JetSDptRaw = -999.9;
            JetSDetaRaw = -999.9;
            JetSDphiRaw = -999.9;
            JetMassPruned = -999.9;
            JetMassTrimmed = -999.9;
            JetTau1 = -999.9;
            JetTau2 = -999.9;
            JetTau3 = -999.9;
            JetTau4 = -999.9;
            JetTau32 = -999.9;
            JetTau21 = -999.9;
            JetSDsubjet0bdisc = -999.9;
            JetSDsubjet1bdisc = -999.9;
            JetSDmaxbdisc = -999.9;
            JetSDmaxbdiscflavHadron = -999.9;
            JetSDmaxbdiscflavParton = -999.9;
            JetSDsubjet0pt = -999.9;
            JetSDsubjet0mass = -999.9;
            JetSDsubjet0eta = -999.9;
            JetSDsubjet0phi = -999.9;
            JetSDsubjet0area = -999.9;
            JetSDsubjet0flavHadron = -999.9;
            JetSDsubjet0flavParton = -999.9;
            JetSDsubjet0matchedgenjetpt = -999.9;
            JetSDsubjet0tau1 = -999.9;
            JetSDsubjet0tau2 = -999.9;
            JetSDsubjet0tau3 = -999.9;
            JetSDsubjet1pt = -999.9;
            JetSDsubjet1mass = -999.9;
            JetSDsubjet1eta = -999.9;
            JetSDsubjet1phi = -999.9;
            JetSDsubjet1area = -999.9;
            JetSDsubjet1flavHadron = -999.9;
            JetSDsubjet1flavParton = -999.9;
            JetSDsubjet1matchedgenjetpt = -999.9;
            JetSDsubjet1tau1 = -999.9;
            JetSDsubjet1tau2 = -999.9;
            JetSDsubjet1tau3 = -999.9;
            //JetPuppiP = -999.9;
            JetPuppiPtRaw = -999.9;
            JetPuppiEtaRaw = -999.9;
            JetPuppiPhiRaw = -999.9;
            JetPuppiMassRaw = -999.9;
            JetPuppiArea = -999.9;
            JetPuppiSDmassRaw = -999.9;
            JetPuppiSDmassSubjetCorr = -999.9;
            JetPuppiSDptRaw = -999.9;
            JetPuppiSDetaRaw = -999.9;
            JetPuppiSDphiRaw = -999.9;
            JetPuppiTau1 = -999.9;
            JetPuppiTau2 = -999.9;
            JetPuppiTau3 = -999.9;
            JetPuppiTau4 = -999.9;
            JetPuppiTau32 = -999.9;
            JetPuppiTau21 = -999.9;
            JetPuppiSDsubjet0bdisc = -999.9;
            JetPuppiSDsubjet1bdisc = -999.9;
            JetPuppiSDmaxbdisc = -999.9;
            JetPuppiSDmaxbdiscflavHadron = -999.9;
            JetPuppiSDmaxbdiscflavParton = -999.9;
            JetPuppiSDsubjet0pt = -999.9;
            JetPuppiSDsubjet0mass = -999.9;
            JetPuppiSDsubjet0eta = -999.9;
            JetPuppiSDsubjet0phi = -999.9;
            JetPuppiSDsubjet0area = -999.9;
            JetPuppiSDsubjet0flavHadron = -999.9;
            JetPuppiSDsubjet0flavParton = -999.9;
            JetPuppiSDsubjet0matchedgenjetpt = -999.9;
            JetPuppiSDsubjet0tau1 = -999.9;
            JetPuppiSDsubjet0tau2 = -999.9;
            JetPuppiSDsubjet0tau3 = -999.9;
            JetPuppiSDsubjet1pt = -999.9;
            JetPuppiSDsubjet1mass = -999.9;
            JetPuppiSDsubjet1eta = -999.9;
            JetPuppiSDsubjet1phi = -999.9;
            JetPuppiSDsubjet1area = -999.9;
            JetPuppiSDsubjet1flavHadron = -999.9;
            JetPuppiSDsubjet1flavParton = -999.9;
            JetPuppiSDsubjet1matchedgenjetpt = -999.9;
            JetPuppiSDsubjet1tau1 = -999.9;
            JetPuppiSDsubjet1tau2 = -999.9;
            JetPuppiSDsubjet1tau3 = -999.9;
            JetPuppiSDECF1 = -999.9;
            JetPuppiSDECF2 = -999.9;
            JetPuppiSDECF3 = -999.9;
            JetPuppiSDECF4 = -999.9;
            JetPuppiSDECF5 = -999.9;
            JetPuppiSDC_2 = -999.9;
            JetPuppiSDD_2 = -999.9;
            JetPuppiSDC_3 = -999.9;
            JetPuppiSDD_3 = -999.9;
            JetPuppiMassPruned = -999.9;
            JetPuppiMassTrimmed = -999.9;
            JetCHF = -999.9;
            JetNHF = -999.9;
            JetCM = -999.9;
            JetNM = -999.9;
            JetNEF = -999.9;
            JetCEF = -999.9;
            JetMF = -999.9;
            JetMult = -999.9;
            JetPuppiCHF = -999.9;
            JetPuppiNHF = -999.9;
            JetPuppiCM = -999.9;
            JetPuppiNM = -999.9;
            JetPuppiNEF = -999.9;
            JetPuppiCEF = -999.9;
            JetPuppiMF = -999.9;
            JetPuppiMult = -999.9;
            JetMassCorrFactor = -999.9;
            JetMassCorrFactorUp = -999.9;
            JetMassCorrFactorDn = -999.9;
            JetCorrFactor = -999.9;
            JetCorrFactorUp = -999.9;
            JetCorrFactorDn = -999.9;
            JetPtSmearFactor = -999.9;
            JetPtSmearFactorUp = -999.9;
            JetPtSmearFactorDn = -999.9;
            JetPuppiMassCorrFactor = -999.9;
            JetPuppiMassCorrFactorUp = -999.9;
            JetPuppiMassCorrFactorDn = -999.9;
            JetPuppiCorrFactor = -999.9;
            JetPuppiCorrFactorUp = -999.9;
            JetPuppiCorrFactorDn = -999.9;
            JetPuppiPtSmearFactor = -999.9;
            JetPuppiPtSmearFactorUp = -999.9;
            JetPuppiPtSmearFactorDn = -999.9;
            JetMatchedGenJetPt = -999.9;
            JetMatchedGenJetMass = -999.9;
            JetPuppiMatchedGenJetPt = -999.9;
            JetPuppiMatchedGenJetMass = -999.9;
            JetGenMatched_TopPt = -999.9;
            JetGenMatched_TopEta = -999.9;
            JetGenMatched_TopPhi = -999.9;
            JetGenMatched_TopMass = -999.9;
            JetGenMatched_bPt = -999.9;
            JetGenMatched_WPt = -999.9;
            JetGenMatched_Wd1Pt = -999.9;
            JetGenMatched_Wd2Pt = -999.9;
            JetGenMatched_Wd1ID = -999.9;
            JetGenMatched_Wd2ID = -999.9;
            JetGenMatched_MaxDeltaRPartonTop = -999.9;
            JetGenMatched_MaxDeltaRWPartonTop = -999.9;
            JetGenMatched_MaxDeltaRWPartonW = -999.9;
            JetGenMatched_DeltaR_t_b = -999.9;
            JetGenMatched_DeltaR_t_W = -999.9;
            JetGenMatched_DeltaR_t_Wd1 = -999.9;
            JetGenMatched_DeltaR_t_Wd2 = -999.9;
            JetGenMatched_DeltaR_W_b1 = -999.9;
            JetGenMatched_DeltaR_W_Wd1 = -999.9;
            JetGenMatched_DeltaR_W_Wd2 = -999.9;
            JetGenMatched_DeltaR_Wd1_Wd2 = -999.9;
            JetGenMatched_DeltaR_Wd1_b = -999.9;
            JetGenMatched_DeltaR_Wd2_b = -999.9;
            JetGenMatched_DeltaR_jet_t = -999.9;
            JetGenMatched_DeltaR_jet_W = -999.9;
            JetGenMatched_DeltaR_jet_b = -999.9;
            JetGenMatched_DeltaR_jet_Wd1 = -999.9;
            JetGenMatched_DeltaR_jet_Wd2 = -999.9;
            JetGenMatched_DeltaR_pup0_b = -999.9;
            JetGenMatched_DeltaR_pup0_Wd1 = -999.9;
            JetGenMatched_DeltaR_pup0_Wd2 = -999.9;
            JetGenMatched_DeltaR_pup1_b = -999.9;
            JetGenMatched_DeltaR_pup1_Wd1 = -999.9;
            JetGenMatched_DeltaR_pup1_Wd2 = -999.9;
            JetGenMatched_partonPt = -999.9;
            JetGenMatched_partonEta = -999.9;
            JetGenMatched_partonPhi = -999.9;
            JetGenMatched_partonMass = -999.9;
            JetGenMatched_partonID = -999.9;
            JetGenMatched_DeltaRjetParton = -999.9;
            HadMETpx = -999.9;
            HadMETpy = -999.9;
            HadMETpt = -999.9;
            HadMETphi = -999.9;
            HadMETsumET = -999.9;
            HadMETgenMET = -999.9;
            HadMETuncorPt = -999.9;
            HadMETshiftedPtJetEnUp = -999.9;
            HadMETshiftedPtJetEnDn = -999.9;
            HadMETshiftedPtElEnUp = -999.9;
            HadMETshiftedPtElEnDn = -999.9;
            HadMETshiftedPtMuEnUp = -999.9;
            HadMETshiftedPtMuEnDn = -999.9;
            HadMETshiftedPtJetResUp = -999.9;
            HadMETshiftedPtJetResDn = -999.9;
            HadMETshiftedPtUnclEnUp = -999.9;
            HadMETshiftedPtUnclEnDn = -999.9;
            HadNvtx = -999.9;
            HadNvtxGood = -999.9;
            HadNPUtrue = -999.9;
            HadRho = -999.9;
            HadEventWeight = -999.9;
            HadPUweight = -999.9;
            HadPUweight_MBup = -999.9;
            HadPUweight_MBdn = -999.9;
            HadGenTTmass = -999.9;
            HTlep = -999.9;
            ST = -999.9;
            ST_CorrDn = -999.9;
            ST_CorrUp = -999.9;
            ST_PtSmearNom = -999.9;
            ST_PtSmearUp = -999.9;
            ST_PtSmearDn = -999.9;
            HadQ2weight_CorrDn = -999.9;
            HadQ2weight_CorrUp = -999.9;
            HadNNPDF3weight_CorrDn = -999.9;
            HadNNPDF3weight_CorrUp = -999.9;
            
            
            //int
            
            parton1id = 0;
            parton2id = 0;
            Wd1_id = 0 ;
            Wd2_id = 0 ;
            bJet_count = 0;
            
            //Int_t
            
            JetNsubjetsSD = -999;
            JetNsubjetsSDPuppi = -999;
            JetGenMatched_TopHadronic = -999;
            HadGenCountHadTop = -999;
            HadRunNum = -999;
            HadLumiBlock = -999;
            HadEventNum = -999;
            HadPassMETFilters = -999;

      }


};
#endif