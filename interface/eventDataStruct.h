#ifndef eventDataStruct_h
#define eventDataStruct_h
#define BADVAL -999.0

#include <vector>
#include <string>

typedef struct eventDataStruct {

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


};
#endif