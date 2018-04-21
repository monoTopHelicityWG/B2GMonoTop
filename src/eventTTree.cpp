#include "Analysis/B2GMonoTop/interface/eventDataStruct.h"
#include "TTree.h"


void bookEventTree(TTree * event_ttree, eventDataStruct * event_ttree_row){

    event_ttree->Branch("mu_preselection"                         , & event_ttree_row->mu_preselection                         ,  "mu_preselection/O"                      );
    event_ttree->Branch("had_preselection"                         , & event_ttree_row->had_preselection                         ,  "had_preselection/O"                      );

    event_ttree->Branch("HLTtriggers"        , "vector<std::string>", &event_ttree_row->HLTtriggers);
    event_ttree->Branch("HLTtriggersPass"        , "vector<bool>", &event_ttree_row->HLTtriggersPass);
    event_ttree->Branch("HLTtriggersPrescales"        , "vector<int>", &event_ttree_row->HLTtriggersPrescales);
    
    event_ttree->Branch("HadTrigPrescales"   , "vector<int>", &event_ttree_row->HadTrigPrescales);
    event_ttree->Branch("HadTrigPass"        , "vector<bool>", &event_ttree_row->HadTrigPass);
    event_ttree->Branch("HadTrigAcceptBits"  , &event_ttree_row->HadTrigAcceptBits);
    
    event_ttree->Branch("HadTrigPrescalesMu"   , "vector<int>", &event_ttree_row->HadTrigPrescalesMu);
    event_ttree->Branch("HadTrigPassMu"        , "vector<bool>", &event_ttree_row->HadTrigPassMu);
    event_ttree->Branch("HadTrigAcceptBitsMu"  , &event_ttree_row->HadTrigAcceptBitsMu);
    
    event_ttree->Branch("HadTrigPrescalesHad"   , "vector<int>", &event_ttree_row->HadTrigPrescalesHad);
    event_ttree->Branch("HadTrigPassHad"        , "vector<bool>", &event_ttree_row->HadTrigPassHad);
    event_ttree->Branch("HadTrigAcceptBitsHad"  , &event_ttree_row->HadTrigAcceptBitsHad);
                         
    event_ttree->Branch("Gen_MET_pT",         & event_ttree_row->Gen_MET_pT          , "Gen_MET_pT/F"        );
    event_ttree->Branch("Gen_MET_phi",        & event_ttree_row->Gen_MET_phi         , "Gen_MET_phi/F"       );
    event_ttree->Branch("Gen_MET_eta",        & event_ttree_row->Gen_MET_eta         , "Gen_MET_eta/F"       );
    
    event_ttree->Branch("Gen_array_t_p4"                   , & event_ttree_row->Gen_array_t_p4                   ,  "Gen_array_t_p4[4]/F"                      );
    event_ttree->Branch("Gen_array_final_t_p4"             , & event_ttree_row->Gen_array_final_t_p4             ,  "Gen_array_final_t_p4[4]/F"                );
    event_ttree->Branch("Gen_array_b_p4"                   , & event_ttree_row->Gen_array_b_p4                   ,  "Gen_array_b_p4[4]/F"                      );
    event_ttree->Branch("Gen_array_W_p4"                   , & event_ttree_row->Gen_array_W_p4                   ,  "Gen_array_W_p4[4]/F"                      );
    event_ttree->Branch("Gen_array_Wd1_p4"                 , & event_ttree_row->Gen_array_Wd1_p4                 ,  "Gen_array_Wd1_p4[4]/F"                    );
    event_ttree->Branch("Gen_array_Wd2_p4"                 , & event_ttree_row->Gen_array_Wd2_p4                 ,  "Gen_array_Wd2_p4[4]/F"                    );
    event_ttree->Branch("Gen_array_hardest_parton_hardScatterOutgoing_p4", & event_ttree_row->Gen_array_hardest_parton_hardScatterOutgoing_p4 ,  "Gen_array_hardest_parton_hardScatterOutgoing_p4[4]/F");
    event_ttree->Branch("Gen_array_second_hardest_parton_hardScatterOutgoing_p4", & event_ttree_row->Gen_array_second_hardest_parton_hardScatterOutgoing_p4 ,  "Gen_array_second_hardest_parton_hardScatterOutgoing_p4[4]/F" );
    event_ttree->Branch("tophadronic"                         , & event_ttree_row->tophadronic                         ,  "tophadronic/O"                      );
    event_ttree->Branch("topleptonic"                         , & event_ttree_row->topleptonic                         ,  "topleptonic/O"                      );
    event_ttree->Branch("parton1id"                           , & event_ttree_row->parton1id                           ,  "parton1id/I"                        );
    event_ttree->Branch("parton2id"                           , & event_ttree_row->parton2id                           ,  "parton2id/I"                        );
    event_ttree->Branch("Wd1_id"                              , & event_ttree_row->Wd1_id                              ,  "Wd1_id/I"                           );
    event_ttree->Branch("Wd2_id"                              , & event_ttree_row->Wd2_id                              ,  "Wd2_id/I"                           ); 
    
    event_ttree->Branch("AK4JetLV_pt","vector<float>",&event_ttree_row->AK4JetLV_pt);
    event_ttree->Branch("AK4JetLV_eta","vector<float>",&event_ttree_row->AK4JetLV_eta);
    event_ttree->Branch("AK4JetLV_phi","vector<float>",&event_ttree_row->AK4JetLV_phi);
    event_ttree->Branch("AK4JetLV_mass","vector<float>",&event_ttree_row->AK4JetLV_mass);
    
    event_ttree->Branch("AK4JetLV_corr", "vector<float>", &event_ttree_row->AK4JetLV_corr);
    event_ttree->Branch("AK4JetLV_corrUp", "vector<float>", &event_ttree_row->AK4JetLV_corrUp);
    event_ttree->Branch("AK4JetLV_corrDn", "vector<float>", &event_ttree_row->AK4JetLV_corrDn);
    event_ttree->Branch("AK4JetLV_SF", "vector<float>", &event_ttree_row->AK4JetLV_SF);
    event_ttree->Branch("AK4JetLV_SFUp", "vector<float>", &event_ttree_row->AK4JetLV_SFUp);
    event_ttree->Branch("AK4JetLV_SFDn", "vector<float>", &event_ttree_row->AK4JetLV_SFDn);
    event_ttree->Branch("AK4JetLV_ptsmear", "vector<float>", &event_ttree_row->AK4JetLV_ptsmear);
    event_ttree->Branch("AK4JetLV_ptsmearUp", "vector<float>", &event_ttree_row->AK4JetLV_ptsmearUp);
    event_ttree->Branch("AK4JetLV_ptsmearDn", "vector<float>", &event_ttree_row->AK4JetLV_ptsmearDn);
    
    event_ttree->Branch("AK8JetLV_corr", "vector<float>", &event_ttree_row->AK8JetLV_corr);
    event_ttree->Branch("AK8JetLV_corrUp", "vector<float>", &event_ttree_row->AK8JetLV_corrUp);
    event_ttree->Branch("AK8JetLV_corrDn", "vector<float>", &event_ttree_row->AK8JetLV_corrDn);
    event_ttree->Branch("AK8JetLV_SF", "vector<float>", &event_ttree_row->AK8JetLV_SF);
    event_ttree->Branch("AK8JetLV_SFUp", "vector<float>", &event_ttree_row->AK8JetLV_SFUp);
    event_ttree->Branch("AK8JetLV_SFDn", "vector<float>", &event_ttree_row->AK8JetLV_SFDn);
    event_ttree->Branch("AK8JetLV_ptsmear", "vector<float>", &event_ttree_row->AK8JetLV_ptsmear);
    event_ttree->Branch("AK8JetLV_ptsmearUp", "vector<float>", &event_ttree_row->AK8JetLV_ptsmearUp);
    event_ttree->Branch("AK8JetLV_ptsmearDn", "vector<float>", &event_ttree_row->AK8JetLV_ptsmearDn);
    
    event_ttree->Branch("AK8JetLV_pt","vector<float>",&event_ttree_row->AK8JetLV_pt);
    event_ttree->Branch("AK8JetLV_eta","vector<float>",&event_ttree_row->AK8JetLV_eta);
    event_ttree->Branch("AK8JetLV_phi","vector<float>",&event_ttree_row->AK8JetLV_phi);
    event_ttree->Branch("AK8JetLV_mass","vector<float>",&event_ttree_row->AK8JetLV_mass);
    
    event_ttree->Branch("AK8SubjetLV_pt","vector<float>",&event_ttree_row->AK8SubjetLV_pt);
    event_ttree->Branch("AK8SubjetLV_eta","vector<float>",&event_ttree_row->AK8SubjetLV_eta);
    event_ttree->Branch("AK8SubjetLV_phi","vector<float>",&event_ttree_row->AK8SubjetLV_phi);
    event_ttree->Branch("AK8SubjetLV_mass","vector<float>",&event_ttree_row->AK8SubjetLV_mass);
    
    event_ttree->Branch("AK4JetBtag_p","vector<float>",&event_ttree_row->AK4JetBtag_p);
    event_ttree->Branch("AK8JetTau1_p","vector<float>",&event_ttree_row->AK8JetTau1_p);
    event_ttree->Branch("AK8JetTau2_p","vector<float>",&event_ttree_row->AK8JetTau2_p);
    event_ttree->Branch("AK8JetTau3_p","vector<float>",&event_ttree_row->AK8JetTau3_p);
    event_ttree->Branch("AK8JetSoftdropMass_p","vector<float>",&event_ttree_row->AK8JetSoftdropMass_p);
    
    event_ttree->Branch("JetPtRaw"                             , & event_ttree_row->JetPtRaw                          ,    "JetPtRaw/F"                               );                                  
    event_ttree->Branch("JetEtaRaw"                            , & event_ttree_row->JetEtaRaw                         ,    "JetEtaRaw/F"                              );                                   
    event_ttree->Branch("JetPhiRaw"                            , & event_ttree_row->JetPhiRaw                         ,    "JetPhiRaw/F"                              );                                   
    event_ttree->Branch("JetMassRaw"                           , & event_ttree_row->JetMassRaw                        ,    "JetMassRaw/F"                             );                                                                                      
    event_ttree->Branch("JetArea"                              , & event_ttree_row->JetArea                           ,    "JetArea/F"                                );                                 
                                        
    event_ttree->Branch("JetSDmassRaw"                         , & event_ttree_row->JetSDmassRaw                      ,    "JetSDmassRaw/F"                           );                                               
    event_ttree->Branch("JetSDmassSubjetCorrL23"                     , & event_ttree_row->JetSDmassSubjetCorrL23                  ,    "JetSDmassSubjetCorrL23/F"                       );                                                                                                     
    event_ttree->Branch("JetSDmassSubjetCorrL123"                    , & event_ttree_row->JetSDmassSubjetCorrL123                 ,    "JetSDmassSubjetCorrL123/F"                      );                                                      
    event_ttree->Branch("JetSDptRaw"                           , & event_ttree_row->JetSDptRaw                        ,    "JetSDptRaw/F"                             );                                                                                                 
    event_ttree->Branch("JetSDetaRaw"                          , & event_ttree_row->JetSDetaRaw                       ,    "JetSDetaRaw/F"                            );                                               
    event_ttree->Branch("JetSDphiRaw"                          , & event_ttree_row->JetSDphiRaw                       ,    "JetSDphiRaw/F"                            );  
    
    event_ttree->Branch("JetMassPruned"                        , & event_ttree_row->JetMassPruned                     ,    "JetMassPruned/F"                          );                                       
    event_ttree->Branch("JetMassTrimmed"                       , & event_ttree_row->JetMassTrimmed                    ,    "JetMassTrimmed/F"                         );                                       
    event_ttree->Branch("JetTau1"                              , & event_ttree_row->JetTau1                           ,    "JetTau1/F"                                );                                 
    event_ttree->Branch("JetTau2"                              , & event_ttree_row->JetTau2                           ,    "JetTau2/F"                                );                                 
    event_ttree->Branch("JetTau3"                              , & event_ttree_row->JetTau3                           ,    "JetTau3/F"                                );                                 
    event_ttree->Branch("JetTau4"                              , & event_ttree_row->JetTau4                           ,    "JetTau4/F"                                );                                 
    event_ttree->Branch("JetTau32"                             , & event_ttree_row->JetTau32                          ,    "JetTau32/F"                               );                                  
    event_ttree->Branch("JetTau21"                             , & event_ttree_row->JetTau21                          ,    "JetTau21/F"                               );                                  
    event_ttree->Branch("JetSDmaxbdisc"                        , & event_ttree_row->JetSDmaxbdisc                     ,    "JetSDmaxbdisc/F"                          );                                       
    event_ttree->Branch("JetSDmaxbdiscflavHadron"              , & event_ttree_row->JetSDmaxbdiscflavHadron           ,    "JetSDmaxbdiscflavHadron/F"                );                                           
    event_ttree->Branch("JetSDmaxbdiscflavParton"              , & event_ttree_row->JetSDmaxbdiscflavParton           ,    "JetSDmaxbdiscflavParton/F"                );  
                                             
    event_ttree->Branch("JetSDsubjet0pt"                       , & event_ttree_row->JetSDsubjet0pt                    ,    "JetSDsubjet0pt/F"                         );    
    event_ttree->Branch("JetSDsubjet0mass"                     , & event_ttree_row->JetSDsubjet0mass                  ,    "JetSDsubjet0mass/F"                       );
    event_ttree->Branch("JetSDsubjet0eta"                      , & event_ttree_row->JetSDsubjet0eta                   ,    "JetSDsubjet0eta/F"                        );
    event_ttree->Branch("JetSDsubjet0phi"                      , & event_ttree_row->JetSDsubjet0phi                   ,    "JetSDsubjet0phi/F"                        );
    event_ttree->Branch("JetSDsubjet0area"                     , & event_ttree_row->JetSDsubjet0area                  ,    "JetSDsubjet0area/F"                       );
    event_ttree->Branch("JetSDsubjet0flavHadron"               , & event_ttree_row->JetSDsubjet0flavHadron            ,    "JetSDsubjet0flavHadron/F"                 );
    event_ttree->Branch("JetSDsubjet0flavParton"               , & event_ttree_row->JetSDsubjet0flavParton            ,    "JetSDsubjet0flavParton/F"                 );
    event_ttree->Branch("JetSDsubjet0matchedgenjetpt"          , & event_ttree_row->JetSDsubjet0matchedgenjetpt       ,    "JetSDsubjet0matchedgenjetpt/F"            ); 
    event_ttree->Branch("JetSDsubjet0tau1"                     , & event_ttree_row->JetSDsubjet0tau1                  ,    "JetSDsubjet0tau1/F"                       );
    event_ttree->Branch("JetSDsubjet0tau2"                     , & event_ttree_row->JetSDsubjet0tau2                  ,    "JetSDsubjet0tau2/F"                       );
    event_ttree->Branch("JetSDsubjet0tau3"                     , & event_ttree_row->JetSDsubjet0tau3                  ,    "JetSDsubjet0tau3/F"                       ); 
    event_ttree->Branch("JetSDsubjet0bdisc"                    , & event_ttree_row->JetSDsubjet0bdisc                 ,    "JetSDsubjet0bdisc/F"                      );                                          
    event_ttree->Branch("JetSDsubjet1pt"                       , & event_ttree_row->JetSDsubjet1pt                    ,    "JetSDsubjet1pt/F"                         );    
    event_ttree->Branch("JetSDsubjet1mass"                     , & event_ttree_row->JetSDsubjet1mass                  ,    "JetSDsubjet1mass/F"                       );
    event_ttree->Branch("JetSDsubjet1eta"                      , & event_ttree_row->JetSDsubjet1eta                   ,    "JetSDsubjet1eta/F"                        );
    event_ttree->Branch("JetSDsubjet1phi"                      , & event_ttree_row->JetSDsubjet1phi                   ,    "JetSDsubjet1phi/F"                        );  
    event_ttree->Branch("JetSDsubjet1area"                     , & event_ttree_row->JetSDsubjet1area                  ,    "JetSDsubjet1area/F"                       );
    event_ttree->Branch("JetSDsubjet1flavHadron"               , & event_ttree_row->JetSDsubjet1flavHadron            ,    "JetSDsubjet1flavHadron/F"                 );
    event_ttree->Branch("JetSDsubjet1flavParton"               , & event_ttree_row->JetSDsubjet1flavParton            ,    "JetSDsubjet1flavParton/F"                 );
    event_ttree->Branch("JetSDsubjet1matchedgenjetpt"          , & event_ttree_row->JetSDsubjet1matchedgenjetpt       ,    "JetSDsubjet1matchedgenjetpt/F"            ); 
    event_ttree->Branch("JetSDsubjet1tau1"                     , & event_ttree_row->JetSDsubjet1tau1                  ,    "JetSDsubjet1tau1/F"                       );
    event_ttree->Branch("JetSDsubjet1tau2"                     , & event_ttree_row->JetSDsubjet1tau2                  ,    "JetSDsubjet1tau2/F"                       );
    event_ttree->Branch("JetSDsubjet1tau3"                     , & event_ttree_row->JetSDsubjet1tau3                  ,    "JetSDsubjet1tau3/F"                       );                                           
    event_ttree->Branch("JetSDsubjet1bdisc"                    , & event_ttree_row->JetSDsubjet1bdisc                 ,    "JetSDsubjet1bdisc/F"                      );                                     
    
    event_ttree->Branch("JetPuppiPtRaw"                           , & event_ttree_row->JetPuppiPtRaw                        ,    "JetPuppiPtRaw/F"                             );                                    
    event_ttree->Branch("JetPuppiEtaRaw"                          , & event_ttree_row->JetPuppiEtaRaw                       ,    "JetPuppiEtaRaw/F"                            );                                     
    event_ttree->Branch("JetPuppiPhiRaw"                          , & event_ttree_row->JetPuppiPhiRaw                       ,    "JetPuppiPhiRaw/F"                            );                                     
    event_ttree->Branch("JetPuppiMassRaw"                         , & event_ttree_row->JetPuppiMassRaw                      ,    "JetPuppiMassRaw/F"                           );                                      
    event_ttree->Branch("JetPuppiArea"                         , & event_ttree_row->JetPuppiArea                      ,    "JetPuppiArea/F"                           );                                      
    
    event_ttree->Branch("JetPuppiMassPruned"                    , & event_ttree_row->JetPuppiMassPruned               ,   "JetPuppiMassPruned/F"                          );
    event_ttree->Branch("JetPuppiMassTrimmed"                   , & event_ttree_row->JetPuppiMassTrimmed              ,   "JetPuppiMassTrimmed/F"                         );
    
    event_ttree->Branch("JetPuppiSDmassRaw"                         , & event_ttree_row->JetPuppiSDmassRaw                    ,    "JetPuppiSDmassRaw/F"                          );
    event_ttree->Branch("JetPuppiSDmassSubjetCorr"                     , & event_ttree_row->JetPuppiSDmassSubjetCorr                ,    "JetPuppiSDmassSubjetCorr/F"                      );
    event_ttree->Branch("JetPuppiSDptRaw"                           , & event_ttree_row->JetPuppiSDptRaw                      ,    "JetPuppiSDptRaw/F"                            );
    event_ttree->Branch("JetPuppiSDetaRaw"                          , & event_ttree_row->JetPuppiSDetaRaw                     ,    "JetPuppiSDetaRaw/F"                           );
    event_ttree->Branch("JetPuppiSDphiRaw"                          , & event_ttree_row->JetPuppiSDphiRaw                     ,    "JetPuppiSDphiRaw/F"                           );
                            
    event_ttree->Branch("JetPuppiTau1"                         , & event_ttree_row->JetPuppiTau1                      ,    "JetPuppiTau1/F"                           );                                      
    event_ttree->Branch("JetPuppiTau2"                         , & event_ttree_row->JetPuppiTau2                      ,    "JetPuppiTau2/F"                           );                                      
    event_ttree->Branch("JetPuppiTau3"                         , & event_ttree_row->JetPuppiTau3                      ,    "JetPuppiTau3/F"                           );                                      
    event_ttree->Branch("JetPuppiTau4"                         , & event_ttree_row->JetPuppiTau4                      ,    "JetPuppiTau4/F"                           );                                      
    event_ttree->Branch("JetPuppiTau32"                        , & event_ttree_row->JetPuppiTau32                     ,    "JetPuppiTau32/F"                          );                                       
    event_ttree->Branch("JetPuppiTau21"                        , & event_ttree_row->JetPuppiTau21                     ,    "JetPuppiTau21/F"                          );                                       
    
    event_ttree->Branch("JetPuppiSDmaxbdisc"                   , & event_ttree_row->JetPuppiSDmaxbdisc                ,    "JetPuppiSDmaxbdisc/F"                     );                                            
    event_ttree->Branch("JetPuppiSDmaxbdiscflavHadron"         , & event_ttree_row->JetPuppiSDmaxbdiscflavHadron      ,    "JetPuppiSDmaxbdiscflavHadron/F"           );                                                
    event_ttree->Branch("JetPuppiSDmaxbdiscflavParton"         , & event_ttree_row->JetPuppiSDmaxbdiscflavParton      ,    "JetPuppiSDmaxbdiscflavParton/F"           );                                                
    event_ttree->Branch("JetPuppiSDsubjet0pt"                  , & event_ttree_row->JetPuppiSDsubjet0pt               ,    "JetPuppiSDsubjet0pt/F"                    );    
    event_ttree->Branch("JetPuppiSDsubjet0mass"                , & event_ttree_row->JetPuppiSDsubjet0mass             ,    "JetPuppiSDsubjet0mass/F"                  );
    event_ttree->Branch("JetPuppiSDsubjet0eta"                 , & event_ttree_row->JetPuppiSDsubjet0eta              ,    "JetPuppiSDsubjet0eta/F"                   );
    event_ttree->Branch("JetPuppiSDsubjet0phi"                 , & event_ttree_row->JetPuppiSDsubjet0phi              ,    "JetPuppiSDsubjet0phi/F"                   );
    event_ttree->Branch("JetPuppiSDsubjet0area"                , & event_ttree_row->JetPuppiSDsubjet0area             ,    "JetPuppiSDsubjet0area/F"                  );
    event_ttree->Branch("JetPuppiSDsubjet0flavHadron"          , & event_ttree_row->JetPuppiSDsubjet0flavHadron       ,    "JetPuppiSDsubjet0flavHadron/F"            );
    event_ttree->Branch("JetPuppiSDsubjet0flavParton"          , & event_ttree_row->JetPuppiSDsubjet0flavParton       ,    "JetPuppiSDsubjet0flavParton/F"            );
    event_ttree->Branch("JetPuppiSDsubjet0matchedgenjetpt"     , & event_ttree_row->JetPuppiSDsubjet0matchedgenjetpt  ,    "JetPuppiSDsubjet0matchedgenjetpt/F"       );
    event_ttree->Branch("JetPuppiSDsubjet0tau1"                , & event_ttree_row->JetPuppiSDsubjet0tau1             ,    "JetPuppiSDsubjet0tau1/F"                  );
    event_ttree->Branch("JetPuppiSDsubjet0tau2"                , & event_ttree_row->JetPuppiSDsubjet0tau2             ,    "JetPuppiSDsubjet0tau2/F"                  );
    event_ttree->Branch("JetPuppiSDsubjet0tau3"                , & event_ttree_row->JetPuppiSDsubjet0tau3             ,    "JetPuppiSDsubjet0tau3/F"                  );
    event_ttree->Branch("JetPuppiSDsubjet0bdisc"               , & event_ttree_row->JetPuppiSDsubjet0bdisc            ,    "JetPuppiSDsubjet0bdisc/F"                 );                                                
    event_ttree->Branch("JetPuppiSDsubjet1pt"                  , & event_ttree_row->JetPuppiSDsubjet1pt               ,    "JetPuppiSDsubjet1pt/F"                    );    
    event_ttree->Branch("JetPuppiSDsubjet1mass"                , & event_ttree_row->JetPuppiSDsubjet1mass             ,    "JetPuppiSDsubjet1mass/F"                  );
    event_ttree->Branch("JetPuppiSDsubjet1eta"                 , & event_ttree_row->JetPuppiSDsubjet1eta              ,    "JetPuppiSDsubjet1eta/F"                   );
    event_ttree->Branch("JetPuppiSDsubjet1phi"                 , & event_ttree_row->JetPuppiSDsubjet1phi              ,    "JetPuppiSDsubjet1phi/F"                   );  
    event_ttree->Branch("JetPuppiSDsubjet1area"                , & event_ttree_row->JetPuppiSDsubjet1area             ,    "JetPuppiSDsubjet1area/F"                  );
    event_ttree->Branch("JetPuppiSDsubjet1flavHadron"          , & event_ttree_row->JetPuppiSDsubjet1flavHadron       ,    "JetPuppiSDsubjet1flavHadron/F"            );
    event_ttree->Branch("JetPuppiSDsubjet1flavParton"          , & event_ttree_row->JetPuppiSDsubjet1flavParton       ,    "JetPuppiSDsubjet1flavParton/F"            );
    event_ttree->Branch("JetPuppiSDsubjet1matchedgenjetpt"     , & event_ttree_row->JetPuppiSDsubjet1matchedgenjetpt  ,    "JetPuppiSDsubjet1matchedgenjetpt/F"       );
    event_ttree->Branch("JetPuppiSDsubjet1tau1"                , & event_ttree_row->JetPuppiSDsubjet1tau1             ,    "JetPuppiSDsubjet1tau1/F"                  );
    event_ttree->Branch("JetPuppiSDsubjet1tau2"                , & event_ttree_row->JetPuppiSDsubjet1tau2             ,    "JetPuppiSDsubjet1tau2/F"                  );
    event_ttree->Branch("JetPuppiSDsubjet1tau3"                , & event_ttree_row->JetPuppiSDsubjet1tau3             ,    "JetPuppiSDsubjet1tau3/F"                  );  
    event_ttree->Branch("JetPuppiSDECF1"                       , & event_ttree_row->JetPuppiSDECF1                    ,    "JetPuppiSDECF1/F"                         );
    event_ttree->Branch("JetPuppiSDECF2"                       , & event_ttree_row->JetPuppiSDECF2                    ,    "JetPuppiSDECF2/F"                         );
    event_ttree->Branch("JetPuppiSDECF3"                       , & event_ttree_row->JetPuppiSDECF3                    ,    "JetPuppiSDECF3/F"                         );
    event_ttree->Branch("JetPuppiSDECF4"                       , & event_ttree_row->JetPuppiSDECF4                    ,    "JetPuppiSDECF4/F"                         );
    event_ttree->Branch("JetPuppiSDECF5"                       , & event_ttree_row->JetPuppiSDECF5                    ,    "JetPuppiSDECF5/F"                         );
    event_ttree->Branch("JetPuppiSDC_2"                       , & event_ttree_row->JetPuppiSDC_2                    ,    "JetPuppiSDC_2/F"                         );
    event_ttree->Branch("JetPuppiSDD_2"                       , & event_ttree_row->JetPuppiSDD_2                    ,    "JetPuppiSDD_2/F"                         );
    event_ttree->Branch("JetPuppiSDC_3"                       , & event_ttree_row->JetPuppiSDC_3                    ,    "JetPuppiSDC_3/F"                         );
    event_ttree->Branch("JetPuppiSDD_3"                       , & event_ttree_row->JetPuppiSDD_3                    ,    "JetPuppiSDD_3/F"                         );
    event_ttree->Branch("JetPuppiSDsubjet1bdisc"               , & event_ttree_row->JetPuppiSDsubjet1bdisc            ,    "JetPuppiSDsubjet1bdisc/F"                 );                                                                                                                        
    
    event_ttree->Branch("JetCHF"                               , & event_ttree_row->JetCHF                            ,    "JetCHF/F"                                 );                                
    event_ttree->Branch("JetNHF"                               , & event_ttree_row->JetNHF                            ,    "JetNHF/F"                                 );                                
    event_ttree->Branch("JetCM"                                , & event_ttree_row->JetCM                             ,    "JetCM/F"                                  );                               
    event_ttree->Branch("JetNM"                                , & event_ttree_row->JetNM                             ,    "JetNM/F"                                  );                               
    event_ttree->Branch("JetNEF"                               , & event_ttree_row->JetNEF                            ,    "JetNEF/F"                                 );                                
    event_ttree->Branch("JetCEF"                               , & event_ttree_row->JetCEF                            ,    "JetCEF/F"                                 );                                
    event_ttree->Branch("JetMF"                                , & event_ttree_row->JetMF                             ,    "JetMF/F"                                  );                               
    event_ttree->Branch("JetMult"                              , & event_ttree_row->JetMult                           ,    "JetMult/F"                                );
    event_ttree->Branch("JetPuppiCHF"                          , & event_ttree_row->JetPuppiCHF                       ,    "JetPuppiCHF/F"                            );                                
    event_ttree->Branch("JetPuppiNHF"                          , & event_ttree_row->JetPuppiNHF                       ,    "JetPuppiNHF/F"                            );                                
    event_ttree->Branch("JetPuppiCM"                           , & event_ttree_row->JetPuppiCM                        ,    "JetPuppiCM/F"                             );                               
    event_ttree->Branch("JetPuppiNM"                           , & event_ttree_row->JetPuppiNM                        ,    "JetPuppiNM/F"                             );                               
    event_ttree->Branch("JetPuppiNEF"                          , & event_ttree_row->JetPuppiNEF                       ,    "JetPuppiNEF/F"                            );                                
    event_ttree->Branch("JetPuppiCEF"                          , & event_ttree_row->JetPuppiCEF                       ,    "JetPuppiCEF/F"                            );                                
    event_ttree->Branch("JetPuppiMF"                           , & event_ttree_row->JetPuppiMF                        ,    "JetPuppiMF/F"                             );                               
    event_ttree->Branch("JetPuppiMult"                         , & event_ttree_row->JetPuppiMult                      ,    "JetPuppiMult/F"                           );                                  
    event_ttree->Branch("JetMassCorrFactor"                    , & event_ttree_row->JetMassCorrFactor                 ,    "JetMassCorrFactor/F"                      );                                           
    event_ttree->Branch("JetMassCorrFactorUp"                  , & event_ttree_row->JetMassCorrFactorUp               ,    "JetMassCorrFactorUp/F"                    );                                             
    event_ttree->Branch("JetMassCorrFactorDn"                  , & event_ttree_row->JetMassCorrFactorDn               ,    "JetMassCorrFactorDn/F"                    );                                             
    event_ttree->Branch("JetCorrFactor"                        , & event_ttree_row->JetCorrFactor                     ,    "JetCorrFactor/F"                          );                                       
    event_ttree->Branch("JetCorrFactorUp"                      , & event_ttree_row->JetCorrFactorUp                   ,    "JetCorrFactorUp/F"                        );                                         
    event_ttree->Branch("JetCorrFactorDn"                      , & event_ttree_row->JetCorrFactorDn                   ,    "JetCorrFactorDn/F"                        );                                         
    event_ttree->Branch("JetPtSmearFactor"                     , & event_ttree_row->JetPtSmearFactor                  ,    "JetPtSmearFactor/F"                       );                                          
    event_ttree->Branch("JetPtSmearFactorUp"                   , & event_ttree_row->JetPtSmearFactorUp                ,    "JetPtSmearFactorUp/F"                     );                                            
    event_ttree->Branch("JetPtSmearFactorDn"                   , & event_ttree_row->JetPtSmearFactorDn                ,    "JetPtSmearFactorDn/F"                     );                                            
    event_ttree->Branch("JetPuppiMassCorrFactor"               , & event_ttree_row->JetPuppiMassCorrFactor            ,    "JetPuppiMassCorrFactor/F"                 );                                                
    event_ttree->Branch("JetPuppiMassCorrFactorUp"             , & event_ttree_row->JetPuppiMassCorrFactorUp          ,    "JetPuppiMassCorrFactorUp/F"               );                                                  
    event_ttree->Branch("JetPuppiMassCorrFactorDn"             , & event_ttree_row->JetPuppiMassCorrFactorDn          ,    "JetPuppiMassCorrFactorDn/F"               );                                                  
    event_ttree->Branch("JetPuppiCorrFactor"                   , & event_ttree_row->JetPuppiCorrFactor                ,    "JetPuppiCorrFactor/F"                     );                                            
    event_ttree->Branch("JetPuppiCorrFactorUp"                 , & event_ttree_row->JetPuppiCorrFactorUp              ,    "JetPuppiCorrFactorUp/F"                   );                                              
    event_ttree->Branch("JetPuppiCorrFactorDn"                 , & event_ttree_row->JetPuppiCorrFactorDn              ,    "JetPuppiCorrFactorDn/F"                   );                                              
    event_ttree->Branch("JetPuppiPtSmearFactor"                , & event_ttree_row->JetPuppiPtSmearFactor             ,    "JetPuppiPtSmearFactor/F"                  );                                               
    event_ttree->Branch("JetPuppiPtSmearFactorUp"              , & event_ttree_row->JetPuppiPtSmearFactorUp           ,    "JetPuppiPtSmearFactorUp/F"                );                                                 
    event_ttree->Branch("JetPuppiPtSmearFactorDn"              , & event_ttree_row->JetPuppiPtSmearFactorDn           ,    "JetPuppiPtSmearFactorDn/F"                );                                                 
    event_ttree->Branch("JetMatchedGenJetPt"                   , & event_ttree_row->JetMatchedGenJetPt                ,    "JetMatchedGenJetPt/F"                     );                                            
    event_ttree->Branch("JetMatchedGenJetMass"                 , & event_ttree_row->JetMatchedGenJetMass              ,    "JetMatchedGenJetMass/F"                   ); 
    event_ttree->Branch("JetPuppiMatchedGenJetPt"              , & event_ttree_row->JetPuppiMatchedGenJetPt           ,    "JetPuppiMatchedGenJetPt/F"                );                                            
    event_ttree->Branch("JetPuppiMatchedGenJetMass"            , & event_ttree_row->JetPuppiMatchedGenJetMass         ,    "JetPuppiMatchedGenJetMass/F"              ); 
                               
    event_ttree->Branch("JetGenMatched_TopHadronic"            , & event_ttree_row->JetGenMatched_TopHadronic         ,    "JetGenMatched_TopHadronic/I"              );      
    event_ttree->Branch("JetGenMatched_TopPt"                  , & event_ttree_row->JetGenMatched_TopPt               ,    "JetGenMatched_TopPt/F"                    );      
    event_ttree->Branch("JetGenMatched_TopEta"                 , & event_ttree_row->JetGenMatched_TopEta              ,    "JetGenMatched_TopEta/F"                   );      
    event_ttree->Branch("JetGenMatched_TopPhi"                 , & event_ttree_row->JetGenMatched_TopPhi              ,    "JetGenMatched_TopPhi/F"                   );      
    event_ttree->Branch("JetGenMatched_TopMass"                , & event_ttree_row->JetGenMatched_TopMass             ,    "JetGenMatched_TopMass/F"                  );      
    event_ttree->Branch("JetGenMatched_bPt"                    , & event_ttree_row->JetGenMatched_bPt                 ,    "JetGenMatched_bPt/F"                      );      
    event_ttree->Branch("JetGenMatched_WPt"                    , & event_ttree_row->JetGenMatched_WPt                 ,    "JetGenMatched_WPt/F"                      );      
    event_ttree->Branch("JetGenMatched_Wd1Pt"                  , & event_ttree_row->JetGenMatched_Wd1Pt               ,    "JetGenMatched_Wd1Pt/F"                    );      
    event_ttree->Branch("JetGenMatched_Wd2Pt"                  , & event_ttree_row->JetGenMatched_Wd2Pt               ,    "JetGenMatched_Wd2Pt/F"                    );      
    event_ttree->Branch("JetGenMatched_Wd1ID"                  , & event_ttree_row->JetGenMatched_Wd1ID               ,    "JetGenMatched_Wd1ID/F"                    );      
    event_ttree->Branch("JetGenMatched_Wd2ID"                  , & event_ttree_row->JetGenMatched_Wd2ID               ,    "JetGenMatched_Wd2ID/F"                    );      
    event_ttree->Branch("JetGenMatched_MaxDeltaRPartonTop"     , & event_ttree_row->JetGenMatched_MaxDeltaRPartonTop  ,    "JetGenMatched_MaxDeltaRPartonTop/F"       );      
    event_ttree->Branch("JetGenMatched_MaxDeltaRWPartonTop"    , & event_ttree_row->JetGenMatched_MaxDeltaRWPartonTop ,    "JetGenMatched_MaxDeltaRWPartonTop/F"      );      
    event_ttree->Branch("JetGenMatched_MaxDeltaRWPartonW"      , & event_ttree_row->JetGenMatched_MaxDeltaRWPartonW   ,    "JetGenMatched_MaxDeltaRWPartonW/F"        );      
    event_ttree->Branch("JetGenMatched_DeltaR_t_b"             , & event_ttree_row->JetGenMatched_DeltaR_t_b          ,    "JetGenMatched_DeltaR_t_b/F"               );      
    event_ttree->Branch("JetGenMatched_DeltaR_t_W"             , & event_ttree_row->JetGenMatched_DeltaR_t_W          ,    "JetGenMatched_DeltaR_t_W/F"               );      
    event_ttree->Branch("JetGenMatched_DeltaR_t_Wd1"           , & event_ttree_row->JetGenMatched_DeltaR_t_Wd1        ,    "JetGenMatched_DeltaR_t_Wd1/F"             );      
    event_ttree->Branch("JetGenMatched_DeltaR_t_Wd2"           , & event_ttree_row->JetGenMatched_DeltaR_t_Wd2        ,    "JetGenMatched_DeltaR_t_Wd2/F"             );      
    event_ttree->Branch("JetGenMatched_DeltaR_W_b1"            , & event_ttree_row->JetGenMatched_DeltaR_W_b1         ,    "JetGenMatched_DeltaR_W_b1/F"              );      
    event_ttree->Branch("JetGenMatched_DeltaR_W_Wd1"           , & event_ttree_row->JetGenMatched_DeltaR_W_Wd1        ,    "JetGenMatched_DeltaR_W_Wd1/F"             );      
    event_ttree->Branch("JetGenMatched_DeltaR_W_Wd2"           , & event_ttree_row->JetGenMatched_DeltaR_W_Wd2        ,    "JetGenMatched_DeltaR_W_Wd2/F"             );      
    event_ttree->Branch("JetGenMatched_DeltaR_Wd1_Wd2"         , & event_ttree_row->JetGenMatched_DeltaR_Wd1_Wd2      ,    "JetGenMatched_DeltaR_Wd1_Wd2/F"           );      
    event_ttree->Branch("JetGenMatched_DeltaR_Wd1_b"           , & event_ttree_row->JetGenMatched_DeltaR_Wd1_b        ,    "JetGenMatched_DeltaR_Wd1_b/F"             );      
    event_ttree->Branch("JetGenMatched_DeltaR_Wd2_b"           , & event_ttree_row->JetGenMatched_DeltaR_Wd2_b        ,    "JetGenMatched_DeltaR_Wd2_b/F"             );      
    event_ttree->Branch("JetGenMatched_DeltaR_jet_t"           , & event_ttree_row->JetGenMatched_DeltaR_jet_t        ,    "JetGenMatched_DeltaR_jet_t/F"             );      
    event_ttree->Branch("JetGenMatched_DeltaR_jet_W"           , & event_ttree_row->JetGenMatched_DeltaR_jet_W        ,    "JetGenMatched_DeltaR_jet_W/F"             );      
    event_ttree->Branch("JetGenMatched_DeltaR_jet_b"           , & event_ttree_row->JetGenMatched_DeltaR_jet_b        ,    "JetGenMatched_DeltaR_jet_b/F"             );      
    event_ttree->Branch("JetGenMatched_DeltaR_jet_Wd1"         , & event_ttree_row->JetGenMatched_DeltaR_jet_Wd1      ,    "JetGenMatched_DeltaR_jet_Wd1/F"           );      
    event_ttree->Branch("JetGenMatched_DeltaR_jet_Wd2"         , & event_ttree_row->JetGenMatched_DeltaR_jet_Wd2      ,    "JetGenMatched_DeltaR_jet_Wd2/F"           );      
    event_ttree->Branch("JetGenMatched_DeltaR_pup0_b"          , & event_ttree_row->JetGenMatched_DeltaR_pup0_b       ,    "JetGenMatched_DeltaR_pup0_b/F"            );      
    event_ttree->Branch("JetGenMatched_DeltaR_pup0_Wd1"        , & event_ttree_row->JetGenMatched_DeltaR_pup0_Wd1     ,    "JetGenMatched_DeltaR_pup0_Wd1/F"          );      
    event_ttree->Branch("JetGenMatched_DeltaR_pup0_Wd2"        , & event_ttree_row->JetGenMatched_DeltaR_pup0_Wd2     ,    "JetGenMatched_DeltaR_pup0_Wd2/F"          );      
    event_ttree->Branch("JetGenMatched_DeltaR_pup1_b"          , & event_ttree_row->JetGenMatched_DeltaR_pup1_b       ,    "JetGenMatched_DeltaR_pup1_b/F"            );      
    event_ttree->Branch("JetGenMatched_DeltaR_pup1_Wd1"        , & event_ttree_row->JetGenMatched_DeltaR_pup1_Wd1     ,    "JetGenMatched_DeltaR_pup1_Wd1/F"          );      
    event_ttree->Branch("JetGenMatched_DeltaR_pup1_Wd2"        , & event_ttree_row->JetGenMatched_DeltaR_pup1_Wd2     ,    "JetGenMatched_DeltaR_pup1_Wd2/F"          );               
    event_ttree->Branch("JetGenMatched_partonPt"               , & event_ttree_row->JetGenMatched_partonPt            ,    "JetGenMatched_partonPt/F"                 );      
    event_ttree->Branch("JetGenMatched_partonEta"              , & event_ttree_row->JetGenMatched_partonEta           ,    "JetGenMatched_partonEta/F"                );      
    event_ttree->Branch("JetGenMatched_partonPhi"              , & event_ttree_row->JetGenMatched_partonPhi           ,    "JetGenMatched_partonPhi/F"                );      
    event_ttree->Branch("JetGenMatched_partonMass"             , & event_ttree_row->JetGenMatched_partonMass          ,    "JetGenMatched_partonMass/F"               );      
    event_ttree->Branch("JetGenMatched_partonID"               , & event_ttree_row->JetGenMatched_partonID            ,    "JetGenMatched_partonID/F"                 );      
    event_ttree->Branch("JetGenMatched_DeltaRjetParton"        , & event_ttree_row->JetGenMatched_DeltaRjetParton     ,    "JetGenMatched_DeltaRjetParton/F"          );      

    event_ttree->Branch("HadMETpx"                        , & event_ttree_row->HadMETpx                     , "HadMETpx/F"                  );
    event_ttree->Branch("HadMETpy"                        , & event_ttree_row->HadMETpy                     , "HadMETpy/F"                  );
    event_ttree->Branch("HadMETpt"                        , & event_ttree_row->HadMETpt                     , "HadMETpt/F"                  );
    event_ttree->Branch("HadMETphi"                       , & event_ttree_row->HadMETphi                    , "HadMETphi/F"                 );
    event_ttree->Branch("HadMETsumET"                     , & event_ttree_row->HadMETsumET                  , "HadMETsumET/F"               );
    event_ttree->Branch("HadMETgenMET"                    , & event_ttree_row->HadMETgenMET                 , "HadMETgenMET/F"              );
    event_ttree->Branch("HadMETuncorPt"                   , & event_ttree_row->HadMETuncorPt                , "HadMETuncorPt/F"             );
    
    event_ttree->Branch("HadMETshiftedPtJetEnUp"      , & event_ttree_row->HadMETshiftedPtJetEnUp   , "HadMETshiftedPtJetEnUp/F"     );
    event_ttree->Branch("HadMETshiftedPtJetEnDn"      , & event_ttree_row->HadMETshiftedPtJetEnDn   , "HadMETshiftedPtJetEnDn/F"     );
    event_ttree->Branch("HadMETshiftedPtElEnUp"       , & event_ttree_row->HadMETshiftedPtElEnUp    , "HadMETshiftedPtElEnUp/F"      );
    event_ttree->Branch("HadMETshiftedPtElEnDn"       , & event_ttree_row->HadMETshiftedPtElEnDn    , "HadMETshiftedPtElEnDn/F"      );
    event_ttree->Branch("HadMETshiftedPtMuEnUp"       , & event_ttree_row->HadMETshiftedPtMuEnUp    , "HadMETshiftedPtMuEnUp/F"      );
    event_ttree->Branch("HadMETshiftedPtMuEnDn"       , & event_ttree_row->HadMETshiftedPtMuEnDn    , "HadMETshiftedPtMuEnDn/F"      );
    event_ttree->Branch("HadMETshiftedPtJetResUp"     , & event_ttree_row->HadMETshiftedPtJetResUp  , "HadMETshiftedPtJetResUp/F"    );
    event_ttree->Branch("HadMETshiftedPtJetResDn"     , & event_ttree_row->HadMETshiftedPtJetResDn  , "HadMETshiftedPtJetResDn/F"    );
    event_ttree->Branch("HadMETshiftedPtUnclEnUp"     , & event_ttree_row->HadMETshiftedPtUnclEnUp  , "HadMETshiftedPtUnclEnUp/F"    );
    event_ttree->Branch("HadMETshiftedPtUnclEnDn"     , & event_ttree_row->HadMETshiftedPtUnclEnDn  , "HadMETshiftedPtUnclEnDn/F"    );
    
    event_ttree->Branch("HadNvtx"                         , & event_ttree_row->HadNvtx                      , "HadNvtx/F"                   );
    event_ttree->Branch("HadNvtxGood"                     , & event_ttree_row->HadNvtxGood                  , "HadNvtxGood/F"               );
    event_ttree->Branch("HadRho"                          , & event_ttree_row->HadRho                       , "HadRho/F"                    );
    event_ttree->Branch("HadEventWeight"                  , & event_ttree_row->HadEventWeight               , "HadEventWeight/F"            );
    event_ttree->Branch("HadPUweight"                     , & event_ttree_row->HadPUweight                  , "HadPUweight/F"            );
    event_ttree->Branch("HadPUweight_MBup"                , & event_ttree_row->HadPUweight_MBup             , "HadPUweight_MBup/F"            );
    event_ttree->Branch("HadPUweight_MBdn"                , & event_ttree_row->HadPUweight_MBdn             , "HadPUweight_MBdn/F"            );
    
    event_ttree->Branch("HadGenTTmass"                    , & event_ttree_row->HadGenTTmass                 , "HadGenTTmass/F"              );
    event_ttree->Branch("HadGenCountHadTop"               , & event_ttree_row->HadGenCountHadTop            , "HadGenCountHadTop/I"         );
      
    event_ttree->Branch("HTlep"                                , & event_ttree_row->HTlep                             , "HTlep/F"                  );
    event_ttree->Branch("ST"                                   , & event_ttree_row->ST                                , "ST/F"                     );
    event_ttree->Branch("ST_CorrDn"                            , & event_ttree_row->ST_CorrDn                         , "ST_CorrDn/F"              );
    event_ttree->Branch("ST_CorrUp"                            , & event_ttree_row->ST_CorrUp                         , "ST_CorrUp/F"              );
    event_ttree->Branch("ST_PtSmearNom"                        , & event_ttree_row->ST_PtSmearNom                     , "ST_PtSmearNom/F"          );
    event_ttree->Branch("ST_PtSmearUp"                         , & event_ttree_row->ST_PtSmearUp                      , "ST_PtSmearUp/F"           );
    event_ttree->Branch("ST_PtSmearDn"                         , & event_ttree_row->ST_PtSmearDn                      , "ST_PtSmearDn/F"           );
      
    event_ttree->Branch("HadQ2weight_CorrDn"              , & event_ttree_row->HadQ2weight_CorrDn           , "HadQ2weight_CorrDn/F"        );
    event_ttree->Branch("HadQ2weight_CorrUp"              , & event_ttree_row->HadQ2weight_CorrUp           , "HadQ2weight_CorrUp/F"        );
    event_ttree->Branch("HadNNPDF3weight_CorrDn"          , & event_ttree_row->HadNNPDF3weight_CorrDn       , "HadNNPDF3weight_CorrDn/F"    );
    event_ttree->Branch("HadNNPDF3weight_CorrUp"          , & event_ttree_row->HadNNPDF3weight_CorrUp       , "HadNNPDF3weight_CorrUp/F"    );
    event_ttree->Branch("HadRunNum"                       , & event_ttree_row->HadRunNum                    , "HadRunNum/I"                 );
    event_ttree->Branch("HadLumiBlock"                    , & event_ttree_row->HadLumiBlock                 , "HadLumiBlock/I"              );
    event_ttree->Branch("HadEventNum"                     , & event_ttree_row->HadEventNum                  , "HadEventNum/I"               );
    event_ttree->Branch("HadPassMETFilters"               , & event_ttree_row->HadPassMETFilters            , "HadPassMETFilters/I"         );
    
    event_ttree->Branch("AK4_uncorr_pt",          & event_ttree_row->AK4_uncorr_pt               , "AK4_uncorr_pt/F");      
    event_ttree->Branch("AK4_corr_pt",          & event_ttree_row->AK4_corr_pt               , "AK4_corr_pt/F");  
    
    event_ttree->Branch("CA12JetPtRaw",          & event_ttree_row->CA12JetPtRaw               , "CA12JetPtRaw/F");                
    event_ttree->Branch("CA12JetEtaRaw",          & event_ttree_row->CA12JetEtaRaw               , "CA12JetEtaRaw/F");                 
    event_ttree->Branch("CA12JetPhiRaw",          & event_ttree_row->CA12JetPhiRaw               , "CA12JetPhiRaw/F");                 
    event_ttree->Branch("CA12JetMassRaw",          & event_ttree_row->CA12JetMassRaw               , "CA12JetMassRaw/F");                
    
    event_ttree->Branch("CA12JetTau1",          & event_ttree_row->CA12JetTau1               , "CA12JetTau1/F");                   
    event_ttree->Branch("CA12JetTau2",          & event_ttree_row->CA12JetTau2               , "CA12JetTau2/F");                   
    event_ttree->Branch("CA12JetTau3",          & event_ttree_row->CA12JetTau3               , "CA12JetTau3/F");                   
    event_ttree->Branch("CA12JetTau4",          & event_ttree_row->CA12JetTau4               , "CA12JetTau4/F");                   
    event_ttree->Branch("CA12JetTau32",          & event_ttree_row->CA12JetTau32               , "CA12JetTau32/F");                  
    event_ttree->Branch("CA12JetTau21",          & event_ttree_row->CA12JetTau21               , "CA12JetTau21/F");                  
    
    event_ttree->Branch("CA12Jetsubjet0bdisc",          & event_ttree_row->CA12Jetsubjet0bdisc               , "CA12Jetsubjet0bdisc/F");           
    event_ttree->Branch("CA12Jetsubjet1bdisc",          & event_ttree_row->CA12Jetsubjet1bdisc               , "CA12Jetsubjet1bdisc/F");           
    event_ttree->Branch("CA12Jetmaxbdisc",          & event_ttree_row->CA12Jetmaxbdisc               , "CA12Jetmaxbdisc/F");               
    
    event_ttree->Branch("CA12Jetsubjet0pt",          & event_ttree_row->CA12Jetsubjet0pt               , "CA12Jetsubjet0pt/F");              
    event_ttree->Branch("CA12Jetsubjet0mass",          & event_ttree_row->CA12Jetsubjet0mass               , "CA12Jetsubjet0mass/F");            
    event_ttree->Branch("CA12Jetsubjet0eta",          & event_ttree_row->CA12Jetsubjet0eta               , "CA12Jetsubjet0eta/F");             
    event_ttree->Branch("CA12Jetsubjet0phi",          & event_ttree_row->CA12Jetsubjet0phi               , "CA12Jetsubjet0phi/F");             
    event_ttree->Branch("CA12Jetsubjet0area",          & event_ttree_row->CA12Jetsubjet0area               , "CA12Jetsubjet0area/F");            
    
    event_ttree->Branch("CA12Jetsubjet1pt",          & event_ttree_row->CA12Jetsubjet1pt               , "CA12Jetsubjet1pt/F");              
    event_ttree->Branch("CA12Jetsubjet1mass",          & event_ttree_row->CA12Jetsubjet1mass               , "CA12Jetsubjet1mass/F");            
    event_ttree->Branch("CA12Jetsubjet1eta",          & event_ttree_row->CA12Jetsubjet1eta               , "CA12Jetsubjet1eta/F");             
    event_ttree->Branch("CA12Jetsubjet1phi",          & event_ttree_row->CA12Jetsubjet1phi               , "CA12Jetsubjet1phi/F");             
    event_ttree->Branch("CA12Jetsubjet1area",          & event_ttree_row->CA12Jetsubjet1area               , "CA12Jetsubjet1area/F");            
    
    event_ttree->Branch("AK4ReconstructedJetPt",          & event_ttree_row->AK4ReconstructedJetPt               , "AK4ReconstructedJetPt/F");         
    event_ttree->Branch("AK4ReconstructedJetEta",          & event_ttree_row->AK4ReconstructedJetEta               , "AK4ReconstructedJetEta/F");        
    event_ttree->Branch("AK4ReconstructedJetPhi",          & event_ttree_row->AK4ReconstructedJetPhi               , "AK4ReconstructedJetPhi/F");        
    event_ttree->Branch("AK4ReconstructedJetMass",          & event_ttree_row->AK4ReconstructedJetMass               , "AK4ReconstructedJetMass/F");       
    
    event_ttree->Branch("AK4bJetPtRaw",      & event_ttree_row->AK4bJetPtRaw ,             "AK4bJetPtRaw/F");       
    event_ttree->Branch("AK4bJetEtaRaw",      & event_ttree_row->AK4bJetEtaRaw ,             "AK4bJetEtaRaw/F");      
    event_ttree->Branch("AK4bJetPhiRaw",      & event_ttree_row->AK4bJetPhiRaw ,             "AK4bJetPhiRaw/F");      
    event_ttree->Branch("AK4bJetMassRaw",      & event_ttree_row->AK4bJetMassRaw ,             "AK4bJetMassRaw/F");     
    
    event_ttree->Branch("AK4bJet_PtSmear",      & event_ttree_row->AK4bJet_PtSmear ,             "AK4bJet_PtSmear/F");      
    event_ttree->Branch("AK4bJet_PtSmearUp",      & event_ttree_row->AK4bJet_PtSmearUp ,             "AK4bJet_PtSmearUp/F");    
    event_ttree->Branch("AK4bJet_PtSmearDn",      & event_ttree_row->AK4bJet_PtSmearDn ,             "AK4bJet_PtSmearDn/F");    
    event_ttree->Branch("AK4bJet_PtUncorr",      & event_ttree_row->AK4bJet_PtUncorr ,             "AK4bJet_PtUncorr/F");   
    event_ttree->Branch("AK4bJet_Corr",      & event_ttree_row->AK4bJet_Corr ,             "AK4bJet_Corr/F");       
    event_ttree->Branch("AK4bJet_CorrUp",      & event_ttree_row->AK4bJet_CorrUp ,             "AK4bJet_CorrUp/F");     
    event_ttree->Branch("AK4bJet_CorrDn",      & event_ttree_row->AK4bJet_CorrDn ,             "AK4bJet_CorrDn/F");     
    event_ttree->Branch("AK4bJet_bDisc",      & event_ttree_row->AK4bJet_bDisc ,             "AK4bJet_bDisc/F");   
    
    event_ttree->Branch("AK4WJetPtRaw",      & event_ttree_row->AK4WJetPtRaw ,             "AK4WJetPtRaw/F");       
    event_ttree->Branch("AK4WJetEtaRaw",      & event_ttree_row->AK4WJetEtaRaw ,             "AK4WJetEtaRaw/F");      
    event_ttree->Branch("AK4WJetPhiRaw",      & event_ttree_row->AK4WJetPhiRaw ,             "AK4WJetPhiRaw/F");      
    event_ttree->Branch("AK4WJetMassRaw",      & event_ttree_row->AK4WJetMassRaw ,             "AK4WJetMassRaw/F");     
    
    event_ttree->Branch("AK4WJet_PtSmear",      & event_ttree_row->AK4WJet_PtSmear ,             "AK4WJet_PtSmear/F");      
    event_ttree->Branch("AK4WJet_PtSmearUp",      & event_ttree_row->AK4WJet_PtSmearUp ,             "AK4WJet_PtSmearUp/F");    
    event_ttree->Branch("AK4WJet_PtSmearDn",      & event_ttree_row->AK4WJet_PtSmearDn ,             "AK4WJet_PtSmearDn/F");    
    event_ttree->Branch("AK4WJet_PtUncorr",      & event_ttree_row->AK4WJet_PtUncorr ,             "AK4WJet_PtUncorr/F");   
    event_ttree->Branch("AK4WJet_Corr",      & event_ttree_row->AK4WJet_Corr ,             "AK4WJet_Corr/F");       
    event_ttree->Branch("AK4WJet_CorrUp",      & event_ttree_row->AK4WJet_CorrUp ,             "AK4WJet_CorrUp/F");     
    event_ttree->Branch("AK4WJet_CorrDn",      & event_ttree_row->AK4WJet_CorrDn ,             "AK4WJet_CorrDn/F");     
    event_ttree->Branch("AK4WJet_bDisc",      & event_ttree_row->AK4WJet_bDisc ,             "AK4WJet_bDisc/F");      
    
    event_ttree->Branch("AK4W2JetPtRaw",      & event_ttree_row->AK4W2JetPtRaw ,             "AK4W2JetPtRaw/F");       
    event_ttree->Branch("AK4W2JetEtaRaw",      & event_ttree_row->AK4W2JetEtaRaw ,             "AK4W2JetEtaRaw/F");      
    event_ttree->Branch("AK4W2JetPhiRaw",      & event_ttree_row->AK4W2JetPhiRaw ,             "AK4W2JetPhiRaw/F");      
    event_ttree->Branch("AK4W2JetMassRaw",      & event_ttree_row->AK4W2JetMassRaw ,             "AK4W2JetMassRaw/F");     
    
    event_ttree->Branch("AK4W2Jet_PtSmear",      & event_ttree_row->AK4W2Jet_PtSmear ,             "AK4W2Jet_PtSmear/F");      
    event_ttree->Branch("AK4W2Jet_PtSmearUp",      & event_ttree_row->AK4W2Jet_PtSmearUp ,             "AK4W2Jet_PtSmearUp/F");    
    event_ttree->Branch("AK4W2Jet_PtSmearDn",      & event_ttree_row->AK4W2Jet_PtSmearDn ,             "AK4W2Jet_PtSmearDn/F");    
    event_ttree->Branch("AK4W2Jet_PtUncorr",      & event_ttree_row->AK4W2Jet_PtUncorr ,             "AK4W2Jet_PtUncorr/F");   
    event_ttree->Branch("AK4W2Jet_Corr",      & event_ttree_row->AK4W2Jet_Corr ,             "AK4W2Jet_Corr/F");       
    event_ttree->Branch("AK4W2Jet_CorrUp",      & event_ttree_row->AK4W2Jet_CorrUp ,             "AK4W2Jet_CorrUp/F");     
    event_ttree->Branch("AK4W2Jet_CorrDn",      & event_ttree_row->AK4W2Jet_CorrDn ,             "AK4W2Jet_CorrDn/F");     
    event_ttree->Branch("AK4W2Jet_bDisc",      & event_ttree_row->AK4W2Jet_bDisc ,             "AK4W2Jet_bDisc/F");      
    
    event_ttree->Branch("MuPhi", "vector<float>",&event_ttree_row->MuPhi);
    event_ttree->Branch("MuPt", "vector<float>",&event_ttree_row->MuPt);
    event_ttree->Branch("MuEta", "vector<float>",&event_ttree_row->MuEta);
    event_ttree->Branch("MuMass", "vector<float>",&event_ttree_row->MuMass);
    event_ttree->Branch("MuIso", "vector<float>",&event_ttree_row->MuIso);
    event_ttree->Branch("MuIsoTrk", "vector<float>",&event_ttree_row->MuIsoTrk);
    event_ttree->Branch("MuMedium", "vector<int>", &event_ttree_row->MuMedium);
    event_ttree->Branch("MuTight", "vector<int>", &event_ttree_row->MuTight);
    
    event_ttree->Branch("Electron_Phi", "vector<float>",&event_ttree_row->Electron_Phi);
    event_ttree->Branch("Electron_Pt", "vector<float>",&event_ttree_row->Electron_Pt);
    event_ttree->Branch("Electron_Eta", "vector<float>",&event_ttree_row->Electron_Eta);
    event_ttree->Branch("Electron_Mass", "vector<float>",&event_ttree_row->Electron_Mass);
    event_ttree->Branch("Elecron_absiso", "vector<float>",&event_ttree_row->Elecron_absiso);
    event_ttree->Branch("Elecron_relIsoWithDBeta", "vector<float>",&event_ttree_row->Elecron_relIsoWithDBeta);
    event_ttree->Branch("Elecron_absiso_EA", "vector<float>",&event_ttree_row->Elecron_absiso_EA);
    event_ttree->Branch("Elecron_relIsoWithEA", "vector<float>",&event_ttree_row->Elecron_relIsoWithEA);
    
    event_ttree->Branch("Electron_iso_passHLTpre", "vector<int>", &event_ttree_row->Electron_iso_passHLTpre);
    event_ttree->Branch("Electron_iso_passLoose", "vector<int>", &event_ttree_row->Electron_iso_passLoose);
    event_ttree->Branch("Electron_iso_passMedium", "vector<int>", &event_ttree_row->Electron_iso_passMedium);
    event_ttree->Branch("Electron_iso_passTight", "vector<int>", &event_ttree_row->Electron_iso_passTight);
    event_ttree->Branch("Electron_iso_passHEEP", "vector<int>", &event_ttree_row->Electron_iso_passHEEP);
    event_ttree->Branch("Electron_noiso_passLoose", "vector<int>", &event_ttree_row->Electron_noiso_passLoose);
    event_ttree->Branch("Electron_noiso_passMedium", "vector<int>", &event_ttree_row->Electron_noiso_passMedium);
    event_ttree->Branch("Electron_noiso_passTight", "vector<int>", &event_ttree_row->Electron_noiso_passTight);
    event_ttree->Branch("Electron_noiso_passHEEP", "vector<int>", &event_ttree_row->Electron_noiso_passHEEP);

}