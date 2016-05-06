TH1F *hPhiAB_blue = new TH1F("hPhiAB_blue", "Phi_AB FOPI Centrality blue", 180, 0., 180.);
TH1F *hPhiAB_gren = new TH1F("hPhiAB_gren", "Phi_AB FOPI Centrality gren", 180, 0., 180.);
TH1F *hPhiAB_redd = new TH1F("hPhiAB_redd", "Phi_AB FOPI Centrality redd", 180, 0., 180.);
//-------------------------------------------------------------------------------------------------
TH1F* hWEPy[11];
TH1F* hOWEPy[11];
TH2F* hSinPsi[6][11];
TH2F* hCosPsi[6][11];
TH2F* hSinFOPI[6][11];
TH2F* hCosFOPI[6][11];
TH1F* hPsiCorr[ 11];
TH1F* hOPsiCorr[11];
for(Int_t i=0;i<11;i++){
    hWEPy[i]  = new TH1F(Form("hWEPy%i",i)   ,Form("hWEP %i META",i)              , 90, -180., 180.); //hWEPy[ i]->Sumw2();
    hOWEPy[i] = new TH1F(Form("hOWEPy%i",i)  ,Form("hOWEP %i MDC",i)              , 90, -180., 180.); //hOWEPy[i]->Sumw2();
    for (Int_t n=0;n<6;n++){
    	hSinPsi[n][i] = new TH2F(Form("hSinPsi%i%i",n,i),Form("hSinPsi %i %i",n,i),15,95.,110.,99, -1.1,1.1);
    	hCosPsi[n][i] = new TH2F(Form("hCosPsi%i%i",n,i),Form("hCosPsi %i %i",n,i),15,95.,110.,99, -1.1,1.1);
    }
    hPsiCorr[ i] = new TH1F(Form("hPsiCorr%i" ,i),Form("hPsiCorr %i after_corr" ,i), 90, -180., 180.);
    hOPsiCorr[i] = new TH1F(Form("hOPsiCorr%i",i),Form("hOPsiCorr %i after_corr",i), 90, -180., 180.);
}
TH1F* hCorFOPI[3];
TH1F* hFOPIWEP[3];
for(Int_t n=0;n<6;n++){
	for(Int_t i=0;i<3;i++){
		hCosFOPI[n][i] = new TH2F(Form("hCosFOPI%i%i",n,i),Form("hCosFOPI %i %i",n,i),15,95.,110.,99, -1.1,1.1);
		hSinFOPI[n][i] = new TH2F(Form("hSinFOPI%i%i",n,i),Form("hSinFOPI %i %i",n,i),15,95.,110.,99, -1.1,1.1);
	}
}
for(Int_t i=0;i<3;i++){
	hCorFOPI[i] = new TH1F(Form("hCorFOPI%i" ,i),Form("hCorFOPI %i after_corr" ,i), 90, -180., 180.);
	hFOPIWEP[i] = new TH1F(Form("hFOPIWEP%i" ,i),Form("hFOPIWEP %i after_corr" ,i), 90, -180., 180.);
}

TH1F* hPCpt = new TH1F("hPCpt","Corrected pt for protons"  ,800,0.,3000.);
TH1F* hpt0  = new TH1F("hpt0" ,"Uncorrected pt for protons",800,0.,3000.);
//-------------------------------------------------------------------------------------------------
//For-Event-Plane-reconstruction-------------------------------------------------------------------------------------------
  TH1F* hFWcellTime    = new TH1F("hFWcellTime",  "FW cell time"            ,  2000,  0.,  200.);
  TH2F* hFW_TimeCell_0 = new TH2F("hFW_TimeCell_0", "FW cell vs. Time", 500, 0., 1000., 300, 0., 300.);
  TH2F* hFW_TimeCell_1 = new TH2F("hFW_TimeCell_1", "FW cell vs. Time", 500, 0., 1000., 300, 0., 300.);
  TH2F* hFW_TimeCell_2 = new TH2F("hFW_TimeCell_2", "FW cell vs. Time", 500, 0., 1000., 300, 0., 300.);
  TH2F* hFW_TimeCell_3 = new TH2F("hFW_TimeCell_3", "FW cell vs. Time", 500, 0., 1000., 300, 0., 300.);

  TH2F* hFW_dEdxCell_0 = new TH2F("hFW_dEdxCell_0", "FW cell vs. dEdx", 700, 0., 1400., 300, 0., 300.);
  TH2F* hFW_dEdxCell_1 = new TH2F("hFW_dEdxCell_1", "FW cell vs. dEdx", 700, 0., 1400., 300, 0., 300.);
  TH2F* hFW_dEdxCell_2 = new TH2F("hFW_dEdxCell_2", "FW cell vs. dEdx", 700, 0., 1400., 300, 0., 300.);
  TH2F* hFW_dEdxCell_3 = new TH2F("hFW_dEdxCell_3", "FW cell vs. dEdx", 700, 0., 1400., 300, 0., 300.);

  TH2F *hFWxyCC        = new TH2F("hFWxyCC"       ,"FW Y vs. X cells centers"              ,500,-1000.,1000.,500,-1000.,1000.);
  TH2F *hFWxyCCsmear   = new TH2F("hFWxyCCsmear"  ,"FW Y vs. X cells smeared"              ,500,-1000.,1000.,500,-1000.,1000.);
  TH2F *hFWxySHsmear   = new TH2F("hFWxySHsmear"  ,"FW Y vs. X cells smeared and shifted"  ,500,-1000.,1000.,500,-1000.,1000.);

  TH2F *hFWxyCCsmearM[8];
  for(Int_t i=0;i<8;i++){
    hFWxyCCsmearM[i]   = new TH2F(Form("hFWxyCCsmearM%i",i)  ,Form("FW Y vs. X cells smeared Mul%i",i),500,-1000.,1000.,500,-1000.,1000.);
  }

  TH2F *hFWxyCCsmearXscanM[8];
  TH2F *hFWxyCCsmearYscanM[8];
  for(Int_t i=0;i<8;i++){
   hFWxyCCsmearXscanM[i] = new TH2F(Form("hFWxyCCsmearXscanM%i",i)  ,Form("FW Y vs. X cells smeared X-scan Mul%i",i) ,500,-1000.,1000.,500,-1000.,1000.);
   hFWxyCCsmearYscanM[i] = new TH2F(Form("hFWxyCCsmearYscanM%i",i)  ,Form("FW Y vs. X cells smeared Y-scan Mul%i",i) ,500,-1000.,1000.,500,-1000.,1000.);
  }

  TH2F *hFWxyCCLsmear  = new TH2F("hFWxyCCLsmear"  ,"FW Y vs. X cells (Large)   smeared"    ,500,-1000.,1000.,500,-1000.,1000.);
  TH2F *hFWxyCCMsmear  = new TH2F("hFWxyCCMsmear"  ,"FW Y vs. X cells (Middle)  smeared"    ,500,-1000.,1000.,500,-1000.,1000.);
  TH2F *hFWxyCCNsmear  = new TH2F("hFWxyCCNsmear"  ,"FW Y vs. X cells (Nearest) smeared"    ,500,-1000.,1000.,500,-1000.,1000.);

  TH1D* h0PhiEP        = new TH1D("h0PhiEP","Nfw>=0: EP phi"                      , 400, -200., 200.);
  TH1D* h0PhiA         = new TH1D("h0PhiA" ,"Nfw>=0: EP phi A 1/2-subevent"       , 400, -200., 200.);
  TH1D* h0PhiB         = new TH1D("h0PhiB" ,"Nfw>=0: EP phi B 1/2-subevent"       , 400, -200., 200.);
  TH1D* h0PhiAB        = new TH1D("h0PhiAB","Nfw>=0: EP phi(A^B) of subevents"    , 400, -200., 200.);

  TH2F *hMfwMmeta      = new TH2F("hMfwMmeta","Mfw vs. Mtof+Mrpc", 100, 0., 400., 100, 0., 100.);
  TH1F *hMfw           = new TH1F("hMfw"     ,"Mfw"              ,                100, 0., 100.);
  TH1F *hMmeta         = new TH1F("hMmeta"   ,"Mmeta"            , 100, 0., 400.               );
  TH1F *hMtof          = new TH1F("hMtof"    ,"Mtof"             , 100, 0., 200.               );
  TH1F *hMrpc          = new TH1F("hMrpc"    ,"Mrpc"             , 100, 0., 200.               );

  TH1F *hRPA           = new TH1F("hRPA" ,"RP angle"                           , 90., -380., 380);

  TF1  *fRPAy[ 11];
  TH1F *hRPAy[ 11];
  TH1F *hRPAyc[11];
  for(Int_t i=0;i<11;i++){
    //---angular-weight-correction-function--[cos1A+cos2A+sin1A +cos3A + sin2A ]-----------------
    fRPAy[i]  = new TF1( Form("fRPAy%i",i) ,  "( 1. +  2.*[0]*TMath::Cos(x[0]*3.1415/180.) + 2.*[1]*TMath::Cos(2.*x[0]*3.1415/180.) + [2]*TMath::Sin(x[0]*3.1415/180.) + 2.*[3]*TMath::Cos(3.*x[0]*3.1415/180.) +  [4]*TMath::Sin(2.*x[0]*3.1415/180.))", -180., 180.);
    //--histograms-for-angular-weighting-correction--(before-correction)-------------------------
    hRPAy[i]  = new TH1F(Form("hRPAy%i",i) ,Form("RP angle within b%i  range"           ,i), 90., -180., 180);
    //--histograms-for-angular-weighting-correction-(after-correction)---------------------------
    hRPAyc[i] = new TH1F(Form("hRPAy%ic",i),Form("RP angle within b%i  range after corr",i), 90., -180., 180);
  }
 
  TH1F *hRPA6Ywc  = new TH1F("hRPA6Ywc","RP angle after corr in 6 b range", 90., -180., 180); //-after-weight-correction-in-6-y-slices--

  //My experiment
  TH1F* hPhiEP[ 6];
  TH1F* hPhiEPrec[6];
  TH1F* hPhiEPc[6];
  TH1F* hQvX[6];
  TH1F* hQvY[6];
  TH1F* hQvXrec[6];
  TH1F* hQvYrec[6];
  for(Int_t n=0;n<6;n++){
    hPhiEP[ n] = new TH1F(Form("hPhiEP%i" ,n) ,Form("EP angle %i harm.",n+1), 400, -200./(n+1), 200./(n+1));                  hPhiEP[ n]  ->Sumw2();
    hPhiEPrec[n]=new TH1F(Form("hPhiEPrec%i",n) ,Form("EP angle after recentering, %i harm.",n+1),400,-200./(n+1),200./(n+1));hPhiEPrec[n]->Sumw2();
    hPhiEPc[n] = new TH1F(Form("hPhiEPc%i",n) ,Form("EP angle after eff. corr. %i harm.",n+1), 400, -200./(n+1), 200./(n+1)); hPhiEPc[n]  ->Sumw2();
    hQvX   [n] = new TH1F(Form("hQvX%i"   ,n) ,Form("Q vector X %i harm.",n+1), 400, -1.1, 1.1);                  hQvX   [n]  ->Sumw2();
    hQvY   [n] = new TH1F(Form("hQvY%i"   ,n) ,Form("Q vector Y %i harm.",n+1), 400, -1.1, 1.1);                  hQvY   [n]  ->Sumw2();
    hQvXrec[n] = new TH1F(Form("hQvXrec%i",n) ,Form("Q vector X after recentering, %i harm.",n+1),400,-1.1,1.1);  hQvXrec[n]  ->Sumw2();
    hQvYrec[n] = new TH1F(Form("hQvYrec%i",n) ,Form("Q vector Y after recentering, %i harm.",n+1),400,-1.1,1.1);  hQvYrec[n]  ->Sumw2();
  }
  TH1F* hPhiEPvect = new TH1F(Form("hPhiEPvect") ,Form("EP angle from vsum"), 760, -380., 380);           hPhiEPvect ->Sumw2();
  TH1F* h0PhiEPvect = new TH1F(Form("h0PhiEPvect") ,Form("EP angle from vsum (nfw>1)"), 760, -380., 380); h0PhiEPvect->Sumw2();
  