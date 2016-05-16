TH1F *hPhiAB[11];
for(Int_t im=0;im<11;im++){
    hPhiAB[im] = new TH1F(Form("hPhiAB%i",im),Form("hPhiAB %i centrality_bins",im), 180,0.,180.);
}
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
TH2F *hFWxyCCsmearM[11];
for(Int_t i=0;i<11;i++){
    hFWxyCCsmearM[i]   = new TH2F(Form("hFWxyCCsmearM%i",i)  ,Form("hFWxyCCsmearM %i",i),500,-1000.,1000.,500,-1000.,1000.);
}
TH2F *hFWxyCCsmearXscanM[11];
TH2F *hFWxyCCsmearYscanM[11];
for(Int_t i=0;i<11;i++){
    hFWxyCCsmearXscanM[i] = new TH2F(Form("hFWxyCCsmearXscanM%i",i)  ,Form("FW Y vs. X cells smeared X-scan Mul%i",i) ,500,-1000.,1000.,500,-1000.,1000.);
    hFWxyCCsmearYscanM[i] = new TH2F(Form("hFWxyCCsmearYscanM%i",i)  ,Form("FW Y vs. X cells smeared Y-scan Mul%i",i) ,500,-1000.,1000.,500,-1000.,1000.);
}
TH1F* hPhiEPvect = new TH1F(Form("hPhiEPvect") ,Form("EP angle from vsum"), 760, -380., 380);           hPhiEPvect ->Sumw2();
TH1F* h0PhiEPvect = new TH1F(Form("h0PhiEPvect") ,Form("EP angle from vsum (nfw>1)"), 760, -380., 380); h0PhiEPvect->Sumw2();
TH1F* hPhiEP[ 6];
TH1F* hPhiEPrec[6];
TH1F* hPhiEPc[6];
TH1F* hQvX[6];
TH1F* hQvY[6];
for(Int_t n=0;n<6;n++){
    hPhiEP[ n] = new TH1F(Form("hPhiEP%i" ,n) ,Form("EP angle %i harm.",n+1), 400, -200./(n+1), 200./(n+1));                  hPhiEP[ n]  ->Sumw2();
    hPhiEPrec[n]=new TH1F(Form("hPhiEPrec%i",n) ,Form("EP angle after recentering, %i harm.",n+1),400,-200./(n+1),200./(n+1));hPhiEPrec[n]->Sumw2();
    hPhiEPc[n] = new TH1F(Form("hPhiEPc%i",n) ,Form("EP angle after eff. corr. %i harm.",n+1), 400, -200./(n+1), 200./(n+1)); hPhiEPc[n]  ->Sumw2();
    hQvX   [n] = new TH1F(Form("hQvX%i"   ,n) ,Form("Q vector X %i harm.",n+1), 400, -1.1, 1.1);                  hQvX   [n]  ->Sumw2();
    hQvY   [n] = new TH1F(Form("hQvY%i"   ,n) ,Form("Q vector Y %i harm.",n+1), 400, -1.1, 1.1);                  hQvY   [n]  ->Sumw2();
}
TH1D* h0PhiEP        = new TH1D("h0PhiEP","Nfw>=2: EP phi"                      , 400, -200., 200.);
TH1D* hPhiCor        = new TH1D("hPhiCor","Nfw>=2: EP phi corr"                 , 400, -200., 200.);
TH1D* h0PhiA         = new TH1D("h0PhiA" ,"Nfw>=2: EP phi A 1/2-subevent"       , 400, -200., 200.);
TH1D* h0PhiB         = new TH1D("h0PhiB" ,"Nfw>=2: EP phi B 1/2-subevent"       , 400, -200., 200.);
TH1D* h0PhiAB        = new TH1D("h0PhiAB","Nfw>=2: EP phi(A^B) of subevents"    , 400, -200., 200.);
TH1F *hRPA           = new TH1F("hRPA" ,"RP angle"                           , 90., -380., 380);

TProfile* hSinPsi[3][6][11];
TProfile* hCosPsi[3][6][11];
TH1F* hPsiEP[  3][ 11];
TH1F* hPsiRcnt[3][ 11];
TH1F* hPsiCorr[3][ 11];
TH1F* hOPsiCorr[3][11];
TProfile* hsumXmean[3][11];
TProfile* hsumYmean[3][11];
TH1F* hQvectX[3][  11];
TH1F* hQvectY[3][  11];
TH1F* hQvXrec[3][  11];
TH1F* hQvYrec[3][  11];
TH2F* hQvRaw[ 3][  11];
TH2F* hQvRec[ 3][  11];
TH1F* hdPsi[  3][  11];
for (Int_t iT=0;iT<3;iT++){
    for(Int_t i=0;i<11;i++){
        hsumXmean[iT][i] = new TProfile(Form("hsumXmean%i%i",iT,i),Form("hsumXmean %i %i",iT,i),15,95.,110.);
        hsumYmean[iT][i] = new TProfile(Form("hsumYmean%i%i",iT,i),Form("hsumYmean %i %i",iT,i),15,95.,110.);
        for (Int_t n=0;n<6;n++){
    	    hSinPsi[iT][n][i] = new TProfile(Form("hSinPsi%i%i%i",iT,n,i),Form("hSinPsi %i %i %i",iT,n,i),15,95.,110.);
            hCosPsi[iT][n][i] = new TProfile(Form("hCosPsi%i%i%i",iT,n,i),Form("hCosPsi %i %i %i",iT,n,i),15,95.,110.);
        }
        hPsiEP[   iT][i] = new TH1F(Form("hPsiEP%i%i"   ,iT,i),Form("hPsiEP %i %i "   ,iT,i), 90, -180., 180.);
        hPsiRcnt[ iT][i] = new TH1F(Form("hPsiRcnt%i%i" ,iT,i),Form("hPsiRcnt %i %i after_recentering" ,iT,i), 90, -180., 180.);
        hPsiCorr[ iT][i] = new TH1F(Form("hPsiCorr%i%i" ,iT,i),Form("hPsiCorr %i %i after_flattering" ,iT,i), 90, -180., 180.);
        hdPsi[    iT][i] = new TH1F(Form("hdPsi%i%i" ,iT,i),Form("hdPsi %i %i after_flattering" ,iT,i), 100, -5., 5.);
        hOPsiCorr[iT][i] = new TH1F(Form("hOPsiCorr%i%i",iT,i),Form("hOPsiCorr %i %i after_corr",iT,i), 90, -180., 180.);
        hQvectX[  iT][i] = new TH1F(Form("hQvectX%i%i",iT,i),Form("hQvectX %i %i FW",iT,i), 600, -2., 2.);
        hQvectY[  iT][i] = new TH1F(Form("hQvectY%i%i",iT,i),Form("hQvectY %i %i FW",iT,i), 600, -2., 2.); 
        hQvXrec[  iT][i] = new TH1F(Form("hQvXrec%i%i",iT,i),Form("hQvXrec %i %i FW",iT,i), 600, -5., 5.);
        hQvYrec[  iT][i] = new TH1F(Form("hQvYrec%i%i",iT,i),Form("hQvYrec %i %i FW",iT,i), 600, -5., 5.);
        hQvRaw[   iT][i] = new TH2F(Form("hQvRaw%i%i",iT,i) ,Form("Qvect %i %i Raw",iT,i), 1000, -2., 2., 1000, -2., 2.);
        hQvRec[   iT][i] = new TH2F(Form("hQvRec%i%i",iT,i) ,Form("Qvect %i %i Rec",iT,i), 1000, -5., 5., 1000, -5., 5.);
    }
}
 
TH1F *hRPA6Ywc  = new TH1F("hRPA6Ywc","RP angle after corr in 6 b range", 90., -180., 180); //-after-weight-correction-in-6-y-slices--
TH2F *hMfwMmeta      = new TH2F("hMfwMmeta","Mfw vs. Mtof+Mrpc", 100, 0., 400., 100, 0., 100.);
TH1F *hMfw           = new TH1F("hMfw"     ,"Mfw"              ,                100, 0., 100.);
TH1F *hMmeta         = new TH1F("hMmeta"   ,"Mmeta"            , 100, 0., 400.               );
TH1F *hMtof          = new TH1F("hMtof"    ,"Mtof"             , 100, 0., 200.               );
TH1F *hMrpc          = new TH1F("hMrpc"    ,"Mrpc"             , 100, 0., 200.               );
//Debug & fixing-----------------------------------------------------------------------------------
TH1F* hPCpt          = new TH1F("hPCpt","Corrected pt for protons"  ,800,0.,3000.);
TH1F* hpt0           = new TH1F("hpt0" ,"Uncorrected pt for protons",800,0.,3000.);
//-------------------------------------------------------------------------------------------------
TH2F* hbetaMomAll    = new TH2F(Form("hbetaMomAll"),"",6000,-3000,3000,600,0,1.5);
TH2F* hBetaMeanP     = new TH2F(Form("hBetaMeanP" ),"",6000,-3000,3000,600,0,1.5);
TH2F* hBetaDownP     = new TH2F(Form("hBetaDownP" ),"",6000,-3000,3000,600,0,1.5);
TH2F* hBetaUpPro     = new TH2F(Form("hBetaUpPro" ),"",6000,-3000,3000,600,0,1.5);
TH2F* hbetaMomPro    = new TH2F(Form("hbetaMomPro"),"",6000,-3000,3000,600,0,1.5);
TH2F* hdEdxMAll      = new TH2F("hdEdxMAll","dEdx (MDC) vs mass",600,-0.1,2.4,600,0.,16.);
TH2F* hdEdxMom       = new TH2F("hdEdxMom" ,"dEdx (MDC) vs mom" ,600,200.,2000.,600,0.,16.);
TH2F* hfdEdxMomLow   = new TH2F("hfdEdxMomLow","dEdx (MDC) Low limit vs mom" ,600,200.,2000.,600,0.,16.);
TF1*  fBetaMeanP     = new TF1("fBetaMeanP",fBeta,600,0,1.5);

TH1F* hcuts          = new TH1F("hcuts"    ,"N_{events} vs cuts",6,0.5,6.5);
TH1F* hPartCuts      = new TH1F("hPartCuts","N_{tracks} vs cuts",4,0.5,4.5);
TH1F* hvertexZ       = new TH1F(Form("hvertexZ") ,"", 1000,-100,100);
TH1F* hvertexX       = new TH1F(Form("hvertexX") ,"", 1000,-100,100);
TH1F* hvertexY       = new TH1F(Form("hvertexY") ,"", 1000,-100,100);
TH2F* hvertexXZ      = new TH2F(Form("hvertexXZ"),"", 300,-80,20,500,-10,10);
TH2F* hvertexXY      = new TH2F(Form("hvertexXY"),"", 500,-10,10,500,-10,10);
TH1F* hvtxCutZ       = new TH1F(Form("hvtxCutZ") ,"", 1000,-100,100);
TH1F* hvtxCutX       = new TH1F(Form("hvtxCutX") ,"", 1000,-100,100);
TH1F* hvtxCutY       = new TH1F(Form("hvtxCutY") ,"", 1000,-100,100);
TH2F* hvtxCutXZ      = new TH2F(Form("hvtxCutXZ"),"", 300,-80,20,500,-10,10);
TH2F* hvtxCutXY      = new TH2F(Form("hvtxCutXY"),"", 500,-10,10,500,-10,10);
TH1F* hTrkMult       = new TH1F(Form("hTrkMult") ,"", 150, 0.,150. );
TH1F* hTrkMultCorr   = new TH1F(Form("hTrkMultCorr") ,"", 150, 0.,150. );
TH2D* hMETAvsTRK     = new TH2D("hMETAvsTRK"    ,"Multiplicity (Tracking):(TOF+RPC)"         , 300, 0., 300, 160, 0., 160.);
TH2D* hMETAvsTRKcorr = new TH2D("hMETAvsTRKcorr","Multiplicity (TrackingCorrected):(TOF+RPC)", 300, 0., 300, 160, 0., 160.);
TH1D* hChi2          = new TH1D("hChi2"   ,"Particle Candidate Chi2 RK"           , 500, -10., 500.);
TH1D* hChi2In        = new TH1D("hChi2In" ,"Particle Candidate Chi2 MDC12"        , 250, -10., 25.);
TH1D* hChi2Out       = new TH1D("hChi2Out","Particle Candidate Chi2 MDC34"        , 250, -10., 25.);
TH1F* hMom           = new TH1F("hMom"       ,"Momentum of proton", 800, 0.,5000.);
TH1F* hMomCorr       = new TH1F("hMomCorr"   ,"Momentum of proton with Energy loss correction 1 method", 800, 0.,5000.);
TH1F* hMomCorrPID    = new TH1F("hMomCorrPID","Momentum of proton with Energy loss correction 2 method", 800, 0.,5000.);
TH1F* hphi           = new TH1F("hphi"   ,"Phi angle from detector",360, 0.,360.);
TH1F* hphiPro        = new TH1F("hphiPro","Phi angle of protons"   ,360, 0.,360.);

TH2F* hNtrMult       = new TH2F("hNtrMult","N_{tracks} vs Multiplicity (MDC)",50,0.,50.,300,0.,300.);
TH2F* hptvsphi       = new TH2F("hptvsphi"   ,"phi vs pt"   ,600,0.,2000.,360,0.,360.);
TH2F* hptvstheta     = new TH2F("hptvstheta" ,"theta vs pt" ,600,0.,2000.,180,0.,180.);
TH2F* hphivstheta    = new TH2F("hphivstheta","phi vs theta",180,0., 180.,360,0.,360.);
TH1F* hFlatDiff[11];
for (Int_t im=0;im<11;im++){
    hFlatDiff[im] = new TH1F(Form("hFlatDiff%i",im),Form("hFlatDiff %i",im),100,-2.,2.);
}
TProfile* hQvsM_X[2];
TProfile* hQvsM_Y[2];
TProfile* hQvFW_X[2];
TProfile* hQvFW_Y[2];
for (Int_t iT=0;iT<2;iT++){
    hQvsM_X[iT] = new TProfile(Form("hQvsM_X%i",iT),Form("hQvsM_X %i",iT),300, 0., 300.);
    hQvsM_Y[iT] = new TProfile(Form("hQvsM_Y%i",iT),Form("hQvsM_Y %i",iT),300, 0., 300.);
    hQvFW_X[iT] = new TProfile(Form("hQvFW_X%i",iT),Form("hQvFW_X %i",iT),100, 0., 100.);
    hQvFW_Y[iT] = new TProfile(Form("hQvFW_Y%i",iT),Form("hQvFW_Y %i",iT),100, 0., 100.);
}
TProfile* CosPsiAB_META[3][2];
TProfile* SinPsiAB_META[3][2];
TProfile* CosPsiAB_FW[  3][2];
TProfile* SinPsiAB_FW[  3][2];
TProfile* CosPsi_META[  3][2][2]; // [TYPE][Harmonic][A/B];
TProfile* SinPsi_META[  3][2][2]; // [TYPE][Harmonic][A/B];
TProfile* CosPsi_FW[    3][2][2];
TProfile* SinPsi_FW[    3][2][2];
for (Int_t iT=0;iT<3;iT++){
    for(Int_t n=0; n<2; n++){
        CosPsiAB_META[iT][n] = new TProfile(Form("CosPsiAB_META%i%i",iT,n),Form("CosPsiAB_META %i %i",iT,n),300,0.,300.);
        SinPsiAB_META[iT][n] = new TProfile(Form("SinPsiAB_META%i%i",iT,n),Form("SinPsiAB_META %i %i",iT,n),300,0.,300.);
        CosPsiAB_FW[  iT][n] = new TProfile(Form("CosPsiAB_FW%i%i"  ,iT,n),Form("CosPsiAB_FW %i %i"  ,iT,n),120,0.,120.);
        SinPsiAB_FW[  iT][n] = new TProfile(Form("SinPsiAB_FW%i%i"  ,iT,n),Form("SinPsiAB_FW %i %i"  ,iT,n),120,0.,120.);
        for (Int_t k=0;k<2;k++){
            CosPsi_META[iT][n][k] = new TProfile(Form("CosPsi_META%i%i%i",iT,n,k),Form("CosPsi_META %i %i %i",iT,n,k),300,0.,300.);
            SinPsi_META[iT][n][k] = new TProfile(Form("SinPsi_META%i%i%i",iT,n,k),Form("SinPsi_META %i %i %i",iT,n,k),300,0.,300.);
            CosPsi_FW[  iT][n][k] = new TProfile(Form("CosPsi_FW%i%i%i"  ,iT,n,k),Form("CosPsi_FW %i %i %i"  ,iT,n,k),120,0.,120.);
            SinPsi_FW[  iT][n][k] = new TProfile(Form("SinPsi_FW%i%i%i"  ,iT,n,k),Form("SinPsi_FW %i %i %i"  ,iT,n,k),120,0.,120.);
        }
    }
}
TH2F*     hMETAvsCent = new TH2F("hMETAvsCent","hMETAvsCent",12 ,1.,12. ,300,0.,300.);
TH2F*     hFWvsCent   = new TH2F("hFWvsCent"  ,"hFWvsCent"  ,12 ,1.,12. ,120,0.,120.);
TH2F*     hMETAvsFW   = new TH2F("hMETAvsFW"  ,"hMETAvsFW"  ,120,0.,120.,300,0.,300.);