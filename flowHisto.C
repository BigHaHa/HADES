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
TF1  *fRPAy[ 11];
TH1F *hRPAy[ 11];
TH1F *hRPAyc[11];
TH1F *hORPAy[ 11];
TH1F *hORPAyc[11];
TH1F* hWEPy[11];
TH1F* hOWEPy[11];
TProfile* hSinPsi[6][11];
TProfile* hCosPsi[6][11];
TH2F* hSinFOPI[6][11];
TH2F* hCosFOPI[6][11];
TH1F* hPsiEP[   11];
TH1F* hPsiRcnt[ 11];
TH1F* hPsiCorr[ 11];
TH1F* hOPsiCorr[11];
TProfile* hsumXmean[11];
TProfile* hsumYmean[11];
TH1F* hQvectX[  11];
TH1F* hQvectY[  11];
TH1F* hQvXrec[  11];
TH1F* hQvYrec[  11];
TH2F* hQvRaw[11];
TH2F* hQvRec[11];
TH1F* hdPsi[11];
for(Int_t i=0;i<11;i++){
    //---angular-weight-correction-function--[cos1A+cos2A+sin1A +cos3A + sin2A ]-----------------
    fRPAy[i]  = new TF1( Form("fRPAy%i",i)  ,  "( 1.+2.*[0]*TMath::Cos(x[0]*3.1415926535/180.)+2.*[1]*TMath::Cos(2.*x[0]*3.1415926535/180.)+2.*[2]*TMath::Cos(3.*x[0]*3.1415926535/180.)+2.*[3]*TMath::Cos(4.*x[0]*3.1415926535/180.) + [4]*TMath::Sin(x[0]*3.1415/180.) + [5]*TMath::Sin(2.*x[0]*3.1415/180.) )", -180., 180.);
    //--histograms-for-angular-weighting-correction--(before-correction)-------------------------
    hRPAy[i]  = new TH1F(Form("hRPAy%i",i)  ,Form("hRPAy %i range"            ,i) , 90, -180., 180);
    //--histograms-for-angular-weighting-correction-(after-correction)---------------------------
    hRPAyc[i] = new TH1F(Form("hRPAc%ic",i) ,Form("hRPAyc %i range_after_corr",i) , 90, -180., 180);
    hORPAy[i] = new TH1F(Form("hORPAy%i",i) ,Form("hORPAy %i range"            ,i), 90, -180., 180.);
    hORPAyc[i]= new TH1F(Form("hORPAc%ic",i),Form("hORPAyc %i range_after_corr",i), 90, -180., 180.);
    hWEPy[i]  = new TH1F(Form("hWEPy%i",i)   ,Form("hWEP %i META",i)              , 90, -180., 180.); //hWEPy[ i]->Sumw2();
    hOWEPy[i] = new TH1F(Form("hOWEPy%i",i)  ,Form("hOWEP %i MDC",i)              , 90, -180., 180.); //hOWEPy[i]->Sumw2();
    hsumXmean[i] = new TProfile(Form("hsumXmean%i",i),Form("hsumXmean %i",i),15,95.,110.);
    hsumYmean[i] = new TProfile(Form("hsumYmean%i",i),Form("hsumYmean %i",i),15,95.,110.);
    for (Int_t n=0;n<6;n++){
    	hSinPsi[n][i] = new TProfile(Form("hSinPsi%i%i",n,i),Form("hSinPsi %i %i",n,i),15,95.,110.);
        hCosPsi[n][i] = new TProfile(Form("hCosPsi%i%i",n,i),Form("hCosPsi %i %i",n,i),15,95.,110.);
    }
    hPsiEP[   i] = new TH1F(Form("hPsiEP%i"   ,i),Form("hPsiEP %i "   ,i), 90, -180., 180.);
    hPsiRcnt[ i] = new TH1F(Form("hPsiRcnt%i" ,i),Form("hPsiRcnt %i after_recentering" ,i), 90, -180., 180.);
    hPsiCorr[ i] = new TH1F(Form("hPsiCorr%i" ,i),Form("hPsiCorr %i after_flattering" ,i), 90, -180., 180.);
    hdPsi[    i] = new TH1F(Form("hdPsi%i" ,i),Form("hdPsi %i after_flattering" ,i), 100, -5., 5.);
    hOPsiCorr[i] = new TH1F(Form("hOPsiCorr%i",i),Form("hOPsiCorr %i after_corr",i), 90, -180., 180.);
    hQvectX[  i] = new TH1F(Form("hQvectX%i",i),Form("hQvectX %i FW",i), 600, -0.05, 0.05);
    hQvectY[  i] = new TH1F(Form("hQvectY%i",i),Form("hQvectY %i FW",i), 600, -0.05, 0.05); 
    hQvXrec[  i] = new TH1F(Form("hQvXrec%i",i),Form("hQvXrec %i FW",i), 600, -5., 5.);
    hQvYrec[  i] = new TH1F(Form("hQvYrec%i",i),Form("hQvYrec %i FW",i), 600, -5., 5.);
    hQvRaw[   i] = new TH2F(Form("hQvRaw%i",i) ,Form("Qvect %i Raw",i), 1000, -0.05, 0.05, 1000, -0.05, 0.05);
    hQvRec[   i] = new TH2F(Form("hQvRec%i",i) ,Form("Qvect %i Rec",i), 1000, -5., 5., 1000, -5., 5.);
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
TProfile* hQvsM_X = new TProfile(Form("hQvsM_X"),Form("hQvsM_X"),300, 0., 300.);
TProfile* hQvsM_Y = new TProfile(Form("hQvsM_Y"),Form("hQvsM_Y"),300, 0., 300.);
TProfile* hQvFW_X = new TProfile(Form("hQvFW_X"),Form("hQvFW_X"),100, 0., 100.);
TProfile* hQvFW_Y = new TProfile(Form("hQvFW_Y"),Form("hQvFW_Y"),100, 0., 100.);
