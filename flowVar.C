//-------------------------Set-constants-and-variables---------------------------
  //----------------------------------------------------------------------------------------------------------------
  const Int_t par              = 0; //0-for-debug-plots, 1-for-PsiEP, 2-for-PsiEP-Recentered, 3-for-PsiEP-Flattened;
  //----------------------------------------------------------------------------------------------------------------
  const Int_t ncent            = 12;
  const Int_t NC               = ncent-1;
  const Int_t NRAPI            = 9;
  const Int_t NPTRN            = 18;
  const Int_t QvRING           = 4;
  const Int_t QvWGHT           = 2;
  const Int_t QvTYPE           = 2;
  const Float_t rad2deg        = 57.2957795130823229;
  const Float_t fTrkMultDay108 =48.78; //-average-EventHeaderTrackingMult-at-the-good-time-of-day-108---
  vector<Float_t> mulVal;
  const Float_t MeanQx[6]      = { 1.4681725259,2.0612748086,0.3247944454,-28.2572955981,-1.3244127418,0.0 };
  const Float_t MeanQy[6]      = {-1.1507586466,1.6045932655,-2.6009353432,0.4215789944,-3.0768050756, 0.0 };
  Float_t X_shM[ncent]         = {  0.,0.,0.,0.,0.,0.,0.,0.,0., 0.,0.,0.  }; 
  Float_t Y_shM[ncent]         = {  0.,0.,0.,0.,0.,0.,0.,0.,0., 0.,0.,0.  };
  const Int_t   Mrang[NC+1]    = {20,30,39,50,59,72,85,103,123,144,163,215};
  const Float_t Yrang[NRAPI+1] = {-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9};
  const Float_t Prang[NPTRN+1] = {200.,300.,400.,500.,600.,700.,800.,900.,1000.,1100.,1200.,1300.,1400.,1500.,1600.,1700.,1800.,1900.,2000.};

  const Int_t   MDCrn[ncent]   = { 0, 0,13,17,22,28,35, 44, 55, 68, 84,140};
  //-------------------------------------------------------------------------------
  
  //----Define-and-preset-variables-----------------
  
  Double_t hpi                 = TMath::Pi()/2.;
  Float_t Ycm                  = 0.740151;  //CM-rapidity for 1.23GeV/c nucleon beam
  Float_t tgtChi2              = 0.0;
  
  Float_t tgtXc                = 0.13; //gen2-be121040{7,8,9,10}
  Float_t tgtYc                = 0.82; //gen2-be121040{7,8,9,10}
  
  Float_t PCr                  = 0.0;
  Float_t PCz                  = 0.0;
  
  Int_t Mtof                   = 0.0;
  Int_t Mrpc                   = 0.0;
  Int_t Mult                   = 0.0;
  
  Int_t   sys                  = 2;   // sys == 0 - Rpc, sys == 1 - Tof
  Float_t mom                  = 0.0; //-later-shall-use: getMomentumPID()---(c)-Alexandr
  Int_t   charge               = 0;
  Float_t metaQa               = 0.0; 
  Float_t beta                 = 0.0;
  Float_t theta                = 0.0;
  Float_t phi                  = 0.0;
  Float_t chi2                 = 0.0;
  Int_t   sec                  = 0.0; // sector is the sixth part of the detector
  Float_t chi2In               = 0.0; // chi2 of inner segment
  Float_t chi2Out              = 0.0; // chi2 of outer segment
  Float_t mass2                = 0.0; // mass2 = mom*mom*(1-beta*beta)/(beta*beta);
  Float_t mass                 = 0.0; // mass  = sqrt(mass2);
  Float_t mdcdEdx              = 0.0;
  Float_t tofdEdx              = 0.0;
  Float_t metaDx               = 0.0; // ?
  Float_t metaDy               = 0.0;
  Float_t PCp     ;
  Float_t PCpz    ;
  Float_t PCpt    ;
  Float_t PCE     ;
  Float_t PCY     ;
  Float_t PCYn    ;
  Float_t PCYo    ;
  Float_t rapidity = 0.0;
  
  //----Boundary-for-nuclei-radius------------------
  Float_t R2   = 12.0;
  
  Int_t        NL              = 0;
  Int_t        NN              = 0;
  Int_t        NM              = 0;
  Float_t      QLx             =-1000.;
  Float_t      QLy             =-1000.;
  Float_t      QMx             =-1000.;
  Float_t      QMy             =-1000.;
  Float_t      QNx             =-1000.;
  Float_t      QNy             =-1000.;
  Float_t      FWdEdxL         = 0.;
  Float_t      FWdEdxM         = 0.;
  Float_t      FWdEdxN         = 0.;
  Float_t      FWdEdxA         = 0.;
  Float_t      FWdEdxB         = 0.;
  
  Float_t      mass2Pip        = HPhysicsConstants::mass(8 )*HPhysicsConstants::mass(8 );
  Float_t      mass2Pim        = HPhysicsConstants::mass(9 )*HPhysicsConstants::mass(9 );
  Float_t      mass2Kap        = HPhysicsConstants::mass(11)*HPhysicsConstants::mass(11);
  Float_t      mass2Pro        = HPhysicsConstants::mass(14)*HPhysicsConstants::mass(14);
  
  TLorentzVector vector;
  Double_t vertexX             = 0.0;
  Double_t vertexY             = 0.0;
  Double_t vertexZ             = 0.0;
  
  //for-HitWall-procedure:---------------------------------
  Int_t   cellNum              = 0;
  Float_t cellCharge           = 0.0;
  Float_t cellTime             = 0.0;
  //-------------------------------------------------------
  //For-Event-plane-reconstruction-------------------------
  static HWallFiredCellsVA cellsVect;
      TRandom3 Random;
      Int_t  multWall=0;
      Float_t wallX=0.0,wallY=0.0,wallZ=0.0;
      Float_t wallXc=0.0, wallYc=0.0; //corrected-reCentered-
      Float_t XfwSmear=0.0, YfwSmear=0.0;
      TVector2 Qvect[QvRING][QvWGHT][QvTYPE];
      TVector2 Qvsum[QvRING][QvWGHT][QvTYPE];
      Bool_t   FWRing[QvRING];
      TVector2 vect(0.,0.);
      TVector2 vsum(0.,0.);
      TVector2 vsumCorr(0.,0.);
      TVector2 vsumCorrA(0.,0.);
      TVector2 vsumCorrB(0.,0.);
      TVector2   eX(1.,0.);
      Float_t dEdxCut=0.0;
      Float_t xyRadius=0.0;
      Double_t  Qx[6],  Qy[6], QxRec[6], QyRec[6]; //-Qnx,y--zero_th-is-not-used-:)---
      Double_t Qax[6], Qay[6]; //-QAx,y--zero_th-is-not-used-:)---
      Double_t Qbx[6], Qby[6]; //-QBx,y--zero_th-is-not-used-:)---
      Float_t phiA             = -1000;
      Float_t phiB             = -1000;
      Float_t phiAB            = -1000;
      Float_t phiCorA          = -1000;
      Float_t phiCorB          = -1000;
      Float_t phiCorAB         = -1000;
      Float_t phiEP[ 6], phiEPrec[6]; 
      Float_t phiEPr[6];
      Float_t dphi, dphiRec, dphiFlt;
      //-weight-for-scalar-product-method-------------------
      Float_t wgh              = 1.0;  //PCpt; //1.0; //-or-could-be-also-equal-to-Pt---
      Int_t nFWhits            = 0;
      Int_t nFWspect           = 0;
      Int_t nFWunderflow       = 0;
      Int_t choiceA            = 1; //
      Int_t choiceB            = 1; // Preparing for A/B subevent method
      Float_t    cellChargeCtime;
      Float_t T1_cut           =   22.0; 
      Float_t T2_cut           =   30.0;
      Float_t Z1_cut_s         =   83.0; 
      Float_t Z1_cut_m         =   84.0; 
      Float_t Z1_cut_l         =   88.0; 
      Float_t X_shift          =  -16.1; 
      Float_t Y_shift          =    7.2;
      Float_t R0_cut           =    0.0; //85.0; 
      Float_t wPhiRPA          =    1.0;
      Int_t   Mcent[ 9]        = {0,20,50,90,130,160,190,220,600};
      Double_t VectPhi_i;
      Int_t multFWcells;
      TVector2 vectA(0.,0.), vectB(0.,0.);
      TVector2 vsumA(0.,0.), vsumB(0.,0.);
      Int_t NA                 = 0;
      Int_t NB                 = 0;
      Float_t levelA, levelB;
      Int_t Mfw                = 0;
      Int_t cellSize           = 0;
      Int_t binX, binY;
      Float_t TrackEff;
      Float_t VectphiEP, VectphiEPr, VectphiCorr, VectphiCorrR;
  //-------------------------------------------------------
  
  //For-looping-over-events--------------------------------
  Int_t          size          = 0;
  //-------------------------------------------------------
  
  Int_t multCand               = 0;
  
  //-------------------------------------------------------
  
  Int_t AbsMinute=0;
  Int_t HR=0, MN=0, SC=0, DAY_NUM;
  TString currFName, currBeName, currFDay, currTime;
  Float_t fAverageMTS[6];
  Float_t fTrkMultScaler=1.0; //-by-default-no-correction-multiplier--
  //-------------------------------------------------------
  
  Int_t CENT                   = -1;
  Int_t RAPI                   = -1;
  Int_t PTRN                   = -1;
  Int_t   nTrkMult, Mtr;
  Float_t nTrkMultCorr, MtrC;
  //Debug & fixing---------------------------------------------------------------------------------
  Float_t pt0,y0;
  Float_t PsiA, PsiB, PsiAB;
  Float_t phiWEP,avSin,avCos,dPsi = 0.0;
  Float_t PsiCorr              = 0.0;
  Float_t PsiCorr2             = 0.0;
  Float_t PsiFOPI              = 0.0;
  Float_t FlatCos[3][6][11][13]; Float_t FlatSin[3][6][11][13];
  Float_t sumXmean[3][11][13];   Float_t sumYmean[3][11][13];
  Float_t sumXsigma[3][11][13];  Float_t sumYsigma[3][11][13];
  Float_t betaMeanP;
  Bool_t pidFlag=kFALSE;
  Float_t Ntrack;
  Float_t wmod=0., wmodA=0.,wmodB=0.;