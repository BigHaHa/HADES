    
  //-------------------------Set-constants-and-variables---------------------------
  const Int_t npt       = 11;
  const Int_t nrpdt     = 11;
  const Int_t nY0       = 10;
  const Int_t npartid   = 3; //defines the number of particles we need to identify;
  const Int_t ncent     = 12;
  const Int_t NC        = 3;
  const Int_t NRAPI      = 9;
  const Int_t NPTRN      = 18;
    
  const Float_t maxpt   = 0.5   *1000;
  const Float_t minpt   = 0.2   *1000;
  const Float_t rad2deg = 57.2957795130823229;
  const Float_t u_p=0.809605;  //2.08335;  //native scale of momentum as seen from CM-system, so u_p=0.809605 (and not 2.08335 which is seen from lab system)
  const Float_t fTrkMultDay108=48.78;
  vector<Float_t> mulVal;
  
  const Double_t binpt  [npt]  = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2};
  
  const Double_t binrpdt[nrpdt]= {0.01,0.21,0.41,0.61,0.81,1.01,1.21,1.41,1.61,1.81,2.01};
  const Float_t  binY0  [nrpdt]= {-0.9,-0.7,-0.5,-0.3,-0.1,0.1,0.3,0.5,0.7,0.9};

  const Float_t MeanQx[6] = { 1.4681725259,2.0612748086,0.3247944454,-28.2572955981,-1.3244127418,0.0 };
  const Float_t MeanQy[6] = {-1.1507586466,1.6045932655,-2.6009353432,0.4215789944,-3.0768050756, 0.0 };

  Float_t X_shM[ncent]  = {  0.,0.,0.,0.,0.,0.,0.,0.,0., 0.,0.,0.  }; 
  Float_t Y_shM[ncent]  = {  0.,0.,0.,0.,0.,0.,0.,0.,0., 0.,0.,0.  };
  const Int_t   Mrang[ncent]  = {20,30,39,50,59,72,85,103,123,144,163,215};


  
  Int_t binsxmom = 6000;
  Int_t binsbeta = 600;
  //-------------------------------------------------------------------------------
  

  //----Define-and-preset-variables-----------------
  
  Double_t hpi   = TMath::Pi()/2.;
  Float_t Ycm    = 0.740151;  //CM-rapidity for 1.23GeV/c nucleon beam
  Float_t Y0     = 0.0;
  Float_t tgtChi2= 0.0;
  
  Float_t tgtXc  = 0.13; //gen2-be121040{7,8,9,10}
  Float_t tgtYc  = 0.82; //gen2-be121040{7,8,9,10}
  
  Float_t PCr    = 0.0;
  Float_t PCz    = 0.0;
  
  Int_t Mtof     = 0.0;
  Int_t Mrpc     = 0.0;
  Int_t Mult     = 0.0;
  
  Int_t   sys      = 2;   // sys == 0 - Rpc, sys == 1 - Tof
  Float_t mom      = 0.0; //-later-shall-use: getMomentumPID()---(c)-Alexandr
  Int_t   charge   = 0;
  Float_t metaQa   = 0.0; 
  Float_t beta     = 0.0;
  Float_t theta    = 0.0;
  Float_t phi      = 0.0;
  Float_t chi2     = 0.0;
  Int_t   sec      = 0.0; // sector is the sixth part of the detector
  Float_t chi2In   = 0.0; // chi2 of inner segment
  Float_t chi2Out  = 0.0; // chi2 of outer segment
  Float_t mass2    = 0.0; // mass2 = mom*mom*(1-beta*beta)/(beta*beta);
  Float_t mass     = 0.0; // mass  = sqrt(mass2);
  Float_t mdcdEdx  = 0.0;
  Float_t tofdEdx  = 0.0;
  Float_t metaDx   = 0.0; // ?
  Float_t metaDy   = 0.0;

  Float_t PCp     ;
  Float_t PCpz    ;
  Float_t PCpt    ;
  Float_t PCE     ;
  Float_t PCY     ;
  Float_t PCYn    ;
  Float_t PCYo    ;

  Float_t P_beta_t;
  Float_t P_u_t   ;
  Float_t P_u_t0  ;
  
  Float_t pz       = 0.0;
  Float_t pt0      = 0.0;
  Float_t px       = 0.0;
  Float_t py       = 0.0;
  
  Float_t rapid    = 0.0;
  Float_t eta      = 0.0;
  
  Float_t betaL    = 0.0;
  Float_t rapidity = 0.0;
  
  //----Boundary-for-nuclei-radius------------------
  Float_t R2   = 12.0;
  
  Int_t        NL     = 0;
  Int_t        NN     = 0;
  Int_t        NM     = 0;
  Float_t      QLx    =-1000.;
  Float_t      QLy    =-1000.;
  Float_t      QMx    =-1000.;
  Float_t      QMy    =-1000.;
  Float_t      QNx    =-1000.;
  Float_t      QNy    =-1000.;
  Float_t      FWdEdxL= 0.;
  Float_t      FWdEdxM= 0.;
  Float_t      FWdEdxN= 0.;
  Float_t      FWdEdxA= 0.;
  Float_t      FWdEdxB= 0.;
  
  Float_t      mass2Pip=HPhysicsConstants::mass(8 )*HPhysicsConstants::mass(8 )/1e6;
  Float_t      mass2Pim=HPhysicsConstants::mass(9 )*HPhysicsConstants::mass(9 )/1e6;
  Float_t      mass2Kap=HPhysicsConstants::mass(11)*HPhysicsConstants::mass(11)/1e6;
  Float_t      mass2Pro=HPhysicsConstants::mass(14)*HPhysicsConstants::mass(14)/1e6;
  
  TLorentzVector vector;
  Float_t      pt = 0.0;
  Float_t      y  = 0.0;  
  
  Long_t       npip = 1;
  Long_t       npim = 1;
  Long_t       nkap = 1;
  Long_t       npro = 1;
  
  Double_t vertexX = 0.0;
  Double_t vertexY = 0.0;
  Double_t vertexZ = 0.0;
  
  //for-HitWall-procedure:---------------------------------
  Int_t   cellNum    = 0;
  Float_t cellCharge = 0.0, cellTime = 0.0, cellPhi = 200.;
  //-------------------------------------------------------

  //For-Event-plane-reconstruction-------------------------
  static HWallFiredCellsVA cellsVect;
      TRandom3 Random;
      Int_t  multWall=0;
      Float_t wallX=0.0,wallY=0.0,wallZ=0.0;
      Float_t wallXc=0.0, wallYc=0.0; //corrected-reCentered-
      Float_t XfwSmear=0.0, YfwSmear=0.0;
      TVector2 vect(0.,0.);
      TVector2 vsum(0.,0.);
      TVector2 vsumCorr(0.,0.);
      TVector2   eX(1.,0.);
      Float_t dEdxCut=0.0;
      Float_t xyRadius=0.0;

      Double_t  Qx[6],  Qy[6], QxRec[6], QyRec[6]; //-Qnx,y--zero_th-is-not-used-:)---
      Double_t Qax[6], Qay[6]; //-QAx,y--zero_th-is-not-used-:)---
      Double_t Qbx[6], Qby[6]; //-QBx,y--zero_th-is-not-used-:)---

      Float_t phiA       = -1000;
      Float_t phiB       = -1000;
      Float_t phiAB      = -1000;
      Float_t phiEP[ 6], phiEPrec[6]; 
      Float_t phiEPr[6];
      Float_t dphi;
      //-weight-for-scalar-product-method-------------------
      Float_t wgh       = 1.0;  //PCpt; //1.0; //-or-could-be-also-equal-to-Pt---

      Int_t nFWhits      = 0;
      Int_t nFWspect     = 0;
      Int_t nFWunderflow = 0;

      Int_t choiceA      = 1; //
      Int_t choiceB      = 1; // Preparing for A/B subevent method

      Float_t    cellChargeCtime;

      Float_t T1_cut     =   22.0; 
      Float_t T2_cut     =   30.0;
      Float_t Z1_cut_s   =   83.0; 
      Float_t Z1_cut_m   =   84.0; 
      Float_t Z1_cut_l   =   88.0; 
      Float_t X_shift    =  -16.1; 
      Float_t Y_shift    =    7.2;
      Float_t R0_cut     =    0.0; //85.0; 
      Float_t wPhiRPA    =    1.0;

      Int_t   Mcent[ 9]  = {0,20,50,90,130,160,190,220,600};

      Double_t VectPhi_i;

      Int_t multFWcells;

      TVector2 vectA(0.,0.), vectB(0.,0.);
      TVector2 vsumA(0.,0.), vsumB(0.,0.);
      Int_t NA           = 0;
      Int_t NB           = 0;
      Float_t levelA, levelB;
      Int_t Mfw          = 0;
      Int_t cellSize     = 0;
      Int_t binX, binY;
      Float_t TrackEff;

      //My experiment
      Float_t Phi_EP[6];
      Float_t v2;
      Float_t v2sum [npt-1];
      Int_t   ncount[npt-1];
      Float_t v2w   [npt-1];
      Float_t v2x   [npt-1];
      Float_t v2e   [npt-1];
      Float_t VectphiEP, VectphiEPr;
  //-------------------------------------------------------
  
  //For-looping-over-events--------------------------------
  Int_t          size    = 0;
  //-------------------------------------------------------
  
  Int_t multCand         = 0;
  
  //-------------------------------------------------------
  
  Int_t AbsMinute=0;
  Int_t HR=0, MN=0, SC=0, DAY_NUM;
  TString currFName, currBeName, currFDay, currTime;
  Float_t fAverageMTS[6];
  Float_t fTrkMultScaler=1.0; //-by-default-no-correction-multiplier--

  //-------------------------------------------------------

  Bool_t FOPI_blue=kFALSE;
  Bool_t FOPI_gren=kFALSE;
  Bool_t FOPI_redd=kFALSE;
  
  Int_t CENT=-1;
  Int_t RAPI=-1;
  Int_t PTRN=-1;

  Int_t   nTrkMult, Mtr;
  Float_t nTrkMultCorr, MtrC;
  