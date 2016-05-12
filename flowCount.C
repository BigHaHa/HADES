#define DST_GEN 8
#include "hades.h"
#include "htool.h"
#include "hphysicsconstants.h"
#include "hrootsource.h"
#include "hiterator.h"
#include "hloop.h"
#include "hdst.h"
#include "haddef.h"
#include "heventheader.h"
#include "hparticledef.h"
#include "hparticlestructs.h"
#include "hparticlecand.h"
#include "hparticlecandsim.h"
#include "hcategorymanager.h"
#include "hparticletracksorter.h"
#include "hparticlet0reco.h"
#include "hparticlevertexfind.h"
#include "hparticleevtinfo.h"
#include "htaskset.h"
#include "henergylosscorrpar.h"
#include "hstart2cal.h"
#include "hstart2hit.h"
#include "hmdcdef.h"
#include "hmdctrackddef.h"
#include "hmdctrackgdef.h"
#include "horadef.h"
#include "horasimdef.h"
#include "horasimdef.h"
#include "hstartdef.h"
#include "richdef.h"
#include "rpcdef.h"
#include "showerdef.h"
#include "simulationdef.h"
#include "tofdef.h"
#include "walldef.h"
#include "hwallhit.h"
#include "htofhit.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TRandom.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TStopwatch.h"
#include "TCutG.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TIterator.h"
#include "htofhit.h"
#include "hrpchit.h"
#include "hrpccluster.h"
#include "hwalleventplane.h"
#include <iostream>
#include <fstream>
#include <map>
#include <stdio.h>
#include <iomanip>
#include <vector>
#include <TH3F.h>
#include <TProfile.h>
#include "/u/parfenov/anaflow/TRandom3.h"          //need because we should not use TRandom::Rndm() 
#include "/u/parfenov/anaflow/EventCountClasses.h" //include all classes and functions we need
using namespace std;
Int_t flowCount(TString inputlist, TString outfile, Int_t nev=-1)
{
    HLoop* loop = new HLoop   (kTRUE    );    // kTRUE : create Hades  (needed to work with standard eventstructure)
           loop ->addMultFiles(inputlist);
    //-----Inicialize-classes-we-need-to-use-----------------
    HMultCorr   cCorr;
    PidParticle p;
    //--------------------------------------------------------------------------------------------------------
    #include "/u/parfenov/anaflow/flowVar.C"     //include all constants and variables we need          
    TH2F*  hProtPtVsY_Eff;
    TFile* foutfile       = new TFile(outfile.Data(), "recreate"); 
    TFile* FileProtEff    = new TFile("EffCorrForProtPtVsY_WithWideBetaPcuts_ShieldGosia.root", "read");
           hProtPtVsY_Eff = (TH2F*) FileProtEff->Get("hProtPtVsY_EffSm");
           foutfile->cd();
    #include "/u/parfenov/anaflow/flowHisto.C"   //include all histograms we need                       
    //if(par!=0){
        FFlow pfy0pt[NC-4][NRAPI][NPTRN-5]; //--flow--9-bins-Yo--12-bins-in-Pt-in-3-centrality-sets--
        FFlow pfR0pt[NC-4][NRAPI][NPTRN-5];
        FFlow pfF0pt[NC-4][NRAPI][NPTRN-5];
        FHist phist[ NC-4][NRAPI][NPTRN-5];
        for (Int_t i=0; i<NC-4;i++){
            for(Int_t j=0; j<NRAPI; j++){
                for(Int_t k=0; k<NPTRN-5; k++){
                    //-proton------------------------------------------
                    pfy0pt[i][j][k].THDeclare("pCent",i,"Yo",j,"Pt",k);
                    pfR0pt[i][j][k].THDeclare("RCent",i,"Yo",j,"Pt",k);
                    pfF0pt[i][j][k].THDeclare("FCent",i,"Yo",j,"Pt",k);
                    phist[ i][j][k].THDeclare("pHist",i,"Yo",j,"Pt",k);
                }
            }
        }
    //}

    fBetaMeanP->SetParameter(0,HPhysicsConstants::mass(14));
    //--------------------------------------------------------------------------------------------------------
    if(!loop->setInput("-*,+HParticleCand,+HParticleEvtInfo,+HStart2Hit,+HStart2Cal,+HWallHit,+HRpcCluster,+HTofHit,+HWallEventPlane")) { exit(1); }
    loop->printCategories();
  
    HParticleEvtInfo* evtinfo;
  
    HParticleCand*   pParticleCand;   HCategory* candCat   = (HCategory*)HCategoryManager::getCategory(catParticleCand);   if(!candCat)  {exit(1);}
    HWallHit*        pWallHit;        HCategory* wallCat   = (HCategory*)HCategoryManager::getCategory(catWallHit);        if(!wallCat)  {exit(1);}
    HWallEventPlane* pWallEventPlane; HCategory* wallEPCat = (HCategory*)HCategoryManager::getCategory(catWallEventPlane); if(!wallEPCat){exit(1);}
    
    //--------------------------CONFIGURATION---------------------------------------------------
    //At begin of the program (outside the event loop)
    HParticleTrackSorter sorter;
    //sorter.setDebug();                          // for debug
    //sorter.setPrintLevel(3);                    // max prints
    //sorter.setIgnoreInnerMDC();                 // do not reject Double_t inner MDC hits
    //sorter.setIgnoreOuterMDC();                 // do not reject Double_t outer MDC hits
    //sorter.setIgnoreMETA();                     // do not reject Double_t META hits
    //sorter.setIgnorePreviousIndex();            // do not reject indices from previous selctions
    sorter.init();                                // get catgegory pointers etc...
    //--------------------------------------------------------------------------------------------

    HEnergyLossCorrPar *momCorPar = new HEnergyLossCorrPar();
    momCorPar->setDefaultPar("apr12");
      
    /////////////////////////////////////////////////////////
    //-----Loop-over-events----------------------------------
    /////////////////////////////////////////////////////////
  
    loop->addMultFiles(inputlist);
    Int_t events = loop->getEntries();
    if(nev < events && nev >= 0 ) events = nev;
  
    for (Int_t i=1; i<events;i++){
    
        if (loop->nextEvent(i)<=0){
            cout << "End recieved with IF exit" << endl;
            break;
        }
        if(par==0) hcuts->Fill(1.);
        if (i%1000 == 0){
            cout << "    event " << i << endl;
        }
    
        if( loop->isNewFile(currFName) ){ 
            printf("--New-file-started--> %s", currFName.Data() ); 
            currFName = currFName( currFName.Last('/')+1,currFName.Length()-currFName.Last('/')-1 );
            currBeName = currFName(0, 13); //-simply-the-file-name-w/o-EBchunk--
            currFDay = currFName( 4, 3 ); //-we-cut-out-day-number-position-starting-from-4th-digit-with-length-of-3digits--
            DAY_NUM = currFDay.Atoi();
            currTime = currFName( 7, 2 ); //-we-cut-out-day-hour--position-starting-from-4th-digit-with-length-of-3digits--
            HR=currTime.Atoi();
            currTime = currFName( 9, 2 ); //-we-cut-out-day-minute-position-starting-from-4th-digit-with-length-of-3digits--
            MN=currTime.Atoi();
            currTime = currFName(11, 2 ); //-we-cut-out-day-second-position-starting-from-4th-digit-with-length-of-3digits--
            SC=currTime.Atoi();
            AbsMinute = MN + HR*60 + DAY_NUM*24*60;
            #include "/u/parfenov/anaflow/FlatSinCos.cc"
            //#include "/u/parfenov/anaflow/FOPICorPar.cc"
            #include "/u/parfenov/anaflow/Recenter.cc"
            #include "/u/parfenov/anaflow/RecenterFW.cc" //recenter FW;
            //#include "/u/parfenov/anaflow/FlatFourier.cc"//for flattening Psi_EP via fit;
            //--Now-we-read-multiplicity-parameters-for-this-beXXXXXXXXXX-file--------------------//
            mulVal = cCorr.getLineValuesAsVectFromCalibFile( currBeName );
            cout<<"mulVal "<< mulVal[0] <<" "<< mulVal[1] <<" "<< mulVal[2] <<" "<< mulVal[3] <<endl;
            fAverageMTS[0]=mulVal[4];  //multiplicity in tracking sector 1
            fAverageMTS[1]=mulVal[5];  //multiplicity in tracking sector 2
            fAverageMTS[2]=mulVal[6];  //...
            fAverageMTS[3]=mulVal[7];
            fAverageMTS[4]=mulVal[8];
            fAverageMTS[5]=mulVal[9];  //multiplicity in tracking sector 6
            printf("fAverageMTS[0,1,2,3,4,5]=[%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,%6.3f,]\n", fAverageMTS[0], fAverageMTS[1], fAverageMTS[2], fAverageMTS[3], fAverageMTS[4], fAverageMTS[5]);
        
            //fTrkMultDay108
            if(mulVal[2]>0.0){
                //-if-correction-is-found-or-if-it-was-not-zero-defined-(in-case-of-not-working-MDCs)---
                fTrkMultScaler = fTrkMultDay108/mulVal[2];
                printf("==>Correction multiplier nTrkMult for this file is: (%f/%f)=%f\n", fTrkMultDay108, mulVal[2], fTrkMultScaler );
            }else{
                //--if-correction-entry-is-not-found-(may-be-not-created-once)------------------
                //--or-if-say-MDCs-were-not-working-and-tracking-multiplicity-was-exectly-zero--
                //-we-keep-correction-as-it-was-in-previous-file-or-if-not-defined-hence-it-would-be-exactly-one=1-or-
                printf("==>Correction multiplier nTrkMult for this file is: (%f/%f)=%f (default or old correction is used)\n", fTrkMultDay108, mulVal[2], fTrkMultScaler);
            }
        }
        HTool::printProgress(i,events,1,"Analyze pairs :");
        evtinfo = HCategoryManager::getObject(evtinfo,catParticleEvtInfo,0);

        if( evtinfo->getSumParticleCandMult()<2 ){if(par==0){ hcuts->Fill(2.); } continue; } //getSumParticleCandMult()>1
        if(!evtinfo->isGoodEvent(
                               Particle::kGoodVertexClust|
                               Particle::kGoodVertexCand|
                               Particle::kNoPileUpSTART|
                               Particle::kGoodSTART|      
                               Particle::kGoodTRIGGER|    
                               Particle::kNoVETO|         
                               Particle::kGoodSTARTVETO|  
                               Particle::kGoodSTARTMETA
                              )){ if(par==0){ hcuts->Fill(3.); } continue ; }// bitwise add flags as suggested by Manuel 

        //------------------------------------------------------------------------
        // clean vectors and index arrays
        sorter.cleanUp();
        //------------------------------------------------------------------------

        sorter.resetFlags(kTRUE,kTRUE,kTRUE,kTRUE);
        Int_t nCandHad        = sorter.fill(selectHadronsQa);
        Int_t nCandHadBest    = sorter.selectBest(HParticleTrackSorter::kIsBestRK,HParticleTrackSorter::kIsHadron);

        Int_t nCandNoLept     = sorter.fill(rejectLeptons);
        Int_t nCandNoLeptBest = sorter.selectBest(HParticleTrackSorter::kIsBestRK,HParticleTrackSorter::kIsLepton);

        //-------Get-vertex-from-combined-fit-of-fitted-inner-segments------------
        HVertex   fPrimVertReco = gHades->getCurrentEvent()->getHeader()->getVertexReco();
        vertexX = fPrimVertReco.getX();
        vertexY = fPrimVertReco.getY();
        vertexZ = fPrimVertReco.getZ();

        if(par==0){
            hvertexXZ->Fill(vertexZ,vertexX);
            hvertexXY->Fill(vertexY,vertexX);
            hvertexX ->Fill(vertexX);
            hvertexY ->Fill(vertexY);
            hvertexZ ->Fill(vertexZ);
        }

        //------Alexandr's-vertex-cuts---------
        tgtChi2 = (Float_t) fPrimVertReco.getChi2();
        if( tgtChi2<2.0 || tgtChi2>40.){ if(par==0){ hcuts->Fill(4.); } continue; }//-skip-events-with-badly-reconstructed-target-------
        if( vertexZ<-59.|| vertexZ>0. ){ if(par==0){ hcuts->Fill(5.); } continue; }//-skip-events-with-primary-vertex-outside-of-target-
        if( (vertexX-tgtXc)*(vertexX-tgtXc)+(vertexY-tgtYc)*(vertexY-tgtYc)>R2 ){ if(par==0){ hcuts->Fill(6.); } continue; }//-skip-events-with-primary-vertex-far-from-target---
        //----end-of-Alexandr's-vertex-cuts----
        if(par==0){ 
            hvtxCutXZ->Fill(vertexZ,vertexX);
            hvtxCutXY->Fill(vertexY,vertexX);
            hvtxCutX ->Fill(vertexX);
            hvtxCutY ->Fill(vertexY);
            hvtxCutZ ->Fill(vertexZ);
        }
                              
        Mtof = evtinfo->getSumTofMult();
        Mrpc = evtinfo->getSumRpcMult();
        Mult = Mtof+Mrpc;
        //---------centrality-classes--------------------------------------------------------------------(begin)--------//
        //-Centrality-selection-for-further-if-selection-statements------------------------//
        nTrkMult = getParticleTrkMult();
        Mtr   = nTrkMult; if(par==0) hTrkMult->Fill(nTrkMult);
        nTrkMultCorr = nTrkMult*fTrkMultScaler; if(par==0) hTrkMultCorr->Fill(nTrkMultCorr);
        if(par==0) hMETAvsTRK    ->Fill(Mtof+Mrpc,nTrkMult    );
        if(par==0) hMETAvsTRKcorr->Fill(Mtof+Mrpc,nTrkMultCorr);

        //--------------------------------------------//
        //--Event plane reconstruction loop-----------//
        //--------------------------------------------// 
        multWall=0;
        cellNum=0;
        cellTime=0.0;
        cellCharge=0.0;
        wallX=0.0,wallY=0.0,wallZ=0.0;
        wallXc=0.0, wallYc=0.0; //corrected-reCentered-
        XfwSmear=0.0, YfwSmear=0.0;
        for (Int_t a=0;a<QvRING;a++){//Ring: 0-full,1-small,2-medium,3-large;
            FWRing[a] = kFALSE;
            for (Int_t b=0;b<QvWGHT;b++){//Weight: 0-wght=1, 1-wght=2;
                for (Int_t c=0;c<QvTYPE;c++){//Type: 0-source from FW, 1-recentered;
                    Qvect[a][b][c].Set(0.,0.);
                    Qvsum[a][b][c].Set(0.,0.);
                }
            }
        }
        vect.Set(0.,0.);
        vsum.Set(0.,0.);
        vsumCorr.Set(0.,0.);
        vsumCorrA.Set(0.,0.);
        vsumCorrB.Set(0.,0.);
        eX  .Set(1.,0.);
        dEdxCut=0.0;
        xyRadius=0.0;
        phiA   = -1000;
        phiB   = -1000;
        phiAB  = -1000;
        phiCorA= -1000;
        phiCorB= -1000;
        PsiA   = -1000;
        PsiB   = -1000;
        cellsVect.Reset();
      
        //-weight-for-scalar-product-method-------------------
        wgh = 1.0;  //PCpt; //1.0; //-or-could-be-also-equal-to-Pt---
        nFWhits   = 0;
        nFWspect  = 0;
        nFWunderflow = 0;
        FWdEdxA=0.0;
        FWdEdxB=0.0;
  
        FWdEdxL=0.0;
        FWdEdxM=0.0;
        FWdEdxN=0.0;
        choiceA = 1; //
        choiceB = 1; // Preparing for A/B subevent method
        wmod    = 0.;
        wmodA   = 0.;
        wmodB   = 0.;
        //-------------------------------------------iteration-WallHits----(begin)--------------------------------------------------//
        for(Int_t j=0; j<(wallCat->getEntries()); j++){
            pWallHit = HCategoryManager::getObject(pWallHit,wallCat,j);
            cellNum    = pWallHit->getCell();
            cellTime   = pWallHit->getTime();
            cellCharge = pWallHit->getCharge();
            cellChargeCtime = cellTime;
            hFW_TimeCell_0->Fill(cellTime  , cellNum);
            hFW_dEdxCell_0->Fill(cellCharge, cellNum);
            if(cellNum>=0   && cellNum<304 ){                   FWRing[0] = kTRUE; }
            if(cellNum< 144                ){ dEdxCut=Z1_cut_s; FWRing[1] = kTRUE; }
            if(cellNum>=144 && cellNum<=207){ dEdxCut=Z1_cut_m; FWRing[2] = kTRUE; }
            if(cellNum> 207                ){ dEdxCut=Z1_cut_l; FWRing[3] = kTRUE; }
            if(cellTime>T1_cut && cellTime<T2_cut && cellCharge>dEdxCut){
           
                nFWspect++;
                hFW_TimeCell_1->Fill(cellTime  , cellNum);
                hFW_dEdxCell_1->Fill(cellCharge, cellNum);
                pWallHit->getXYZLab(wallX,wallY,wallZ);
                hFWxyCC->Fill(wallX,wallY); //cell centers
                
                //-recentering--FW-----------------------------------------------------
                for (Int_t im=1;im<12;im++){
                    if( (Mtof+Mrpc)>=Mrang[im-1] && (Mtof+Mrpc)<Mrang[im]){ wallXc = wallX - X_shM[im]; wallYc = wallY - Y_shM[im]; }
                }
                //-here-we-fill-d2N/(dxdy)-inX-band--and--y-band--for-auto-recentering-
                //-this-makes-realistic-distribution-on-FW-face------------------------
                if (cellNum>=0   && cellNum<=143) { cellSize = 40;  }
                if (cellNum>=144 && cellNum<=207) { cellSize = 80;  }
                if (cellNum>=210 && cellNum<=301) { cellSize = 160; }
                XfwSmear = wallX + ( Random.Rndm(1) - 0.5 )*cellSize;
                YfwSmear = wallY + ( Random.Rndm(1) - 0.5 )*cellSize;
           
                for (Int_t im=0;im<11;im++){                    
                    if( (Mtof+Mrpc)>=  Mrang[im] && (Mtof+Mrpc)< Mrang[im+1] ){
                        hFW_TimeCell_0->Fill(cellTime  , cellNum);
                        hFW_dEdxCell_0->Fill(cellCharge, cellNum);
                    }
                }
                //-this-cut-was-forgotten-for-Feb2014-HADES-CM-report--//
                //-I-reintroduce-it-14.03.2014----(c)-Sadovsky---------//
                if(  wallXc*wallXc + wallYc*wallYc >= R0_cut*R0_cut  /*50.*50.*/ ){
                    hFW_TimeCell_2->Fill(cellTime  , cellNum);
                    hFW_dEdxCell_2->Fill(cellCharge, cellNum);
                    //-spectators-selected--
                    multWall++;
                    /*for (Int_t a=0;a<QvRING;a++){
                       for (Int_t b=0;b<QvWGHT;b++){
                            if (FWRing[a]){
                                Qvect[a][b][0].Set(wallX, wallY);
                                Qvect[a][b][0] = Qvect[a][b][0].Unit();
                                if(b == 1) Qvect[a][b][0] *= cellCharge;
                                Qvsum[a][b][0] = Qvsum[a][b][0] + Qvect[a][b][0];
                            }
                        }
                    }*/
                    vect.Set(wallX, wallY);
                    vect = vect.Unit();
                    vect *= cellCharge;
                    wmod += cellCharge;
                    vsum += vect;
                    cellsVect.SetCellVect(vect, cellCharge); //-note-[vect]-is-a-unit-vector-this-is-for-further-subevent-(A/B)-split--
                    hFWxyCCsmear->Fill(XfwSmear        ,YfwSmear        ); //cells smeared but not shifted
                    hFWxySHsmear->Fill(XfwSmear-X_shift,YfwSmear-Y_shift); //cells smeared and shifted as a whole
                    //-center-of-gravity-study-as-a-course-of-centrality-----------------------------
                    for (Int_t im=0;im<11;im++){
                        if( (Mtof+Mrpc)>=Mrang[im] && (Mtof+Mrpc)<Mrang[im+1]){ 
                            hFWxyCCsmearM[im]->Fill(XfwSmear,YfwSmear); 
                            hFWxyCCsmearXscanM[im]->Fill(XfwSmear - X_shM[im],YfwSmear - Y_shM[im]); 
                        }
                    }
                }//-endif-R0_cut--//
            }//end-if-cells-cut---
        }//end-for-HitWall-loop---
        vsum      /= wmod;
        for (Int_t im=0;im<11;im++){
            if((Mtof+Mrpc)>=Mrang[im] && (Mtof+Mrpc)<Mrang[im+1]){
                vsumCorr = cellsVect.Recenter(vsum,sumXmean[0][im][DAY_NUM-96],sumYmean[0][im][DAY_NUM-96],sumXsigma[0][im][DAY_NUM-96],sumYsigma[0][im][DAY_NUM-96]);
            }
        }
        VectphiEP =  vsum.DeltaPhi(eX)    *rad2deg;
        VectphiEPr=  vsum.DeltaPhi(eX); //radians
        VectphiCorr= vsumCorr.DeltaPhi(eX)*rad2deg;
        VectphiCorrR= vsumCorr.DeltaPhi(eX);//radians
        hPhiEPvect->Fill(VectphiEP);
        //-now-we-go-over-spectators-and-make-calculations-for-A/B-subevents
        NA=0;
        NB=0;
        multFWcells = cellsVect.GetNumbOfCells();
        Float_t choice;
        vectA.Set(0.,0.);
        vectB.Set(0.,0.);
        vsumA.Set(0.,0.);
        vsumB.Set(0.,0.);
        for(Int_t ic=0; ic<multFWcells; ic++){
            choice = Random.Rndm(1);
            //-this-is-(A/B)-subevents-split--(begin)-(c)-Sadovsky-----------------------------
            //-I-can-be-proud-that-my-algorithm-has-more-flat-distribution--:)-(c)-Sadovsky----
            levelA = (multFWcells/2.-choiceA+1.)/Float_t(multFWcells);
            levelB = (multFWcells/2.-choiceB+1.)/Float_t(multFWcells);
            if(choice < levelA/(levelA+levelB)){
                vectA = cellsVect.GetCellVector(ic);
                vsumA += vectA;
                choiceA++;
                NA++;
                vsumA += vectA;
                choiceA++;
                NA++;
                FWdEdxA=FWdEdxA + cellsVect.GetCellCharge(ic); //-total-Eloss-from-all-spectators---
                wmodA += cellsVect.GetCellCharge(ic);
            }else{
                vectB = cellsVect.GetCellVector(ic);
                vsumB += vectB;
                choiceB++;
                NB++;
                vectB = cellsVect.GetCellVector(ic);
                vsumB += vectB;
                choiceB++;
                NB++;
                VectPhi_i = vectB.DeltaPhi(eX); //-Phi-of-i_th-spectator-particle-from-subevent-B--
                FWdEdxB=FWdEdxB + cellsVect.GetCellCharge(ic); //-total-Eloss-from-all-spectators---
                wmodB += cellsVect.GetCellCharge(ic);
            }//-HWallHit-second-loop-for-reaction-plane-resolution--( end )---//
            //-this-is-(A/B)-subevents-split--( end )------------------------------------------
        }//-endfor-multFWcells-loop--for-A/B-method--
        //-calculating-eventplane-angles---------------
        vsumA /= wmodA;
        vsumB /= wmodB;
        for (Int_t im=0;im<11;im++){
            if((Mtof+Mrpc)>=Mrang[im] && (Mtof+Mrpc)<Mrang[im+1]){
                vsumCorrA = cellsVect.Recenter(vsumA,sumXmean[1][im][DAY_NUM-96],sumYmean[1][im][DAY_NUM-96],sumXsigma[1][im][DAY_NUM-96],sumYsigma[1][im][DAY_NUM-96]);
                vsumCorrB = cellsVect.Recenter(vsumB,sumXmean[2][im][DAY_NUM-96],sumYmean[2][im][DAY_NUM-96],sumXsigma[2][im][DAY_NUM-96],sumYsigma[2][im][DAY_NUM-96]);
            }
        }
        phiA    = vsumA.DeltaPhi(eX)            *rad2deg;
        phiB    = vsumB.DeltaPhi(eX)            *rad2deg;
        phiAB   = vsumA.DeltaPhi(vsumB)         *rad2deg;
        phiCorAB= vsumCorrA.DeltaPhi(vsumCorrB) *rad2deg;
        phiCorA = vsumCorrA.DeltaPhi(eX)        *rad2deg;
        phiCorB = vsumCorrB.DeltaPhi(eX)        *rad2deg;
        PsiAB   = phiCorA-phiCorB;
        if (PsiAB> 180.){ PsiAB-=360; }
        if (PsiAB<-180.){ PsiAB+=360; }
        Mfw     = NA+NB;
        if (Mfw > 3 && NA>1 && NB>1){
            for (Int_t im=0;im<11;im++){
                if((Mtof+Mrpc)>=Mrang[im] && (Mtof+Mrpc)<Mrang[im+1]){
                    hQvectX[0][im]->Fill(vsum.X());          hQvectY[0][im]->Fill(vsum.Y());
                    hQvXrec[0][im]->Fill(vsumCorr.X());      hQvYrec[0][im]->Fill(vsumCorr.Y());
                    hQvRaw[0][ im]->Fill(vsum.X(),vsum.Y()); hQvRec[0][ im]->Fill(vsumCorr.X(),vsumCorr.Y());
                }
            }
            for (Int_t n=1;n<=6;n++){
                for (Int_t im=0;im<11;im++){
                    if ((Mtof+Mrpc)>=Mrang[im] && (Mtof+Mrpc)<Mrang[im+1]){
                        hSinPsi[0][n-1][im]->Fill(DAY_NUM,sin(n*VectphiCorrR),1);
                        hCosPsi[0][n-1][im]->Fill(DAY_NUM,cos(n*VectphiCorrR),1);
                        hSinPsi[1][n-1][im]->Fill(DAY_NUM,sin(n*phiCorA*hpi/90.),1);
                        hCosPsi[1][n-1][im]->Fill(DAY_NUM,cos(n*phiCorA*hpi/90.),1);
                        hSinPsi[2][n-1][im]->Fill(DAY_NUM,sin(n*phiCorB*hpi/90.),1);
                        hCosPsi[2][n-1][im]->Fill(DAY_NUM,cos(n*phiCorB*hpi/90.),1);
                    }
                }
            }
            for (Int_t im=0;im<11;im++){
                if((Mtof+Mrpc)>=Mrang[im] && (Mtof+Mrpc)<Mrang[im+1]){
                    dPsi = 0.;
                    for (Int_t n=0;n<6;n++){
                        cellsVect.SetFlatt(n,FlatSin[0][n][im][DAY_NUM-96],FlatCos[0][n][im][DAY_NUM-96]);
                    }
                    PsiCorr = cellsVect.Flattening(VectphiCorr);
                    //hdPsi[im]->Fill(dPsi*rad2deg);
                    //PsiCorr = atan2(sin(VectphiCorrR+dPsi),cos(VectphiCorrR+dPsi))*rad2deg;
                    //PsiCorr2 = VectphiCorr + dPsi;
                    //hFlatDiff[im]->Fill(PsiCorr/PsiCorr2);
                    if (PsiCorr > 180. ) PsiCorr-=180.;
                    if (PsiCorr < -180.) PsiCorr+=180.;
                    hPsiEP[0][im]->Fill(VectphiEP); hPsiRcnt[0][im]->Fill(VectphiCorr); hPsiCorr[0][im]->Fill(PsiCorr);

                    hsumXmean[0][im]->Fill(DAY_NUM,vsum.X(),1); 
                    hsumYmean[0][im]->Fill(DAY_NUM,vsum.Y(),1); 

                    hQvsM_X[0]->Fill(Mtof+Mrpc,vsum.X()); 
                    hQvsM_Y[0]->Fill(Mtof+Mrpc,vsum.Y()); 
                    hQvFW_X[0]->Fill(nFWspect,vsum.X()); 
                    hQvFW_Y[0]->Fill(nFWspect,vsum.Y()); 
                    hQvsM_X[1]->Fill(Mtof+Mrpc,vsumCorr.X()); 
                    hQvsM_Y[1]->Fill(Mtof+Mrpc,vsumCorr.Y()); 
                    hQvFW_X[1]->Fill(nFWspect,vsumCorr.X()); 
                    hQvFW_Y[1]->Fill(nFWspect,vsumCorr.Y()); 
                }
            }
            for (Int_t im=0;im<11;im++){
                if((Mtof+Mrpc)>=Mrang[im] && (Mtof+Mrpc)<Mrang[im+1]){
                    for (Int_t n=0;n<6;n++){
                        cellsVect.SetFlatt(n,FlatSin[1][n][im][DAY_NUM-96],FlatCos[1][n][im][DAY_NUM-96]);
                    }
                    PsiA = cellsVect.Flattening(phiCorA);
                    if (PsiA > 180. ) PsiA-=180.;
                    if (PsiA < -180.) PsiA+=180.;
                    hPsiEP[1][im]->Fill(phiA); hPsiRcnt[1][im]->Fill(phiCorA); hPsiCorr[1][im]->Fill(PsiA);

                    hsumXmean[1][im]->Fill(DAY_NUM,vsumA.X(),1); 
                    hsumYmean[1][im]->Fill(DAY_NUM,vsumA.Y(),1); 
                }
            }
            for (Int_t im=0;im<11;im++){
                if((Mtof+Mrpc)>=Mrang[im] && (Mtof+Mrpc)<Mrang[im+1]){
                    for (Int_t n=0;n<6;n++){
                        cellsVect.SetFlatt(n,FlatSin[2][n][im][DAY_NUM-96],FlatCos[2][n][im][DAY_NUM-96]);
                    }
                    PsiB = cellsVect.Flattening(phiCorB);
                    if (PsiB > 180. ) PsiB-=180.;
                    if (PsiB < -180.) PsiB+=180.;
                    hPsiEP[2][im]->Fill(phiA); hPsiRcnt[2][im]->Fill(phiCorA); hPsiCorr[2][im]->Fill(PsiB);

                    hsumXmean[2][im]->Fill(DAY_NUM,vsumB.X(),1); 
                    hsumYmean[2][im]->Fill(DAY_NUM,vsumB.Y(),1); 
                }
            }

            for(Int_t im=0;im<11;im++){
                if((Mtof+Mrpc)>=Mrang[im] && (Mtof+Mrpc)<Mrang[im+1]){
                    for(Int_t n=0;n<2;n++){
                        CosPsiAB_META[0][n]->Fill(Mtof+Mrpc,cos((n+1)*phiAB   ),1);
                        CosPsiAB_META[1][n]->Fill(Mtof+Mrpc,cos((n+1)*phiCorAB),1);
                        CosPsiAB_META[2][n]->Fill(Mtof+Mrpc,cos((n+1)*PsiAB   ),1);
                        SinPsiAB_META[0][n]->Fill(Mtof+Mrpc,sin((n+1)*phiAB   ),1);
                        SinPsiAB_META[1][n]->Fill(Mtof+Mrpc,sin((n+1)*phiCorAB),1);
                        SinPsiAB_META[2][n]->Fill(Mtof+Mrpc,sin((n+1)*PsiAB   ),1);

                        CosPsiAB_FW[  0][n]->Fill(nFWspect ,cos((n+1)*phiAB   ),1);
                        CosPsiAB_FW[  1][n]->Fill(nFWspect ,cos((n+1)*phiCorAB),1);
                        CosPsiAB_FW[  2][n]->Fill(nFWspect ,cos((n+1)*PsiAB   ),1);
                        SinPsiAB_FW[  0][n]->Fill(nFWspect ,sin((n+1)*phiAB   ),1);
                        SinPsiAB_FW[  1][n]->Fill(nFWspect ,sin((n+1)*phiCorAB),1);
                        SinPsiAB_FW[  2][n]->Fill(nFWspect ,sin((n+1)*PsiAB   ),1);

                        CosPsi_META[0][n][0]->Fill(Mtof+Mrpc,cos((n+1)*phiA   ),1);
                        CosPsi_META[1][n][0]->Fill(Mtof+Mrpc,cos((n+1)*phiCorA),1);
                        CosPsi_META[2][n][0]->Fill(Mtof+Mrpc,cos((n+1)*PsiA   ),1);
                        SinPsi_META[0][n][0]->Fill(Mtof+Mrpc,sin((n+1)*phiA   ),1);
                        SinPsi_META[1][n][0]->Fill(Mtof+Mrpc,sin((n+1)*phiCorA),1);
                        SinPsi_META[2][n][0]->Fill(Mtof+Mrpc,sin((n+1)*PsiA   ),1);

                        CosPsi_FW[  0][n][0]->Fill(nFWspect ,cos((n+1)*phiA   ),1);
                        CosPsi_FW[  1][n][0]->Fill(nFWspect ,cos((n+1)*phiCorA),1);
                        CosPsi_FW[  2][n][0]->Fill(nFWspect ,cos((n+1)*PsiA   ),1);
                        SinPsi_FW[  0][n][0]->Fill(nFWspect ,sin((n+1)*phiA   ),1);
                        SinPsi_FW[  1][n][0]->Fill(nFWspect ,sin((n+1)*phiCorA),1);
                        SinPsi_FW[  2][n][0]->Fill(nFWspect ,sin((n+1)*PsiA   ),1);

                        CosPsi_META[0][n][1]->Fill(Mtof+Mrpc,cos((n+1)*phiB   ),1);
                        CosPsi_META[1][n][1]->Fill(Mtof+Mrpc,cos((n+1)*phiCorB),1);
                        CosPsi_META[2][n][1]->Fill(Mtof+Mrpc,cos((n+1)*PsiB   ),1);
                        SinPsi_META[0][n][1]->Fill(Mtof+Mrpc,sin((n+1)*phiB   ),1);
                        SinPsi_META[1][n][1]->Fill(Mtof+Mrpc,sin((n+1)*phiCorB),1);
                        SinPsi_META[2][n][1]->Fill(Mtof+Mrpc,sin((n+1)*PsiB   ),1);

                        CosPsi_FW[  0][n][1]->Fill(nFWspect ,cos((n+1)*phiB   ),1);
                        CosPsi_FW[  1][n][1]->Fill(nFWspect ,cos((n+1)*phiCorB),1);
                        CosPsi_FW[  2][n][1]->Fill(nFWspect ,cos((n+1)*PsiB   ),1);
                        SinPsi_FW[  0][n][1]->Fill(nFWspect ,sin((n+1)*phiB   ),1);
                        SinPsi_FW[  1][n][1]->Fill(nFWspect ,sin((n+1)*phiCorB),1);
                        SinPsi_FW[  2][n][1]->Fill(nFWspect ,sin((n+1)*PsiB   ),1);
                    }
                    hMETAvsCent->Fill(im+1,Mtof+Mrpc    );
                    hFWvsCent  ->Fill(im+1,nFWspect     );
                    hMETAvsFW  ->Fill(nFWspect,Mtof+Mrpc);
                }
            }

            h0PhiEPvect->Fill(VectphiEP);
            h0PhiEP->Fill(VectphiEP);
            hPhiCor->Fill(VectphiCorr);
            h0PhiA ->Fill(phiA);
            h0PhiB ->Fill(phiB);
            h0PhiAB->Fill(phiAB);
        }
        wPhiRPA = 1.;
        //-extended-version-introduced-at-09.07.15--which-corresponds-to-META-mult-equivalent-to-numTracking-Glauber--//
        if(NA>1 && NB>1 && Mfw>3){
            hRPA->Fill(VectphiEP);
            for(Int_t im=0;im<11;im++){
                if( (Mtof+Mrpc)>=Mrang[im] && (Mtof+Mrpc)<Mrang[im+1]){ hPhiAB[im]->Fill(fabs(phiAB),wPhiRPA); }
            }
            hRPA6Ywc->Fill(VectphiCorr);
        }
        //-here-we-make-a-control-plot-for-pileup-events--//
        hMfwMmeta->Fill(Mtof+Mrpc, Mfw );
        hMfw     ->Fill(           Mfw );
        hMmeta   ->Fill(Mtof+Mrpc      );
        hMtof    ->Fill(           Mtof);
        hMrpc    ->Fill(           Mrpc);
        //-------------------------------------------iteration-WallHits----( end )--------------------------------------------------//
        
        ////////////////////////////////////////////
        //----getting-ready-to-loop-over-tracks---//
        ////////////////////////////////////////////
        size = candCat->getEntries();
        //note-new-event----------------------------
        //if(par!=0){
            for (Int_t icent=0;icent<NC-4;icent++){
                for (Int_t irpi=0;irpi<NRAPI;irpi++){
                    for (Int_t ipt=0;ipt<NPTRN-5;ipt++){
                        pfy0pt[icent][irpi][ipt].NewEvt();
                        pfR0pt[icent][irpi][ipt].NewEvt();
                        pfF0pt[icent][irpi][ipt].NewEvt();
                        phist[ icent][irpi][ipt].NewEvt();
                    }
                }
            }
        //}
        Ntrack = 0;
        /////////////////////////////////////////////////////
        //--------Now-lets-try-to-get-tracks-----------------
        /////////////////////////////////////////////////////
        for(Int_t j=0;j<size;j++){
            hPartCuts->Fill(1.);
            pParticleCand = HCategoryManager::getObject(pParticleCand,candCat,j);
            if(!pParticleCand->isFlagBit(Particle::kIsUsed)){ continue; }
            hPartCuts->Fill(2.);
            //pidFlag==kFALSE;
            
            sys      = pParticleCand->getSystem();           // sys == 0 - Rpc, sys == 1 - Tof
            mom      = pParticleCand->getMomentum();         //-later-shall-use: getMomentumPID()---(c)-Alexandr
            charge   = pParticleCand->getCharge();
            metaQa   = pParticleCand->getMetaMatchQuality(); 
            beta     = pParticleCand->getBeta();
            theta    = pParticleCand->getTheta();
            phi      = pParticleCand->getPhi();
            chi2     = pParticleCand->getChi2();
            sec      = pParticleCand->getSector();           // sector is the sixth part of the detector
            chi2In   = pParticleCand->getInnerSegmentChi2(); // chi2 of inner segment (mdc 1,2)
            chi2Out  = pParticleCand->getOuterSegmentChi2(); // chi2 of outer segment (mdc 3,4)
            mass2    = pParticleCand->getMass2();            // mass2 = mom*mom*(1-beta*beta)/(beta*beta);
            mass     = pParticleCand->getMass();             // mass  = sqrt(mass2);
            mdcdEdx  = pParticleCand->getMdcdEdx();
            tofdEdx  = pParticleCand->getTofdEdx();
            metaDx   = pParticleCand->getRkMetaDx();         
            metaDy   = pParticleCand->getRkMetaDy();         
            PCr      = pParticleCand->getR();
            PCz      = pParticleCand->getZ();
            hphi    -> Fill(phi);
            dphi     = phi - VectphiEP;
            if (dphi> 180.){dphi = dphi - 360.; }
            if (dphi<-180.){dphi = dphi + 360.; }
            dphiRec  = phi - VectphiCorr;
            if (dphiRec> 180.){dphiRec = dphiRec - 360.; }
            if (dphiRec<-180.){dphiRec = dphiRec + 360.; }
            dphiFlt  = phi - PsiCorr;
            if (dphiFlt> 180.){dphiFlt = dphiFlt - 360.; }
            if (dphiFlt<-180.){dphiFlt = dphiFlt + 360.; }
            // apply track QA cuts
            if (metaQa>2){ continue; }
            hPartCuts->Fill(3.);
            if (mass2 <0){ continue; }
            hPartCuts->Fill(4.);
            //Proton selection
            betaMeanP = fBetaMeanP->Eval(mom);
            hbetaMomAll ->Fill(mom*charge,beta);
            if(charge>0) hBetaMeanP  ->Fill(mom*charge,betaMeanP);
            if(charge>0) hBetaDownP  ->Fill(mom*charge,betaMeanP-0.05);
            if(charge>0) hBetaUpPro  ->Fill(mom*charge,betaMeanP+0.05);
            hdEdxMAll   ->Fill(charge*mass/1000 ,mdcdEdx);
            hdEdxMom    ->Fill(mom*charge,mdcdEdx);
            hfdEdxMomLow->Fill(mom*charge,fdEdxVsMomLowLimit(mom));
            hChi2In     ->Fill(chi2In);
            hChi2Out    ->Fill(chi2Out);
            hChi2       ->Fill(chi2);
            pt0       = mom*sin(theta*hpi/90.); //uncorrected pt
            y0        = (0.5*log((1+beta*cos(theta*hpi/90.))/(1-beta*cos(theta*hpi/90.)))-Ycm)/Ycm; //uncorrected Yo
            for(Int_t i=0;i<11-4;i++){
                for(Int_t j=0;j<9;j++){
                    for(Int_t k=0;k<18-5;k++){
                        if( (Mtof+Mrpc)>=Mrang[i+4] && (Mtof+Mrpc)<Mrang[i+5] && y0>Yrang[j] && y0<Yrang[j+1] && pt0>Prang[k] && pt0<Prang[k+1]){
                            phist[i][j][k].FillMass(mass2/1e6);
                        }
                    }
                }
            }
            if (p.fPID(14, mom,beta,charge)){ hbetaMomPro->Fill(mom*charge,beta); /*pidFlag == kTRUE;*/ }
            if ( mass>600 && mass<1250 && NA>0 && NB>0 && charge>0 && metaQa<5 && chi2<200. && chi2In>0.1 && chi2In<12. && chi2Out>0.1 && chi2Out<12. && sec>=0 && sec<=5 && mdcdEdx>fdEdxVsMomLowLimit(mom)){
                hphiPro ->Fill(phi);
                hMom    ->Fill(mom*charge);
                hMomCorr->Fill(charge*momCorPar->getCorrMom(14,mom,theta));
                PCp       = pParticleCand->getCorrectedMomentumPID(14);
                hMomCorrPID->Fill(charge*PCp);
                PCpz      = PCp*cos(theta*hpi/90.);
                PCE       = sqrt(mass2Pro + PCp*PCp        );
                PCY       = 0.5*(log((PCE+PCpz)/(PCE-PCpz))); //corrected rapidity for protons
                PCYn      = PCY/(2.*Ycm); //normalized to projectile rapidity Y/Yproj
                PCYo      = (PCY-Ycm)/Ycm; //normalized to projectile rapidity (Y(cm)/Yproj(cm))
                PCpt      = PCp*sin(theta*hpi/90.); //corrected pt
                hPCpt->Fill(PCpt);
                hpt0 ->Fill( pt0);
                //---get-track-efficiency-------//
                binX      =   hProtPtVsY_Eff->GetXaxis()->FindBin((Double_t) PCY );  //-Pt:Y-correction--
                binY      =   hProtPtVsY_Eff->GetYaxis()->FindBin((Double_t) PCpt);  //-Pt:Y-correction--
                TrackEff  =   hProtPtVsY_Eff->GetBinContent(binX,binY);
                //------------------------------//

                CENT=-1;
                if( (Mtof+Mrpc)>=163 && (Mtof+Mrpc)<215 ){ CENT = 10; }
                if( (Mtof+Mrpc)>=144 && (Mtof+Mrpc)<163 ){ CENT =  9; }
                if( (Mtof+Mrpc)>=123 && (Mtof+Mrpc)<144 ){ CENT =  8; }
                if( (Mtof+Mrpc)>=103 && (Mtof+Mrpc)<123 ){ CENT =  7; }
                if( (Mtof+Mrpc)>= 85 && (Mtof+Mrpc)<103 ){ CENT =  6; }
                if( (Mtof+Mrpc)>= 72 && (Mtof+Mrpc)< 85 ){ CENT =  5; }
                if( (Mtof+Mrpc)>= 59 && (Mtof+Mrpc)< 72 ){ CENT =  4; }
                if( (Mtof+Mrpc)>= 50 && (Mtof+Mrpc)< 59 ){ CENT =  3; }
                if( (Mtof+Mrpc)>= 39 && (Mtof+Mrpc)< 50 ){ CENT =  2; }
                if( (Mtof+Mrpc)>= 30 && (Mtof+Mrpc)< 39 ){ CENT =  1; }
                if( (Mtof+Mrpc)>= 20 && (Mtof+Mrpc)< 30 ){ CENT =  0; }
                RAPI=-1;
                if( PCYo>-0.9 && PCYo<-0.7 ){ RAPI=0; }
                if( PCYo>-0.7 && PCYo<-0.5 ){ RAPI=1; }
                if( PCYo>-0.5 && PCYo<-0.3 ){ RAPI=2; }
                if( PCYo>-0.3 && PCYo<-0.1 ){ RAPI=3; }
                if( PCYo>-0.1 && PCYo< 0.1 ){ RAPI=4; }
                if( PCYo> 0.1 && PCYo< 0.3 ){ RAPI=5; }
                if( PCYo> 0.3 && PCYo< 0.5 ){ RAPI=6; }
                if( PCYo> 0.5 && PCYo< 0.7 ){ RAPI=7; }
                if( PCYo> 0.7 && PCYo< 0.9 ){ RAPI=8; }
                PTRN=-1;
                if( PCpt> 200 && PCpt< 300 ){ PTRN= 0; }
                if( PCpt> 300 && PCpt< 400 ){ PTRN= 1; }
                if( PCpt> 400 && PCpt< 500 ){ PTRN= 2; }
                if( PCpt> 500 && PCpt< 600 ){ PTRN= 3; }
                if( PCpt> 600 && PCpt< 700 ){ PTRN= 4; }
                if( PCpt> 700 && PCpt< 800 ){ PTRN= 5; }
                if( PCpt> 800 && PCpt< 900 ){ PTRN= 6; }
                if( PCpt> 900 && PCpt<1000 ){ PTRN= 7; }
                if( PCpt>1000 && PCpt<1100 ){ PTRN= 8; }
                if( PCpt>1100 && PCpt<1200 ){ PTRN= 9; }
                if( PCpt>1200 && PCpt<1300 ){ PTRN=10; }
                if( PCpt>1300 && PCpt<1400 ){ PTRN=11; }
                if( PCpt>1400 && PCpt<1500 ){ PTRN=12; }
                if( PCpt>1500 && PCpt<1600 ){ PTRN=13; }
                if( PCpt>1600 && PCpt<1700 ){ PTRN=14; }
                if( PCpt>1700 && PCpt<1800 ){ PTRN=15; }
                if( PCpt>1800 && PCpt<1900 ){ PTRN=16; }
                if( PCpt>1900 && PCpt<2000 ){ PTRN=17; }
                if( CENT>=4 && RAPI>=0 && PTRN>=0 && PTRN<=12 ){
                    pfy0pt[ CENT-4][RAPI][PTRN].Fill(dphi   , fabs(phiAB), VectphiEP  , 1.*TrackEff);
                    pfR0pt[ CENT-4][RAPI][PTRN].Fill(dphiRec, fabs(phiAB), VectphiCorr, 1.*TrackEff);
                    pfF0pt[ CENT-4][RAPI][PTRN].Fill(dphiFlt, fabs(phiAB), PsiCorr    , 1.*TrackEff); 
                    phist[  CENT-4][RAPI][PTRN].Fill(mass2/1e6,phi); 
                    hptvsphi->Fill(PCpt,phi);
                    hptvstheta->Fill(PCpt,theta);
                    hphivstheta->Fill(theta,phi);
                }
            }//end-proton-selection
            //------Now-make-PID----------end----------
            Ntrack++;
        }//j-----Loop-over-entries-----
        hNtrMult->Fill(Ntrack,Mtof+Mrpc);
    }//i-----Loop-over-events----------

    sorter.finalize();

    foutfile->Write();
    foutfile->Close();
  
    return 0;  
}