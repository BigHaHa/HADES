// . /cvmfs/hades.gsi.de/install/5.34.01/hydra2-4.9f/defall.sh
// root -l -b -q EventCount.C++(\"/hera/hades/dst/apr12/gen8/108/root/be1210812550001.hld_dst_apr12_1.root\","histo.root")
//
//or
//
//root -l
//.x EventCount.C++(\"/hera/hades/dst/apr12/gen8/108/root/be1210812550001.hld_dst_apr12_1.root\","histo.root")

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
#include <iostream>
#include <fstream>
#include <map>
#include <stdio.h>
#include <iomanip>
#include <vector>
#include <TH3F.h>
#include "/u/parfenov/anaflow/TRandom3.h" //need because we should not use TRandom::Rndm() 
#include "/u/parfenov/anaflow/EventCountClasses.h" //include all classes and functions we need

using namespace std;

Int_t EventCount(TString inputlist, TString outfile, Int_t nev=-1)
{  

  HLoop* loop = new HLoop   (kTRUE    );  // kTRUE : create Hades  (needed to work with standard eventstructure)
         loop ->addMultFiles(inputlist);

//-----Inicialize-classes-we-need-to-use-----------------
  PidParticle p;
  dEdx        d;
  HMultCorr   cCorr;
//-------------------------------------------------------

//--------------------------------------------------------------------------------------------------------
#include "/u/parfenov/anaflow/EventCountVar.h"     //include all constants and variables we need          
#include "/u/parfenov/anaflow/outsigma.C"          //Declaration of the 3-sigma cuts for our pt & rapidity
TH2F*  hProtPtVsY_Eff;
TFile* foutfile       = new TFile(outfile.Data(), "recreate"); 
TFile* FileProtEff    = new TFile("EffCorrForProtPtVsY_WithWideBetaPcuts_ShieldGosia.root", "read");
       hProtPtVsY_Eff = (TH2F*) FileProtEff->Get("hProtPtVsY_EffSm");
       foutfile->cd();

#include "/u/parfenov/anaflow/EventCountHisto.h"   //include all histograms we need                       

FFlow pfy0pt[NC][NRAPI][NPTRN]; //--flow--9-bins-Yo--12-bins-in-Pt-in-3-centrality-sets--
  for (Int_t i=0; i<NC;i++){
   for(Int_t j=0; j<NRAPI; j++){
     for(Int_t k=0; k<NPTRN; k++){
       //-proton------------------------------------------
       pfy0pt[i][j][k].THDeclare("pCent",i,"Yo",j,"Pt",k);
     }
   }
  }
//--------------------------------------------------------------------------------------------------------

  if(!loop->setInput("-*,+HParticleCand,+HParticleEvtInfo,+HStart2Hit,+HStart2Cal,+HWallHit,+HRpcCluster,+HTofHit,+HWallEventPlane")) { exit(1); }
  loop->printCategories();
  
  TIterator* iterWallHit = 0;
  if (loop->getCategory("HWallHit")  ) iterWallHit   = loop->getCategory("HWallHit")->MakeIterator();
  
  HParticleEvtInfo* evtinfo;
  
  HParticleCand* pParticleCand; HCategory* candCat = (HCategory*)HCategoryManager::getCategory(catParticleCand); if(!candCat){exit(1);}
  HWallHit*      pWallHit;      HCategory* wallCat = (HCategory*)HCategoryManager::getCategory(catWallHit);      if(!wallCat){exit(1);}

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
    if (i%1000 == 0){
      cout << "    event " << i << endl;
    }
    
    if( loop->isNewFile(currFName) ){ 
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

        #include "/u/parfenov/anaflow/anaKaons_params_gen7_Auto_RC_RW_updated_gen108RminEqZero_11MultRecent.cc" //-for recentering---

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
    
      //-------Get-vertex-from-combined-fit-of-fitted-inner-segments------------
      HVertex            fPrimVertReco = gHades->getCurrentEvent()->getHeader()->getVertexReco();
      vertexX = fPrimVertReco.getX();
      vertexY = fPrimVertReco.getY();
      vertexZ = fPrimVertReco.getZ();
      
      hvertexXZ->Fill(vertexZ,vertexX);
      hvertexXY->Fill(vertexY,vertexX);
      hvertexX ->Fill(vertexX);
      hvertexY ->Fill(vertexY);
      hvertexZ ->Fill(vertexZ);
      
      //------Alexandr's-vertex-cuts---------
      
      tgtChi2 = (Float_t) fPrimVertReco.getChi2();

      hTargChi2->Fill(tgtChi2);


      if( tgtChi2<2.0 || tgtChi2>40.) continue; //-skip-events-with-badly-reconstructed-target-------
      if( vertexZ<-59.|| vertexZ>0. ) continue; //-skip-events-with-primary-vertex-outside-of-target-

      hTargX_EVT->Fill(vertexX);
      hTargY_EVT->Fill(vertexY);
      hTargZ_EVT->Fill(vertexZ);

      if( (vertexX-tgtXc)*(vertexX-tgtXc)+(vertexY-tgtYc)*(vertexY-tgtYc)>R2 ) continue; //-skip-events-with-primary-vertex-far-from-target---
      hCutSteps->Fill(5.5);

      hTargXYsel->Fill(vertexX,vertexY);
      hTargXZsel->Fill(vertexZ,vertexX);
      //---------end-of-Alexandr's-vertex-cut--------------
      
      //------HWall-iteration------------------------------
      if (iterWallHit){
    iterWallHit->Reset();
    HWallHit* pWallHit;
        
    while ( (pWallHit = (HWallHit*)iterWallHit->Next()) != 0 ){
         
         cellNum    = pWallHit->getCell();
             cellTime   = pWallHit->getTime();
             cellCharge = pWallHit->getCharge();
             cellPhi    = pWallHit->getPhi();
         
         hWallHitTIME->Fill(cellNum-0.1, cellTime  );
         hWallHitCHRG->Fill(cellNum-0.1, cellCharge);
         
         if (cellCharge > d.dEdxCut(cellNum)/10. && cellTime > 16. && cellTime < 44.){
           hWallHitTIMEcut->Fill(cellNum-0.1, cellTime  );
           hWallHitCHRGcut->Fill(cellNum-0.1, cellCharge);
         }//end-cut
         
    }//------end-HWall-iteration-----------------------
      }//endif

    if( evtinfo->getSumParticleCandMult()<2 ) continue;  //getSumParticleCandMult()>1
      if(!evtinfo->isGoodEvent(
                               Particle::kGoodVertexClust|
                               Particle::kGoodVertexCand|
                               Particle::kNoPileUpSTART|
                               Particle::kGoodSTART|      
                               Particle::kGoodTRIGGER|    
                               Particle::kNoVETO|         
                               Particle::kGoodSTARTVETO|  
                               Particle::kGoodSTARTMETA
                              )) continue ; // bitwise add flags as suggested by Manuel 
                              
      Mtof = evtinfo->getSumTofMult();
      Mrpc = evtinfo->getSumRpcMult();
      Mult = Mtof+Mrpc;
      
      hmultTofRpc->Fill(Mult);                                   // Fill( name , weight )
      
      
      //---------centrality-classes--------------------------------------------------------------------(begin)--------//
      //-Centrality-selection-for-further-if-selection-statements------------------------//
      Int_t   nTrkMult = getParticleTrkMult();
      Mtr   = nTrkMult;
      Float_t nTrkMultCorr = nTrkMult*fTrkMultScaler;
      MtrC  = nTrkMultCorr; //-after-correction-to-average-day108--

      if( nTrkMultCorr>33 && nTrkMultCorr<=49 ){ FOPI_blue=kTRUE; }else{ FOPI_blue=kFALSE; }        // Berusz's estimate
      if( nTrkMultCorr>49 && nTrkMultCorr<=82 ){ FOPI_gren=kTRUE; }else{ FOPI_gren=kFALSE; }        // based on Glauber model
      if( nTrkMultCorr>82                     ){ FOPI_redd=kTRUE; }else{ FOPI_redd=kFALSE; }        // (lucky me, seems very similar to my naive guess!)--Alexandr-(c)--

      if( FOPI_blue ){ hPhiAB_blue->Fill(fabs(phiAB), wPhiRPA); }
      if( FOPI_gren ){ hPhiAB_gren->Fill(fabs(phiAB), wPhiRPA); }
      if( FOPI_redd ){ hPhiAB_redd->Fill(fabs(phiAB), wPhiRPA); }

      multCand = evtinfo->getSumSelectedParticleCandMult();

      hmultSumSPCand->Fill(multCand);
     
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
      vect.Set(0.,0.);
      vsum.Set(0.,0.);
      eX  .Set(1.,0.);
      dEdxCut=0.0;
      xyRadius=0.0;
      
      QLx=0.0, QLy=0.0, NL=0; 
      QMx=0.0, QMy=0.0, NM=0;
      QNx=0.0, QNy=0.0, NN=0;

      for(Int_t n=0; n<6; n++){
         Qx[ n]    =0.0; Qy[ n]   =0.0;
         QxRec[n]  =0.0; QyRec[n] =0.0;
         Qax[n]    =0.0; Qay[n]   =0.0;
         Qbx[n]    =0.0; Qby[n]   =0.0;
         phiEP[ n] = -1000.;
         phiEPr[n] = 0.;
      }

      phiA   = -1000;
      phiB   = -1000;
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

      //-------------------------------------------iteration-WallHits----(begin)--------------------------------------------------//
      for(Int_t j=0; j<(wallCat->getEntries()); j++){
        pWallHit = HCategoryManager::getObject(pWallHit,wallCat,j);

         cellNum    = pWallHit->getCell();
         cellTime   = pWallHit->getTime();
         cellCharge = pWallHit->getCharge();

         cellChargeCtime = cellTime;
         hFW_TimeCell_0->Fill(cellTime  , cellNum);
         hFW_dEdxCell_0->Fill(cellCharge, cellNum);

         if(cellNum< 144                ){ dEdxCut=Z1_cut_s; }
         if(cellNum>=144 && cellNum<=207){ dEdxCut=Z1_cut_m; }
         if(cellNum> 207                ){ dEdxCut=Z1_cut_l; }

         if(cellTime>T1_cut && cellTime<T2_cut && cellCharge>dEdxCut){
           
           nFWspect++;

           hFW_TimeCell_1->Fill(cellTime  , cellNum);
           hFW_dEdxCell_1->Fill(cellCharge, cellNum);

           pWallHit->getXYZLab(wallX,wallY,wallZ);
           hFWxyCC->Fill(wallX,wallY); //cell centers
           
           for (Int_t im=1;im<12;im++){
            if( (Mtof+Mrpc)>=Mrang[im-1] && (Mtof+Mrpc)<Mrang[im]){ wallXc = wallX - X_shM[12-im]; wallYc = wallY - Y_shM[12-im]; }
           }

           //-here-we-fill-d2N/(dxdy)-inX-band--and--y-band--for-auto-recentering-
           //-this-makes-realistic-distribution-on-FW-face------------------------
           if (cellNum<=143)                 { cellSize = 40;  }
           if (cellNum>=144 && cellNum<=207) { cellSize = 80;  }
           if (cellNum>=210 && cellNum<=301) { cellSize = 160; }
           XfwSmear = wallXc + ( Random.Rndm(1) - 0.5 )*cellSize;
           YfwSmear = wallYc + ( Random.Rndm(1) - 0.5 )*cellSize;

           
           for (Int_t im=0;im<8;im++){
            //-X-beam-center-scan---
            if( (Mtof+Mrpc)>=  Mcent[im] && (Mtof+Mrpc)< Mcent[im+1] ){ hFWxyCCsmearXscanM[im]->Fill(XfwSmear,YfwSmear); }
            //-Y-beam-center-scan---
            if( (Mtof+Mrpc)>=  Mcent[im] && (Mtof+Mrpc)< Mcent[im+1] ){ hFWxyCCsmearYscanM[im]->Fill(XfwSmear,YfwSmear); }
           }


           //-this-cut-was-forgotten-for-Feb2014-HADES-CM-report--//
           //-I-reintroduce-it-14.03.2014----(c)-Sadovsky---------//
           if(  wallXc*wallXc + wallYc*wallYc >= R0_cut*R0_cut  /*50.*50.*/ ){
             hFW_TimeCell_2->Fill(cellTime  , cellNum);
             hFW_dEdxCell_2->Fill(cellCharge, cellNum);


             //-spectators-selected--
             multWall++;

             vect.Set(wallXc, wallYc);
             vect =  vect.Unit();
             vsum += vect;

             //-Qnx-method--for-n=1,2,3,4,5------
             VectPhi_i = TMath::ATan2(YfwSmear,XfwSmear); //-Phi-of-i_th-spectator-particle---
             for(Int_t n=1; n<6; n++){//-loop-over-harmonics-(n)--
               Qx[n] = Qx[n] + wgh*cos(double(n)*VectPhi_i);
               Qy[n] = Qy[n] + wgh*sin(double(n)*VectPhi_i);
             }

             cellsVect.SetCellVect(vect, cellCharge); //-note-[vect]-is-a-unit-vector-this-is-for-further-subevent-(A/B)-split--


             hFWxyCCsmear->Fill(XfwSmear+X_shift,YfwSmear+Y_shift); //cells smeared but not shifted
             hFWxySHsmear->Fill(XfwSmear        ,YfwSmear        ); //cells smeared and shifted as a whole

             //-center-of-gravity-study-as-a-course-of-centrality-----------------------------
             for (Int_t im=0;im<8;im++){
              if( (Mtof+Mrpc)>=Mcent[im] && (Mtof+Mrpc)<Mcent[im+1]){ hFWxyCCsmearM[im]->Fill(XfwSmear,YfwSmear); }
             }
             
           }//-endif-R0_cut--//
           
        }//end-if-cells-cut---
     }//end-for-HitWall-loop---

     VectphiEP =  vsum.DeltaPhi(eX)   *rad2deg;
     VectphiEPr=  vsum.DeltaPhi(eX); //radians
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
           VectPhi_i = vectA.DeltaPhi(eX); //-Phi-of-i_th-spectator-particle-from-subevent-A--
           
           for(Int_t n=1; n<6; n++){//-loop-over-harmonics-(n)--
              Qax[n] = Qax[n] + wgh*cos(double(n)*VectPhi_i);
              Qay[n] = Qay[n] + wgh*sin(double(n)*VectPhi_i);
           }

           FWdEdxA=FWdEdxA + cellsVect.GetCellCharge(ic); //-total-Eloss-from-all-spectators---
           
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
           
           for(Int_t n=1; n<6; n++){//-loop-over-harmonics-(n)--
              Qbx[n] = Qbx[n] + wgh*cos(double(n)*VectPhi_i);
              Qby[n] = Qby[n] + wgh*sin(double(n)*VectPhi_i);
           }

           FWdEdxB=FWdEdxB + cellsVect.GetCellCharge(ic); //-total-Eloss-from-all-spectators---
        }//-HWallHit-second-loop-for-reaction-plane-resolution--( end )---//
        //-this-is-(A/B)-subevents-split--( end )------------------------------------------
      }//-endfor-multFWcells-loop--for-A/B-method--

      //-calculating-eventplane-angles---------------
      phiA   = vsumA.DeltaPhi(eX)       *rad2deg;
      phiB   = vsumB.DeltaPhi(eX)       *rad2deg;
      phiAB  = vsumA.DeltaPhi(vsumB)    *rad2deg;
      Mfw = NA+NB;

      if (Mfw > 1){
        for (Int_t n=1;n<6;n++){
            QxRec[n]     = Qx[n] - MeanQx[n];
            QyRec[n]     = Qy[n] - MeanQy[n];
            phiEP[n ]    =  TMath::ATan2(Qy[n],Qx[n])/n       *rad2deg;
            phiEPr[n]    =  TMath::ATan2(Qy[n],Qx[n])/n               ;//radians
            phiEPrec[n]  =  TMath::ATan2(QyRec[n],QxRec[n])/n *rad2deg;
            hPhiEP[n-1]  -> Fill(phiEP[n]);
            hPhiEPrec[n-1]->Fill(phiEPrec[n]);
            hQvX  [n-1]  -> Fill(Qx[n]);
            hQvY  [n-1]  -> Fill(Qy[n]);
            hQvXrec[n-1] -> Fill(QxRec[n]);
            hQvYrec[n-1] -> Fill(QyRec[n]);
        }
        h0PhiEPvect->Fill(VectphiEP);
      }

      if(Mfw>=0){
        h0PhiEP->Fill(VectphiEP);
        h0PhiA ->Fill(phiA);
        h0PhiB ->Fill(phiB);
        h0PhiAB->Fill(phiAB);
      }
      //-extended-version-introduced-at-09.07.15--which-corresponds-to-META-mult-equivalent-to-numTracking-Glauber--//
      if(NA>0 && NB>0){
        hRPA->Fill(VectphiEP);
        for(Int_t im=0;im<11;im++){
            if( (Mtof+Mrpc)>=Mrang[im] && (Mtof+Mrpc)<Mrang[im+1] ){ hRPAy[im] ->Fill(VectphiEP); wPhiRPA = 1./fRPAy[im] ->Eval(VectphiEP); hRPAyc[im] ->Fill(VectphiEP,wPhiRPA); }
        }
        hRPA6Ywc->Fill(VectphiEP,wPhiRPA);
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
      Bool_t pimFlag[size];    // pi-
      Bool_t pipFlag[size];    // pi+
      Bool_t proFlag[size];    // proton
      Bool_t kapFlag[size];    // K+ meson        
      
      //note-new-event----------------------------
      for (Int_t icent=0;icent<NC;icent++){
        for (Int_t irpi=0;irpi<NRAPI;irpi++){
            for (Int_t ipt=0;ipt<NPTRN;ipt++){
                pfy0pt[icent][irpi][ipt].NewEvt();
            }
        }
      }   
      
      /////////////////////////////////////////////////////
      //--------Now-lets-try-to-get-tracks-----------------
      /////////////////////////////////////////////////////
      for(Int_t j=0;j<size;j++){
    
    pimFlag[j] = kFALSE;
    pipFlag[j] = kFALSE;
    proFlag[j] = kFALSE;
    kapFlag[j] = kFALSE;
        
    hEntries->Fill(size);
    
    pParticleCand = HCategoryManager::getObject(pParticleCand,candCat,j);
    if(!pParticleCand->isFlagBit(Particle::kIsUsed)) continue;
    
    
    sys      = pParticleCand->getSystem();           // sys == 0 - Rpc, sys == 1 - Tof
    mom      = pParticleCand->getMomentum();         //-later-shall-use: getMomentumPID()---(c)-Alexandr
    charge   = pParticleCand->getCharge();
    metaQa   = pParticleCand->getMetaMatchQuality(); 
    beta     = pParticleCand->getBeta();
    theta    = pParticleCand->getTheta();
    phi      = pParticleCand->getPhi();
    chi2     = pParticleCand->getChi2();
    sec      = pParticleCand->getSector();           // sector is the sixth part of the detector
    chi2In   = pParticleCand->getInnerSegmentChi2(); // chi2 of inner segment
    chi2Out  = pParticleCand->getOuterSegmentChi2(); // chi2 of outer segment
    mass2    = pParticleCand->getMass2();            // mass2 = mom*mom*(1-beta*beta)/(beta*beta);
    mass     = pParticleCand->getMass();             // mass  = sqrt(mass2);
    mdcdEdx  = pParticleCand->getMdcdEdx();
    tofdEdx  = pParticleCand->getTofdEdx();
    metaDx   = pParticleCand->getRkMetaDx();         
    metaDy   = pParticleCand->getRkMetaDy();         
    
    PCr              = pParticleCand->getR();
        PCz              = pParticleCand->getZ();
    
    pz       = mom*cos(theta*hpi/90.)/1000;
    pt0      = mom*sin(theta*hpi/90.)/1000;
    px       = pt0*cos(phi*hpi/90.);
    py       = pt0*sin(phi*hpi/90.);
    
    rapid    = 0.5*TMath::Log((1+beta)/(1-beta));
    eta      =    -TMath::Log(tan(theta/2));
    
    betaL    = beta*cos(theta*hpi/90.);
    rapidity = 0.5*TMath::Log((1+betaL)/(1-betaL));
    Y0       = (rapidity - Ycm)/Ycm;
    dphi     = phi - VectphiEP;
    if (dphi>180.){dphi = dphi - 360.; }
    if (dphi<-180.){dphi = dphi + 360.; }

    hPhi->Fill(phi);
    
        // apply track QA cuts
        if (metaQa>2) continue;
    if (mass2 <0) continue;
    
    hRapidity->Fill(rapidity);
    
    hM2fromM->Fill(mass);
    
    hCharge->Fill(charge);
    
     //-Tetiana's-plots-----------
           hChi2       ->Fill(chi2);
           hChi2In     ->Fill(chi2In);
           hChi2Out    ->Fill(chi2Out);
           //Fill Theta vs Phi Plots
           hthetaPhiAll->Fill(phi,theta);
           //Fill Beta vs Momentum
           hbetaMomAll ->Fill(mom*charge,beta);
           hmommassAll ->Fill(charge*mass,mom);
       
       hsys    ->Fill(sys);
       hmom    ->Fill(mom);
       hmetaQa ->Fill(metaQa);
       hsec    ->Fill(sec);
       hmdcdEdx->Fill(mdcdEdx);
       htofdEdx->Fill(tofdEdx);
       hmetaDx ->Fill(metaDx);
       hmetaDy ->Fill(metaDy);
       hPCr    ->Fill(PCr);
       hPCz    ->Fill(PCz);
       
       hpx     ->Fill(px);
       hpy     ->Fill(py);
       hpz     ->Fill(pz);
       hpt     ->Fill(pt);
       
       hdEdxMAll->Fill(mass/1000,mdcdEdx);
       
       
       //------Now-make-PID-------begin-----------
       
       pimFlag[j] = p.fPID(9, mom,beta,charge);
       pipFlag[j] = p.fPID(8, mom,beta,charge);
       kapFlag[j] = p.fPID(11,mom,beta,charge);
       proFlag[j] = p.fPID(14,mom,beta,charge);
       
       //fill-particle-numeration-integral-----------------------------------------------------
                                      hmassAll->Fill(charge*mass2/1000000);
       if (pipFlag[j] == kTRUE && mdcdEdx < 5) {
         hPIDPip ->Fill(charge*mass2/1000000); 
         pParticleCand->calc4vectorProperties(HPhysicsConstants::mass( 8));
         vector = (*pParticleCand);
         pt     = vector.Pt();
         y      = vector.Rapidity();
         hyPIDPip->Fill(y);
         hRPIDPip->Fill(rapidity);
         npip++;
      } // pi+    selection (Pip);
       if (pimFlag[j] == kTRUE && mdcdEdx < 5) {
         hPIDPim ->Fill(charge*mass2/1000000);
         pParticleCand->calc4vectorProperties(HPhysicsConstants::mass( 9));
         vector = (*pParticleCand);
         pt     = vector.Pt();
         y      = vector.Rapidity();
         hyPIDPim->Fill(y);
         hRPIDPim->Fill(rapidity);
         npim++;
      } // pi-    selection (Pim);
       if (kapFlag[j] == kTRUE && mdcdEdx < 5) { 
         hPIDKap ->Fill(charge*mass2/1000000);
         pParticleCand->calc4vectorProperties(HPhysicsConstants::mass(11));
         vector = (*pParticleCand);
         pt     = vector.Pt();
         y      = vector.Rapidity();
         hyPIDKap->Fill(y);
         hRPIDKap->Fill(rapidity);
         nkap++;
      } // K+     selection (Kap);
       if (proFlag[j] == kTRUE && mdcdEdx < 5) { 
         hPIDPro ->Fill(charge*mass2/1000000);
         pParticleCand->calc4vectorProperties(HPhysicsConstants::mass(14));
         vector = (*pParticleCand);
         pt     = vector.Pt();
         y      = vector.Rapidity();
         hyPIDPro->Fill(y);
         hRPIDPro->Fill(rapidity);
         npro++;
      } // proton selection (Pro);
       //--------------------------------------------------------------------------------------
           
       for (Int_t ipt=0; ipt<npt-1;ipt++){
        
         if (pt0>binpt[ipt] && pt0<binpt[ipt+1]){
                                    hM2Y   [ipt]->Fill(charge*mass2/1000000,rapidity);
         if (pipFlag[j] == kTRUE) { hM2YPip[ipt]->Fill(charge*mass2/1000000,rapidity); }
         if (pimFlag[j] == kTRUE) { hM2YPim[ipt]->Fill(charge*mass2/1000000,rapidity); }
         if (kapFlag[j] == kTRUE) { hM2YKap[ipt]->Fill(charge*mass2/1000000,rapidity); }
         if (proFlag[j] == kTRUE) { hM2YPro[ipt]->Fill(charge*mass2/1000000,rapidity); }
         }
         
         for (Int_t irpdt=0; irpdt<nrpdt-1;irpdt++){
           if (pt0>binpt[ipt] && pt0<binpt[ipt+1] && rapidity>binrpdt[irpdt] && rapidity<binrpdt[irpdt+1]){
         if (                                                                                                       mdcdEdx < 5) { hM2   [ipt][irpdt]->Fill( charge*mass2/1000000);}
             
         if (pipFlag[j] == kTRUE && mass2/1e6<=mass2Pip+sPip[ipt][irpdt] && mass2/1e6>=mass2Pip-sPip[ipt][irpdt] && mdcdEdx < 5) { hM2Pip[ipt][irpdt]->Fill( charge*mass2/1000000);}
         if (pimFlag[j] == kTRUE && mass2/1e6<=mass2Pim+sPim[ipt][irpdt] && mass2/1e6>=mass2Pim-sPim[ipt][irpdt] && mdcdEdx < 5) { hM2Pim[ipt][irpdt]->Fill( charge*mass2/1000000);}
         if (kapFlag[j] == kTRUE && mass2/1e6<=mass2Kap+sKap[ipt][irpdt] && mass2/1e6>=mass2Kap-sKap[ipt][irpdt] && mdcdEdx < 5) { hM2Kap[ipt][irpdt]->Fill( charge*mass2/1000000);}
         if (proFlag[j] == kTRUE && mass2/1e6<=mass2Pro+sPro[ipt][irpdt] && mass2/1e6>=mass2Pro-sPro[ipt][irpdt] && mdcdEdx < 5) { hM2P  [ipt][irpdt]->Fill( charge*mass2/1000000);} 
         
         hdEdxMass[ipt][irpdt]->Fill(charge*mass/1000,mdcdEdx);
         
           }//endif
         }//irpdt
       }//ipt

        for (Int_t ipt=0; ipt<npt-1;ipt++){
        
         if (pt0>binpt[ipt] && pt0<binpt[ipt+1] && (rapidity-Ycm)>=-Ycm && (rapidity-Ycm)<=Ycm){
          v2            = cos(2*((phi-180-VectphiEP)*hpi/90.));
          v2sum [ipt]  += v2;
          ncount[ipt]++;
          hv2phi[ipt]->Fill((phi-180-VectphiEP),v2,wgh);

          hNphi [ipt]->Fill((phi-180-VectphiEP),wgh);
         }
        }
        
        if (mass>600 && mass<1250 && NA>0 && NB>0 && charge>0 && metaQa<5 && chi2<200. && chi2In>0.1 && chi2In<12. && chi2Out>0.1 && chi2Out<12. && sec>=0 && sec<=5 && mdcdEdx>fdEdxVsMomLowLimit(mom)){

            PCp     = pParticleCand->getCorrectedMomentumPID(14);
            PCpz      = PCp*cos(theta*hpi/90.);

            PCE     = sqrt(938.272*938.272 + PCp*PCp );
            PCY     = 0.5*(log((PCE+PCpz)/(PCE-PCpz))); //corrected rapidity for protons
            PCYn    = PCY/(2.*Ycm); //normalized to projectile rapidity Y/Yproj
            PCYo    = (PCY-Ycm)/Ycm; //normalized to projectile rapidity (Y(cm)/Yproj(cm))
            PCpt    = PCp*sin(theta*hpi/90.); //corrected pt

            //---get-track-efficiency-------//
            binX     =   hProtPtVsY_Eff->GetXaxis()->FindBin((Double_t) rapidity );  //-Pt:Y-correction--
            binY     =   hProtPtVsY_Eff->GetYaxis()->FindBin((Double_t) PCpt );  //-Pt:Y-correction--
            TrackEff =   hProtPtVsY_Eff->GetBinContent(binX,binY);
            //------------------------------//
            for (Int_t n=1; n<6; n++){
                hPhiEPc[n-1]->Fill(phiEP[n],wgh*TrackEff);
            }

            CENT=-1;
            if( FOPI_redd ){ CENT=0; }
            if( FOPI_gren ){ CENT=1; }
            if( FOPI_blue ){ CENT=2; }
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
            if( CENT>=0 && RAPI>=0 && PTRN>=0 ){
                   pfy0pt[CENT][RAPI][PTRN].Fill(dphi, fabs(phiAB), VectphiEP, wgh*TrackEff);
                 }
        }
        
       
       //------Now-make-PID----------end----------
       
      }//j-----Loop-over-entries-----

  }//i-----Loop-over-events----------

  for (Int_t i=0;i<npt-1;i++){
        v2x[i] = binpt[i] + 0.5*(binpt[i+1]-binpt[i]);
        v2w[i] = hv2phi[i]->GetMean(2);
        hv2   -> Fill(v2x[i],v2w[i]);
  }

  Float_t yPID [4] = {hyPIDPim->GetMean(), hyPIDPip->GetMean(), hyPIDKap->GetMean(), hyPIDPro->GetMean()};
  Float_t ybet [4] = {hRPIDPim->GetMean(), hRPIDPip->GetMean(), hRPIDKap->GetMean(), hRPIDPro->GetMean()};
  Float_t dyPID[4] = {hyPIDPim->GetRMS() , hyPIDPip->GetRMS() , hyPIDKap->GetRMS() , hyPIDPro->GetRMS() };
  Float_t dybet[4] = {hRPIDPim->GetRMS() , hRPIDPip->GetRMS() , hRPIDKap->GetRMS() , hRPIDPro->GetRMS() };
  Float_t dy   [4];
  Float_t MeanY = 0.0;
  for (Int_t i=0;i<4;i++){
    dy[i] = 0.5*sqrt(ybet[i]*ybet[i]*dyPID[i]*dyPID[i]+yPID[i]*yPID[i]*dybet[i]*dybet[i])/(ybet[i]*ybet[i]);
    MeanY+=yPID[i]/ybet[i];
  }
  MeanY = MeanY/4;
  
  for (Int_t i=0;i<4;i++){
    hyvsPID ->Fill(i+0.5,yPID[i]/ybet[i]);
    hnormPID->Fill(i+0.5,MeanY);
  }
  
  foutfile->Write();
  foutfile->Close();
  
  return 0;
  
}