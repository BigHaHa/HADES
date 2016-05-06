#include "hades.h"
#include "htool.h"
#include "hphysicsconstants.h"
#include "hrootsource.h"
#include "hiterator.h"
#include "hloop.h"
#include "hdst.h"

#include "haddef.h"
#include "heventheader.h"
//#include "hgeantdef.h"

#include "hparticledef.h"
#include "hparticlestructs.h"
#include "hparticlecand.h"
#include "hparticlecandsim.h"
#include "hcategorymanager.h"
#include "hparticletracksorter.h"
#include "hparticlet0reco.h"  // new
#include "hparticlevertexfind.h"
#include "hparticleevtinfo.h"
#include "htaskset.h"
//#include "hgeantkine.h"


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

#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TCutG.h"
#include "TF1.h"
#include "TLegend.h"
#include "TLegendEntry.h"
//#include "GeomFunct.h"

//_____________________________appended via mhimmelreich
#include "TProfile.h"
#include "hntuple.h"
#include "hparticletool.h"
//#include "TArray.h"
#include "TVector2.h"
#include "TVector3.h"
#include "hwallhit.h"
#include "hcategorymanager.h"

#include <iostream>
#include <map>
#include <stdio.h>

using namespace std;

const static Float_t mParticle = 938.0F;                 //protonmass
//HPhysicsConstants::mass(14)



Bool_t isTriggerPT3(){ return gHades->getCurrentEvent()->getHeader()->isTBit(13); }
Bool_t isTriggerPT2(){ return gHades->getCurrentEvent()->getHeader()->isTBit(12); }


Bool_t selectStartCorr() {
    HCategory *fCatStartHit = HCategoryManager::getCategory(catStart2Hit,1,"catStart2Hit");
    if(!((HStart2Hit*)fCatStartHit->getObject(0)))  return kFALSE;    //object not available (3%)
    if(((HStart2Hit*)fCatStartHit->getObject(0))->getCorrFlag() == -1) return kFALSE;  //No start time found (5%)
    //    if(((HStart2Hit*)fCatStartHit->getObject(0))->getTime() >  2) return kFALSE;
    //    if(((HStart2Hit*)fCatStartHit->getObject(0))->getTime() < -2) return kFALSE;
    return kTRUE; // one ore two hits in Start (92%)
}


Bool_t selectStart()
{
    HCategory *fCatStartCal = HCategoryManager::getCategory(catStart2Cal,1,"catStart2Cal");
    HCategory *fCatStartHit = HCategoryManager::getCategory(catStart2Hit,1,"catStart2Hit");

    Int_t nStartHits = fCatStartHit->getEntries();
    if(nStartHits != 1) {
	return kFALSE;
    }
    //    HStart2Hit *start = (HStart2Hit*)fCatStartHit->getObject(0);
    //    Int_t mult_start  = start->getMultiplicity();
    //    if(mult_start > 1) return kFALSE;

    Int_t nStartCals = fCatStartCal->getEntries();
    Int_t mult_start = 0;
    for(Int_t n = 0; n < nStartCals; n++) {
	HStart2Cal *start_cal = (HStart2Cal*)fCatStartCal->getObject(n);
	if(start_cal->getModule() == 0
	   && TMath::Abs(start_cal->getTime(1)) < 10.
	  ) {
	    mult_start++;
	}
    }

    if(mult_start > 1) return kFALSE;
    return kTRUE;
}


// calculates velocity as a function of mass and momentum
Double_t fBeta(Double_t* x_val, Double_t* par)
{
    Double_t Momentum = TMath::Abs(x_val[0]);
    Double_t Mass     = par[0];
    Double_t Beta     = 0.0;
    Double_t poverm   = 0.0;

    if(Mass > 0.0)
    {
	poverm = Momentum/Mass;
	Beta   = poverm*1/(sqrt(poverm*poverm+1.0));
    }
    return Beta;
}


static Bool_t selectHadronsQa(HParticleCand* pcand)
{
    Bool_t test = kFALSE;

    if(pcand->isFakeRejected()) return kFALSE;

    if( pcand->isFlagAND(4,
			 Particle::kIsAcceptedHitInnerMDC,
			 Particle::kIsAcceptedHitOuterMDC,
			 Particle::kIsAcceptedHitMETA,
			 Particle::kIsAcceptedRK
			)
       &&
       pcand->getChi2() < 400
       &&
       pcand->getMetaMatchQuality() < 2
      ) test = kTRUE;

    return test;
}


static Bool_t rejectLeptons(HParticleCand* pcand){ return kFALSE; }

Int_t getCentralityClassFixedCuts10()
{
    // Return centrality class 10% of total cross section
    HCategory* fParticleEvtInfoCat =  (HCategory*)HCategoryManager::getCategory(catParticleEvtInfo,kTRUE,"catParticleEvtInfo");

    if(!fParticleEvtInfoCat)
    { cout<<"No catParticleEvtInfo in input!"<<endl;
	return kTRUE;
    }

    HParticleEvtInfo *event_info =(HParticleEvtInfo*)fParticleEvtInfoCat->getObject(0);
    if(!event_info) {std::cout<<"No Event Info"<<std::endl;
	return kTRUE;}

    Int_t buffer = event_info->getSumTofMultCut() + event_info->getSumRpcMultHitCut();
    Int_t CentrClasses[10] = {250,160,121,88,60,35,0,0,0,0};
    // 10% gen8 TOFRPC(timecut)

    Int_t i = 0;
    while ((buffer < CentrClasses[i]) && ( CentrClasses[i] > 0) && (i< 10)) { i++; }
    return i;
}


Double_t DiffAngle(Double_t PhiRPlane, Double_t phi)
{
    return ( PhiRPlane == phi) ? 0 : ((PhiRPlane > phi) ? (PhiRPlane - phi) :(phi-PhiRPlane));
}


Bool_t selectTracking(HParticleCand* cand){
    Float_t metaQa   = cand->getMetaMatchQuality();
    Float_t chi2     = cand->getChi2();
    Float_t inchi2   = cand->getInnerSegmentChi2();
    Float_t outchi2  = cand->getOuterSegmentChi2();
    Float_t beta     = cand->getBeta();

    //RKchi2 <50

    //if( metaQa < 4 && chi2 < 50  && inchi2 < 50 && outchi2 < 50){
    //if( metaQa < 4 && chi2 < 250 && chi2>0 && inchi2 < 50 && inchi2 > 0 && outchi2 < 50 && outchi2 > 0 && beta>0 && beta<1.2)
    
    if( chi2 < 50 && chi2>0 && inchi2 < 50 && inchi2 > 0 && outchi2 < 50 && outchi2 > 0 && beta>0 && beta<1.2)
    {
	return  kTRUE;
    }
    else return kFALSE;
}


//=============================================================================================================
//=============================================================================================================
Int_t flow_Analysis(TString inputlist, TString outfile, Int_t nev=-1)
{

    HLoop* loop = new HLoop(kTRUE);  // kTRUE : create Hades  (needed to work with standard eventstructure)

    Bool_t ret = kFALSE;
    if (inputlist.Contains(",")) ret = loop->addMultFiles(inputlist);
    else if(inputlist.Contains(".list")) ret = loop->addFilesList(inputlist);
    else ret = loop->addFiles(inputlist);


    if(!loop->setInput("-*,+HParticleCand,+HParticleEvtInfo,+HStart2Hit,+HStart2Cal,+HWallHit")) { exit(1); }
    loop->printCategories();                                                                                                    // Alles was ein H vorne dran hat--> besser nachfragen (J)

    HParticleCand* cand;
    HParticleEvtInfo* evtinfo;

    HCategory* candCat    = (HCategory*)HCategoryManager::getCategory(catParticleCand);
    if(!candCat) { exit(1); }

    //____________________________________________________________________________________________________MUSS ausserhalb der EVENT-LOOP
    HCategory* fCatWallHit = (HCategory*)HCategoryManager::getCategory(catWallHit,kTRUE,"catWallHit");   //static HCategory * 	getCategory (Short_t num, Int_t warn=0, TString name="")
    
    if(!fCatWallHit) {cout<<"No catWallHit in input!"<<endl;}
    HIterator *iterWallHit = 0;
    iterWallHit = (HIterator*)fCatWallHit->MakeIterator("native");


    Int_t binsxmom = 6000;
    Int_t binsbeta = 600;

    Int_t binMinPhi = -4;
    Int_t binMaxPhi = +4;
    Int_t binsPhi = 800;


    // All tracks
    //TH1F* hChi2    = new TH1F("hChi2",    "", 1000,-100,900);
    //TH1F* hChi2In  = new TH1F("hChi2In",  "", 1000,-100,900);
    //TH1F* hChi2Out = new TH1F("hChi2Out", "", 1000,-100,900);

    // Multiplicity
    TH1F* hmultTofRpc = new TH1F("hmultTofRpc","",250,0,250);
    hmultTofRpc->SetXTitle("MULT ToF + RPC");
    hmultTofRpc->SetYTitle("counts");

    // Vertex
    TH1F* hvertexZ   = new TH1F("hvertexZ", "", 1000,-100,100);
    TH1F* hvertexX   = new TH1F("hvertexX", "", 1000,-100,100);
    TH1F* hvertexY   = new TH1F("hvertexY", "", 1000,-100,100);
    TH2F* hvertexXZ  = new TH2F("hvertexXZ","",300,-80,20,500,-10,10);
    TH2F* hvertexXY  = new TH2F("hvertexXY","",500,-10,10,500,-10,10);

    //_____________________________________________________________________________________________________
    TH2F* hmetaMatchQa[6];
    TH2F* hthetaPhi[6];
    //TH1F* hmass[6];
    TH2F* hmommass[6];
    TH2F* hbetaMom[6];

    for(Int_t s=0; s<6; s++){                                        //___________________________________________SECTION
        //hmetaMatchQa[s] = new TH2F(Form("hmetaMatchQa_sec%i",s),"",binsxmom,-3000,3000,80,-2,18);
	hthetaPhi[s] = new TH2F(Form("hthetaPhi_sec%i",s),"",360,0,360,45,0,90);
	//hmass[s] = new TH1F(Form("hmass_sec%i",s),"",12000,-6000,6000);
	hmommass[s] = new TH2F(Form("hmommass_sec%i",s),"",800,700,2000,400,0,2000);
	//hbetaMom[s] = new TH2F(Form("hbetaMom_sec%i",s),"",binsxmom,-3000,3000,binsbeta,-1.5,1.5);
    }

    //____________________________________________________________________________________________________________
    TH2F* hmommassMult[6];
    TH1F* hphiPMult[6];
    TH2F* hptyP[6];

    for(Int_t m = 0; m < 6; m ++){                                    //__________________________________________MULTIPLICITY
	hmommassMult[m] = new TH2F(Form("hmommassMult_mult%i",m),"",800,-1000,3000,400,0,2000);
	hphiPMult[m]   = new TH1F(Form("hphiPMult_mult%i",m),"",360,-4,4);
        hptyP[m]   = new TH2F(Form("hptyP_mult%i",m)  ,"",40,-1.,3.,75,0,1500);
    }


    //____________________________________________________________________________________________________________
    TH1F* hthetaPMult[6][6];

    TH2F* hthetaPhiPMult[6][6];
    for(Int_t s=0; s<6; s++){
	for(Int_t m=0; m<6; m++){                                   //____________________________________________SECTION_MULTIPLICITY
	    hthetaPMult[s][m]   = new TH1F(Form("hthetaPMult_sec%i_mult%i",s,m),"",360,0,180);
	    //hthetaPhiPMult[s][m]   = new TH2F(Form("hthetaPhiPMult_sec%i_mult%i",s,m),"",360,0,360,45,0,90);
	}
    }

    //___________________________________________________________________________________Protons TWO Analysis     ALL  and HALF (Upper/Lower)
    TH1F * hmassALL = new TH1F("hmassALL", "",4000,0,2000);                                                       //
    TH1F * hmassHALF = new TH1F("hmassHALF", "", 4000,0,2000);                                                    //
                                                                                                                      //
    TH2F* hbetaMomP   = new TH2F("hbetaMomP","hbetaMomP",binsxmom, 0,3500,binsbeta,0.2,1.2);                          //
    TH2F* hbetaMomP_HALF   = new TH2F("hbetaMomP_HALF","hbetaMomP_HALF",binsxmom, 0,3500,binsbeta,0.2,1.2);           //
                                                                                                                      //
    TH1F* hthetaP    = new TH1F("hthetaP","hthetaP",360,0,180);                                                       //
    TH1F* hthetaP_HALF     = new TH1F("hthetaP_HALF","hthetaP_HALF",360,0,180);                                       //
                                                                                                                      //
    TH1F* hphiP     = new TH1F("hphiP","",binsPhi,binMinPhi,binMaxPhi);                                               //
    TH1F* hphiP_HALF     = new TH1F("hphiP_HALF","",binsPhi,binMinPhi,binMaxPhi);                                     //
                                                                                                                      //
    //___________________________________________________________________________________Protons TWO Analysis     ALL  and HALF (Upper/Lower)

    TH1F* hPhi     = new TH1F("hPhi","",800,-4,4);
    TH1F* hRPlanephi     = new TH1F("hRPlanephi","",800,-4,4);

    TH2F * hBetaMean = new TH2F("hBetaMean","",binsxmom,-3000,3000,binsbeta,0,1.5);

    TH1F * hdNdy_vs_y = new TH1F("hdNdy_vs_y","hdNdy_vs_y",400, 0, 2);
    TH1F * hdNdy_vs_y_corr = new TH1F("hdNdy_vs_y_corr","hdNdy_vs_y_corr",400, -1, 1.5);

    TH2F * hV1vsPtP = new TH2F ("hV1vsPtP", "hV1vsPtP", 30, 0, 3000, 2000, -2, 2);
    TH2F * hV2vsPtP = new TH2F ("hV2vsPtP", "hV2vsPtP",  30, 0, 3000, 2000, -2, 2);

                                                                                                                                                              
    //=================================================================================================30.11.2015
    //=======================================================================================================V1

    /*TH2F* hV1vsPtPmult[5][15];                                                                             //[mult] und [y_Bins]
    for(Int_t m = 1; m < 6; m++){
	for(Int_t y = 0; y < 15; y++){
	    hV1vsPtPmult[m][y] = new TH2F(Form("hV1vsPtPmult%i_yBin%i", m, y) ,"", 30, 0, 2500, 2000, -2, 2);
	}
    }
    */

    //=======================================================================================================V2
    /*TH2F* hV2vsPtPmult[5][15];                                                                             //[mult] und [y_Bins]
    for(Int_t m = 1; m < 6; m++){
	for(Int_t y = 0; y < 15; y++){
	    hV2vsPtPmult[m][y] = new TH2F(Form("hV2vsPtPmult%i_yBin%i", m, y) ,"", 30, 0, 2500, 2000, -2, 2);
	}
    }
    */

    //TProfile *profV1vsPt[5][15];
    //TProfile *profV2vsPt[5][15];


    //=====================================================================================================11.12.2015
    TProfile *hv1ProtcandPtY[15];
    TProfile *hv2ProtcandPtY[15];
    //TProfile *hv3ProtcandPtY[15];

    TProfile *hv1ProtcandPtY_HALF[15];
    TProfile *hv2ProtcandPtY_HALF[15];

    Float_t minX = 0;
    Float_t maxX = 2000;
    Float_t binsX=(maxX-minX)/50.;                                                                            // binning p_t


    Float_t minY=-1;
    Float_t maxY= 1;

    for(Int_t k = 0; k < 15; k++){
	hv1ProtcandPtY[k] = new TProfile(Form("hv1ProtcandPtY%i",k), "", binsX, minX, maxX, minY, maxY);
	hv2ProtcandPtY[k] = new TProfile(Form("hv2ProtcandPtY%i",k), "", binsX, minX, maxX, minY, maxY);
	//hv3ProtcandPtY[k] = new TProfile(Form("hv3ProtcandPtY%i",k), "", binsX, minX, maxX, minY, maxY);
        hv1ProtcandPtY_HALF[k] = new TProfile(Form("hv1ProtcandPtY%i_HALF",k), "", binsX, minX, maxX, minY, maxY);
	hv2ProtcandPtY_HALF[k] = new TProfile(Form("hv2ProtcandPtY%i_HALF",k), "", binsX, minX, maxX, minY, maxY);
    }



    //=================================================================================================01.12.2015
    //TH1F * hQVecTrackPhi = new TH1F("hQVecTrackPhi", "hQVecTrackPhi", 200, 0, 400);     //for QVecTrackPhi in Deg


    TH2F * hQVecFWvsTrackPhi = new TH2F("hQVecFWvsTrackPhi", "hQVecFWvsTrackPhi", 360, -4.,4., 360, -4.,4.);         //for QVecTrackPhi in Rad
 

    TH1F * hQVecTrackPhi = new TH1F("hQVecTrackPhi", "hQVecTrackPhi", 360, -4.,4.);         //for QVecTrackPhi in Rad
    TH1F * hQVecFWPhi = new TH1F("hQVecFWPhi", "hQVecFWPhi", 360, -4., 4.);         //for QVecFWPhi in Rad
    TH1F * hQVecDeltaPhiFWvsTrack = new TH1F("hQVecDeltaPhiFWvsTrack", "hQVecDeltaPhiFWvsTrack", 360, -4., 4.);         //Delta between hQVecTrackPhi and hQVecFWPhi   SOPHISTICATED
    TH1F * hQVecDeltaPhiFWvsTrackAdvanced  = new TH1F("hQVecDeltaPhiFWvsTrackAdvanced", "hQVecDeltaPhiFWvsTrackAdvanced", 360, -4., 4.); 
    TH1F * hDiff = new TH1F("hDiff", "hDiff", 360, -4., 4.);


    TF1 *fBetaMeanP = new TF1("fBetaMeanP",fBeta,binsbeta,0,1.5);
    fBetaMeanP->SetParameter(0,HPhysicsConstants::mass(14));


    //==============================================================================================================
    TFile * cutfile_betamom = new TFile("/hera/hades/user/tscheib/apr12/ID_Cuts/BetaMomIDCuts_PionsProtons_gen8_DATA_RK400_PionConstMom.root");
    if (!cutfile_betamom) {cout<<"ID cut root-file not found "<<endl; return 0;}                                  //prevents segmentation violation message


    static TCutG* betamom_cut_tof=(TCutG*)cutfile_betamom->Get("BetaCutProton_TOF_1.0");                          //CUT#1         //CUT#2 _2.0       //CUT#3 _3.5
    static TCutG* betamom_cut_rpc=(TCutG*)cutfile_betamom->Get("BetaCutProton_RPC_1.0");

    if (!betamom_cut_tof) {cout<<"betamom_cut_tof in ID cut root-file not found "<<endl; return 0;}               //prevents segmentation violation message
    if (!betamom_cut_rpc) {cout<<"betamom_cut_rpc in ID cut root-file not found "<<endl; return 0;}               //prevents segmentation violation message

    TFile* out = new TFile(outfile.Data(),"RECREATE");
                                                                                                                   // Standardcuts --> systematisch variieren, drei pt Bereiche waehlen
    TVector2 QVectorTracks;
    TVector2 QVectorFW;


    //--------------------------CONFIGURATION----------------------------------------------------
    //At begin of the program (outside the event loop)
    HParticleTrackSorter sorter;
    //sorter.setDebug();                                            // for debug
    //sorter.setPrintLevel(3);                                      // max prints
    //sorter.setIgnoreInnerMDC();                                   // do not reject Double_t inner MDC hits
    //sorter.setIgnoreOuterMDC();                                   // do not reject Double_t outer MDC hits
    //sorter.setIgnoreMETA();                                       // do not reject Double_t META hits
    //sorter.setIgnorePreviousIndex();                              // do not reject indices from previous selctions
    sorter.init();                                                  // get catgegory pointers etc...
    //--------------------------------------------------------------------------------------------

    TStopwatch timer;
    timer.Reset();
    timer.Start();

    HEnergyLossCorrPar eLoss;
    eLoss.setDefaultPar("apr12");

    Int_t entries = loop->getEntries();
    if(nev < entries && nev >= 0 ) entries = nev;

    //_________________________________________________________________________________ Event LOOP
    //____________________________________________________________________________________________
    for(Int_t i = 1; i < entries; i++)
    {
	//----------break if last event is reached-------------
	//if(!gHades->eventLoop(1)) break;

	if(loop->nextEvent(i) <= 0) {
	    cout<<" end recieved " << endl;
	    break;
	} // last event reached

	//t0reco.execute();    // do T0 rco

	HTool::printProgress(i,entries,1,"Analyze pairs :");
	//----------------------------------------------------------------------------looping data

	// track multiplicity / event
	evtinfo = HCategoryManager::getObject(evtinfo,catParticleEvtInfo,0);
	Int_t mult = 0;
	Int_t multTof = 0;
	Int_t multRpc = 0;
	Int_t multSelectedTracks = evtinfo->getSumSelectedParticleCandMult();

	multTof = evtinfo->getSumTofMultCut();
	multRpc = evtinfo->getSumRpcMultHitCut();

	Int_t multbin = 0;
	mult = multTof+multRpc;
	hmultTofRpc->Fill(mult);

	multbin =getCentralityClassFixedCuts10();
	if (multbin>5) multbin=0;

	//----------------------------------------------------------------------- RP determination
	//----------------------------------------------------------------------------------------

	Double_t RPlanePhi = evtinfo->getRPlanePhi();
	if (RPlanePhi == -1000.) continue;
	Double_t PsiA  = evtinfo->getPhiA();  //in deg
	Double_t PsiB  = evtinfo->getPhiB();  //in deg
	Double_t PsiAB = evtinfo->getPhiAB();  //in deg


	if(RPlanePhi   ==-1000 || PsiA  ==-1000 || PsiB  ==-1000 || PsiAB ==-1000) {
	    continue;
	}

	//	if(RPlanePhi<0) RPlanePhi += 360;

	Double_t  EventPlaneFWPsi =0;

	if(RPlanePhi  !=-1000) EventPlaneFWPsi=TMath::DegToRad()*RPlanePhi;
	hRPlanephi->Fill(EventPlaneFWPsi);

	//------------------------------------------------------------------------EVENT SELECTION
	//---------------------------------------------------------------------------------------
	if(!evtinfo->isGoodEvent(Particle::kGoodTRIGGER) ) continue;
	if(!evtinfo->isGoodEvent(Particle::kGoodSTART)) continue;
        if(!evtinfo->isGoodEvent(Particle::kNoPileUpSTART)) continue;
        if(!evtinfo->isGoodEvent(Particle::kGoodVertexClust)) continue;
	if(!evtinfo->isGoodEvent(Particle::kGoodVertexCand)) continue;
	if(!evtinfo->isGoodEvent(Particle::kNoVETO)) continue;
	if(!evtinfo->isGoodEvent(Particle::kGoodSTARTVETO)) continue;
	if(!evtinfo->isGoodEvent(Particle::kGoodSTARTMETA)) continue;

	//------------------------------------------------------------------------
	// clean vectors and index arrays
	sorter.cleanUp();


	//------------------------------------------------------------------------
	sorter.resetFlags(kTRUE,kTRUE,kTRUE,kTRUE);
	Int_t nCandHad     = sorter.fill(selectHadronsQa);
	Int_t nCandHadBest = sorter.selectBest(HParticleTrackSorter::kIsBestRKRKMETA,HParticleTrackSorter::kIsHadron);

	Int_t nCandNoLept     = sorter.fill(rejectLeptons);
	Int_t nCandNoLeptBest = sorter.selectBest(HParticleTrackSorter::kIsBestRKRKMETA,HParticleTrackSorter::kIsLepton);

	// Get vertex from combined fit of fitted inner segments
	Double_t vertexX = gHades->getCurrentEvent()->getHeader()->getVertexReco().getX();
	Double_t vertexY = gHades->getCurrentEvent()->getHeader()->getVertexReco().getY();
	Double_t vertexZ = gHades->getCurrentEvent()->getHeader()->getVertexReco().getZ();

	hvertexXZ->Fill(vertexZ,vertexX);
	hvertexXY->Fill(vertexY,vertexX);
	hvertexX->Fill(vertexX);
	hvertexY->Fill(vertexY);
	hvertexZ->Fill(vertexZ);


	//----------------looping data-------------------------
	Int_t size = candCat->getEntries();

	Bool_t pFlag[size];

	Double_t x0 = 0;
	Double_t y0 = 0;
	TVector2 QVectorTracks(x0, y0);


	//-----------------------------------------------------------------------------FW QVector
	//---------------------------------------------------------------------------------------

	//HWallHit
	HWallHit *wall = 0;
	iterWallHit->Reset();

	Float_t mult_wall = 0;
	Float_t betaFW    = 0;
	Float_t wallX     = 0;
	Float_t wallY     = 0;
	Float_t wallZ     = 0;

	TVector2 QVectorFW(0., 0.);          // von 0 bis 2pi
	TVector2 vect(0.,0.);
	TVector2 vsum(0.,0.);


	//cout<<">>>  fCatWallHit in macro adress "<< fCatWallHit << endl;
	//cout<<">>>  fCatWallHit->getEntries() in macro "<< fCatWallHit->getEntries() << endl;
	//cout << (HWallHit *)iterWallHit->Next() << "   WallHit???"<< endl;

	while ((wall = (HWallHit *)iterWallHit->Next()) != 0)
	{
            //cout << "DEBUG, in loop" << endl;

	    betaFW = wall->getDistance()/wall->getTime()/299.792458;

	    if(wall->getCell() < 144 && wall->getCharge() > 80)                                      //small
	    {                               
		if(betaFW > 0.84 && betaFW < 1.0)
		{
		    wall->getXYZLab(wallX,wallY,wallZ);
		    vect.Set(wallX, wallY);
		    vsum += vect;
		}
	    }

	    else if(wall->getCell() >= 144 && wall->getCell() < 208 && wall->getCharge() > 85)        //medium
	    {
		if(betaFW > 0.85 && betaFW < 1.0) {
		    wall->getXYZLab(wallX,wallY,wallZ);
		    vect.Set(wallX, wallY);
		    vsum += vect;
		}
	    }

	    else if(wall->getCell() >= 208 && wall->getCharge() > 86)                                //large
	    {
		if(betaFW > 0.8 && betaFW < 1.0)
		{
		    wall->getXYZLab(wallX,wallY,wallZ);
		    vect.Set(wallX, wallY);
		    vsum += vect;
		}
	    }
	    QVectorFW = vsum;                // von 0 bis 2pi
	}



	//___________________________________________________________________Particle LOOP
	//_______________________________________________________________________________

	for(Int_t j = 0; j < size; j ++)
	{
	    pFlag[j]   = kFALSE;

	    Float_t betaMeanP = 0;

	    // Track candidates
	    cand = HCategoryManager::getObject(cand,candCat,j);

	    //_______________________________________________________________________________
            // Track cuts                                                                    
	    if(!cand->isFlagBit(Particle::kIsUsed)) continue;

            cand->setMomentum(1.0065*eLoss.getCorrMom(14,cand->getMomentum(),cand->getTheta())); //better precision in eloss, factor comes from mag field

	    cand->calc4vectorProperties();
	    if(!cand->isFlagBit(Particle::kIsUsed)) continue;
	    //if(cand->getPID()!=14) continue;
	    if(!selectTracking(cand)) continue;      // apply track QA cuts

            HGeomVector fEventVertex = HParticleTool::getGlobalVertex(Particle::kVertexParticle);
	    HGeomVector base, dir;

	    HParticleTool::calcSegVector(cand->getZ(), cand->getR(), (TMath::DegToRad())*cand->getPhi(), (TMath::DegToRad())*cand->getTheta(), base, dir);
	    Float_t dist = HParticleTool::calculateMinimumDistanceStraightToPoint(base, dir, fEventVertex );

	    //Add constraint that particle is primary
	    if(dist>10) continue;


	    Int_t sys        = cand->getSystem();
	    Int_t sysU       = cand->getSystemUsed();
	    Float_t mom      = cand->getMomentum();      
	    Int_t charge     = cand->getCharge();
	    //Float_t metaQa   = cand->getMetaMatchQuality();
	    Float_t beta     = cand->getBeta();
	    Float_t theta    = cand->getTheta();
                                                           

      	    // ------------------------------------------------------------------------------
	    Float_t phi      = cand->getPhi();

	    //convert to radians
	    phi = TMath::DegToRad()*phi;

	    //Define phi range same as Psi range,
	    if (phi>=TMath::Pi()) phi-= 2*TMath::Pi();
	    if (phi<-TMath::Pi()) phi+= 2*TMath::Pi();            


	    Float_t angle = (phi-EventPlaneFWPsi);

	    //Define angle between Psi and phi in the same range as phi
	    if (angle>=TMath::Pi()) angle-= 2*TMath::Pi();
	    if (angle<-TMath::Pi()) angle+= 2*TMath::Pi();

	    hPhi->Fill(angle);


	    // ------------------------------------------------------------------------------
	    //Float_t chi2rk   = cand->getChi2();                                //rk == RungeKutta
	    Int_t sec        = cand->getSector();
	    //Float_t chi2In   = cand->getInnerSegmentChi2();
	    //Float_t chi2Out  = cand->getOuterSegmentChi2();
	    Float_t mass2    = cand->getMass2();  //mass2 = mom*mom*(1-beta*beta)/(beta*beta);
	    Float_t mass     = TMath::Sqrt(cand->getMass2());
	    Float_t mdcdEdx  = cand->getMdcdEdx();
	    Float_t tofdEdx  = cand->getTofdEdx();

	    TLorentzVector tlv1 = (*cand);
	    Float_t mt       = tlv1.Mt();
	    Float_t pt       = tlv1.Pt();

	    Float_t y        = tlv1.Rapidity();
	    Float_t y_corr   = (tlv1.Rapidity()-0.74);       //um Vorwaertsboost "zurueckzusetzen"/ Kurve nach hinten zu schieben

	    


	    //----------------------------------------------------------
	    // allow only system 0 + 1
	    if(sys==2){
		if(cand->isRpcClstUsed() || cand->isShowerUsed() ) sys = 0;
		if(cand->isTofHitUsed()  || cand->isTofClstUsed()) sys = 1;
	    }
	    if(sys==2) {cout<<"still system 2!!! skipped"<<endl;
		continue;
	    }
	    if(sys==-1) {
		cout<<"system -1!!! skipped"<<endl;
		continue;
	    }
	    //----------------------------------------------------------

	    // hChi2->Fill(chi2rk);
	    // hChi2In->Fill(chi2In);
	    // hChi2Out->Fill(chi2Out);

	    // hmetaMatchQa[sec]->Fill(mom*charge,metaQa);
	    // hthetaPhi[sec]->Fill(phi,theta);

	    // hmass[sec]->Fill(charge*TMath::Sqrt(mass2));
            hmassALL->Fill(mass);


	    //hbetaMom[sec]->Fill(mom*charge,beta);
	    hmommass[sec]->Fill(charge*TMath::Sqrt(mass2),mom);
	    hmommassMult[multbin]->Fill(charge*TMath::Sqrt(mass2),mom);

	    //apply track QA cuts
            /*
	    if(chi2rk>10 || chi2rk < 0 || metaQa > 4 )  continue;
	    if(chi2In<0)  continue;
	    if(chi2Out<0) continue;
	    if(metaQa>2)  continue;
	    */

      	    //=========================================================================== Select PROTONS
	    fBetaMeanP->SetParameter(0,HPhysicsConstants::mass(14));
	    betaMeanP = fBetaMeanP->Eval(mom);
	    hBetaMean->Fill(mom*charge,betaMeanP);

	    Int_t ybins    = 15.;
	    Double_t ymin  = 0.09;
	    Double_t ymax  = 1.59;

	    // _____________________________________________________theta-Verteilung   //reicht obiger CUT 3.5  		//(beta<(betaMeanP+0.05) && (beta>(betaMeanP-0.05))) &&
	    if(1 == charge){
		if(sys==0 || sys==1){
		    // _______________________________________________________________________________________________________

		    if((y >  ymin) &&  (y <  ymax)){
			Float_t step_y = TMath::Abs((ymax-ymin)/ybins);
			Int_t yBin = (y - ymin)/step_y;

			/*if(                                                                                                //Attention, this used cut of Timo cuts PROTONS. If Pions are of interest, it should be commented
			   (
			    (sysU==1 && betamom_cut_tof->IsInside(charge*mom,beta))   ||    (sysU==0 && betamom_cut_rpc->IsInside(charge*mom,beta))
			   )
			  )*/
                        if( mass > mParticle*0.95 && mass < mParticle*1.05)                                                 //MASS CUT
			{
			    pFlag[j] = kTRUE;
			    hdNdy_vs_y->Fill(y);
			    hdNdy_vs_y_corr->Fill(y_corr);

			    hthetaP->Fill(theta);
			    hbetaMomP->Fill(mom*charge,beta);
			    hphiP->Fill(angle);

			    //Double_t trackv1 = cos(angle);
			    //Double_t trackv2 = cos(2*angle);

			    //if(multbin==0) continue;

			    //hV1vsPtPmult[multbin][yBin]->Fill(pt,cos(angle));
			    //hV2vsPtPmult[multbin][yBin]->Fill(pt,cos(2*angle));

			    //ALL
			    hv1ProtcandPtY[yBin]->Fill(pt,cos(angle));
			    hv2ProtcandPtY[yBin]->Fill(pt,cos(2*angle));              // was ist der Fehler in TProfile?

			}
			//____________________________________________________________________________________________________________MASS_CUT
			//if( mass > mParticle*1.2 && mass < mParticle*1.5)                                                // Upper Sample  MASS_CUT
			if( mass > mParticle*0.5 && mass < mParticle*0.8)                                         // Lower Sample  MASS_CUT

			    //if(beta > betaMeanP)                                 // Upper Sample_____________________________________betaMomCuts
			    //if(beta<betaMeanP)                                 // Lower Sample
			{
			    hmassHALF->Fill(mass);

			    hthetaP_HALF->Fill(theta);
			    hbetaMomP_HALF->Fill(mom*charge,beta);
			    hphiP_HALF->Fill(angle);
			    hv1ProtcandPtY_HALF[yBin]->Fill(pt,cos(angle));     // was ist der Fehler in TProfile?
			    hv2ProtcandPtY_HALF[yBin]->Fill(pt,cos(2*angle));
			}
		    }
		}
	    } else pFlag[j] = kFALSE;


	    if(pFlag[j] == kTRUE){
		hptyP[multbin]->Fill(y,pt);
		//hthetaPhiPMult[sec][multbin]->Fill(phi,theta);
		hthetaPMult[sec][multbin]->Fill(theta);
		hphiPMult[multbin]->Fill(angle);

		TVector2 QVectTemp;
		QVectTemp.Set(tlv1.X(), tlv1.Y());

		QVectorTracks +=QVectTemp;
		//QVectorTracks.Print();
		//cout <<endl;
	    }

	    /*if(pFlag[j] == kTRUE){                                                     //"Gleiches Ergebnis" wie mit TVec2,  nur mit Array

	    Double_t TofX = 0.;
	    Double_t TofY = 0.;

	    Int_t size = candCat->getEntries();
	    Double_t QVecTofContainer[size][size];

	    for(Int_t w = 0; w < size ; w ++)                                                          //Int_t size = candCat->getEntries();
	    {
	    TofX += tlv1.X();                                                                       //TLorentzVector tlv1 = (*cand);
	    TofY += tlv1.Y();
	    //cout << w << " w" <<endl;
	    QVecTofContainer[w][w] = (TofX, TofY);
	    }
	    //cout << TofX << " TofX" <<endl;
	    //cout << TofY << " TofY" <<endl;
	    hTofXvsTofY->Fill(TofX, TofY);

	    }  */

	}
        //______________________________________________________________________________________end__Particle LOOP
	//________________________________________________________________________________________________________


	Double_t QVecTrackPhi = TVector2::Phi_mpi_pi(QVectorTracks.Phi())  ;
	hQVecTrackPhi->Fill(QVecTrackPhi);

	Double_t QVectFWPhi = TVector2::Phi_mpi_pi(QVectorFW.Phi())  ;
	hQVecFWPhi->Fill(QVectFWPhi);

        hQVecFWvsTrackPhi->Fill(QVectFWPhi, QVecTrackPhi);

        //____________________________________________________________________________________________________

	Double_t QVecDeltaPhiFWvsTrack = TVector2::Phi_mpi_pi(QVectFWPhi-QVecTrackPhi);                                         //NOT advanced

	//Define QVecDeltaPhiFWvsTrack in the same range as QVecTrackPhi and QVectFWPhi
	//if (QVecDeltaPhiFWvsTrack>=TMath::Pi()) QVecDeltaPhiFWvsTrack-= 2*TMath::Pi();
	//if (QVecDeltaPhiFWvsTrack<-TMath::Pi()) QVecDeltaPhiFWvsTrack+= 2*TMath::Pi();
	hQVecDeltaPhiFWvsTrack->Fill(QVecDeltaPhiFWvsTrack);

	/*
	TVector2 Test1(1., 0.);
	TVector2 Test2(0., 1.);
        Double_t Test1Phi = Test1.Phi_mpi_pi(Test1.Phi())  ;
        Double_t Test2Phi = Test2.Phi_mpi_pi(Test2.Phi())  ;
	Double_t TestSimple = (Test1Phi - Test2Phi);
        hQVecDeltaPhiFWvsTrack->Fill(TestSimple);
        */

	Double_t QVecDeltaPhiFWvsTrackAdvanced = QVectorFW.DeltaPhi(QVectorTracks);                            //advanced
        hQVecDeltaPhiFWvsTrackAdvanced->Fill(QVecDeltaPhiFWvsTrackAdvanced);


	Double_t Difference = QVecDeltaPhiFWvsTrack - QVecDeltaPhiFWvsTrackAdvanced;
        if (Difference>=TMath::Pi()) Difference-= 2*TMath::Pi();
	if (Difference<-TMath::Pi()) Difference+= 2*TMath::Pi();
	hDiff->Fill(Difference);
	/*
	TVector2 Test3(1., 0.);
	TVector2 Test4(0., 1.);
	Double_t TestAdvanced = Test3.DeltaPhi(Test4);
        hQVecDeltaPhiFWvsTrackAdvanced->Fill(TestAdvanced);
        */

        //TVector3 Test3D1(1., 0. , 0.);
	//TVector3 Test3D2(0. ,1. , 0.);
        //Double_t Test3DPhi = Test3D1.Angle(Test3D2);                               // Falls man mal TVec3 hat
	

    } // end event-loop


    sorter.finalize();

//==============================================================================
   /* for(Int_t yBin = 0; yBin< 15; yBin++){
	for(Int_t m = 1; m < 6; m++){
	    profV1vsPt[m][yBin] =  hV1vsPtPmult[m][yBin]->ProfileX();
	    profV2vsPt[m][yBin] =  hV2vsPtPmult[m][yBin]->ProfileX();
	}
    }*/


    hV1vsPtP->Write();
    hV2vsPtP->Write();

    /*
    for(Int_t i = 0; i<15; i++){
	for(Int_t m = 1; m<6; m++){
	    hV1vsPtPmult[m][i]->Write();
	    hV2vsPtPmult[m][i]->Write();
	    profV1vsPt[m][i]->Write();
	    profV2vsPt[m][i]->Write();
	}
    }*/
    //TDirectory * DirMomMass = TDirectoryfile(DirMomMass, "DirMomMass");
    //TFile * fOpen
    for(Int_t s=0; s<6; s++){
	//hmetaMatchQa[s]->Write();
	//hthetaPhi[s]->Write();
	hmommass[s]->Write();
	//hmass[s]->Write();
	//hbetaMom[s]->Write();
    }


    timer.Stop();
    out->cd();

    //________________________________________________________________________
    hdNdy_vs_y->Write();
    hdNdy_vs_y_corr->Write();

    // hChi2->Write();
    // hChi2In->Write();
    // hChi2Out->Write();

    for(Int_t m=0; m<6; m++)
    {
	hphiPMult[m]->Write();
        //hmass[m]->Write();
    }

    for(Int_t i=0; i<15; i++)
    {
	//ALL
	hv1ProtcandPtY[i]->Write();
	hv2ProtcandPtY[i]->Write();
        //hv3ProtcandPtYm0[i]->Write();
	//
	hv1ProtcandPtY_HALF[i]->Write();
	hv2ProtcandPtY_HALF[i]->Write();
    }

    
    hQVecTrackPhi->Write();
    hQVecFWPhi->Write();
    hQVecDeltaPhiFWvsTrack->Write();
    hQVecDeltaPhiFWvsTrackAdvanced->Write();
    hDiff->Write();

    hQVecFWvsTrackPhi->Write();

    //___________________________________________________________________Protons TWO Analysis     ALL and HALF (Upper/Lower)

    hmassALL->Write();
    hmassHALF->Write();

    hBetaMean->Write();

    hthetaP->Write();
    hthetaP_HALF->Write();
                                                                                               
    hbetaMomP->Write();
    hbetaMomP_HALF->Write();

    hPhi->Write();
    hphiP->Write();
    hphiP_HALF->Write();

    hmultTofRpc->Write();

    hRPlanephi->Write();



    //____________________________________________________________________________//Canvas

    out->Save();
    out->Close();

    //delete myHades;
    cout<<"####################################################"<<endl;

    return 0 ;
}
