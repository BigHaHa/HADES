{
	#include "ReadHistFunc.C"
	#include "FuncPlot.C"

	FileHist f;
    TPlot    p;

    //---------------------------------------------------------------------------------------------------------------------------------------------------------//
    
	TH2F* hFW[11];
	TH1F* hPhiPR[3];
	for(Int_t iC=0; iC<11; iC++){
		hFW[iC] = f.ReadTH2F("hFWxyCCsmearM",iC);
	}
    hPhiPR[0] = f.ReadTH1F("hPhiPRpCent",4,0,3);
    hPhiPR[1] = f.ReadTH1F("hPhiPRRCent",4,0,3);
    hPhiPR[2] = f.ReadTH1F("hPhiPRFCent",4,0,3);

    for (Int_t i=0;i<3;i++){ hPhiPR[i]->Rebin(10); }

    TCanvas* canv;

    canv = p.makePlot(hPhiPR[0],"Raw",hPhiPR[1],"Recentered","#varphi-#Psi_{EP}, deg","#frac{dN}{d(#varphi-#Psi_{EP})}","R");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Doc/phiPRR%i%i%i.png",4,0,3),"recreate");

    canv = p.makePlot(hPhiPR[0],"Raw",hPhiPR[2],"Flattened","#varphi-#Psi_{EP}, deg","#frac{dN}{d(#varphi-#Psi_{EP})}","R");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Doc/phiPRF%i%i%i.png",4,0,3),"recreate");
    
    //---------------------------------------------------------------------------------------------------------------------------------------------------------//
    
    
    TProfile* QxMult[2][2];
    TProfile* QyMult[2][2];

    for (Int_t iT=0;iT<2;iT++){
        QxMult[iT][0] = f.ReadTProfile("hQvsM_X",iT,0);
        QxMult[iT][1] = f.ReadTProfile("hQvFW_X",iT,0);
        QyMult[iT][0] = f.ReadTProfile("hQvsM_Y",iT,0);
        QyMult[iT][1] = f.ReadTProfile("hQvFW_Y",iT,0);
    }
    
    canv = p.makePlot(QxMult[0][0],"Raw",QxMult[1][0],"Recentered","TOF+RPC","<Q_{x}>",0.,300.,-0.09,0.05,"S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Doc/QxMETA.png"),"recreate");
    canv = p.makePlot(QyMult[0][0],"Raw",QyMult[1][0],"Recentered","TOF+RPC","<Q_{y}>","S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Doc/QyMETA.png"),"recreate");

    canv = p.makePlot(QxMult[0][1],"Raw",QxMult[1][1],"Recentered","FW","<Q_{x}>","S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Doc/QxFW.png"),"recreate");
    canv = p.makePlot(QyMult[0][1],"Raw",QyMult[1][1],"Recentered","FW","<Q_{y}>","S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Doc/QyFW.png"),"recreate");
    
    //---------------------------------------------------------------------------------------------------------------------------------------------------------//
    
    
    TProfile* CosPsiAB[2][3][2]; //[iM=META/FW][iT=Raw/Recentered/Flattened][n=1,2 harmonics]
    TProfile* SinPsiAB[2][3][2];

    for(Int_t iT=0;iT<3;iT++){
        for(Int_t n=0;n<2;n++){
            CosPsiAB[0][iT][n] = f.ReadTProfile("CosPsiAB_META",iT,n);
            CosPsiAB[1][iT][n] = f.ReadTProfile("CosPsiAB_FW"  ,iT,n);
            SinPsiAB[0][iT][n] = f.ReadTProfile("SinPsiAB_META",iT,n);
            SinPsiAB[1][iT][n] = f.ReadTProfile("SinPsiAB_FW"  ,iT,n);
        }
    }

    for (Int_t iM=0;iM<2;iM++){
        for(Int_t iT=0;iT<3;iT++){
            for(Int_t n=0;n<2;n++){
                if(iM==0) CosPsiAB[iM][iT][n]->Rebin(10);
                if(iM==0) SinPsiAB[iM][iT][n]->Rebin(10);
                if(iM==1) CosPsiAB[iM][iT][n]->Rebin(5);
                if(iM==1) SinPsiAB[iM][iT][n]->Rebin(5);
            }
        }
    }

    canv = p.makePlot(CosPsiAB[0][0][0],"Raw",CosPsiAB[0][1][0],"Recentered",CosPsiAB[0][2][0],"Flattened","TOF+RPC","<cos(#varphi_{A}-#varphi_{B})>","S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Resolution/CosPsiAB_META.png"),"recreate");
    canv = p.makePlot(SinPsiAB[0][0][0],"Raw",SinPsiAB[0][1][0],"Recentered",SinPsiAB[0][2][0],"Flattened","TOF+RPC","<sin(#varphi_{A}-#varphi_{B})>","S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Resolution/SinPsiAB_META.png"),"recreate");
    canv = p.makePlot(CosPsiAB[1][0][0],"Raw",CosPsiAB[1][1][0],"Recentered",CosPsiAB[1][2][0],"Flattened","FW"     ,"<cos(#varphi_{A}-#varphi_{B})>","S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Resolution/CosPsiAB_FW.png"),"recreate");
    canv = p.makePlot(SinPsiAB[1][0][0],"Raw",SinPsiAB[1][1][0],"Recentered",SinPsiAB[1][2][0],"Flattened","FW"     ,"<sin(#varphi_{A}-#varphi_{B})>","S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Resolution/SinPsiAB_FW.png"),"recreate");

    canv = p.makePlot(CosPsiAB[0][0][1],"Raw",CosPsiAB[0][1][1],"Recentered",CosPsiAB[0][2][1],"Flattened","TOF+RPC","<cos(2(#varphi_{A}-#varphi_{B}))>","S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Resolution/CosPsi2AB_META.png"),"recreate");
    canv = p.makePlot(SinPsiAB[0][0][1],"Raw",SinPsiAB[0][1][1],"Recentered",SinPsiAB[0][2][1],"Flattened","TOF+RPC","<sin(2(#varphi_{A}-#varphi_{B}))>","S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Resolution/SinPsi2AB_META.png"),"recreate");
    canv = p.makePlot(CosPsiAB[1][0][1],"Raw",CosPsiAB[1][1][1],"Recentered",CosPsiAB[1][2][1],"Flattened","FW"     ,"<cos(2(#varphi_{A}-#varphi_{B}))>","S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Resolution/CosPsi2AB_FW.png"),"recreate");
    canv = p.makePlot(SinPsiAB[1][0][1],"Raw",SinPsiAB[1][1][1],"Recentered",SinPsiAB[1][2][1],"Flattened","FW"     ,"<sin(2(#varphi_{A}-#varphi_{B}))>","S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Resolution/SinPsi2AB_FW.png"),"recreate");
    
    //---------------------------------------------------------------------------------------------------------------------------------------------------------//

    
    TProfile* CosPsi[2][3][2][2]; //[iM=META/FW][TYPE][n][k=A/B]
    TProfile* SinPsi[2][3][2][2];

    for (Int_t iT=0;iT<3;iT++){
        for (Int_t n=0;n<2;n++){
            for (Int_t k=0;k<2;k++){
                CosPsi[0][iT][n][k] = f.ReadTProfile("CosPsi_META",iT,n,k);
                CosPsi[1][iT][n][k] = f.ReadTProfile("CosPsi_FW"  ,iT,n,k);
                SinPsi[0][iT][n][k] = f.ReadTProfile("SinPsi_META",iT,n,k);
                SinPsi[1][iT][n][k] = f.ReadTProfile("SinPsi_FW"  ,iT,n,k);
            }
        }
    }
    for (Int_t iM=0;iM<2;iM++){
        for (Int_t iT=0;iT<3;iT++){
            for (Int_t n=0;n<2;n++){
                for (Int_t k=0;k<2;k++){
                    if(iM==0) CosPsi[iM][iT][n][k]->Rebin(10);
                    if(iM==0) SinPsi[iM][iT][n][k]->Rebin(10);
                    if(iM==1) CosPsi[iM][iT][n][k]->Rebin(5 );
                    if(iM==1) SinPsi[iM][iT][n][k]->Rebin(5 );
                }
            }
        }
    }

    canv = p.makePlot(CosPsi[0][0][0][0],"Raw",CosPsi[0][1][0][0],"Recentered",CosPsi[0][2][0][0],"Flattened","TOF+RPC","<cos(#varphi_{A})>","S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Resolution/CosPsiA_META.png"),"recreate");
    canv = p.makePlot(CosPsi[0][0][0][1],"Raw",CosPsi[0][1][0][1],"Recentered",CosPsi[0][2][0][1],"Flattened","TOF+RPC","<cos(#varphi_{B})>","S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Resolution/CosPsiB_META.png"),"recreate");
    canv = p.makePlot(SinPsi[0][0][0][0],"Raw",SinPsi[0][1][0][0],"Recentered",SinPsi[0][2][0][0],"Flattened","TOF+RPC","<sin(#varphi_{A})>","S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Resolution/SinPsiA_META.png"),"recreate");
    canv = p.makePlot(SinPsi[0][0][0][1],"Raw",SinPsi[0][1][0][1],"Recentered",SinPsi[0][2][0][1],"Flattened","TOF+RPC","<sin(#varphi_{B})>","S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Resolution/SinPsiB_META.png"),"recreate");

    canv = p.makePlot(CosPsi[1][0][0][0],"Raw",CosPsi[1][1][0][0],"Recentered",CosPsi[1][2][0][0],"Flattened","FW"     ,"<cos(#varphi_{A})>","S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Resolution/CosPsiA_FW.png"),"recreate");
    canv = p.makePlot(CosPsi[1][0][0][1],"Raw",CosPsi[1][1][0][1],"Recentered",CosPsi[1][2][0][1],"Flattened","FW"     ,"<cos(#varphi_{B})>","S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Resolution/CosPsiB_FW.png"),"recreate");
    canv = p.makePlot(SinPsi[1][0][0][0],"Raw",SinPsi[1][1][0][0],"Recentered",SinPsi[1][2][0][0],"Flattened","FW"     ,"<sin(#varphi_{A})>","S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Resolution/SinPsiA_FW.png"),"recreate");
    canv = p.makePlot(SinPsi[1][0][0][1],"Raw",SinPsi[1][1][0][1],"Recentered",SinPsi[1][2][0][1],"Flattened","FW"     ,"<sin(#varphi_{B})>","S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Resolution/SinPsiB_FW.png"),"recreate");


    canv = p.makePlot(CosPsi[0][0][1][0],"Raw",CosPsi[0][1][1][0],"Recentered",CosPsi[0][2][1][0],"Flattened","TOF+RPC","<cos(2#varphi_{A})>","S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Resolution/Cos2PsiA_META.png"),"recreate");
    canv = p.makePlot(CosPsi[0][0][1][1],"Raw",CosPsi[0][1][1][1],"Recentered",CosPsi[0][2][1][1],"Flattened","TOF+RPC","<cos(2#varphi_{B})>","S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Resolution/Cos2PsiB_META.png"),"recreate");
    canv = p.makePlot(SinPsi[0][0][1][0],"Raw",SinPsi[0][1][1][0],"Recentered",SinPsi[0][2][1][0],"Flattened","TOF+RPC","<sin(2#varphi_{A})>","S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Resolution/Sin2PsiA_META.png"),"recreate");
    canv = p.makePlot(SinPsi[0][0][1][1],"Raw",SinPsi[0][1][1][1],"Recentered",SinPsi[0][2][1][1],"Flattened","TOF+RPC","<sin(2#varphi_{B})>","S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Resolution/Sin2PsiB_META.png"),"recreate");

    canv = p.makePlot(CosPsi[1][0][1][0],"Raw",CosPsi[1][1][1][0],"Recentered",CosPsi[1][2][1][0],"Flattened","FW"     ,"<cos(2#varphi_{A})>","S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Resolution/Cos2PsiA_FW.png"),"recreate");
    canv = p.makePlot(CosPsi[1][0][1][1],"Raw",CosPsi[1][1][1][1],"Recentered",CosPsi[1][2][1][1],"Flattened","FW"     ,"<cos(2#varphi_{B})>","S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Resolution/Cos2PsiB_FW.png"),"recreate");
    canv = p.makePlot(SinPsi[1][0][1][0],"Raw",SinPsi[1][1][1][0],"Recentered",SinPsi[1][2][1][0],"Flattened","FW"     ,"<sin(2#varphi_{A})>","S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Resolution/Sin2PsiA_FW.png"),"recreate");
    canv = p.makePlot(SinPsi[1][0][1][1],"Raw",SinPsi[1][1][1][1],"Recentered",SinPsi[1][2][1][1],"Flattened","FW"     ,"<sin(2#varphi_{B})>","S");
    canv->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/Resolution/Sin2PsiB_FW.png"),"recreate");
    
    //---------------------------------------------------------------------------------------------------------------------------------------------------------//
}