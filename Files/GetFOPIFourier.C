{
	gROOT->Reset(); 
	TString currFName = gROOT->GetFile()->GetName();
    cout <<  gROOT->GetFile()->GetName() << endl;
    TFile *currFile = TFile::Open( currFName.Data(),"READONLY" );
    currFile->Print();
    //----------------------------------------------------------------------------------------------------//  
    gStyle->SetPalette(1);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadBorderSize(0);
    gStyle->SetPadRightMargin(0.16);
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetPadBottomMargin(0.08);
    gStyle->SetLabelOffset(0.015,"xyz");
    gStyle->SetLabelSize(0.042,"xyz");
    gStyle->SetLabelFont(62,"xyz") ;
    //-----------------------//------------------//
    //gStyle->SetTitleX(0.0) //title X location -//
    gStyle->SetTitleW(0.9);  //width-of-title ---//
    gStyle->SetTitleFontSize(0.1); //shifts-title-a-bit-down-
    gStyle->SetTitleFont(42,"t") ;
    gROOT->ForceStyle(); //-needed-to-apply-global-settings-(otherwice-histograms-keep-their-own-settings)--
    //----------------------------------------------------------------------------------------------------//
    TIter next( gDirectory->GetListOfKeys() ); //-dir-content--

    Int_t oncePRsum=0;
    const Int_t color[3] = {2,8,9};

    TH1F *hEPyc[ 3];
    TH1F *hOEPyc[3];
    TH1F *hEPyr[ 3];
    TH1F *hCnCent[3][6];
    TH1F *hCnCent[3][6];
    for(Int_t i=0;i<3;i++){
        for(Int_t n=0;n<6;n++){
            hCnCent[i][n] = new TH1F(Form("hCnCent%i%i",i,n),Form("C_{%i} vs Centrality %i",n+1,i+1),3,1.,4.);
            hCnCent[i][n]->SetLineColor(color[i]);
        }
    }

    TF1 flows("flows", "[6]*( 1.+2.*[0]*TMath::Cos(x[0]*3.1415926535/180.)+2.*[1]*TMath::Cos(2.*x[0]*3.1415926535/180.)+2.*[2]*TMath::Cos(3.*x[0]*3.1415926535/180.)+2.*[3]*TMath::Cos(4.*x[0]*3.1415926535/180.) + [4]*TMath::Cos(5.*x[0]*3.1415/180.) + [5]*TMath::Cos(6.*x[0]*3.1415/180.) + [7]*TMath::Sin(x[0]*3.1415/180.) + [8]*TMath::Sin(2.*x[0]*3.1415/180.) )", -180., 180.);
    
    FILE* output = fopen("Fourier.txt","wt");
    //--reading-flow-histograms-from-ROOT-file-------------(begin)----
    TKey *key; 
    while( (key = (TKey*)next()) ){ 
    	//printf("Class=%s Name=%s Title=\"%s\"\n", key->GetClassName(), key->GetName(), key->GetTitle());
    	TObject *obj = key->ReadObj();
        if ( obj->IsA()->InheritsFrom( TH1::Class() ) ) {
        	TH1 *h = (TH1*)obj;
        	Char_t ch1[1024], ch2[1024];
            Int_t iC;
            sscanf(key->GetTitle(), "%s %i %s "    , ch1, &iC, ch2);
            TString s(key->GetName());
            if( s.Contains("hFOPIWEP") ){
            	printf("[EP]-EP--------%s[%i]--->\n",ch1,iC);
                hEPyc[iC] = new TH1F;
                hEPyc[iC] = (TH1F*)obj;
            }
            if( s.Contains("hRecFOPI") ){
            	printf("[EP]-Recent----%s[%i]--->\n",ch1,iC);
                hEPyr[iC] = new TH1F;
                hEPyr[iC] = (TH1F*)obj;
            }
            if( s.Contains("hCorFOPI") ){
            	printf("[EP]-Corrected-%s[%i]--->\n",ch1,iC);
                hOEPyc[iC] = new TH1F;
                hOEPyc[iC] = (TH1F*)obj;
            }
        }
    }

    cout << "reading complete" << endl;

    Float_t vi[6];
    Float_t ei[6];
    Int_t    div = 3, k = 1, l = 0, dim = 1;
    TText txt;
    Char_t ds[100];
    TCanvas* c[dim];
    for(Int_t i=0;i<dim;i++){
    	c[i] = new TCanvas(Form("c[%i]",i),Form("PhiEP %i",i+1),1000,300);
    	c[i]-> Divide(3,1);
    }

    for (Int_t im=0;im<3;im++){
    	if(l > dim){ break; }
    	c[l]->cd(k);
    	Float_t summ  = hEPyc[im]->GetSumOfWeights(); 
        Int_t phibins = hEPyc[im]->GetXaxis()->GetNbins();
        hEPyc[im]->SetMinimum(0.);
        //hOEPyc[im]->SetMinimum(0.);
        //if (im == 1) hEPyc[im]->GetYaxis()->SetRangeUser(0.,15e3);
        //if (im == 3) hEPyc[im]->GetYaxis()->SetRangeUser(0.,6e4);
        //if (im  > 3) hEPyc[im]->GetYaxis()->SetRangeUser(0.,1e5);
        hEPyc[ im]->SetLineColor(2);
        hEPyr[ im]->SetLineColor(8);
        hOEPyc[im]->SetLineColor(9);
        hEPyc[im]->Draw();
        hOEPyc[im]->Draw("same");
        hEPyr[im]->Draw("same");
        fprintf(output,"-------------EP----------------\n");
        if( summ > 0 ){ 
        	flows.SetLineColor(2);
        	hEPyc[im]->Fit("flows");
            for(Int_t i=0; i<6; i++){
                vi[i] = flows->GetParameter(i);
                ei[i] = flows->GetParError( i);
                if(i==0) { fprintf(output,"[*]-N-Cent-[%i]\n", im); }
                fprintf(output,"Fitted: v%i=%+5.6f +/- %5.6f;\n",i+1,vi[i],ei[i]);
                sprintf(ds, "c%i=%+5.6f +/- %5.6f",i+1, vi[i], ei[i]);
                //txt.DrawTextNDC(0.25, 0.36+0.05  , "META");
                //txt.DrawTextNDC(0.15, 0.36-0.05*i, ds);
                hCnCent[0][i]->SetBinContent(im+1,vi[i]);
                hCnCent[0][i]->SetBinError(  im+1,ei[i]);
            }
        }
        fprintf(output,"----------Recentering----------\n");
        if( hEPyr[im]->GetSumOfWeights() > 0 ){
            flows.SetLineColor(8);
            hEPyr[im]->Fit("flows");
            for(Int_t i=0; i<6; i++){
                vi[i] = flows->GetParameter(i);
                ei[i] = flows->GetParError( i);
                if(i==0) { fprintf(output,"[*]-N-Cent-[%i]\n", im); }
                fprintf(output,"Fitted: v%i=%+5.4f +/- %5.4f;\n",i+1,vi[i],ei[i]);
                sprintf(ds, "v%i=%+5.4f +/- %5.4f",i+1, vi[i], ei[i]);
                //txt.DrawTextNDC(0.62, 0.36+0.05 , "MDC");
                //txt.DrawTextNDC(0.5, 0.36-0.05*i, ds);
                hCnCent[1][i]->SetBinContent(im+1,vi[i]);
                hCnCent[1][i]->SetBinError(  im+1,ei[i]);
            }
        }
        fprintf(output,"----------Flattering-----------\n");
        if( hOEPyc[im]->GetSumOfWeights() > 0 ){
            flows.SetLineColor(9);
            hOEPyc[im]->Fit("flows");
            for(Int_t i=0; i<6; i++){
                vi[i] = flows->GetParameter(i);
                ei[i] = flows->GetParError( i);
                if(i==0) { fprintf(output,"[*]-N-Cent-[%i]\n", im); }
                fprintf(output,"Fitted: v%i=%+5.4f +/- %5.4f;\n",i+1,vi[i],ei[i]);
                sprintf(ds, "v%i=%+5.4f +/- %5.4f",i+1, vi[i], ei[i]);
                //txt.DrawTextNDC(0.62, 0.36+0.05 , "MDC");
                //txt.DrawTextNDC(0.5, 0.36-0.05*i, ds);
                hCnCent[2][i]->SetBinContent(im+1,vi[i]);
                hCnCent[2][i]->SetBinError(  im+1,ei[i]);
            }
        }

        c[0]->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/PhiEP/PhiEPFOPI%2.0d.png",1),"recreate");
        if(k >= div){ k=0; l++; }
       
        k++;
    }
    TCanvas* canvCent = new TCanvas("canvCent","Coeff. vs Centrality bins",1000,700);
    TLatex* Ltxt = new TLatex();
    TLine*  zero = new TLine();
    canvCent->Divide(3,2);
    for(Int_t i=0;i<3;i++){
        for(Int_t n=0;n<6;n++){
    	    canvCent->cd(n+1);
    	    zero->SetLineWidth(2);
    	    zero->SetLineStyle(2);
            zero->DrawLine(1.,0.,4.,0.);
    	    if(n==0) hCnCent[i][n]->GetYaxis()->SetRangeUser(-0.08,0.02);
    	    if(n==1) hCnCent[i][n]->GetYaxis()->SetRangeUser(-0.001,0.01);
    	    if(n==2) hCnCent[i][n]->GetYaxis()->SetRangeUser(-0.01,0.002);
            if(n==3) hCnCent[i][n]->GetYaxis()->SetRangeUser(-0.0025,0.006);
    	    if(n==4) hCnCent[i][n]->GetYaxis()->SetRangeUser(-0.02,0.005);
            if(n==5) hCnCent[i][n]->GetYaxis()->SetRangeUser(-0.005,0.015);
    	    if(i==0) hCnCent[i][n]->Draw("e");
    	    if(i!=0) hCnCent[i][n]->Draw("esame");
    	    //Ltxt->SetTextSize(0.02);
    	    //Ltxt->DrawLatexNDC(0.1,0.15,"#frac{dN}{d#Psi}=1+2c_{1}cos(#Psi)+2c_{2}cos(2#Psi)+2c_{3}cos(3#Psi)+2c_{4}cos(4#Psi)+c_{5}sin(#Psi)+c_{6}sin(2#Psi)");
        }
    }
    canvCent->SaveAs(Form("/home/peter/Documents/WorkLocal/WorkFiles/HADES/pionflow/anaflow/Results/PhiEP/CoeffFOPI.png"),"recreate");
}