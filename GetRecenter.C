{
	gROOT->Reset(); 
	TString currFName = gROOT->GetFile()->GetName();
    cout <<  gROOT->GetFile()->GetName() << endl;
    TFile *currFile = TFile::Open( currFName.Data(),"READONLY" );
    TIter next( gDirectory->GetListOfKeys() ); //-dir-content--
    const Int_t NHARM = 6;
    const Int_t NMULT = 11;
    const Int_t NDAY  = 13;

    FILE* output = fopen("Recenter.cc","wt");
    TKey *key; 
    TProfile* hsumX[3][NMULT];
    TProfile* hsumY[3][NMULT];
    while( (key = (TKey*)next()) ){ 
    	//printf("Class=%s Name=%s Title=\"%s\"\n", key->GetClassName(), key->GetName(), key->GetTitle());
    	TObject *obj = key->ReadObj();
        if ( obj->IsA()->InheritsFrom( TProfile::Class() ) ) {
        	TProfile *h = (TProfile*)obj;
        	Char_t ch1[1024], ch2[1024];
            Int_t iC1, iC2; // nMultiplicity;
            sscanf(key->GetTitle(), "%s %i %i" , ch1, &iC1, &iC2);
            TString s(key->GetName());
            if( s.Contains("hsumXmean") ){
            	printf("[EP]-%s[%i][%i]--->\n",ch1,iC1,iC2);
                hsumX[iC1][iC2] = new TProfile;
                hsumX[iC1][iC2] = (TProfile*)obj;
            }
            if( s.Contains("hsumYmean") ){
            	printf("[EP]-%s[%i][%i]--->\n",ch1,iC1,iC2);
                hsumY[iC1][iC2] = new TProfile;
                hsumY[iC1][iC2] = (TProfile*)obj;
            }
        }
    }

    cout << "reading complete" << endl;
    fprintf(output,"//X and Y shift in terms of normalized TVector2;\n");
    for (Int_t iT=0;iT<3;iT++){
        for(Int_t imult=0;imult<NMULT;imult++){
            for(Int_t iday=0;iday<NDAY;iday++){
                fprintf(output,"sumXmean[%i][%2i][%2i] = %2.10f; sumYmean[%i][%2i][%2i] = %2.10f; sumXsigma[%i][%2i][%2i] = %2.10f; sumYsigma[%i][%2i][%2i] = %2.10f;\n",iT,imult,iday,hsumX[iT][imult]->GetBinContent(iday+2),iT,imult,iday,hsumY[iT][imult]->GetBinContent(iday+2),iT,imult,iday,hsumX[iT][imult]->GetBinError(iday+2)*sqrt(hsumX[iT][imult]->GetEntries()),iT,imult,iday,hsumY[iT][imult]->GetBinError(iday+2)*sqrt(hsumY[iT][imult]->GetEntries()));
    	    }
        }
    }
    fclose(output);


    TObject* object;
    TProfile* QvX[2][2];
    TProfile* QvY[2][2];

    FILE* output1 = fopen("RecenterMETA.cc","wt");
    for (Int_t iT=0;iT<3;iT++){
        key = gDirectory->GetKey(Form("hQvsM_X%i%i",0,iT));
        if( key = gDirectory->GetKey(Form("hQvsM_X%i%i",0,iT)) ) cout << key->GetName() << endl;
        object = key->ReadObj();
        QvX[0][iT] = (TProfile*)object;

        key = gDirectory->GetKey(Form("hQvsM_Y%i%i",0,iT));
        if( key = gDirectory->GetKey(Form("hQvsM_Y%i%i",0,iT)) ) cout << key->GetName() << endl;
        object = key->ReadObj();
        QvY[0][iT] = (TProfile*)object;

        for (Int_t i=0;i<215;i++){
            fprintf(output1,"Qxmean[%i][%i] = %2.10f; Qymean[%i][%i] = %2.10f; Qxsigm[%i][%i] = %2.10f; Qysigm[%i][%i] = %2.10f;\n",i,iT,QvX[0][iT]->GetBinContent(i+1),i,iT,QvY[0][iT]->GetBinContent(i+1),i,iT,QvX[0][iT]->GetBinError(i+1),i,iT,QvY[0][iT]->GetBinError(i+1));
        }
    }
    cout << "reading complite 2" << endl;
    fclose(output1);

    FILE* output2 = fopen("RecenterFW.cc","wt");
    for (Int_t iT=0;iT<3;iT++){
        key = gDirectory->GetKey(Form("hQvFW_X%i%i",0,iT));
        if( key = gDirectory->GetKey(Form("hQvFW_X%i%i",0,iT)) ) cout << key->GetName() << endl;
        object = key->ReadObj();
        QvX[1][iT] = (TProfile*)object;

        key = gDirectory->GetKey(Form("hQvFW_Y%i%i",0,iT));
        if( key = gDirectory->GetKey(Form("hQvFW_Y%i%i",0,iT)) ) cout << key->GetName() << endl;
        object = key->ReadObj();
        QvY[1][iT] = (TProfile*)object;

        for (Int_t i=0;i<100;i++){
            fprintf(output2,"mQxFW[%i][%i] = %2.10f; mQyFW[%i][%i] = %2.10f; sQxFW[%i][%i] = %2.10f; sQyFW[%i][%i] = %2.10f;\n",i,iT,QvX[1][iT]->GetBinContent(i+1),i,iT,QvY[1][iT]->GetBinContent(i+1),i,iT,QvX[1][iT]->GetBinError(i+1),i,iT,QvY[1][iT]->GetBinError(i+1));
        }
    }
    cout << "reading complite 3" << endl;
    fclose(output2);
}