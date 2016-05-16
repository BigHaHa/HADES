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
    TProfile* QvX[2];
    TProfile* QvY[2];
    FILE* output1 = fopen("RecenterMETA.cc","wt");
    key = gDirectory->GetKey(Form("hQvsM_X%i",0));
    if( key = gDirectory->GetKey(Form("hQvsM_X%i",0)) ) cout << key->GetName() << endl;
    object = key->ReadObj();
    QvX[0] = (TProfile*)object;

    key = gDirectory->GetKey(Form("hQvsM_Y%i",0));
    if( key = gDirectory->GetKey(Form("hQvsM_Y%i",0)) ) cout << key->GetName() << endl;
    object = key->ReadObj();
    QvY[0] = (TProfile*)object;

    cout << "reading complite 2" << endl;

    for (Int_t i=0;i<215;i++){
        fprintf(output1,"Qxmean[%i] = %2.10f; Qymean[%i] = %2.10f; Qxsigm[%i] = %2.10f; Qysigm[%i] = %2.10f;\n",i,QvX[0]->GetBinContent(i+1),i,QvY[0]->GetBinContent(i+1),i,QvX[0]->GetBinError(i+1),i,QvY[0]->GetBinError(i+1));
    }
    fclose(output1);

    FILE* output2 = fopen("RecenterFW.cc","wt");
    key = gDirectory->GetKey(Form("hQvFW_X%i",0));
    if( key = gDirectory->GetKey(Form("hQvFW_X%i",0)) ) cout << key->GetName() << endl;
    object = key->ReadObj();
    QvX[1] = (TProfile*)object;

    key = gDirectory->GetKey(Form("hQvFW_Y%i",0));
    if( key = gDirectory->GetKey(Form("hQvFW_Y%i",0)) ) cout << key->GetName() << endl;
    object = key->ReadObj();
    QvY[1] = (TProfile*)object;

    cout << "reading complite 3" << endl;

    for (Int_t i=0;i<100;i++){
        fprintf(output2,"mQxFW[%i] = %2.10f; mQyFW[%i] = %2.10f; sQxFW[%i] = %2.10f; sQyFW[%i] = %2.10f;\n",i,QvX[1]->GetBinContent(i+1),i,QvY[1]->GetBinContent(i+1),i,QvX[1]->GetBinError(i+1),i,QvY[1]->GetBinError(i+1));
    }
    fclose(output2);
}