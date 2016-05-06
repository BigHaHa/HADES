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
    TH2F* hshift[NMULT];
    while( (key = (TKey*)next()) ){ 
        //printf("Class=%s Name=%s Title=\"%s\"\n", key->GetClassName(), key->GetName(), key->GetTitle());
        TObject *obj = key->ReadObj();
        if ( obj->IsA()->InheritsFrom( TH2::Class() ) ) {
            TH2 *h = (TH2*)obj;
            Char_t ch1[1024], ch2[1024];
            Int_t iC1; // nMultiplicity;
            sscanf(key->GetTitle(), "%s %i" , ch1, &iC1);
            TString s(key->GetName());
            if( s.Contains("hFWxyCCsmearM") ){
                printf("[EP]-%s[%i]--->\n",ch1,iC1);
                hshift[iC1] = new TH2F;
                hshift[iC1] = (TH2F*)obj;
            }
        }
    }
    cout << "reading complete" << endl;
    fprintf(output,"//X and Y shift in terms of FW hits;\n");
    for(Int_t imult=0;imult<NMULT;imult++){
        fprintf(output,"FWmeanX[%2i] = %2.6f; FWmeanY[%2i] = %2.6f;\n",imult,hshift[imult]->GetMean(1),imult,hshift[imult]->GetMean(2),imult);
    }
    fclose(output);
}