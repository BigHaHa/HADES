#include <cstring>
class FileHist{
  public:
 	FileHist(){
        ReadFileName();
 	}

    void  ReadFileName(){
        currFName = gROOT->GetFile()->GetName();
        cout <<  gROOT->GetFile()->GetName() << endl;
        currFile = TFile::Open( currFName.Data(),"READONLY" );
    }

    TH1F* ReadTH1F(TString histname){
        key = gDirectory->GetKey(Form(histname.Data()));
        obj = key->ReadObj();
        if(strcmp("TH1F",key->GetClassName())==0){//if ( obj->IsA()->InheritsFrom( TH1::Class() ) ){
            cout << key->GetClassName() << " " <<key->GetName() << endl;
            return (TH1F*)obj;
        }
        else{
    	    cout << "Wrong file type (TH1F) TH1F::ReadHist(TString histogram_name)" << endl;
    	    return 0;
        }
    }

    TH2F* ReadTH2F(TString histname){
        key = gDirectory->GetKey(Form(histname.Data()));
        obj = key->ReadObj();
        if(strcmp("TH2F",key->GetClassName())==0){//if ( obj->IsA()->InheritsFrom( TH2::Class() ) ){
            cout << key->GetClassName() << " " <<key->GetName() << endl;
            return (TH2F*)obj;
        }
        else{
    	    cout << "Wrong file type (TH2F) TH2F::ReadHist(TString histogram_name)" << endl;
    	    return 0;
        }
    }

    TProfile* ReadTProfile(TString histname){
        key = gDirectory->GetKey(Form(histname.Data()));
        obj = key->ReadObj();
        if(strcmp("TProfile",key->GetClassName())==0){//if ( obj->IsA()->InheritsFrom( TProfile::Class() ) ){
            cout << key->GetClassName() << " " <<key->GetName() << endl;
            return (TProfile*)obj;
        }
        else{
            cout << "Wrong file type (TH2F) TH2F::ReadHist(TString histogram_name)" << endl;
            return 0;
        }
    }

    TH1F* ReadTH1F(TString histname, Int_t iC){
        key = gDirectory->GetKey(Form("%s%i",histname.Data(),iC));
        obj = key->ReadObj();
        if(strcmp("TH1F",key->GetClassName())==0){//if ( obj->IsA()->InheritsFrom( TH1::Class() ) ){
            cout << Form("%2i",iC) << "||" << key->GetClassName() << " " << key->GetName() << endl;
            return (TH1F*)obj;
        }
        else{
    	    cout << "Wrong file type (TH1F) TH1F::ReadHist(TString histogram_name, Int_t iC)" << endl;
    	    return 0;
        }
    }

    TH2F* ReadTH2F(TString histname, Int_t iC){
        key = gDirectory->GetKey(Form("%s%i",histname.Data(),iC));
        obj = key->ReadObj();
        if(strcmp("TH2F",key->GetClassName())==0){//if ( obj->IsA()->InheritsFrom( TH2::Class() ) ){
            cout << Form("%2i",iC) << "||" << key->GetClassName() << " " <<key->GetName() << endl;
            return (TH2F*)obj;
        }
        else{
    	    cout << "Wrong file type (TH2F) TH2F::ReadHist(TString histogram_name, Int_t iC)" << endl;
    	    return 0;
        }
    }


    TProfile* ReadTProfile(TString histname, Int_t iC){
        key = gDirectory->GetKey(Form("%s%i",histname.Data(),iC));
        obj = key->ReadObj();
        if(strcmp("TProfile",key->GetClassName())==0){//if ( obj->IsA()->InheritsFrom( TProfile::Class() ) ){
            cout << Form("%2i",iC) << "||" << key->GetClassName() << " " <<key->GetName() << endl;
            return (TProfile*)obj;
        }
        else{
            cout << "Wrong file type (TH2F) TH2F::ReadHist(TString histogram_name, Int_t iC)" << endl;
            return 0;
        }
    }

    TH1F* ReadTH1F(TString histname, Int_t iC, Int_t iY){
        key = gDirectory->GetKey(Form("%s%i%i",histname.Data(),iC,iY));
        obj = key->ReadObj();
        if(strcmp("TH1F",key->GetClassName())==0){//if ( obj->IsA()->InheritsFrom( TH1::Class() ) ){
            cout << Form("%2i|%2i",iC,iY) << "||" << key->GetClassName() << " " <<key->GetName() << endl;
            return (TH1F*)obj;
        }
        else{
    	    cout << "Wrong file type (TH1F) TH1F::ReadHist(TString histogram_name, Int_t iC, Int_t iY)" << endl;
    	    return 0;
        }
    }

    TH2F* ReadTH2F(TString histname, Int_t iC, Int_t iY){
        key = gDirectory->GetKey(Form("%s%i%i",histname.Data(),iC,iY));
        obj = key->ReadObj();
        if(strcmp("TH2F",key->GetClassName())==0){//if ( obj->IsA()->InheritsFrom( TH2::Class() ) ){
            cout << Form("%2i|%2i",iC,iY) << "||" << key->GetClassName() << " " <<key->GetName() << endl;
            return (TH2F*)obj;
        }
        else{
    	    cout << "Wrong file type (TH2F) TH2F::ReadHist(TString histogram_name, Int_t iC, Int_t iY)" << endl;
    	    return 0;
        }
    }

    TProfile* ReadTProfile(TString histname, Int_t iC, Int_t iY){
        key = gDirectory->GetKey(Form("%s%i%i",histname.Data(),iC,iY));
        obj = key->ReadObj();
        if(strcmp("TProfile",key->GetClassName())==0){//if ( obj->IsA()->InheritsFrom( TProfile::Class() ) ){
            cout << Form("%2i|%2i",iC,iY) << "||" << key->GetClassName() << " " <<key->GetName() << endl;
            return (TProfile*)obj;
        }
        else{
            cout << "Wrong file type (TH2F) TH2F::ReadHist(TString histogram_name, Int_t iC, Int_t iY)" << endl;
            return 0;
        }
    }

    TH1F* ReadTH1F(TString histname, Int_t iC, Int_t iY, Int_t iP){
        key = gDirectory->GetKey(Form("%s_%i_Yo_%i_Pt_%i",histname.Data(),iC,iY,iP));
        obj = key->ReadObj();
        if(strcmp("TH1F",key->GetClassName())==0){//if ( obj->IsA()->InheritsFrom( TH1::Class() ) ){
            cout << Form("%2i|%2i|%2i",iC,iY,iP) << "||" << key->GetClassName() << " " <<key->GetName() << endl;
            return (TH1F*)obj;
        }
        else{
    	    cout << "Wrong file type (TH1F) TH1F::ReadHist(TString histogram_name, Int_t iC, Int_t iY, Int_t iP)" << endl;
    	    return 0;
        }
    }

    TH2F* ReadTH2F(TString histname, Int_t iC, Int_t iY, Int_t iP){
        key = gDirectory->GetKey(Form("%s_%i_Yo_%i_Pt_%i",histname.Data(),iC,iY,iP));
        obj = key->ReadObj();
        if(strcmp("TH2F",key->GetClassName())==0){//if ( obj->IsA()->InheritsFrom( TH2::Class() ) ){
            cout << Form("%2i|%2i|%2i",iC,iY,iP) << "||" << key->GetClassName() << " " <<key->GetName() << endl;
            return (TH2F*)obj;
        }
        else{
    	    cout << "Wrong file type (TH2F) TH2F::ReadHist(TString histogram_name, Int_t iC, Int_t iY, Int_t iP)" << endl;
    	    return 0;
        }
    }

    TProfile* ReadTProfile(TString histname, Int_t iC, Int_t iY, Int_t iP){
        key = gDirectory->GetKey(Form("%s%i%i%i",histname.Data(),iC,iY,iP));
        obj = key->ReadObj();
        if(strcmp("TProfile",key->GetClassName())==0){//if ( obj->IsA()->InheritsFrom( TProfile::Class() ) ){
            cout << Form("%2i|%2i|%2i",iC,iY,iP) << "||" << key->GetClassName() << " " <<key->GetName() << endl;
            return (TProfile*)obj;
        }
        else{
            cout << "Wrong file type (TH2F) TH2F::ReadHist(TString histogram_name, Int_t iC, Int_t iY, Int_t iP)" << endl;
            return 0;
        }
    }

    TKey* GetHistKey(){
 	    return key;
    }

  private:
 	TObject* obj;
    TKey* key;
 	TFile* currFile;
 	TString currFName;
};