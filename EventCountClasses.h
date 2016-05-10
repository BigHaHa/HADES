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

    if(pcand->isFakeRejected()!=0) return kFALSE;

    if(pcand->isFlagAND(4,
                        Particle::kIsAcceptedHitInnerMDC,
                        Particle::kIsAcceptedHitOuterMDC,
                        Particle::kIsAcceptedHitMETA,
                        Particle::kIsAcceptedRK)
       && pcand->getChi2() < 400
      )        test = kTRUE;

    if(!HParticleTool::isGoodMetaCell(pcand,4,kTRUE)) return kFALSE;
    if(pcand->getBeta()<0) return kFALSE;

    return test;
}

static Bool_t rejectLeptons(HParticleCand* pcand){ return kFALSE; }

TF1 fdEdxVsMomLowLimit("fdEdxVsMomLowLimit","0.1+1200./(x-60)" , 0., 4000.);  //-lower-limit-for-dE/dx-to-select-protons-from-pions-optimized-for-gen5

Int_t getParticleTrkMult()
{
     //-using-Berusz's-approach-
     HCategory* fParticleEvtInfoCat =  (HCategory*)HCategoryManager::getCategory(catParticleEvtInfo,kTRUE,"catParticleEvtInfo");
     if(!fParticleEvtInfoCat) { cout<<"No catParticleEvtInfo in input!"<<endl; return -1;}
     HParticleEvtInfo *event_info = (HParticleEvtInfo*)fParticleEvtInfoCat->getObject(0);
     if(!event_info) {std::cout<<"No Event INFo"<<std::endl;return -1;}

     Int_t num = event_info->getSumSelectedParticleCandMult();
     return num;
}


class HMultCorr{
  public: 
    HMultCorr(){
      ReadOnce();
    }
  
    void ReadOnce(){
      //FMulCal=fopen("filenames_mult_calib_day096_day107.txt","r");
      FMulCal=fopen("filenames_mult_calib_allDays096_125_ver1.txt","r");
      if(FMulCal){}else{ printf("\n\n===Can not open muliplicity calibration file===\n\n"); exit(0); }
      char ch[1024], ch1[1024], ch2[1024];
      TString str;
      Int_t absMinute;
      Float_t tm1, tm2, tm3, tm4, mS1, mS2, mS3, mS4, mS5, mS6;
      //while (fgets (ch, 1024, FMulCal) != NULL){
      //  printf("%s", ch);
      //}
      //fscanf(inFILE1, "%f",&delta);  gaCalibWall[i] = gaCalibWall[i] - delta;
      Bool_t problem=0;
      Int_t i=0;
      TString fileName, charMinute;
      while( (fscanf(FMulCal,"%s %s %i %f %f %f %f %s %f %f %f %f %f %f"                    , ch, ch1, &absMinute, &tm1, &tm2, &tm3, &tm4, ch2, &mS1, &mS2, &mS3, &mS4, &mS5, &mS6 ) != EOF) ){
        if(i<50){
          printf("%s %s %i %6.3f %6.3f %6.3f %6.3f %s %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n", ch, ch1,  absMinute,  tm1,  tm2,  tm3,  tm4, ch2,  mS1,  mS2,  mS3,  mS4,  mS5,  mS6 );
        }

        if(absMinute<=0) problem=1;
        str=ch1; if(!(str.Contains("==ABS=MINUTE==") || str.Contains("==XXX=MINUTE==")) ){ printf("\nCalibration file does not contain \'==ABS=MINUTE==\' separator\n"); problem=1; }
        str=ch2; if(! str.Contains("==MULT=6SECT==")                                    ){ printf("\nCalibration file does not contain \'==MULT=6SECT==\' separator\n"); problem=1; }
        if(problem){ printf("\n\n==Error detected while reading muliplicity calibration file===\n\n"); exit(0); }

        fileName=ch;

        vector<Float_t> vals;
        vals.push_back(tm1) ;
        vals.push_back(tm2) ;
        vals.push_back(tm3) ;
        vals.push_back(tm4) ;
        vals.push_back(mS1) ;
        vals.push_back(mS3) ;
        vals.push_back(mS3) ;
        vals.push_back(mS4) ;
        vals.push_back(mS5) ;
        vals.push_back(mS6) ;
        BeToMultMap[fileName] = vals;

        i++;
      }

      printf("===============/ First 200 lines printed /==============\n");
      printf("==============/ Alltogether %i lines read /=============\n\n", i);

      fclose(FMulCal);
    }

    void printWholeCalibFile(TString filename, map<TString,vector<Float_t> > &my)
    {
        //-Inclusion-from-Jochen Markert---j.markert@gsi.de-----------//
        map<TString,vector<Float_t> >::iterator iter = my.find(filename);
        if (iter == my.end()) {
            // file not in list
        printf("printWholeCalibFile() : File [%s] not found in map!",filename.Data());
        } else {
              cout<< setw(30)<< filename
                  << " "     << iter->second[0]
                  << " "     << iter->second[1]
                  << " "     << iter->second[2]
                  << " "     << iter->second[3]
                  << " "     << iter->second[4]
                  << " "     << iter->second[5]
                  << " "     << iter->second[6]
                  << " "     << iter->second[7]
                  << " "     << iter->second[8]
                  << " "     << iter->second[9]
                  <<endl;
        }
    }

    //vector<Float_t>& getLineValuesAsVectFromCalibFile(TString filename, map<TString,vector<Float_t> > &my)
    //{
    //    //-Inclusion-from-Jochen Markert---j.markert@gsi.de-----------//
    //    map<TString,vector<Float_t> >::iterator iter = my.find(filename);
    //    if (iter == my.end()) {
    //       // file not in list
    //       printf("getLineValuesAsVectFromCalibFile() : File [%s] not found in map!",filename.Data());
    //    } else {
    //       return iter->second;
    //    }
    //}

    vector<Float_t> getLineValuesAsVectFromCalibFile(TString filename)
    {
        //-Inclusion-from-Jochen Markert---j.markert@gsi.de-----------//
        map<TString,vector<Float_t> >::iterator iter = BeToMultMap.find(filename);
        if (iter == BeToMultMap.end()) {
           // file not in list
           printf("getLineValuesAsVectFromCalibFile() : File [%s] not found in map!",filename.Data());
           //-in-this-case-we-build-a-zero-vector-and-return-it--
           vector<Float_t> valZero;
           valZero.push_back(0.0) ;
           valZero.push_back(0.0) ;
           valZero.push_back(0.0) ;
           valZero.push_back(0.0) ;
           valZero.push_back(0.0) ;
           valZero.push_back(0.0) ;
           valZero.push_back(0.0) ;
           valZero.push_back(0.0) ;
           valZero.push_back(0.0) ;
           valZero.push_back(0.0) ;
           return valZero;
           
        } else {
           return iter->second;
        }
    }


    void print(map<TString, vector<Float_t> > &my)
    {
        //-Inclusion-from-Jochen Markert---j.markert@gsi.de-----------//
        Int_t ct=0;
        for(map< TString, vector<Float_t> >::iterator iter = my.begin(); iter != my.end(); ++iter ) {
              ct++;
              cout<< setw(5)      << ct
                  << " "<<setw(30)<< iter->first
                  << " "          << iter->second[0]
                  << " "          << iter->second[1]
                  << " "          << iter->second[2]
                  << " "          << iter->second[3]
                  << " "          << iter->second[4]
                  << " "          << iter->second[5]
                  << " "          << iter->second[6]
                  << " "          << iter->second[7]
                  << " "          << iter->second[8]
                  << " "          << iter->second[9]
                  <<endl;
        }
    }

    
  protected:
    FILE *FMulCal;
    map<TString, vector<Float_t> > BeToMultMap;

    //-how-to-use------------------------------------------------------------------
    ////print(BeToMultMap);
    //vector<Float_t>& result = getLineValuesAsVectFromCalibFile("12096220603",BeToMultMap);
    //cout<<"result "<< result[0] <<" "<< result[1] <<" "<< result[2] <<" "<< result[3] <<endl;
    ////printWholeCalibFile("12096220603",BeToMultMap); 
};

class PidParticle{
private:
  TF1*    fBetaMean;
  Int_t   pidcharge;
  Float_t betaMean;
  
public:
  PidParticle(){
    InitFunc();
  }
  void InitFunc(){
    fBetaMean = new TF1("fBetaMean",fBeta,600,0,1.5); // Need for PID function!
  }
// PID as boolean function of a charge, mass and momentum
 Bool_t fPID(Int_t pid, Float_t mom_input, Float_t beta_input, Int_t charge_input)
 {
  fBetaMean  ->SetParameter(0,HPhysicsConstants::mass(pid));
  pidcharge  = HPhysicsConstants::charge(pid);
  betaMean   = fBetaMean->Eval(mom_input);
  
  if(beta_input<(betaMean+0.05) && (beta_input>(betaMean-0.05)) && charge_input==pidcharge)
  {
    return kTRUE;
  }
  else 
  { 
    return kFALSE;
  }
 }
};
 
class dEdx{
private:
  Int_t ADCth[303];
public:
  Int_t dEdxCut(Int_t cellNum_input){
  
  ADCth[ 0]=661;
  ADCth[ 1]=447;
  ADCth[ 2]=535;
  ADCth[ 3]=482;
  ADCth[ 4]=485;
  ADCth[ 5]=585;
  ADCth[ 6]=573;
  ADCth[ 7]=667;
  ADCth[ 8]=489;
  ADCth[ 9]=736;
  ADCth[10]=749;
  ADCth[11]=485;
  ADCth[12]=435;
  ADCth[13]=585;
  ADCth[14]=464;
  ADCth[15]=548;
  ADCth[16]=433;
  ADCth[17]=504;
  ADCth[18]=516;
  ADCth[19]=724;
  ADCth[20]=780;
  ADCth[21]=667;
  ADCth[22]=464;
  ADCth[23]=535;
  ADCth[24]=0  ;
  ADCth[25]=491;
  ADCth[26]=478;
  ADCth[27]=535;
  ADCth[28]=661;
  ADCth[29]=592;
  ADCth[30]=639;
  ADCth[31]=742;
  ADCth[32]=454;
  ADCth[33]=686;
  ADCth[34]=667;
  ADCth[35]=686;
  ADCth[36]=818;
  ADCth[37]=617;
  ADCth[38]=736;
  ADCth[39]=736;
  ADCth[40]=655;
  ADCth[41]=435;
  ADCth[42]=573;
  ADCth[43]=705;
  ADCth[44]=818;
  ADCth[45]=592;
  ADCth[46]=736;
  ADCth[47]=447;
  ADCth[48]=742;
  ADCth[49]=460;
  ADCth[50]=454;
  ADCth[51]=591;
  ADCth[52]=611;
  ADCth[53]=504;
  ADCth[54]=648;
  ADCth[55]=812;
  ADCth[56]=523;
  ADCth[57]=504;
  ADCth[58]=749;
  ADCth[59]=661;
  ADCth[60]=454;
  ADCth[61]=780;
  ADCth[62]=661;
  ADCth[63]=525;
  ADCth[64]=742;
  ADCth[65]=0  ;
  ADCth[66]=0  ;
  ADCth[67]=761;
  ADCth[68]=705;
  ADCth[69]=617;
  ADCth[70]=567;
  ADCth[71]=592;
  ADCth[72]=648;
  ADCth[73]=761;
  ADCth[74]=542;
  ADCth[75]=434;
  ADCth[76]=554;
  ADCth[77]=0  ;
  ADCth[78]=0  ;
  ADCth[79]=441;
  ADCth[80]=535;
  ADCth[81]=692;
  ADCth[82]=768;
  ADCth[83]=548;
  ADCth[84]=448;
  ADCth[85]=579;
  ADCth[86]=529;
  ADCth[87]=495;
  ADCth[88]=523;
  ADCth[89]=535;
  ADCth[90]=724;
  ADCth[91]=818;
  ADCth[92]=686;
  ADCth[93]=542;
  ADCth[94]=465;
  ADCth[95]=554;
  ADCth[96]=655;
  ADCth[97]=629;
  ADCth[98]=592;
  ADCth[99]=504;
  
  ADCth[100]=818;
  ADCth[101]=0  ;
  ADCth[102]=812;
  ADCth[103]=529;
  ADCth[104]=742;
  ADCth[105]=485;
  ADCth[106]=661;
  ADCth[107]=673;
  ADCth[108]=623;
  ADCth[109]=598;
  ADCth[110]=661;
  ADCth[111]=579;
  ADCth[112]=598;
  ADCth[113]=510;
  ADCth[114]=611;
  ADCth[115]=680;
  ADCth[116]=617;
  ADCth[117]=658;
  ADCth[118]=592;
  ADCth[119]=579;
  ADCth[120]=505;
  ADCth[121]=592;
  ADCth[122]=529;
  ADCth[123]=730;
  ADCth[124]=610;
  ADCth[125]=573;
  ADCth[126]=579;
  ADCth[127]=654;
  ADCth[128]=617;
  ADCth[129]=573;
  ADCth[130]=629;
  ADCth[131]=699;
  ADCth[132]=830;
  ADCth[133]=768;
  ADCth[134]=711;
  ADCth[135]=655;
  ADCth[136]=636;
  ADCth[137]=749;
  ADCth[138]=711;
  ADCth[139]=667;
  ADCth[140]=576;
  ADCth[141]=730;
  ADCth[142]=642;
  ADCth[143]=585;
  ADCth[144]=521;
  ADCth[145]=585;
  ADCth[146]=623;
  ADCth[147]=598;
  ADCth[148]=642;
  ADCth[149]=585;
  ADCth[150]=711;
  ADCth[151]=560;
  ADCth[152]=560;
  ADCth[153]=567;
  ADCth[154]=699;
  ADCth[155]=636;
  ADCth[156]=554;
  ADCth[157]=636;
  ADCth[158]=567;
  ADCth[159]=655;
  ADCth[160]=585;
  ADCth[161]=504;
  ADCth[162]=558;
  ADCth[163]=573;
  ADCth[164]=636;
  ADCth[165]=611;
  ADCth[166]=655;
  ADCth[167]=642;
  ADCth[168]=654;
  ADCth[169]=585;
  ADCth[170]=648;
  ADCth[171]=535;
  ADCth[172]=548;
  ADCth[173]=768;
  ADCth[174]=0  ;
  ADCth[175]=560;
  ADCth[176]=692;
  ADCth[177]=730;
  ADCth[178]=554;
  ADCth[179]=705;
  ADCth[180]=529;
  ADCth[181]=642;
  ADCth[182]=648;
  ADCth[183]=573;
  ADCth[184]=573;
  ADCth[185]=567;
  ADCth[186]=441;
  ADCth[187]=548;
  ADCth[188]=523;
  ADCth[189]=542;
  ADCth[190]=686;
  ADCth[191]=724;
  ADCth[192]=642;
  ADCth[193]=648;
  ADCth[194]=692;
  ADCth[195]=667;
  ADCth[196]=705;
  ADCth[197]=636;
  ADCth[198]=592;
  ADCth[199]=629;

  ADCth[200]=611;
  ADCth[201]=567;
  ADCth[202]=573;
  ADCth[203]=742;
  ADCth[204]=673;
  ADCth[205]=680;
  ADCth[206]=636;
  ADCth[207]=655;
  ADCth[208]=0  ;
  ADCth[209]=0  ;
  ADCth[210]=567;
  ADCth[211]=724;
  ADCth[212]=524;
  ADCth[213]=542;
  ADCth[214]=768;
  ADCth[215]=573;
  ADCth[216]=629;
  ADCth[217]=0  ;
  ADCth[218]=0  ;
  ADCth[219]=0  ;
  ADCth[220]=554;
  ADCth[221]=585;
  ADCth[222]=570;
  ADCth[223]=620;
  ADCth[224]=485;
  ADCth[225]=598;
  ADCth[226]=667;
  ADCth[227]=516;
  ADCth[228]=535;
  ADCth[229]=0  ;
  ADCth[230]=548;
  ADCth[231]=579;
  ADCth[232]=516;
  ADCth[233]=560;
  ADCth[234]=560;
  ADCth[235]=411;
  ADCth[236]=523;
  ADCth[237]=617;
  ADCth[238]=555;
  ADCth[239]=611;
  ADCth[240]=573;
  ADCth[241]=592;
  ADCth[242]=542;
  ADCth[243]=585;
  ADCth[244]=560;
  ADCth[245]=454;
  ADCth[246]=573;
  ADCth[247]=573;
  ADCth[248]=943;
  ADCth[249]=579;
  ADCth[250]=573;
  ADCth[251]=598;
  ADCth[252]=542;
  ADCth[253]=548;
  ADCth[254]=554;
  ADCth[255]=479;
  ADCth[256]=347;
  ADCth[257]=812;
  ADCth[258]=717;
  ADCth[259]=686;
  ADCth[260]=642;
  ADCth[261]=554;
  ADCth[262]=560;
  ADCth[263]=592;
  ADCth[264]=535;
  ADCth[265]=585;
  ADCth[266]=554;
  ADCth[267]=611;
  ADCth[268]=648;
  ADCth[269]=542;
  ADCth[270]=623;
  ADCth[271]=542;
  ADCth[272]=598;
  ADCth[273]=592;
  ADCth[274]=447;
  ADCth[275]=629;
  ADCth[276]=485;
  ADCth[277]=496;
  ADCth[278]=479;
  ADCth[279]=472;
  ADCth[280]=452;
  ADCth[281]=460;
  ADCth[282]=0  ;
  ADCth[283]=466;
  ADCth[284]=412;
  ADCth[285]=434;
  ADCth[286]=447;
  ADCth[287]=611;
  ADCth[288]=535;
  ADCth[289]=648;
  ADCth[290]=472;
  ADCth[291]=554;
  ADCth[292]=0  ;
  ADCth[293]=0  ;
  ADCth[294]=0  ;
  ADCth[295]=0  ;
  ADCth[296]=686;
  ADCth[297]=523;
  ADCth[298]=686;
  ADCth[299]=660;

  ADCth[300]=490;
  ADCth[302]=555;
  
   return ADCth[cellNum_input];
  
  }
};

class HWallFiredCellsVA{
  //-Fired Cells Vector Array---------------------------------------------------------------------//
  //----------------------------------------------------------------------------------------------//
  //- this class was introduced to avoid two loops over events and hence two identical cuts      -//
  //- in each of sample (while in reality it triggers small mistakes, such that cut in one sample-//
  //- is updated, but in the second loop the same update of the cut could be forgotten)          -//
  //- so it is better to enable array where identical selection from the first loop is stored    -//
  //----------------------------------------------------------------------------------------------//
  public:
    HWallFiredCellsVA(){
      Reset();
    }

    void Reset(){
      nCells=0;
    }

    void SetCellVect(TVector2 vc){
      if(nCells<304){
        vCells[nCells] = vc;
        nCells++;
      }else{
        printf("Error from class HWallFiredCellsVA: attempt to ADD to many cells, can't be so many!!\n");
      }
    }

    void SetCellVect(TVector2 vc, Float_t ch){
      if(nCells<304){
        vCells[nCells] = vc;
        eCells[nCells] = ch;
        nCells++;
      }else{
        printf("Error from class HWallFiredCellsVA: attempt to ADD to many cells, can't be so many!!\n");
      }
    }

    TVector2 GetCellVector(Int_t n){
      if(n>=0 && n<nCells){ 
        return vCells[n]; 
      }else{
        TVector2 emptyVect; emptyVect.Set(0.0,0.0);
        printf("Error from class HWallFiredCellsVA: attempt to GET cell vector which was not defined, n=%i > nCells=%i!!\n", n, nCells);
        return emptyVect;
      }
    }

    Float_t GetCellCharge(Int_t n){
      if(n>=0 && n<nCells){ 
        return eCells[n]; 
      }else{
        printf("Error from class HWallFiredCellsVA: attempt to GET cell dE/dx which was not defined, n=%i > nCells=%i!!\n", n, nCells);
        return 0.0;
      }
    }

    void Print(){
      printf("=== HWallFiredCellsVA::Print(//begin//) ===\n");
      for(Int_t i=0; i<nCells; i++){
        printf("==> [%i] ==  dE/dx=%6.3f  ",i, eCells[i]); vCells[i].Print();
      }
      printf("=== HWallFiredCellsVA::Print(// end //) ===\n");
    }

    Int_t GetNumbOfCells(){return nCells;}

    TVector2 Recenter(Int_t cellNum, Float_t Qx, Float_t Qy, Float_t Sx, Float_t Sy){
      v.Set((vCells[cellNum].X()-Qx)/Sx,(vCells[cellNum].Y()-Qy)/Sy);
      return v;
    }

    void SetFlatt(Int_t n, Float_t sin_in, Float_t cos_in){
      avSin[n] = sin_in;
      avCos[n] = cos_in;
    }

    Float_t Flattening(Float_t Psi){
      dPsi = 0;
      Psi  *= TMath::Pi()/180.;
      for (Int_t n=0;n<6;n++){
        dPsi += 2*(-avSin[n]*cos((n+1)*Psi) + avCos[n]*sin((n+1)*Psi))/(n+1);
      }
      return atan2(sin(Psi+dPsi),cos(Psi+dPsi))*57.2957795130823229;
    }

    Float_t GetdPsi(){
      return dPsi;
    }
    
  protected:
    Int_t    nCells;
    TVector2 vCells[304]; //-vector-------//
    Float_t  eCells[304]; //-energy-dE/dx-//
    Float_t avSin[6];
    Float_t avCos[6];
    TVector2 v;
    Float_t dPsi;
};

class FFlow{
  public:
    FFlow(){
      num=0;
    }

    void THDeclare(TString sDescr1, Int_t i1){
      printf("%s %i\n", sDescr1.Data(), i1);
      Char_t chN[2048], chT[2048];
      sprintf(chN, "hPhiPR%s_%i", sDescr1.Data(), i1 );                //(P-EP), 90, -180.0, 180.0);
      sprintf(chT, "hPhiPR%s %i", sDescr1.Data(), i1 ); hA = new TH1F(chN, chT, 360, -180.0, 180.0); hA->Sumw2(); hA->SetMinimum(0);

      sprintf(chN, "hPhiAB%s_%i", sDescr1.Data(), i1 );                   //(A^B), 180, 0.0, 180.0
      sprintf(chT, "hPhiAB%s %i", sDescr1.Data(), i1 ); hB = new TH1F(chN, chT, 180,    0.0, 180.0); hB->Sumw2(); hB->SetMinimum(0);

      sprintf(chN, "hPhiEP%s_%i", sDescr1.Data(), i1 );                  //PhiEP 90,  -180.0, 180.0
      sprintf(chT, "hPhiEP%s %i", sDescr1.Data(), i1 ); hC = new TH1F(chN, chT,  90,  -180.0, 180.0); hC->Sumw2(); hC->SetMinimum(0); //-self-test-eventplane-distribution--

      sprintf(chN, "hMulTr%s_%i", sDescr1.Data(), i1 );
      sprintf(chT, "hMulTr%s %i", sDescr1.Data(), i1 ); hT = new TH1F(chN, chT,  50,    0.0,  50.0);                                 //-self-test-track-multiplicity-------
    }

    void THDeclare(TString sDescr1, Int_t i1, TString sDescr2, Int_t i2){
      printf("%s %i %s %i\n", sDescr1.Data(), i1, sDescr2.Data(), i2);
      Char_t chN[2048], chT[2048];
      sprintf(chN, "hPhiPR%s_%i_%s_%i", sDescr1.Data(), i1, sDescr2.Data(), i2);
      sprintf(chT, "hPhiPR%s %i %s %i", sDescr1.Data(), i1, sDescr2.Data(), i2); hA = new TH1F(chN, chT, 360, -180.0, 180.0); hA->Sumw2(); hA->SetMinimum(0);

      sprintf(chN, "hPhiAB%s_%i_%s_%i", sDescr1.Data(), i1, sDescr2.Data(), i2);
      sprintf(chT, "hPhiAB%s %i %s %i", sDescr1.Data(), i1, sDescr2.Data(), i2); hB = new TH1F(chN, chT, 180,    0.0, 180.0); hB->Sumw2(); hB->SetMinimum(0);

      sprintf(chN, "hPhiEP%s_%i_%s_%i", sDescr1.Data(), i1, sDescr2.Data(), i2);
      sprintf(chT, "hPhiEP%s %i %s %i", sDescr1.Data(), i1, sDescr2.Data(), i2); hC = new TH1F(chN, chT, 90,  -180.0, 180.0); hC->Sumw2(); hC->SetMinimum(0); //-self-test-eventplane-distribution--

      sprintf(chN, "hMulTr%s_%i_%s_%i", sDescr1.Data(), i1, sDescr2.Data(), i2);
      sprintf(chT, "hMulTr%s %i %s %i", sDescr1.Data(), i1, sDescr2.Data(), i2); hT = new TH1F(chN, chT,  50,    0.0,  50.0);                                 //-self-test-track-multiplicity-------
    }

    void THDeclare(TString sDescr1, Int_t i1, TString sDescr2, Int_t i2, TString sDescr3, Int_t i3){
      //printf("%s %i %s %i %s %i\n", sDescr1.Data(), i1, sDescr2.Data(), i2, sDescr3.Data(), i3);
      Char_t chN[2048], chT[2048];
      sprintf(chN, "hPhiPR%s_%i_%s_%i_%s_%i", sDescr1.Data(), i1, sDescr2.Data(), i2, sDescr3.Data(), i3);
      sprintf(chT, "hPhiPR%s %i %s %i %s %i", sDescr1.Data(), i1, sDescr2.Data(), i2, sDescr3.Data(), i3); hA = new TH1F(chN, chT, 360, -180.0, 180.0); hA->Sumw2(); hA->SetMinimum(0);

      sprintf(chN, "hPhiAB%s_%i_%s_%i_%s_%i", sDescr1.Data(), i1, sDescr2.Data(), i2, sDescr3.Data(), i3);
      sprintf(chT, "hPhiAB%s %i %s %i %s %i", sDescr1.Data(), i1, sDescr2.Data(), i2, sDescr3.Data(), i3); hB = new TH1F(chN, chT, 180,    0.0, 180.0); hB->Sumw2(); hB->SetMinimum(0);

      sprintf(chN, "hPhiEP%s_%i_%s_%i_%s_%i", sDescr1.Data(), i1, sDescr2.Data(), i2, sDescr3.Data(), i3);
      sprintf(chT, "hPhiEP%s %i %s %i %s %i", sDescr1.Data(), i1, sDescr2.Data(), i2, sDescr3.Data(), i3); hC = new TH1F(chN, chT, 90,  -180.0, 180.0); hC->Sumw2(); hC->SetMinimum(0); //-self-test-eventplane-distribution--

      sprintf(chN, "hMulTr%s_%i_%s_%i_%s_%i", sDescr1.Data(), i1, sDescr2.Data(), i2, sDescr3.Data(), i3);
      sprintf(chT, "hMulTr%s %i %s %i %s %i", sDescr1.Data(), i1, sDescr2.Data(), i2, sDescr3.Data(), i3); hT = new TH1F(chN, chT,  50,    0.0,  50.0);                                 //-self-test-track-multiplicity-------
    }

    void NewEvt(){
      hT->Fill(num);
      num=0;
    }

    void Fill(Float_t pr, Float_t ab, Float_t ep, Float_t weight=1.0){
      hA->Fill(pr,weight);
      hB->Fill(ab,weight);
      hC->Fill(ep,weight);
      num++;
    }

  private:
    TString str;
    TH1F *hA;
    TH1F *hB;
    TH1F *hC;
    Int_t num;
    TH1F *hT;
};

class FHist{
  public:
    FHist(){
      num=0;
    }

    void THDeclare(TString sDescr1, Int_t i1, TString sDescr2, Int_t i2, TString sDescr3, Int_t i3){
      //printf("%s %i %s %i %s %i\n", sDescr1.Data(), i1, sDescr2.Data(), i2, sDescr3.Data(), i3);
      Char_t chN[2048], chT[2048];
      sprintf(chN, "hRawM2%s_%i_%s_%i_%s_%i", sDescr1.Data(), i1, sDescr2.Data(), i2, sDescr3.Data(), i3);
      sprintf(chT, "hRawM2%s %i %s %i %s %i", sDescr1.Data(), i1, sDescr2.Data(), i2, sDescr3.Data(), i3); hA = new TH1F(chN, chT, 600,-0.1,2.4);

      sprintf(chN, "hProM2%s_%i_%s_%i_%s_%i", sDescr1.Data(), i1, sDescr2.Data(), i2, sDescr3.Data(), i3);
      sprintf(chT, "hProM2%s %i %s %i %s %i", sDescr1.Data(), i1, sDescr2.Data(), i2, sDescr3.Data(), i3); hB = new TH1F(chN, chT, 600,-0.1,2.4);

      sprintf(chN, "hPhiTr%s_%i_%s_%i_%s_%i", sDescr1.Data(), i1, sDescr2.Data(), i2, sDescr3.Data(), i3);
      sprintf(chT, "hPhiTr%s %i %s %i %s %i", sDescr1.Data(), i1, sDescr2.Data(), i2, sDescr3.Data(), i3); hC = new TH1F(chN, chT, 90,  0.0, 360.0); //-self-test-eventplane-distribution--

      sprintf(chN, "hNumTr%s_%i_%s_%i_%s_%i", sDescr1.Data(), i1, sDescr2.Data(), i2, sDescr3.Data(), i3);
      sprintf(chT, "hNumTr%s %i %s %i %s %i", sDescr1.Data(), i1, sDescr2.Data(), i2, sDescr3.Data(), i3); hT = new TH1F(chN, chT,  1,    0.0,  2.0);                                 //-self-test-track-multiplicity-------
    }

    void NewEvt(){
      num=0;
    }

    void FillMass(Float_t m2){
      hA->Fill(m2);
    }

    void Fill(Float_t m2, Float_t phi){
      hB->Fill(m2);
      hC->Fill(phi);
      hT->Fill(1.);
      num++;
    }

  private:
    TString str;
    TH1F *hA;
    TH1F *hB;
    TH1F *hC;
    Int_t num;
    TH1F *hT;
};