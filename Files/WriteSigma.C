////////////////////////////////
//Parfenov Peter 
//29.12.2015
//
//root -l -b -q WriteSigma.C++()
////////////////////////////////

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
//#include "GeomFunct.h"

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
using namespace std;
 
void WriteSigma(){
  
  const Int_t npt       = 11;
  const Int_t nrpdt     = 11;
  Float_t sPim[npt-1][nrpdt-1];
  Float_t sPip[npt-1][nrpdt-1];
  Float_t sKap[npt-1][nrpdt-1];
  Float_t sPro[npt-1][nrpdt-1];
  
  FILE*  sigmafile = fopen("/u/parfenov/anaflow/3SigmaOut.list","rt");
    
  for (Int_t ipt=0; ipt<npt-1;ipt++){
    for (Int_t irpdt=0; irpdt<nrpdt-1;irpdt++){
  
     fscanf(sigmafile,"%f%f%f%f",&sPim[ipt][irpdt],&sPip[ipt][irpdt],&sKap[ipt][irpdt],&sPro[ipt][irpdt]);
     
    }
  }
  
      
  fclose(sigmafile);
  
  FILE* outsigma   = fopen("/u/parfenov/anaflow/outsigma.C","wt");
  fprintf(outsigma,"//const Int_t npt            = %2i;\n",npt  );
  fprintf(outsigma,"//const Int_t nrpdt          = %2i;\n",nrpdt);
  fprintf(outsigma,"\nFloat_t sPim[npt-1][nrpdt-1] = {");
  for (Int_t ipt=0; ipt<npt-1;ipt++){
    if(ipt ==0){ fprintf(outsigma,"{");                                }
    if(ipt > 0){ fprintf(outsigma,"                                {");}
    for (Int_t irpdt=0; irpdt<nrpdt-1;irpdt++){
      if(irpdt>0 && irpdt<nrpdt-1){fprintf(outsigma,",");}
      fprintf(outsigma,"%3.5f",sPim[ipt][irpdt]);
    }
  if(ipt<npt-2){fprintf(outsigma,"},\n");}
  else{         fprintf(outsigma,"}" );}
  }
  fprintf(outsigma,"};\n");
  
  
  fprintf(outsigma,"\nFloat_t sPip[npt-1][nrpdt-1] = {");
  for (Int_t ipt=0; ipt<npt-1;ipt++){
    if(ipt ==0){ fprintf(outsigma,"{");                                }
    if(ipt > 0){ fprintf(outsigma,"                                {");}
    for (Int_t irpdt=0; irpdt<nrpdt-1;irpdt++){
      if(irpdt>0 && irpdt<nrpdt-1){fprintf(outsigma,",");}
      fprintf(outsigma,"%3.5f",sPip[ipt][irpdt]);
    }
  if(ipt<npt-2){fprintf(outsigma,"},\n");}
  else{         fprintf(outsigma,"}" );}
  }
  fprintf(outsigma,"};\n");
  
  fprintf(outsigma,"\nFloat_t sKap[npt-1][nrpdt-1] = {");
  for (Int_t ipt=0; ipt<npt-1;ipt++){
    if(ipt ==0){ fprintf(outsigma,"{");                                }
    if(ipt > 0){ fprintf(outsigma,"                                {");}
    for (Int_t irpdt=0; irpdt<nrpdt-1;irpdt++){
      if(irpdt>0 && irpdt<nrpdt-1){fprintf(outsigma,",");}
      fprintf(outsigma,"%3.5f",sKap[ipt][irpdt]);
    }
  if(ipt<npt-2){fprintf(outsigma,"},\n");}
  else{         fprintf(outsigma,"}" );}
  }
  fprintf(outsigma,"};\n");
  
  fprintf(outsigma,"\nFloat_t sPro[npt-1][nrpdt-1] = {");
  for (Int_t ipt=0; ipt<npt-1;ipt++){
    if(ipt ==0){ fprintf(outsigma,"{");                                }
    if(ipt > 0){ fprintf(outsigma,"                                {");}
    for (Int_t irpdt=0; irpdt<nrpdt-1;irpdt++){
      if(irpdt>0 && irpdt<nrpdt-1){fprintf(outsigma,",");}
      fprintf(outsigma,"%3.5f",sPro[ipt][irpdt]);
    }
  if(ipt<npt-2){fprintf(outsigma,"},\n");}
  else{         fprintf(outsigma,"}" );}
  }
  fprintf(outsigma,"};\n");
  fclose(outsigma);
  
}