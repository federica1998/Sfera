#include <iostream>
#include <stdlib.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"




int main( int argc, char* argv[] ) {

  if( argc!= 2 ) {
    std::cout << "USAGE: ./measFromPShape [rootFileName]" << std::endl;
    exit(1);
  }

  std::string fileName(argv[1]);
  
  TFile* file = TFile::Open(fileName.c_str());
  TTree* tree = (TTree*)file->Get("tree");

  std::cout << "-> Opened file " << fileName.c_str() << std::endl;
  std::cout << "-> Will check the measurements from the digitizer" << std::endl;

  TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
  c1->cd();
  
  int ev;
  int nch;
  float base[128];
  float vcharge[128];
  float vamp[128];
  float pshape[128][1024];

  tree->SetBranchAddress( "ev" , &ev     );
  tree->SetBranchAddress( "nch"   , &nch );
  tree->SetBranchAddress( "base" , &base );
  tree->SetBranchAddress( "vamp" , &vamp );
  tree->SetBranchAddress( "vcharge" , &vcharge );
  tree->SetBranchAddress( "pshape", &pshape );

  int k = 0;
  int NBASE = 100;
  int NEV = tree->GetEntries();
  
  tree->GetEntry(3);
  std::cout << " Analizing " << nch << " channels and "<< NEV << " events..." << std::endl;
  
  float check_base[nch * NEV], dig_base[nch * NEV];
  float check_vamp[nch * NEV], dig_vamp[nch * NEV];
  float check_vcharge[nch * NEV], dig_vcharge[nch * NEV];

  float sumB, sumC, maxV;

  
  for( unsigned e =0; e < NEV; ++e ) { // sum over events
    
    std::cout << "Analizing the event " << e + 1 <<  " of " << NEV << std::endl;

    k = e ;  
    tree->GetEntry(e);

    for (unsigned ch = 0 ; ch < nch; ch++){
      
      sumB = 0.;
      maxV = 0.;
      sumC = 0.;

      for( unsigned i=0; i < 1024; ++i ){ 

        //baseline
        if ( i <= NBASE){ 
        
         sumB += pshape[ch][i] ;

        }

        //amplitude
        if (abs(pshape[ch][i]) > abs(maxV)){
          maxV = pshape[ch][i];
        }

        //charge
         sumC += pshape[ch][i];

      } //all points (1024)
      
      
      sumB /= NBASE;
      check_base[k*nch + ch] = sumB;
      dig_base[k*nch + ch] = base[ch];
      //std::cout << " my base = " << sumB << "\t dig = " << base[ch] << std:: endl;
      check_vcharge[k*nch + ch] = sumC ;
      //- sumB * 1024;
      dig_vcharge[k*nch + ch] = vcharge[ch];
      //std::cout << " my charge = " <<  check_vcharge[k*nch + ch] << "\t dig = " << vcharge[ch] << std:: endl;
      check_vamp[k*nch + ch] = maxV;
      dig_vamp[k*nch + ch] = vamp[ch];
      //std::cout << " my amp = " <<  check_vamp[k*nch + ch] << "\t dig = " << vamp[ch] << std:: endl;


  }// for all channels
} // for events

  TGraph* Gbase = new TGraph(NEV * nch, dig_base, check_base);
  Gbase->SetTitle("CONFRONTO DELLE BASELINE");
  Gbase->GetXaxis()->SetTitle("base (digitzer) [V]");
  Gbase->GetYaxis()->SetTitle("base [V]");
  TGraph* Gamp = new TGraph(NEV * nch, dig_vamp, check_vamp);
  Gamp->Fit("pol1");
  Gamp->SetTitle("CONFRONTO DELLE AMPIEZZE");
  Gamp->GetXaxis()->SetTitle("amp (digitzer)");
  Gamp->GetYaxis()->SetTitle("amp [V] ");
  TGraph* Gcharge = new TGraph(NEV * nch, dig_vcharge, check_vcharge);
  Gcharge->Fit("pol1");
  Gcharge->SetTitle("CONFRONTO DELLE CHARCHE");
  Gcharge->GetXaxis()->SetTitle("charge (digitzer) [nC]");
  Gcharge->GetYaxis()->SetTitle("charge [?] ");
  
  size_t pos = 0;
  std::string prefix;
  if((pos = fileName.find(".")) != std::string::npos) {
    prefix = fileName.substr(0, pos);
  }

  std::string plotsDir(Form("check_plots/%s", prefix.c_str()));
  system( Form("mkdir -p %s", plotsDir.c_str()) );

  Gbase->Draw("A*");
  c1->SaveAs(Form("%s/baseline.pdf",plotsDir.c_str()));
  Gamp->Draw("A*");
  c1->SaveAs(Form("%s/amplitude.pdf",plotsDir.c_str()));  
  Gcharge->Draw("A*");
  c1->SaveAs(Form("%s/charge.pdf",plotsDir.c_str()));


  return 0;

}
