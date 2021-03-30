#include <iostream>
#include <stdlib.h>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TAxis.h"


float purify (std::vector<float> vec){
  
  //std::cout << "... purify ..." << std::endl;

  int N = vec.size();
  float mean = 0., sd = 0.;
  float N_S = 5.;
  
  for (unsigned int i = 0; i < N ; i++ ){
    mean += vec[i];
    sd += vec[i] * vec[i];
  }

  mean /= N;
  sd =  sd/N - mean * mean;

  //std::cout << " mean = " << mean << "\t sd = " << sd << std::endl;

  float new_mean = mean * N;
  int Nereased = 0;
  
  for (unsigned int j = 0; j < N ; j++ ){

    if ( abs(vec[j] - mean )/sd > N_S){
      new_mean -= vec[j];
      Nereased++;
      //std::cout << "removed : " << vec[j] << std::endl;
    }
  }
  
  new_mean /= (N - Nereased);
  //std::cout << " new_mean = " << new_mean << std::endl;
  
  return new_mean ;
}



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

  int NBASE = 70;
  int NEV = tree->GetEntries();
  
  tree->GetEntry(0);
  std::cout << " Analizing " << nch << " channels and "<< NEV << " events..." << std::endl;

  std::vector<float> vcheck_base;
  float check_vamp, check_vcharge, check_base;
  float K = 0.0202904;
  int NBINS = 100;
  

  TH1F* hBase = new TH1F("baseline", "", NBINS ,-0.004 , 0.004);
  TH1F* hVamp = new TH1F("ampiezza", "", NBINS ,-0.01 , 0.01);
  TH1F* hCharge = new TH1F("carica", "", NBINS ,-0.5, 0.5);

  for( unsigned e =0; e < NEV; ++e ) { // sum over events
    
    if ( (e+1 % 100) == 0) std::cout << "Analizing the event " << e + 1 <<  " of " << NEV << std::endl;
    //if (e > 1 ) break;
   
    tree->GetEntry(e);

    for (unsigned ch = 0 ; ch < nch; ch++){
      
      if (ch == 11 ) continue;

      check_vamp = 0.;
      check_vcharge = 0.;

      for( unsigned i=0; i < 1024; ++i ){ 

        //baseline
        if ( i <= NBASE){ 
          vcheck_base.push_back(pshape[ch][i]);
        }

        //amplitude
        if (abs(pshape[ch][i]) > abs(check_vamp)){
          check_vamp = pshape[ch][i];
        }

        //charge
        check_vcharge += pshape[ch][i];

      } //all points (1024)
      
      
      check_base = purify(vcheck_base);
      hBase->Fill(check_base - base[ch]);
      vcheck_base.clear();
      //std::cout << " my base = " << sumB << "\t dig = " << base[ch] << std:: endl;
      
      check_vcharge -=  base[ch] * 1024;
      hCharge->Fill( (check_vcharge / K - vcharge[ch])/vcharge[ch]);
      //std::cout << " my charge = " <<  check_vcharge[k*nch + ch] << "\t dig = " << vcharge[ch] << std:: endl;
     
      check_vamp -= check_base;
      hVamp->Fill(check_vamp - vamp[ch]);
      //std::cout << " my amp = " <<  check_vamp[k*nch + ch] << "\t dig = " << vamp[ch] << std:: endl;
    

    }// for all channels
  } // for events

  TAxis* xaxis = hBase->GetXaxis();
  TAxis* yaxis = hBase->GetYaxis();
  float dim = 0.02;

  //c1->SetLogy();

  xaxis->SetTitle("baseline ([calcolata] - [digitizer])");
  yaxis->SetTitle("counts");
  xaxis->SetLabelSize(dim);
  yaxis->SetLabelSize(dim);
  //hBase->Fit("gaus");

  xaxis = hVamp->GetXaxis();
  yaxis = hVamp->GetYaxis();
  xaxis->SetTitle("amplitude ( [calcolata] - [digitizer] )");
  yaxis->SetTitle("counts");
  xaxis->SetLabelSize(dim);
  yaxis->SetLabelSize(dim);
  //hVamp->Fit("gaus");
  
  xaxis = hCharge->GetXaxis();
  yaxis = hCharge->GetYaxis();
  xaxis->SetTitle("carica ([calcolata] - [digitizer])/ [digitizer]");
  yaxis->SetTitle("counts");
  xaxis->SetLabelSize(dim);
  yaxis->SetLabelSize(dim);
  std::cout << "K : " << hCharge->GetMean() + 1. << std::endl;
  //hCharge->Fit("gaus");
  

  size_t pos = 0;
  std::string prefix;
  if((pos = fileName.find(".")) != std::string::npos) {
    prefix = fileName.substr(0, pos);
  }

  std::string plotsDir(Form("check_plots/%s", prefix.c_str()));
  system( Form("mkdir -p %s", plotsDir.c_str()) );

  
  hBase->Draw();
  c1->SaveAs(Form("%s/baseline_diff.pdf",plotsDir.c_str()));
  hVamp->Draw();
  c1->SaveAs(Form("%s/amplitudes_diff.pdf",plotsDir.c_str()));
  hCharge->Draw();
  c1->SaveAs(Form("%s/charge_diff.pdf",plotsDir.c_str()));


  return 0;

}
