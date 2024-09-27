#include <iostream>
#include <stdlib.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"



int main( int argc, char* argv[] ) {

  if( argc!= 3 ) {
    std::cout << "USAGE: ./checkPulseShape [rootFileName] [event] " << std::endl;
    exit(1);
  }

  std::string fileName(argv[1]);
  int event  (atoi(argv[2]));

  TFile* file = TFile::Open(fileName.c_str());
  TTree* tree = (TTree*)file->Get("tree");

  std::cout << "-> Opened file " << fileName.c_str() << std::endl;
  std::cout << "-> Will check pulse shape of event: " << event << std::endl;

  TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
  c1->cd();

  int ev;
  int nch;
  float vamp[128];
  float pshape[128][1024];

  tree->SetBranchAddress( "ev" , &ev     );
  tree->SetBranchAddress( "nch"   , &nch    );
  tree->SetBranchAddress( "vamp"   , &vamp    );
  tree->SetBranchAddress( "pshape", &pshape );

  TH1D* h0 = new TH1D("h0", "", 1024, 0., 1024. );
  TH1D* h1 = new TH1D("h1", "", 1024, 0., 1024. );
  TH1D* h2 = new TH1D("h2", "", 1024, 0., 1024. );
  int nentries = tree->GetEntries();

  int which_channel;

  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    tree->GetEntry(iEntry);
    
    if( ev!=event ) continue;
    
    if (ev == event){

      std::cout << vamp[0] << "\t" << vamp[1] << "\t" << vamp[2]<< std::endl;

      if ( (abs(vamp[0]) > abs(vamp[1])) && (abs(vamp[0]) > abs(vamp[2])) ) which_channel = 0;
        
      if ( (abs(vamp[1]) > abs(vamp[0])) && (abs(vamp[1]) > abs(vamp[2])) ) which_channel = 1;
        
      if ( (abs(vamp[2]) > abs(vamp[0])) && (abs(vamp[2]) > abs(vamp[1])) ) which_channel = 2;


    }

    
    for( unsigned i=0; i<1024; ++i ) {
    
      h0->SetBinContent( i+1, pshape[0][i] );
      h1->SetBinContent( i+1, pshape[1][i] );
      h2->SetBinContent( i+1, pshape[2][i] );
    }



  } // for entries

  if ( which_channel == 0 ){

    std::cout << "*** EVENTO SUL CHANNEL 0 ***" << std::endl;

    h0->Draw();
    h1->Draw("SAME");
    h2->Draw("SAME");

  }
  if (which_channel == 1 ) {
    
    std::cout << "*** EVENTO SUL CHANNEL 1 ***" << std::endl;

    h1->Draw();
    h0->Draw("SAME");
    h2->Draw("SAME");
  }
  
  if(which_channel == 2){

    std::cout << "*** EVENTO SUL CHANNEL 2 ***" << std::endl;

    h2->Draw();
    h1->Draw("SAME");
    h0->Draw("SAME");

  }
  
  
  h1->SetLineColor(2);
  h2->SetLineColor(3);

  TLegend* legend = new TLegend(0.8, 0.8, 0.99, 0.99);
  legend->AddEntry(h0, "ch 0");
  legend->AddEntry(h1, "ch 1");
  legend->AddEntry(h2, "ch 2");
  legend->Draw();

  size_t pos = 0;
  std::string prefix;
  if((pos = fileName.find(".")) != std::string::npos) {
    prefix = fileName.substr(0, pos);
  }

  std::string plotsDir(Form("compared_plots/%s", prefix.c_str()));
  system( Form("mkdir -p %s", plotsDir.c_str()) );

  c1->SaveAs(Form("%s/pulseShapes_ev%d.pdf",plotsDir.c_str(),event));

  return 0;

}
