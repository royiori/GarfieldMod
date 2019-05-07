//------------------------------------------------------------------//
// Author: Qian LIU
//------------------------------------------------------------------//
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include <TSystem.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoPcon.h>
#include <TGeoHalfSpace.h>
#include <TGeoMatrix.h>
#include <TGeoCompositeShape.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>

using namespace std;

char Path[100] = "./simdata";
TString RootFile = TString(Path)+TString("/result.root");

int main() 
{
    TFile *f = TFile::Open(RootFile,"READ");
    if (!f) { return -1; }

    TTree *t; 
    f->GetObject("tree",t);

    int ne, ni;
    vector<double> *e1hit = 0;
    vector<double> *e2hit = 0;
    vector<double> *i1hit = 0;
    vector<double> *i2hit = 0;

    TBranch *be1hit = 0;
    TBranch *be2hit = 0;
    TBranch *bi1hit = 0;
    TBranch *bi2hit = 0;

    t->SetBranchAddress("e1hit",&e1hit, &be1hit);
    t->SetBranchAddress("e2hit",&e2hit, &be2hit);
    t->SetBranchAddress("i1hit",&i1hit, &bi1hit);
    t->SetBranchAddress("i2hit",&i2hit, &bi2hit);
    t->SetBranchAddress("ne", &ne);
    t->SetBranchAddress("ni", &ni);

    for (Int_t i = 0; i < t->GetEntries(); i++) {
      t->GetEntry(i);
      cout<<ne<<" "<<ni<<endl;

      for (UInt_t j = 0; j < e1hit->size(); j+=5) {
          //cout<<e1hit->at(j)<<endl;
      }
    }

    t->ResetBranchAddresses();

    return 0;
}
