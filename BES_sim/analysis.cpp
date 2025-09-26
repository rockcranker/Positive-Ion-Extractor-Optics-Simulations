#include <random>
#include <fstream>
#include <iomanip>
#include <limits>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <cmath>
#include "TROOT.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TLine.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TF1.h"
#define _USE_MATH_DEFINES
#include "math.h"
#include <TGraph.h>
#include <TVector.h>

using namespace std;

int main( int argc, char **argv )
{
  TFile *file = TFile::Open("results.root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Error opening file!" << std::endl;
    return 1;
  }

  // Retrieve the histogram
  TH1D *h1 = dynamic_cast<TH1D*>(file->Get("h1"));  
  if (!h1) {
    std::cerr << "Histogram not found!" << std::endl;
    file->Close();
    return 1;
  }

  size_t nbins = h1->GetNbinsX();
  double xmin = h1->GetBinLowEdge(1);
  double xmax = h1->GetBinLowEdge(nbins+1);

  double theta = 75 * M_PI / 180;
  double c = 299792458;
  double ha = 656.297;

  double E = 10e3; // ev
  
  double m1 = 938.27e6; // ev
  double m2 = 2 * m1;
  double m3 = 3 * m1;

  double b1 = sqrt(1 - (m1/(E+m1)) * (m1/(E+m1)));
  double b2 = sqrt(1 - (m2/(E+m2)) * (m2/(E+m2)));
  double b3 = sqrt(1 - (m3/(E+m3)) * (m3/(E+m3)));

  double h1c = ha * (1 - b1 * cos(theta))/(sqrt(1 - b1 * b1));
  double h2c = ha * (1 - b2 * cos(theta))/(sqrt(1 - b2 * b2));
  double h3c = ha * (1 - b3 * cos(theta))/(sqrt(1 - b3 * b3));
  double edge_offset = 0.2; //nm
  double e12 = (h1c + h2c)/2;
  double e23 = (h2c + h3c)/2;
  double e1 = h1c - edge_offset;
  double e3 = h3c + edge_offset;
  double overlap_offset_1 = 0.03;
  double overlap_offset_2 = 0.01;
  

    
  TF1 *f1 = new TF1("f1","[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[1])/[4])^2)",e1,e12 - overlap_offset_1);
  TF1 *f2 = new TF1("f2","[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[1])/[4])^2)",e12 + overlap_offset_1, e23 - overlap_offset_2);
  TF1 *f3 = new TF1("f3","[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[1])/[4])^2)", e23 + overlap_offset_2, e3);
  f1->SetParameter(0,200);
  f1->SetParLimits(0,0,500000);
  f1->SetParameter(1,h1c);
  f1->SetParLimits(1,400,800);
  f1->SetParameter(2,0.02);
  f1->SetParLimits(2,0,1);
  f1->SetParameter(3,100);
  f1->SetParLimits(3,0,500000);
  f1->SetParameter(4,0.04);
  f1->SetParLimits(4,0,1);
  
  f2->SetParameter(0,200);
  f2->SetParLimits(0,0,500000);
  f2->SetParameter(1,h2c);
  f2->SetParLimits(1,400,800);
  f2->SetParameter(2,0.02);
  f2->SetParLimits(2,0.001,1);
  f2->SetParameter(3,100);
  f2->SetParLimits(3,0,500000);
  f2->SetParameter(4,0.03);
  f2->SetParLimits(4,0,1);
  
  f3->SetParameter(0,200);
  f3->SetParLimits(0,0,500000);
  f3->SetParameter(1,h3c);
  f3->SetParLimits(1,400,800);
  f3->SetParameter(2,0.02);
  f3->SetParLimits(2,0,1);
  f3->SetParameter(3,100);
  f3->SetParLimits(3,0,500000);
  f3->SetParameter(4,0.03);
  f3->SetParLimits(4,0,1);
  
  f1->SetNpx(1000);
  h1->Fit(f1, "LR");
  
  f2->SetNpx(1000);
  h1->Fit(f2, "LR+");
  
  f3->SetNpx(1000);
  h1->Fit(f3, "LR+");

  

  double std1 = f1->GetParameter(2);
  double w1c = sqrt(2) * std1;
  cout<<"1_e width core 1 = " << w1c <<endl;
  double mid1 = f1->GetMaximumX();
  double max1 = f1->GetMaximum();
  double xl1 = f1->GetX(max1/M_E,e1, mid1);
  double xr1 = f1->GetX(max1/M_E, mid1, e12);
  double w1 = (xr1 - xl1)/2;
  cout<<"1_e width 1 = " << w1 <<endl;

  double std2 = f2->GetParameter(2);
  double w2c = sqrt(2) * std2;
  cout<<"1_e width  core 2 = " << w2c <<endl;
  double mid2 = f2->GetMaximumX();
  double max2 = f2->GetMaximum();
  double xl2 = f2->GetX(max2/M_E,e12, mid2);
  double xr2 = f2->GetX(max2/M_E, mid2, e23);
  double w2 = (xr2 - xl2)/2;
  cout<<"1_e width 2 = " << w2 <<endl;

  double std3 = f3->GetParameter(2);
  double w3c = sqrt(2) * std3;
  cout<<"1_e width core 3 = " << w3c <<endl;
  double mid3 = f3->GetMaximumX();
  double max3 = f3->GetMaximum();
  double xl3 = f3->GetX(max3/M_E, e23, mid3);
  double xr3 = f3->GetX(max3/M_E, mid3, e3);
  double w3 = (xr3 - xl3)/2;
  cout<<"1_e width 3 = " << w3 <<endl;

  TCanvas *c1 = new TCanvas("c1", "BES Results", 2000, 1000);
  h1->SetLineWidth(2);
  h1->GetXaxis()->CenterTitle(true);
  h1->GetYaxis()->CenterTitle(true);
  h1->SetTitle("BES results; Wavelength (nm); Counts");
  h1->Draw("HIST 9");
  f1->Draw("Same");
  f2->Draw("Same");
  f3->Draw("Same");
  c1->SaveAs("h_fit.png");

  ofstream resultfile( "hist_results.txt" );
  resultfile << "1_e width core 1 = " << w1c <<endl;
  resultfile << "1_e width 1 = " << w1 <<endl;
  resultfile << "1_e width core 2 = " << w2c <<endl;
  resultfile << "1_e width 2 = " << w2 <<endl;
  resultfile << "1_e width core 3 = " << w3c <<endl;
  resultfile << "1_e width 3 = " << w3 <<endl;

  resultfile << "\nPeak 1 parameters:\n";
  resultfile << "A1 = \t" << f1->GetParameter(0) << "\t" << f1->GetParError(0) << "\n";
  resultfile << "m = \t" << f1->GetParameter(1) << "\t" << f1->GetParError(1) << "\n";
  resultfile << "std1 = \t" << f1->GetParameter(2) << "\t" << f1->GetParError(2) << "\n";
  resultfile << "A2 = \t" << f1->GetParameter(3) << "\t" << f1->GetParError(3) << "\n";
  resultfile << "std2 = \t" << f1->GetParameter(4) << "\t" << f1->GetParError(4) << "\n";

  resultfile << "\nPeak 2 parameters:\n";
  resultfile << "A1 = \t" << f2->GetParameter(0) << "\t" << f2->GetParError(0) << "\n";
  resultfile << "m = \t" << f2->GetParameter(1) << "\t" << f2->GetParError(1) << "\n";
  resultfile << "std1 = \t" << f2->GetParameter(2) << "\t" << f2->GetParError(2) << "\n";
  resultfile << "A2 = \t" << f2->GetParameter(3) << "\t" << f2->GetParError(3) << "\n";
  resultfile << "std2 = \t" << f2->GetParameter(4) << "\t" << f2->GetParError(4) << "\n";

  resultfile << "\nPeak 3 parameters:\n";
  resultfile << "A1 = \t" << f3->GetParameter(0) << "\t" << f3->GetParError(0) << "\n";
  resultfile << "m = \t" << f3->GetParameter(1) << "\t" << f3->GetParError(1) << "\n";
  resultfile << "std1 = \t" << f3->GetParameter(2) << "\t" << f3->GetParError(2) << "\n";
  resultfile << "A2 = \t" << f3->GetParameter(3) << "\t" << f3->GetParError(3) << "\n";
  resultfile << "std2 = \t" << f3->GetParameter(4) << "\t" << f3->GetParError(4) << "\n";

  resultfile.close();
  
  return( 0 );
}
