#include <random>
#include <fstream>
#include <iomanip>
#include <limits>
#include "meshvectorfield.hpp"
#include "dxf_solid.hpp"
#include "mydxffile.hpp"
#include "gtkplotter.hpp"
#include "geomplotter.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "trajectorydiagnostics.hpp"
#include "particledatabase.hpp"
#include "particlediagplotter.hpp"

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

// Intersection of two lines, x and z directions only, y = 0
Vec3D intPointZX(Vec3D p1, double m1, Vec3D p2, double m2){
  double z1 = p1[2];
  double x1 = p1[0];

  double z2 = p2[2];
  double x2 = p2[0];

  double z = (m1 * z1 - m2 * z2 + x2 - x1)/(m1-m2);
  double x = m1 * (z - z1) + x1;

  return Vec3D(x, 0 , z);
}

// Angle between two vectors
double vec2ang(Vec3D v1, Vec3D v2){
  return acos((v1 * v2) /(v1.norm2() * v2.norm2()));
}

// Random uniform variable in bounds
double getRandom(double min, double max) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(min, max);
    return dis(gen);
}

// Random gaussian variable
double getRandomGauss(double mean, double sigma){
  static std::random_device rd;
  static std::mt19937 gen(rd());
  std::normal_distribution<> dis(mean, sigma);
  return dis(gen);
}

// Main simulation code
void simu( int argc, char **argv )
{
  // Read input particle, geometry data
  std::ifstream is_geom( "geom.dat" );
  Geometry geom( is_geom );
  is_geom.close();
  
  std::ifstream is_pdb( "pdb.dat" );
  ParticleDataBase3D pdb( is_pdb, geom);
  is_pdb.close();
  
  size_t npart = 0;
  size_t pmax = pdb.size();
  
  vector<Particle3D> pvec;
  vector<double> px;
  vector<double> py;
  vector<double> pz;
  vector<double> pvx;
  vector<double> pvy;
  vector<double> pvz;
  vector<double> pm;
  vector<double> pI;

  std::ifstream inputSet;
  inputSet.open("parameters.txt");
  std::string dirName;
  std::string line;
  double lSim;
  if ( inputSet.is_open() ) {
    inputSet >> dirName;
    getline(inputSet, line);
    getline(inputSet, line);
    lSim = stod(line);
  }
  inputSet.close();

  // Extract all particles exiting the beam simulation
  for( size_t k = 0; k < pmax; k++ ) {

    Particle3D &pp = pdb.particle( k );

    // Skip ions not at the end
    if( pp(PARTICLE_Z) < lSim - 0.0005 ){
      continue;
    }
    
    pvec.push_back(pp); // macroparticle object
    px.push_back(pp(1)); // x
    pvx.push_back(pp(2)); // v_x
    py.push_back(pp(3)); // y 
    pvy.push_back(pp(4)); // v_y
    pz.push_back(pp(5)); // z
    pvz.push_back(pp(6)); //v_z
    pm.push_back(pp.m());  // mass
    pI.push_back(pp.IQ()); // current
    npart++;
  }

  vector<double> pE;
  vector<double> xAng;
  vector<double> yAng;

  // Interaction rate parameters
  double xs = 5e-22;
  double uptime = 1; //s
  double pressure = 0.03; //pa
  double kb = 1.38e-23;
  double numdens = pressure/(kb * 293.15);
  double XS = xs * numdens;
  
  // Telescope position setup
  double acceptance = 1.7e-3;
  double optic_rad = 5e-3;
  double beam_lower_bound = -0.02;
  double xT = 0.099;
  double xTl = xT - beam_lower_bound;
  double theta = 75 * M_PI/180;
  double del = 0.01;
  double zTl = (xTl)/tan(theta);
  double zT = lSim + del + zTl;
  Vec3D vecT = Vec3D(xT,0,zT);
  double Lmax = sqrt(xTl * xTl + zTl * zTl);

  Vec3D vLOS = Vec3D(xTl, 0, zTl);
  vLOS.normalize();
  Vec3D off1 = Vec3D(0,0,1) - (Vec3D(0,0,1)*vLOS)*vLOS;
  Vec3D off2 = cross(vLOS, off1);
  double umin = cos(acceptance/2);
  double u;
  double moff;
  double aoff;
  Vec3D off;
  Vec3D em;
  Vec3D plane_int;
  Vec3D plane_vec;
  double d_int;

  // simulation bundles parameters
  int n_seg = 40;
  double n_bund = 1;

  // Stop particles
  double BES_width = acceptance * Lmax + 2 * optic_rad;

  // Starting plane
  double xp1 = xT + optic_rad * cos(theta);
  double zp1 = zT - optic_rad * sin(theta);
  double yp1 = 0;
  Vec3D Vp1 = Vec3D(xp1,yp1,zp1);
  double mp1 = tan(theta - acceptance/2);

  // Ending plane
  double xp2 = xT - optic_rad * cos(theta);
  double zp2 = zT + optic_rad * sin(theta);
  double yp2 = 0;
  Vec3D Vp2 = Vec3D(xp2,yp2,zp2);
  double mp2 = tan(theta + acceptance/2);

  double mx;
  double my;
  Vec3D pstart;
  Vec3D pvel;
  Vec3D int1;
  Vec3D int2;
  Vec3D P1;
  Vec3D P2;
  Vec3D track;
  double Ltrack;
  double Lstep;
  double colls;
  double bundle;

  double vmag;

  // Wavelength parameters
  double alpha;
  double beta;
  double ha = 656.297; //nm
  double qe = 1.602e-19;
  double h;
  double fpart;

  // Output detections to text file
  bool text_out = true;
  int events_per_out = 1000;
  ofstream fileBES( "BES_data.txt" );
  if(text_out){
    fileBES << setw(12) << "w";
    fileBES << setw(12) << "alpha";
    fileBES << setw(12) << "beta";
    fileBES << setw(12) << "h";
    fileBES << setw(12) << "m ";
    fileBES << "\n";
  }

  vector<Vec3D> P1el;
  vector<Vec3D> P2el;
  vector<Vec3D> draw_sample;

  vector<double> H;
  vector<double> W;

  Vec3D travel;
  double track_fraction;
  Vec3D seg;
  Vec3D seg2tel;

  unsigned long long detections = 0;
  int particles_stopped = 0;
  int draw_count = 0;

  // Find particle trajectories, simulate interactions
  for (size_t j = 0; j < npart; j++){
    if(j % 100000 == 0){
      cout << "particle " << j << "\n";
      cout << detections << "\n";	
    }

    // find intersection with planes to define simulated track
    pstart = Vec3D(px[j], py[j], pz[j]);
    pvel = Vec3D(pvx[j], pvy[j], pvz[j]);
    vmag = pvel.norm2();
    my = pvel[1]/pvel[2];
    mx = pvel[0]/pvel[2];

    int1 = intPointZX(pstart, mx, Vp1, mp1);
    P1 = Vec3D(int1[0], pstart[1] + my * (int1[2] - pstart[2]), int1[2]);
    int2 = intPointZX(pstart, mx, Vp2, mp2);
    P2 = Vec3D(int2[0], pstart[1] + my * (int2[2] - pstart[2]), int2[2]);
    track = P2 - P1;
    Ltrack = track.norm2();

    // Discard particles based on y and z
    if (abs(P1[1]) >= BES_width){
      continue;
    }
    if (P1[2] <= lSim){
      cout << "Discarded particle based on z < lsim \n";
      continue;
    }
    
    P1el.push_back(P1);
    P2el.push_back(P2);
    particles_stopped++;

    // divide track into segments, emissions from segment into bundle
    fpart = pI[j] / qe;
    Lstep = Ltrack/n_seg;
    colls = Lstep*XS*fpart* uptime;
    colls = colls * (acceptance * acceptance)/16;
    bundle = colls/n_bund;

    // simulate emissions
    for(int k = 0; k < n_seg; k++){
      track_fraction = (k + 0.5)/n_seg;
      travel = track_fraction * track;
      seg = P1 + travel;
      seg2tel = vecT - seg;

      for(int l = 0; l < n_bund; l++){
	// sample direction
	u = getRandom(umin, 1);
	moff = tan(acos(u));
	aoff = getRandom(0, 2*M_PI);
	off = sin(aoff) * moff * off1 + cos(aoff) * moff * off2;
	em = vLOS + off;
	em.normalize();

	// find intersection
	plane_int = seg + (((seg2tel)*vLOS)/(em*vLOS))*em;
	plane_vec = plane_int - vecT;
	d_int = plane_vec.norm2();

	// record detection if in range
	if (d_int <= optic_rad){
	  // calculate wavelength
	  alpha = vec2ang(pvel, em);
	  beta = vmag/(299792458);
	  // ripple effect
	  beta = beta * (1 + getRandomGauss(0, 0.01)); // 1% ripple
	  h = ha * (1 - beta*cos(alpha))/(sqrt(1 - beta * beta));
	  // random offset
	  h += getRandomGauss(0, 0.025); // 25 picometers error
	  
	  W.push_back(bundle);
	  H.push_back(h);

	  // visualization
	  if(detections % 999 == 0){
	    draw_sample.push_back(seg);
	    draw_count++;
	  }

	  // text output
	  if(text_out && (detections % events_per_out == 0)){
	      fileBES << setw(12) << bundle;
	      fileBES << setw(12) << alpha;
	      fileBES << setw(12) << beta;
	      fileBES << setw(12) << h;
	      fileBES << setw(12) << pm[j];
	      fileBES << "\n";
	   }

	  detections++;
	}
      }
    }
    
  }
  
  // Construct unshifted peak
  // Physical constants
  constexpr double k_B = 1.380649e-23;      // Boltzmann constant (J/K)
  constexpr double T_b = 298.0;               // Temperature in K
  constexpr double m_2 = 3.35e-27;            // Mass of H2 molecule in kg
  double sigma = sqrt(k_B * T_b / m_2);

  double unshift_frac = 0.35;
  int n_unshift = (int)detections * unshift_frac;
  
  double vx;
  double vy;
  double vz;
  
  vector<double> HU;
  vector<double> WU;
  
  for(int m = 0; m < n_unshift; m++){
	  vx = getRandomGauss(0, sigma);
	  vy = getRandomGauss(0, sigma);
	  vz = getRandomGauss(0, sigma);
	  vmag = sqrt( vx * vx + vy * vy  + vz * vz);
	  
	  alpha = atan((sqrt(vx * vx + vy * vy))/vz);
	  beta = vmag/(299792458);
	  h = ha * (1 - beta*cos(alpha))/(sqrt(1 - beta * beta));
	  h += getRandomGauss(0, 0.025); // 25 picometers error
	  
	  WU.push_back(W[m]);
	  HU.push_back(h);
  }
  
  fileBES.close();
  cout << "Total particles stopped: " << particles_stopped << "\n";
  cout << "Total events detected: " << detections << "\n";
  cout << "Drawing particles" << draw_count << "\n";
  
  
  // Add to parameter file
  ofstream param("parameters.txt", std::ios::app);
  param << "BES angle: " << theta << "\n";
  param << "Acceptance angle: " << acceptance << "\n";
  param << "Optical diameter: " << 2*optic_rad << "\n";
  param.close();
  
  // Histogram setup
  cout << "Creating hist \n";
  double max_h = *std::max_element(HU.begin(), HU.end());
  double min_h = *std::min_element(H.begin(), H.end());
  size_t nbins = 1000; //  should be set as instrumental bin width

  cout << "Creating file \n";
  TFile fout("results.root","recreate");

  auto h1 = new TH1D("h1", "BES Results", nbins, min_h,max_h);
  h1->FillN(detections,H.data(),W.data());
  h1->FillN(n_unshift, HU.data(), WU.data());

  int b_max = h1->GetMaximumBin();
  TLine *line_ha = new TLine(ha, 0, ha, h1->GetBinContent(b_max));
  line_ha->SetLineColor(kRed);
  line_ha->SetLineWidth(4);
  h1->SetLineWidth(4);

  TCanvas *ch = new TCanvas("ch", "Canvas BES Hist", 10000, 5000);
  h1->SetStats(0);
  h1->Draw("hist");
  line_ha->Draw();
  cout << "Creating graph \n";
  ch->SaveAs("BES_out.png");
  fout.WriteObject(h1, "h1");
  fout.Close();


  // Write histogram data to file
  cout << "Hist data file \n";
  ofstream histdata("histdata.txt");
  double low_edge, high_edge, bin_counts;
  for(size_t bin = 1; bin < nbins; bin++){
    low_edge = h1->GetBinLowEdge(bin);
    high_edge = h1->GetBinLowEdge(bin+1);
    bin_counts = h1->GetBinContent(bin);
    histdata << low_edge << "\t" << high_edge << "\t" << bin_counts << "\n";
  }
  histdata.close();
  
  
  // visualize stopped particle trajectories, detections, in z-x and z-y
  bool BES_sample_viz = false;
  if (BES_sample_viz){
	  cout << "visualizing sample \n";
    int nBES = P1el.size();
    double xarr[2*nBES];
    double yarr[2*nBES];
    double zarr[2*nBES];
    double collx[draw_count];
    double colly[draw_count];
    double collz[draw_count];
    for(int i = 0; i < nBES; i++){
      xarr[i] = P1el[i][0];
      yarr[i] = P1el[i][1];
      zarr[i] = P1el[i][2];
    }
    for(int i = 0; i < nBES; i++){
      xarr[nBES+i] = P2el[i][0];
      yarr[nBES+i] = P2el[i][1];
      zarr[nBES+i] = P2el[i][2];
    }
    for(int i = 0; i < draw_count; i++){
      collx[i] = draw_sample[i][0];
      colly[i] = draw_sample[i][1];
      collz[i] = draw_sample[i][2];
    }
    TGraph *g = new TGraph(2*nBES, zarr, xarr);
    TCanvas *c1 = new TCanvas("c1", "z-x", 800, 400);
    c1->SetFixedAspectRatio();
    g->Draw("AP");
    c1->SaveAs("z_x.png");
    TGraph *g2 = new TGraph(nBES, zarr, yarr);
    TCanvas *c2 = new TCanvas("c2", "z-y", 800, 400);
    c2->SetFixedAspectRatio();
    g2->Draw("AP");
    c2->SaveAs("z_y.png");
    TGraph *g3 = new TGraph(draw_count, collz, collx);
    TCanvas *c3 = new TCanvas("c3", "detections zx", 1920, 1080);
    c3->SetFixedAspectRatio();
    g3->Draw("AP");
    c3->SaveAs("detections zx.png");
    TGraph *g4 = new TGraph(draw_count, collz, colly);
    TCanvas *c4 = new TCanvas("c4", "detections zy", 1920, 1080);
    c4->SetFixedAspectRatio();
    g4->Draw("AP");
    c4->SaveAs("detections zy.png");
  }

  
}


int main( int argc, char **argv )
{
  ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
  ibsimu.set_thread_count( 4 );
  simu( argc, argv );
  return( 0 );
}
