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
#include <sstream>

#include "TROOT.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"

#include <sys/stat.h>  // mkdir
#include <unistd.h>    // access
#include <cerrno>
#include <cstring>


using namespace std;

tuple<double, double, double, double, double> root_hist_2xGauss( vector<double> var,
								 int bins,
								 string name,
								 string label,
								 string root_file,
								 bool to_png,
								 double p1g,
								 double p2g,
								 double p3g,
								 double p4g,
								 double p5g)
{
  double max = *std::max_element(var.begin(), var.end());
  double min = *std::min_element(var.begin(), var.end());

  TFile fout(root_file.c_str(),"UPDATE");
  std::vector<double> w(var.size(),1); // weights vector
  auto h1 = new TH1D(name.c_str(), label.c_str(), bins, min,max);
  h1->FillN(var.size(),var.data(),w.data());
  TF1 *f1 = new TF1("f1","[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[1])/[4])^2)", min, max); //
  f1->SetParameter(0,p1g);
  f1->SetParLimits(0,0,100000);
  f1->SetParameter(1,p2g);
  f1->SetParameter(2,p3g);
  f1->SetParLimits(2,0.001,0.01);
  f1->SetParameter(3,p4g);
  f1->SetParLimits(3,0,100000);
  f1->SetParameter(4,p5g);
  f1->SetParLimits(4,0.001,0.01);
  f1->SetNpx(1000);
  h1->Fit(f1, "Q");

  if (to_png){
    // Plot
    TCanvas *c1 = new TCanvas("c1", "Canvas 1", 2000, 1000);
    h1->Draw();
    c1->SaveAs((name + ".png").c_str());
    fout.WriteObject(h1, name.c_str());
    fout.WriteObject(f1, (name + "_fit").c_str());
    fout.Close();
  }

  return make_tuple(f1->GetParameter(0), f1->GetParameter(1), f1->GetParameter(2), f1->GetParameter(3), f1->GetParameter(4));
}

tuple<double, double, double, double, double> root_hist_Gauss( vector<double> var,
								 int bins,
								 string name,
								 string label,
								 string root_file,
								 bool to_png,
								 double p1g,
								 double p2g,
								 double p3g
								 )
{
  double max = *std::max_element(var.begin(), var.end());
  double min = *std::min_element(var.begin(), var.end());

  TFile fout(root_file.c_str(),"UPDATE");
  std::vector<double> w(var.size(),1); // weights vector
  auto h1 = new TH1D(name.c_str(), label.c_str(), bins, min,max);
  h1->FillN(var.size(),var.data(),w.data());
  TF1 *f1 = new TF1("f1","gaus", min, max); //
  //f1->SetParameter(0,p1g);
  //f1->SetParLimits(0,0,100000);
  //f1->SetParameter(1,p2g);
  //f1->SetParameter(2,p3g);
  //f1->SetParLimits(2,0.0001,0.01);
  f1->SetNpx(1000);
  h1->Fit(f1, "LQ");

  if (to_png){
    // Plot
    TCanvas *c1 = new TCanvas("c1", "Canvas 1", 2000, 1000);
    h1->Draw();
    c1->SaveAs((name + ".png").c_str());
    fout.WriteObject(h1, name.c_str());
    fout.WriteObject(f1, (name + "_fit").c_str());
    fout.Close();
  }

  return make_tuple(f1->GetParameter(0), f1->GetParameter(1), f1->GetParameter(2), 0, 0);
}

double width(double p0, double p1, double p2, double p3, double p4){
  double min = p1 - 3 * (p2 + p4);
  double max = p1 + 3 * (p2 + p4);
  TF1 *f1 = new TF1("f1","[0]*exp(-0.5*((x-[1])/[2])^2) + [3]*exp(-0.5*((x-[1])/[4])^2)", min, max); //
  f1->SetParameter(0,p0);
  f1->SetParameter(1,p1);
  f1->SetParameter(2,p2);
  f1->SetParameter(3,p3);
  f1->SetParameter(4,p4);
  f1->SetNpx(1000);

  double fmax1 = f1->GetMaximum();
  double fmax1_x = f1->GetMaximumX();
  double w_l1 = f1->GetX(fmax1/M_E,min, fmax1_x);
  double w_r1 = f1->GetX(fmax1/M_E, fmax1_x, max);
  return (w_r1-w_l1)/2;
}

void simu( int argc, char **argv )
{
  // Settings
  const bool remove_traj = false; // Remove trajectories (reduces filesize for large simulations )

  // Load input files and reconstruct simulation domain
  std::ifstream is_geom( "geom.dat" );
  Geometry geom( is_geom );
  is_geom.close();

  std::ifstream is_epot( "epot.dat" );
  EpotField epot( is_epot, geom );
  is_epot.close();

  std::ifstream is_pdb( "pdb.dat" );
  ParticleDataBase3D pdb( is_pdb, geom);
  is_pdb.close();
  if(remove_traj){
    pdb.clear_trajectories();
  }

  EpotEfield efield( epot );
  field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				   FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				   FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
  efield.set_extrapolation( efldextrpl );

  MeshScalarField tdens( geom );
  pdb.build_trajectory_density_field( tdens );

  // Plot XZ and XY views
  GeomPlotter geomplotter( geom );
  geomplotter.set_size( 1280, 1280 );
  geomplotter.set_epot( &epot );
  geomplotter.set_particle_database( &pdb );
  geomplotter.set_view( VIEW_ZX);
  geomplotter.plot_png( "plot_zx.png" );
  geomplotter.set_view( VIEW_ZY);
  geomplotter.plot_png( "plot_zy.png" );

  // Write output file containing all particles
  ofstream fileOut( "particles_out.txt" );
  // column headers
  fileOut << setw(12) << "Current"<< " ";
  fileOut << setw(12) << "Mass"<< " ";
  fileOut << setw(12) << "Mass #"<< " ";
  fileOut << setw(12) << "Sim time"<< " ";
  fileOut << setw(12) << "x"<< " ";
  fileOut << setw(12) << "vx"<< " ";
  fileOut << setw(12) << "y"<< " ";
  fileOut << setw(12) << "vy"<< " ";
  fileOut << setw(12) << "z"<< " ";
  fileOut << setw(12) << "vz\n";

  // particle data setup
  size_t npart = 0;
  size_t pmax = pdb.size();
  double beam_I = 0;
  double extracted_I = 0;
  
  vector<Particle3D> pvec;
  vector<double> px;
  vector<double> py;
  vector<double> pz;
  vector<double> pvx;
  vector<double> pvy;
  vector<double> pvz;
  vector<double> pm;
  vector<int> pmn;

  // Read settings from parameter file
  std::ifstream inputSet;
  inputSet.open("parameters.txt");
  std::string dirName;
  std::string line;
  double lSim;
  double a;
  if ( inputSet.is_open() ) {
    inputSet >> dirName;
    getline(inputSet, line);
    getline(inputSet, line);
    lSim = stod(line);
    getline(inputSet, line);
    a = stod(line);
  }
  inputSet.close();

  // mass separation setup
  double m1 = 1.66054e-27;
  double masses[3] = {m1, 2*m1, 3*m1};
  int min_i = 0;
  double min_del = 1;
  double del = 0;

  // New reduced pdb that only stores relevant particles
  ParticleDataBase3D pdbrel( geom );
  bool pmirror[6] = { false, false, false, false, false, false };
  pdbrel.set_mirror( pmirror );
  ParticleDataBase3D pdbdep1( geom );
  pdbdep1.set_mirror( pmirror );
  ParticleDataBase3D pdbdep2( geom );
  pdbdep2.set_mirror( pmirror );

  // Count particles
  for( size_t k = 0; k < pmax; k++ ) {

    Particle3D &pp = pdb.particle( k );

    // count particles that are actually extractred, add to new pdb
    if( pp(PARTICLE_Z) > 0.004){
      pdbrel.add_particle(pp);
      extracted_I += pp.IQ();
    }

    // count particles on first grid
    if(pp(PARTICLE_Z) > 0.003 + a - 0.0005 && pp(PARTICLE_Z) < 0.003 + a + 0.003 + 0.0005){
      pdbdep1.add_particle(pp);
    }

    // count particles on second grid
    if(pp(PARTICLE_Z) > 0.003 + a + 0.005 - 0.0005 && pp(PARTICLE_Z) < 0.003 + a + 0.008 + 0.0005){
      pdbdep2.add_particle(pp);
    }

    // Skip ions not at the end
    if( pp(PARTICLE_Z) < lSim - 0.001){
      continue;
    }
    
    // Extract particle I, m, coordinates
    // 3D has 7 coordinates
    pvec.push_back(pp);
    px.push_back(pp(1));
    pvx.push_back(pp(2));
    py.push_back(pp(3));
    pvy.push_back(pp(4));
    pz.push_back(pp(5));
    pvz.push_back(pp(6));
    pm.push_back(pp.m());

    beam_I += pp.IQ();

    // determine particle mass number based on mass
    min_i = 0;
    min_del = 1;
    for(size_t l = 0; l < 3; l++){
      del = abs(masses[l] - pp.m());
      if (del < min_del){
	      min_del = del;
	      min_i = l;
      }
    }

    pmn.push_back(min_i + 1);

    // write particle data
    fileOut << setw(12) << pp.IQ() << " ";
    fileOut << setw(12) << pp.m() << " ";
    fileOut << setw(12) << min_i + 1 << " ";
    for( size_t j = 0; j < 7; j ++ )
      fileOut << setw(12) << pp(j) << " ";
    fileOut << "\n";

    npart++;
  }
  fileOut.close();
  std::cout << "particles: " << npart << "\n";

  // overwrite original pdb file
  pdbrel.save( "pdb.dat" );
  

  // Divergence setup
  vector<double> pE;
  vector<double> xAng;
  vector<double> yAng;
  vector<double> xA1;
  vector<double> xA2;
  vector<double> xA3;
  vector<double> yA1;
  vector<double> yA2;
  vector<double> yA3;

  size_t N1 = 0;
  size_t N2 = 0;
  size_t N3 = 0;
  
  double E_temp;
  double Ax;
  double Ay;
  
  // Calculate travel angles
  for(size_t i = 0; i < pvec.size(); i++){
    E_temp = 0.5* pm[i] * (pvx[i] * pvx[i] + pvy[i] * pvy[i] + pvz[i] * pvz[i]);
    pE.push_back(E_temp);

    Ax = atan2(pvx[i], pvz[i]);
    Ay = atan2(pvy[i], pvz[i]);
    xAng.push_back(Ax);
    yAng.push_back(Ay);

    // Separate by mass number
    switch (pmn[i]){
    case 1:
      xA1.push_back(Ax);
      yA1.push_back(Ay);
      N1++;
      break;
    case 2:
      xA2.push_back(Ax);
      yA2.push_back(Ay);
      N2++;
      break;
    case 3:
      xA3.push_back(Ax);
      yA3.push_back(Ay);
      N3++;
      break;
    }
  }

  // Create divergence histograms
  size_t nbins = 300;
  TFile foutx("output.root","RECREATE");
  foutx.Close();

  // Plot and measure width of all histograms
  // all particles
  auto x_div = root_hist_2xGauss( xAng, nbins, "hdx", "X Divergence All", "output.root", true, 200, 0, 0.001, 10, 0.001);
  double wx = width(get<0>(x_div), get<1>(x_div), get<2>(x_div), get<3>(x_div), get<4>(x_div));
  auto y_div = root_hist_Gauss( yAng, nbins, "hdy", "Y Divergence All", "output.root", true, 200, 0, 0.001);
  double wy = width(get<0>(y_div), get<1>(y_div), get<2>(y_div), get<3>(y_div), get<4>(y_div));

  // m = 1
  auto x_div_1 = root_hist_2xGauss( xA1, nbins, "hdx1", "X Divergence 1", "output.root", true, 300, 0, 0.001, 50, 0.001);
  double wx1 = width(get<0>(x_div_1), get<1>(x_div_1), get<2>(x_div_1), get<3>(x_div_1), get<4>(x_div_1));
  auto y_div_1 = root_hist_2xGauss( yA1, nbins, "hdy1", "Y Divergence 1", "output.root", true, 50, 0, 0.001, 1, 0.002);
  double wy1 = width(get<0>(y_div_1), get<1>(y_div_1), get<2>(y_div_1), get<3>(y_div_1), get<4>(y_div_1));

  // m = 2
  auto x_div_2 = root_hist_2xGauss( xA2, nbins, "hdx2", "X Divergence 2", "output.root", true, 300, 0, 0.001, 50, 0.001);
  double wx2 = width(get<0>(x_div_2), get<1>(x_div_2), get<2>(x_div_2), get<3>(x_div_2), get<4>(x_div_2));
  auto y_div_2 = root_hist_2xGauss( yA2, nbins, "hdy2", "Y Divergence 2", "output.root", true, 50, 0, 0.001, 1, 0.002);
  double wy2 = width(get<0>(y_div_2), get<1>(y_div_2), get<2>(y_div_2), get<3>(y_div_2), get<4>(y_div_2));

  // m = 3
  auto x_div_3 = root_hist_2xGauss( xA3, nbins, "hdx3", "X Divergence 3", "output.root", true, 300, 0, 0.001, 50, 0.001);
  double wx3 = width(get<0>(x_div_3), get<1>(x_div_3), get<2>(x_div_3), get<3>(x_div_3), get<4>(x_div_3));
  auto y_div_3 = root_hist_2xGauss( yA3, nbins, "hdy3", "Y Divergence 3", "output.root", true, 50, 0, 0.001, 1, 0.002);
  double wy3 = width(get<0>(y_div_3), get<1>(y_div_3), get<2>(y_div_3), get<3>(y_div_3), get<4>(y_div_3));

  // Particle deposition tracking
  double v;
  double gamma;
  double E;
  double IQ;
  double qe = 1.602e-19;
  double c = 299792459;

  // grid 1
  double dep1_x[pdbdep1.size()];
  double dep1_y[pdbdep1.size()];
  double dep1_z[pdbdep1.size()];
  double dep1_I = 0;
  double dep1_P = 0;
  
  // Location, current, and energy of each deposited particle
  for(size_t k = 0; k < pdbdep1.size(); k++){
    Particle3D &pp = pdbdep1.particle( k );
    dep1_x[k] = pp(1);
    dep1_y[k] = pp(3);
    dep1_z[k] = pp(5);
    
    v = sqrt(pp(2) * pp(2) + pp(4) * pp(4) + pp(6) * pp(6));
    gamma = 1 /sqrt(1 - v * v /(c * c));
    E = (gamma - 1) * pp.m() * c * c;
    IQ = pp.IQ();
    dep1_I += IQ;
    dep1_P += (IQ/qe) * E;
  }
  int size = sizeof(dep1_x) / sizeof(dep1_x[0]);
  double dep1_minx = *min_element(dep1_x, dep1_x + size);
  double dep1_maxx = *max_element(dep1_x, dep1_x + size);
  double dep1_miny = *min_element(dep1_y, dep1_y + size);
  double dep1_maxy = *max_element(dep1_y, dep1_y + size);

  double size_x = geom.size(0) * geom.h();
  double size_y = geom.size(1) * geom.h();
  
  // Plot deposition profile
  TGraph *g = new TGraph(pdbdep1.size(), dep1_x, dep1_y);
  TCanvas *c2 = new TCanvas("c2", "Canvas 2", 10000, 5000);
  g->Draw("AP");
  g->GetXaxis()->SetLimits(-0.5*size_x, 0.5*size_x);
  g->SetMinimum(-0.5*size_y);
  g->SetMaximum(0.5*size_y);
  c2->SetFixedAspectRatio();
  c2->SaveAs("dep1.png");

  // grid 2
  double dep2_x[pdbdep2.size()];
  double dep2_y[pdbdep2.size()];
  double dep2_z[pdbdep2.size()];
  double dep2_I = 0;
  double dep2_P = 0;
  
  for(size_t k = 0; k < pdbdep2.size(); k++){
    Particle3D &pp = pdbdep2.particle( k );
    dep2_x[k] = pp(1);
    dep2_y[k] = pp(3);
    dep2_z[k] = pp(5);
    
    v = sqrt(pp(2) * pp(2) + pp(4) * pp(4) + pp(6) * pp(6));
    gamma = 1 /sqrt(1 - v * v /(c * c));
    E = (gamma - 1) * pp.m() * c * c;
    IQ = pp.IQ();
    dep2_I += IQ;
    dep2_P += (IQ/qe) * E;
  }
  size = sizeof(dep2_x) / sizeof(dep2_x[0]);
  double dep2_minx = *min_element(dep2_x, dep2_x + size);
  double dep2_maxx = *max_element(dep2_x, dep2_x + size);
  double dep2_miny = *min_element(dep2_y, dep2_y + size);
  double dep2_maxy = *max_element(dep2_y, dep2_y + size);

  // GG deposition profile
  TGraph *g2 = new TGraph(pdbdep2.size(), dep2_x, dep2_y);
  c2->cd();
  c2->Clear();
  g2->Draw("AP");
  g2->GetXaxis()->SetLimits(-0.5*size_x, 0.5*size_x);
  g2->SetMinimum(-0.5*size_y);
  g2->SetMaximum(0.5*size_y);
  c2->Modified();
  c2->Update();
  c2->SaveAs("dep2.png");

  // Output beam profile
  TGraph *g3 = new TGraph(px.size(), px.data(), py.data());
  c2->cd();
  c2->Clear();
  g3->Draw("AP");
  g3->GetXaxis()->SetLimits(-0.5*size_x, 0.5*size_x);
  g3->SetMinimum(-0.5*size_y);
  g3->SetMaximum(0.5*size_y);
  g3->SetMarkerSize(1);
  c2->Modified();
  c2->Update();
  c2->SaveAs("beam_out.png");
  c2->SaveAs("beam_out.eps", "eps");

  // Write results file
  ofstream file3( "results.txt" );
  file3 << "N beam: " << pdbrel.size() << "\n";
  file3 << "N transmitted: " << npart << "\n";
  file3 << "H1: " << N1 << "\tH2: " << N2 << "\tH3: " << N3 << "\n";
  file3 << "Transmitted Fraction: " << (double)npart/(double)pdbrel.size() << "\n";
  file3 << "X Divergence 1/e: " << wx << "\n";
  file3 << "Y Divergence 1/e: " << wy << "\n";
  file3 << "X Divergence Mean: " << get<1>(x_div) << "\n";
  file3 << "---\n Parameters:\n";
  file3 << "\t A1 \t\t m \t\t sig \t\t A2 \t\t sig2 \n";
  file3 << "X:\t " << get<0>(x_div) << "\t\t " << get<1>(x_div)<< "\t\t " << get<2>(x_div)<< "\t\t " << get<3>(x_div)<< "\t\t " << get<4>(x_div) <<"\n";
  file3 << "Y:\t " << get<0>(y_div) << "\t\t " << get<1>(y_div)<< "\t\t " << get<2>(y_div)<< "\t\t " << get<3>(y_div)<< "\t\t " << get<4>(y_div) <<"\n";
  file3 << "---\n\n Current and Power:\n";
  file3 << "Extracted Current: " << extracted_I << "\n";
  file3 << "Current Deposition on Grid 1: " << dep1_I << "\n";
  file3 << "Current Deposition on Grid 2: " << dep2_I << "\n";
  file3 << "Final Beam Current: " << beam_I << "\n";
  file3 << "Grid 1 Power: " << dep1_P << "\n";
  file3 << "Between x values: " << dep1_minx << "\t" << dep1_maxx <<"\n";
  file3 << "Between y values: " << dep1_miny << "\t" << dep1_maxy <<"\n";
  file3 << "Grid 2 Power: " << dep2_P << "\n";
  file3 << "Between x values: " << dep2_minx << "\t" << dep2_maxx <<"\n";
  file3 << "Between y values: " << dep2_miny << "\t" << dep2_maxy <<"\n";
  file3 << "\n---\n H1:\n";
  file3 << "X Divergence 1/e: " << wx1 << "\n";
  file3 << "Y Divergence 1/e: " << wy1 << "\n";
  file3 << "X Divergence Mean: " << get<1>(x_div_1) << "\n";
  file3 << "---\n H2:\n";
  file3 << "X Divergence 1/e: " << wx2 << "\n";
  file3 << "Y Divergence 1/e: " << wy2 << "\n";
  file3 << "X Divergence Mean: " << get<1>(x_div_2) << "\n";
  file3 << "---\n H3:\n";
  file3 << "X Divergence 1/e: " << wx3 << "\n";
  file3 << "Y Divergence 1/e: " << wy3 << "\n";
  file3 << "X Divergence Mean: " << get<1>(x_div_3) << "\n";

  
  file3.close();

  
  // Create folder and move files
  std::string folderName = dirName; // from parameters.txt file
  // Files to look for
  std::vector<std::string> filesToMove = {"geom.dat", "pdb.dat", "epot.dat", "angles.txt", "parameters.txt", "particles_out.txt",
					  "plot_zx.png", "plot_zy.png", "output.root", "results.txt",
					  "hdx.png", "hdx1.png", "hdx2.png", "hdx3.png",
					  "hdy.png", "hdy1.png", "hdy2.png", "hdy3.png",
					  "dep1.png", "dep2.png", "beam_out.png" 
  };

  bool new_folder = true; // whether to make new folder

  if(new_folder){
    // Create folder
    if (mkdir(folderName.c_str(), 0755) != 0) {
      if (errno != EEXIST) {
        std::cerr << "Failed to create directory: " << strerror(errno) << std::endl;
      }
    }
    std::cout << "Created folder: " << folderName << std::endl;

    folderName.append("/");
    std::string destination = "";

    // Move files
    for (const auto& file : filesToMove) {
      destination = folderName;
      destination.append(file);
      if (std::rename(file.c_str(), destination.c_str()) != 0) {
        std::cerr << "Failed to move file: " << strerror(errno) << std::endl;
      }
    }
  }

  
}

int main( int argc, char **argv )
{
  ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
  ibsimu.set_thread_count( 4 );
  simu( argc, argv );
  return( 0 );
}