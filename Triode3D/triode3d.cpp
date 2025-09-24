#include "epot_bicgstabsolver.hpp"
#include "particledatabase.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "meshvectorfield.hpp"
#include "ibsimu.hpp"
#include "error.hpp"
#include "particlediagplotter.hpp"
#include "geomplotter.hpp"
#include "config.h"
#include "gtkplotter.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <sstream>

using namespace std;

// Plasma parameters
const double Te = 5.; // Electron temperature
const double Ti = 0.1; // Ion temperature, can be eparated into tranvere and longitudinal if needed
const double Up = 5.0;  // Plasma potential

// Extractor geometry, in meters
const double a = 0.009; // Spacing PG to EG
const double b = 0.002; // Spacing EG to GG
const double t1 = 0.002; // PG thickness
const double t2 = 0.003; // EG and GG thickness

const double wSlit = 0.036; // Slit width
const double hSlit = 0.004; // Slit height

// Extractor configuration
const double Vext = -10e3; // Total extraction voltage
const double dVgg = 300.0; // Diffference between EG and GG voltage
const double J = 330.0; // Extracted current

// B field components, uniform, Tesla
double Bx = 0.0;
double By = 0.0;
double Bz = 0.0;

// Space charge neutralization past extractor
const bool SC_neutralization = true;

// SImulation domain setup
const double lAcc = 0.001 + t1 + a + t2 + b + t2 + 0.001;
const double lSim = lAcc; // simulation length, z
const double wSim = 0.048; // Simulation width, x
const double hSim = 0.014; // Simulation height, y

const double lMesh = 0.0004; // Cubic mesh size, m
const int ppc = 45; // Particles per mesh cell in beam
const double H1_frac = (double)1/3; // Fraction of H1+
const double H2_frac = (double)1/3; // Fraction of H1+
const double H3_frac = (double)1/3; // Fraction of H1+

// Beam size, can be smaller than Sim size, but must be larger than slit
const double wBeam = 0.048;
const double hBeam = 0.014;

// Cell size and beam calculations
const int nCellW = (int)(wSim/lMesh);
const int nCellH = (int)(hSim/lMesh);
const int nCellL = (int)(lSim/lMesh);
const int NBeam = (int)(wBeam/lMesh) * (int)(hBeam/lMesh) * ppc; // Expected particles in beam
const int N1 = (int)NBeam*H1_frac;
const int N2 = (int)NBeam*H2_frac;
const int N3 = (int)NBeam*H3_frac;
const int N = N1 + N2 + N3; // Actual particle number
const int NSlit = (int)(wSlit/lMesh) * (int)(hSlit/lMesh) * ppc; // Expected particles exracted from slit


// Simulation iteration cycles
const int N_cycles = 10;

// PG geometry w/ Pierce angle 
bool solid1( double x, double y, double z )
{
  return ( (z <= 0.001 + t1 && z >= 0.001) &&
	   ((x <= -0.02 || x >= 0.02 || y <= -0.00296 || y >= 0.00296) ||
	   ((z >= 0.001 && z <= 0.0024) && !(std::abs(x) <= 0.018 || std::abs(x) <= 5*z + 0.008)) ||
	    ((z >= 0.001 && z <= 0.0024) && !(std::abs(y) <= 0.002 ||  abs(y) <= 2.4*z - 0.0028))
	    ));
}

// EG geometry with rounded slit
bool solid2( double x, double y, double z )
{
  double w1 = 0.04;
  double h1 = 0.003;
  double r1 = (x - (w1-h1)/2)*(x - (w1-h1)/2) + y * y;
  double r2 = (x + (w1-h1)/2)*(x + (w1-h1)/2) + y * y;

  return  ( z <= a + t1 + t2 + 0.001 && z >= a + t1 + 0.001 &&
	    (y >= h1/2 || y <= -h1/2 || x >= (w1-h1)/2 || x <= -(w1-h1)/2) &&
	    r1 >= (h1/2)*(h1/2) &&   r2 >= (h1/2)*(h1/2)
	    );
}

// GG geometry with rounded slit
bool solid3( double x, double y, double z )
{
  double w1 = 0.04;
  double h1 = 0.0036;
  double r1 = (x - (w1-h1)/2)*(x - (w1-h1)/2) + y * y;
  double r2 = (x + (w1-h1)/2)*(x + (w1-h1)/2) + y * y;
  return  ( z <= a + t1 + t2 + t2 + b + 0.001 && z >= a + 0.001 + t1 + t2 + b &&
	    (y >= h1/2 || y <= -h1/2 || x >= (w1-h1)/2 || x <= -(w1-h1)/2) &&
	    r1 >= (h1/2)*(h1/2) &&   r2 >= (h1/2)*(h1/2)
	    );
}

void simu( int *argc, char ***argv )
{
  // Geometry set-up
  Geometry geom( MODE_3D, Int3D(nCellW + 1, nCellH + 1, nCellL + 1), Vec3D(-wSim/2,-hSim/2,0), lMesh); // Simulation domain
  // Grids
  Solid *s1 = new FuncSolid( solid1 );
  geom.set_solid( 7, s1 );
  Solid *s2 = new FuncSolid( solid2 );
  geom.set_solid( 8, s2 );
  Solid *s3 = new FuncSolid( solid3 );
  geom.set_solid( 9, s3 );
  // BC
  geom.set_boundary( 1, Bound(BOUND_NEUMANN,    0.0) ); // negative x
  geom.set_boundary( 2, Bound(BOUND_NEUMANN,    0.0) ); // positive x
  geom.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) ); // negative y
  geom.set_boundary( 4, Bound(BOUND_NEUMANN,    0.0) ); // positive y
  geom.set_boundary( 5, Bound(BOUND_NEUMANN,    0.0) ); // negative z
  geom.set_boundary( 6, Bound(BOUND_DIRICHLET, Vext) ); // positive z
  geom.set_boundary( 7, Bound(BOUND_DIRICHLET,  Up)  ); // PG
  geom.set_boundary( 8, Bound(BOUND_DIRICHLET, Vext - dVgg) ); // EG
  geom.set_boundary( 9, Bound(BOUND_DIRICHLET, Vext) ); // GG
  geom.build_mesh();
    
  EpotBiCGSTABSolver solver( geom );
  
  // Inital plasma region
  InitialPlasma initp( AXIS_Z, 3e-3 );
  solver.set_initial_plasma( Up, &initp );

  // E field
  EpotField epot( geom );
  MeshScalarField scharge( geom );
  EpotEfield efield( epot );
  field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
            FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
            FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
  efield.set_extrapolation( efldextrpl );

  // B field, uniform
  bool fout[3] = { true, true, true };
  MeshVectorField bfield( geom, fout );
  for( uint32_t i = 0; i < bfield.size(0); i++ ) {
    for( uint32_t j = 0; j < bfield.size(1); j++ ) {
      for ( uint32_t k = 0; k < bfield.size(2); k++ ) {
        bfield.set( i, j, k, Vec3D( Bx, By, Bz ) );
      }
    }
  }

  // Particle setup
  ParticleDataBase3D pdb( geom );
  bool pmirror[6] = { false, false, false, false, false, false };
  pdb.set_mirror( pmirror );

  // Iterative solver
  for( size_t i = 0; i < N_cycles; i++ ) {
    if( i == 1 ) {
        double rhoe = pdb.get_rhosum();
        solver.set_pexp_plasma( -rhoe, Te, Up );
    }


    // Space Charge Suppression past accelerator
    if (SC_neutralization){
      for( uint32_t k = 0; k < scharge.size(0); k++ ) {
        for( uint32_t l = 0; l < scharge.size(1); l++ ) {
          for ( uint32_t m = 0; m < scharge.size(2); m++ ) {
            if ((m+1) * lMesh > lAcc){
              scharge( k, l, m) = 0;
            }
          }
        }
      }
    }

    solver.solve( epot, scharge );
    efield.recalculate();

    // Particle beam
    pdb.clear();
    // H1+
    pdb.add_rectangular_beam_with_energy( N1, J * H1_frac, 1.0, 1.0, 
                  5.0, Ti, Ti, 
                  Vec3D(0,0,0), Vec3D(1,0,0), Vec3D(0,1,0), 
                  wBeam/2, hBeam/2);
    // H2+
    pdb.add_rectangular_beam_with_energy( N2, J * H2_frac, 1.0, 2.0, 
                  5.0, Ti, Ti, // Initial energy and temperature
                  Vec3D(0,0,0), Vec3D(1,0,0), Vec3D(0,1,0), 
                  wBeam/2, hBeam/2);
    // H3+
    pdb.add_rectangular_beam_with_energy( N3, J * H3_frac, 1.0, 3.0, 
                  5.0, Ti, Ti, 
                  Vec3D(0,0,0), Vec3D(1,0,0), Vec3D(0,1,0), 
                  wBeam/2, hBeam/2 );
    pdb.iterate_trajectories( scharge, efield, bfield );

    // Optional plotting of intermediate results
    /*
    if (i == 8){
      GeomPlotter geomplotter1( geom );
      geomplotter1.set_size( 1280, 1280 );
      geomplotter1.set_epot( &epot );
      geomplotter1.set_particle_database( &pdb );
      geomplotter1.set_view( VIEW_ZY);
      geomplotter1.plot_png( "plot_zy_intermediate.png" );
    }
    */
  }
  
  // Output files
  geom.save( "geom.dat" );
  epot.save( "epot.dat" );
  pdb.save( "pdb.dat" );

  // Write output file containing parameters
  
  std::ostringstream oss; // Create an output string stream object
  
  std::ofstream fileOut( "parameters.txt" );
  
  // convert spacing to string
  std::string astr;
  oss << (a*1000);    
  astr =  oss.str();
  astr.erase(std::remove(astr.begin(), astr.end(), '.'), astr.end());
  astr.erase(astr.find_last_not_of("0") + 1, std::string::npos);
  oss.str("");   // Clear the content
  oss.clear();   // Clear the state flags
  
  // convert temperature to string
  std::string Tstr;
  oss<<(Ti);
  Tstr = oss.str();
  std::replace(Tstr.begin(), Tstr.end(), '.','-');
  Tstr.erase(Tstr.find_last_not_of("0") + 1, std::string::npos);

  // Write parameters to file
  fileOut << "J" << J << "a" << astr << "N" << NSlit << "T" << Tstr << "\n"; // Name string with relevant parameters
  fileOut << lSim << "\n"; // Simulation length, read by postp and BES_sim programs
  // any relevant parameters output for future reference
  fileOut << a << "\n";
  fileOut << "a: " << a << ",   ";
  fileOut << "b: " << b << ",   ";
  fileOut << "t1: " << t1 << ",   ";
  fileOut << "t2: " << t2 << ",   ";
  fileOut << "\n";
  fileOut << "Voltage: " << Vext << " ";
  fileOut << "\n";
  fileOut << "Current: " << J << " ";
  fileOut << "\n";
  fileOut << "N simulated: " << NBeam << "\n";
  fileOut << "N extracted estimate: " << NSlit << "\n";
  fileOut << "Mesh size: " << lMesh << "\n";
  fileOut << "Particles per cell: " << ppc << "\n";
  fileOut << "Particle ratios: \nH1: " << H1_frac << "\t\tH2: " << H2_frac << "\t\tH3: " << H2_frac << "\n";
  fileOut.close();
  
  // Trajectory plots
  GeomPlotter geomplotter( geom );
  geomplotter.set_size( 1280, 1280 );
  geomplotter.set_epot( &epot );
  geomplotter.set_particle_database( &pdb );
  geomplotter.set_view( VIEW_ZX);
  geomplotter.plot_png( "plot_zx.png" );
  geomplotter.set_view( VIEW_ZY);
  geomplotter.plot_png( "plot_zy.png" );

  // Particle profile cross section in center of grids
  //geomplotter.set_view( VIEW_XY, (int)(0.0135/lMesh));
  //geomplotter.plot_png( "grid1xs.png" );
  //geomplotter.set_view( VIEW_XY, (int)(0.0185/lMesh));
  //geomplotter.plot_png( "grid2xs.png" );
  
  // Interactive plotter
  if( false ) { // change to true to enable
      MeshScalarField tdens( geom );
      pdb.build_trajectory_density_field( tdens );
      GTKPlotter plotter( argc, argv );
      plotter.set_geometry( &geom );
      plotter.set_epot( &epot );
      plotter.set_bfield( &bfield );
      plotter.set_efield( &efield );
      plotter.set_scharge( &scharge );
      plotter.set_trajdens( &tdens );
      plotter.set_particledatabase( &pdb );
      plotter.new_geometry_plot_window();
      plotter.new_particle_plot_window( AXIS_Z, 0.03, PARTICLE_DIAG_PLOT_SCATTER, DIAG_VZ, DIAG_VX);
      plotter.run();	
  }
  
}

// main simulation code
int main( int argc, char **argv )
{
  try {
      ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
      ibsimu.set_thread_count( 4 );
      simu( &argc, &argv );
  } catch ( Error e ) {
      e.print_error_message( ibsimu.message( 0 ) );
      exit( 1 );
  }

  return( 0 );
}
