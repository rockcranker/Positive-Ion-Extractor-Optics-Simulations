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


const double Te = 1.0;
const double Up = 5.0;
const double spc = 0.009;

const double Vext = -10e3;
const double dVgg = 300.0;

const double J = 330.0;
const int N = 30000;

bool solid1( double x, double y, double z )
{
  return ( (z <= 0.003 && z >= 0.001) &&
	   ((x <= -0.02 || x >= 0.02 || y <= -0.00296 || y >= 0.00296) ||
	   ((z >= 0.001 && z <= 0.0024) && !(abs(x) <= 0.018 || abs(x) <= 5*z + 0.008)) ||
	    ((z >= 0.001 && z <= 0.0024) && !(abs(y) <= 0.002 ||  abs(y) <= 2.4*z - 0.0028))
	    ));
}


bool solid2( double x, double y, double z )
{
    return  ( z <= spc + 0.006 && z >= spc + 0.003 &&
              (y >= 0.0015 || y <= -0.0015 || x >= 0.02 || x <= -0.02)
	      );
}

bool solid3( double x, double y, double z )
{
    return  ( z <= spc + 0.011 && z >= spc + 0.008 &&
              (y >= 0.0018 || y <= -0.0018 || x >= 0.02 || x <= -0.02)
	      );
}


void simu( int *argc, char ***argv )
{
    Geometry geom( MODE_3D, Int3D(201,41,121), Vec3D(-0.025,-0.005,0), 0.00025 );
    Solid *s1 = new FuncSolid( solid1 );
    geom.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    geom.set_solid( 8, s2 );
    Solid *s3 = new FuncSolid( solid3 );
    geom.set_solid( 9, s3 );
    geom.set_boundary( 1, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 2, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 4, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 5, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 6, Bound(BOUND_DIRICHLET,    Vext) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,  Up)  );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, Vext - dVgg) );
    geom.set_boundary( 9, Bound(BOUND_DIRICHLET, Vext) );
    geom.build_mesh();

    
    EpotBiCGSTABSolver solver( geom );
    InitialPlasma initp( AXIS_Z, 3e-3 );
    solver.set_initial_plasma( Up, &initp );

    EpotField epot( geom );
    MeshScalarField scharge( geom );
    MeshVectorField bfield;
    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    ParticleDataBase3D pdb( geom );
    bool pmirror[6] = { false, false, false, false, false, false };
    pdb.set_mirror( pmirror );

    for( size_t i = 0; i < 15; i++ ) {

      if( i == 1 ) {
	double rhoe = pdb.get_rhosum();
	solver.set_pexp_plasma( -rhoe, Te, Up );
      }

      solver.solve( epot, scharge );
      efield.recalculate();

      pdb.clear();
      pdb.add_rectangular_beam_with_energy( N, J, 1.0, 1.0, 
					    5.0, 0.5, 0.5, 
					    Vec3D(0,0,0), Vec3D(1,0,0), Vec3D(0,1,0), 
					    0.025, 0.005 );
      0.025, 0.005 );
    pdb.iterate_trajectories( scharge, efield, bfield );
  }
    
  geom.save( "geom.dat" );
  epot.save( "epot.dat" );
  pdb.save( "pdb.dat" );
    
  GeomPlotter geomplotter( geom );
  geomplotter.set_size( 1280, 1280 );
  geomplotter.set_epot( &epot );
  geomplotter.set_particle_database( &pdb );
  geomplotter.set_view( VIEW_ZX);
  geomplotter.plot_png( "plot_zx.png" );
  geomplotter.set_view( VIEW_ZY);
  geomplotter.plot_png( "plot_zy.png" );

}

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
