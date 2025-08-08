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

#ifdef GTK3
#include "gtkplotter.hpp"
#endif


const double Te = 5.0;
const double Up = 5.0;


bool solid1( double x, double y, double z )
{
  return (( x <= 0.003 && x >= 0.001 && (y >= 0.02 || y <= -0.02)) ||
	     ( x <= 0.0024 && x >= 0.001 && ((y >= 0.018 && y >= 5*x + 0.008) || (y <= -0.018 && y <= -5*x - 0.008))));
}


bool solid2( double x, double y, double z )
{
    return  ( x <= 0.015 && x >= 0.012 &&
             (y >= 0.02 || y <= -0.02) );
}

bool solid3( double x, double y, double z )
{
    return  ( x <= 0.020 && x >= 0.018 &&
             (y >= 0.02 || y <= -0.02) );
}


void simu( int *argc, char ***argv )
{
    Geometry geom( MODE_2D, Int3D(83,101,1), Vec3D(0,-0.025,0), 0.0005 );
    Solid *s1 = new FuncSolid( solid1 );
    geom.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    geom.set_solid( 8, s2 );
    Solid *s3 = new FuncSolid( solid3 );
    geom.set_solid( 9, s3 );
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET,  5) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET,  -10e3) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 4, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,  5)  );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, -10e3) );
    geom.set_boundary( 9, Bound(BOUND_DIRICHLET, -10e3) );
    geom.build_mesh();

    EpotBiCGSTABSolver solver( geom );
    InitialPlasma initp( AXIS_X, 2e-3 );
    solver.set_initial_plasma( 5.0, &initp );

    EpotField epot( geom );
    MeshScalarField scharge( geom );
    MeshVectorField bfield;
    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    ParticleDataBase2D pdb( geom );
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
	pdb.add_2d_beam_with_energy( 10000, 330.0, 1.0, 1.0, 
				     5.0, 0.0, 0.5, 
				     0.0, -0.025, 
				     0.0, 0.025 );
	pdb.iterate_trajectories( scharge, efield, bfield );
    }
    GeomPlotter geomplotter( geom );
    geomplotter.set_size( 750, 750 );
    geomplotter.set_epot( &epot );
    geomplotter.set_particle_database( &pdb );
    geomplotter.plot_png( "plot1.png" );

    
    

#ifdef GTK3
    GTKPlotter plotter( argc, argv );
    plotter.set_geometry( &geom );
    plotter.set_epot( &epot );
    plotter.set_scharge( &scharge );
    plotter.set_particledatabase( &pdb );
    plotter.new_geometry_plot_window();
    plotter.run();
#endif

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
