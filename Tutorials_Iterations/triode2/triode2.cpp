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


const double Te = 1.0;
const double Up = 5.0;
const double spc = 0.009;



bool solid1( double x, double y, double z )
{
  return (( x <= 0.003 && x >= 0.001 && (y >= 0.00296 || y <= -0.00296)) ||
	     ( x <= 0.0024 && x >= 0.001 && ((y >= 0.002 && y >= 2.4*x - 0.0028) || (y <= -0.002 && y <= -2.4*x + 0.0028))));
}


bool solid2( double x, double y, double z )
{
    return  ( x <= spc + 0.006 && x >= spc + 0.003 &&
             (y >= 0.0015 || y <= -0.0015) );
}

bool solid3( double x, double y, double z )
{
    return  ( x <= spc + 0.011 && x >= spc + 0.008 &&
             (y >= 0.0018 || y <= -0.0018) );
}


void simu( int *argc, char ***argv )
{
    Geometry geom( MODE_2D, Int3D(251,101,1), Vec3D(0,-0.005,0), 0.0001 );
    Solid *s1 = new FuncSolid( solid1 );
    geom.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    geom.set_solid( 8, s2 );
    Solid *s3 = new FuncSolid( solid3 );
    geom.set_solid( 9, s3 );
    geom.set_boundary( 1, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET, -8.0e3) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 4, Bound(BOUND_NEUMANN,    0.0) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,  Up)  );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, -8.3e3) );
    geom.set_boundary( 9, Bound(BOUND_DIRICHLET, -8.0e3) );
    geom.build_mesh();

    EpotBiCGSTABSolver solver( geom );
    InitialPlasma initp( AXIS_X, 3e-3 );
    solver.set_initial_plasma( Up, &initp );

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
	pdb.add_2d_beam_with_energy( 30000, 330.0, 1.0, 1.0, 
				     5.0, 0.5, 0.5, 
				     0.0, -0.025, 
				     0.0, 0.025 );
	pdb.iterate_trajectories( scharge, efield, bfield );
    }

    
    GeomPlotter geomplotter( geom );
    geomplotter.set_size( 1280, 750 );
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
