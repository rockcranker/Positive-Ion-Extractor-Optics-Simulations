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
#include <palette.hpp>



using namespace std;

void simu( int argc, char **argv )
{
    std::ifstream is_geom( argv[1] );
    if( !is_geom.good() )
	throw( Error( ERROR_LOCATION, (string)"couldn\'t open file \'" + argv[1] + "\'" ) );
    Geometry geom( is_geom );
    is_geom.close();
    geom.build_surface();

    std::ifstream is_epot( argv[2] );
    if( !is_epot.good() )
	throw( Error( ERROR_LOCATION, (string)"couldn\'t open file \'" + argv[2] + "\'" ) );
    EpotField epot( is_epot, geom );
    is_epot.close();

    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    MeshScalarField tdens( geom );
    
    GeomPlotter geomplotter( geom );
    geomplotter.set_size( 1920, 1920 );
    geomplotter.set_epot( &epot );
    geomplotter.set_efield( &efield);
    vector<double> epot_lines = {0.0,-40};
    geomplotter.set_eqlines_manual(epot_lines);
    geomplotter.set_fieldgraph_plot(FIELD_EFIELD);

    geomplotter.set_view( VIEW_ZX);
    geomplotter.plot_png( "plot_zx.png" );
    geomplotter.set_view( VIEW_ZY);
    geomplotter.plot_png( "plot_zy.png" );

    FieldGraph* fgp = geomplotter.fieldgraph();
    
    Palette& p0 = fgp->palette();
    Palette p1 = p0;
    p1.clear();
    p1.push_back(Vec3D(0,0,1), 0.0);
    p1.push_back(Vec3D(0,0,0), 0.000001);
    p1.push_back(Vec3D(1,0,0), 0.01);
    p1.push_back(Vec3D(1,.5,0), 0.2);
    p1.push_back(Vec3D(1,1,0), 0.6);
    p1.push_back(Vec3D(1,1,1), 1.0);
    p1.set_steps(0);
    p1.debug_print (std::cout);

    /*
    Palette::Entry e1(Vec3D(0,0,0), 0.0);
    Palette::Entry e2(Vec3D(1,1,1), 1.0);
    vector<Palette::Entry> pev = {e1, e2};
    Palette p1(pev);
    p1.set_steps(0);
    */
    
    fgp->set_palette(p1);

    const double width1 = 2000;
    const double height1 = 1000;
    Frame frame;
    SolidGraph sldgrf(geom);
    sldgrf.set_view( VIEW_ZY, 101);
    frame.set_geometry( width1, height1, 0, 0 );
    frame.enable_colormap_legend( true );
    frame.set_font_size( 20 );

    MeshColormap& g3p = *fgp;
    Graph& gp = g3p;
    frame.add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, &gp );
    frame.add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, &sldgrf );

    cairo_surface_t *surface;
    surface = cairo_image_surface_create( CAIRO_FORMAT_ARGB32, width1, height1 );
    cairo_t *cairo = cairo_create( surface );
    frame.draw( cairo );
    cairo_surface_write_to_png( surface, "testfield.png");
    cairo_destroy( cairo );
    cairo_surface_destroy( surface );
    
}


int main( int argc, char **argv )
{
    if( argc <= 2 ) {
	cerr << "Usage: analysis geom epot pdb <bfield>\n";
	exit( 1 );
    }

    try {
	ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
	ibsimu.set_thread_count( 4 );
	simu( argc, argv );
    } catch( Error e ) {
	e.print_error_message( ibsimu.message( 0 ) );
        exit( 1 );
    }

    return( 0 );
}
