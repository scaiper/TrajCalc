#include "plot.h"

#include <core_functions/propagate_lagrangian.h>
#include <core_functions/ic2par.h>

void plot_planet(std::ostream& os, const kep_toolbox::planet& pl, const kep_toolbox::epoch& t0, unsigned N , double units)
{
    const double PI = acos(-1.0);
    
	os << std::setprecision(std::numeric_limits<double>::digits10 + 2);
	
    kep_toolbox::array3D r, v;
    pl.get_eph(t0, r, v);
    kep_toolbox::array6D E;
    kep_toolbox::ic2par(r, v, pl.get_mu_central_body(), E);
    
    double a = E[0];
    double dt = 2 * PI * sqrt(a*a*a/pl.get_mu_central_body()) / N;
    
    for (unsigned i = 0; i <= N; ++i)
    {
        if (i != 0) kep_toolbox::propagate_lagrangian(r, v, dt, pl.get_mu_central_body());
        os << r[0] / units << "\t" << r[1] / units << "\t" << r[2] / units << std::endl;
    }
    os << "e" << std::endl;
}

void plot_lambert(std::ostream& os, const kep_toolbox::lambert_problem& lp, bool draw_dst, size_t sol, unsigned N, double units)
{
	os << std::setprecision(std::numeric_limits<double>::digits10 + 2);

    kep_toolbox::array3D r = lp.get_r1();
    kep_toolbox::array3D v = lp.get_v1()[0];
    double T = lp.get_tof();
    double dt = T / N;
    os << r[0] / units << "\t" << r[1] / units << "\t" << r[2] / units << std::endl;
    os << "e" << std::endl;
    for (unsigned i = 0; i <= N; ++i)
    {
        if (i != 0) kep_toolbox::propagate_lagrangian(r, v, dt, lp.get_mu());
        os << r[0] / units << "\t" << r[1] / units << "\t" << r[2] / units << std::endl;
    }
    os << "e" << std::endl;
    
    if (draw_dst)
    {
        os << r[0] / units << "\t" << r[1] / units << "\t" << r[2] / units << std::endl;
        os << "e" << std::endl;
    }
}
