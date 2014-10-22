#include "mga.h"
#include "fb_calc.h"
#include "plot.h"

#include <core_functions/array3D_operations.h>

#include <boost/range/algorithm.hpp>

#include <map>

using namespace kep_toolbox;

void mul(array3D& out, const array3D& v1, double a)
{
    out[0] = a * v1[0];
    out[1] = a * v1[1];
    out[2] = a * v1[2];
}

MGA::MGA(std::vector< planet_ptr > _seq,
        const epoch& t_begin, const epoch& t_end,
        double tof_min, double tof_max,
        bool zeroFB, bool optimizeToF):
    base(
		_seq.size(),
		0,
		optimizeToF ? 2 : 1,
		4,
		4,
		{0.1, 0.0, 0.0, 0.0}
	),
    seq(std::move(_seq)),
    zeroFB(zeroFB),
    optimizeToF(optimizeToF),
    arrivalC(ArrivalConstraint::IGNORE),
    common_mu(seq.front()->get_mu_central_body()),
	maxToF(tof_max)
{
    for (auto&& pl: seq)
    {
        if (pl->get_mu_central_body() != common_mu) {
            //throw
        }
    }
    
    pagmo::decision_vector lb, ub;
    lb.push_back(t_begin.mjd());
    ub.push_back(t_end.mjd());
    lb.insert(lb.end(), seq.size() - 1, 1.0e-5);
    ub.insert(ub.end(), seq.size() - 1, tof_max);
	
    set_bounds(lb, ub);
}

void MGA::setArrivalIgnore()
{
    arrivalC = ArrivalConstraint::IGNORE;
	insertR = boost::none;
}

void MGA::setArrivalInsert(boost::optional<double> R)
{
    arrivalC = ArrivalConstraint::INSERT_LO;
	insertR = R;
}

void MGA::setArrivalAerobrake(double max_Venc)
{
    arrivalC = ArrivalConstraint::AEROBRAKE;
    ab_maxV = max_Venc;
	insertR = boost::none;
}

void MGA::setMaxDV(boost::optional<double> maxDV)
{
	this->maxDV = maxDV;
}

void MGA::objfun_impl(pagmo::fitness_vector& f, const pagmo::decision_vector& x) const
{
    do_calc(x, &f);
}

void MGA::calc_objfun_constr(pagmo::fitness_vector& f, pagmo::constraint_vector& c, const pagmo::decision_vector& x) const
{
	do_calc(x, &f, nullptr, nullptr, &c);
}

pagmo::problem::base_ptr MGA::clone() const
{
    return pagmo::problem::base_ptr(new MGA(*this));
}

std::string MGA::human_readable_extra() const
{
    std::string s;
    for (auto&& pl: seq)
    {
        s += s.empty() ? "" : ", ";
        s += pl->get_name();
    }
    return "\n\tSequence: " + s;
}

void MGA::pretty(std::ostream& os, const pagmo::decision_vector& x) const
{
    do_calc(x, nullptr, &os);
}

void MGA::plot(std::ostream& os, const pagmo::decision_vector& x) const
{
    do_calc(x, nullptr, nullptr, &os);
}

void MGA::compute_constraints_impl(pagmo::constraint_vector& c, const pagmo::decision_vector& x) const
{
	std::fill(c.begin(), c.end(), 0.0);
	do_calc(x, nullptr, nullptr, nullptr, &c);
}

void MGA::do_calc(const pagmo::decision_vector& x, pagmo::fitness_vector* f, std::ostream* os, std::ostream* plot, pagmo::constraint_vector* c) const
{
    std::vector<double> T(seq.size() - 1);

    for (size_t i = 0; i < T.size(); ++i)
    {
        T[i] = x[1 + i] / ASTRO_DAY2YEAR;
    }
    
    std::vector< kep_toolbox::epoch > t_P;
    std::vector< kep_toolbox::array3D > r_P, v_P;
    std::vector< double > DV, RC, VC;
    
    epoch departureEpoch(x[0], epoch::type::MJD);
    
    double sumT = 0.0;
    for (size_t i = 0; i < seq.size(); ++i)
    {
        if (i != 0) sumT += T[i - 1];
        t_P.emplace_back(departureEpoch.mjd() + sumT, epoch::type::MJD);
        array3D r, v;
        seq[i]->get_eph(t_P.back(), r, v);
        r_P.push_back(std::move(r));
        v_P.push_back(std::move(v));
    }
    
    if (plot)
    {
        std::map< std::string , planet_ptr > planets;
        for (auto&& pl: seq)
        {
            planets[pl->get_name()] = pl;
        }
        std::vector< std::pair< double , planet_ptr > > sorted_planets;
        for (auto&& pl: planets)
        {
            array3D r, v;
            pl.second->get_eph(departureEpoch, r, v);
            sorted_planets.emplace_back(norm(r), pl.second);
        }
        boost::sort(sorted_planets);
        *plot << "set view 0, 0" << std::endl;
        *plot << "set view equal xyz" << std::endl;
        *plot << "set grid" << std::endl;
        *plot << "set xyplane at 0" << std::endl;
        *plot << "unset border" << std::endl;
        *plot << "set key Left" << std::endl;
        *plot << "splot ";
        *plot << "'-' with points notitle lt rgb 'yellow' pt 7 ps 4, ";
        for (auto&& pl: sorted_planets)
        {
            *plot << "'-' with lines notitle lt -1, ";
        }
        std::vector< std::string > lc{"red", "blue", "dark-magenta", "orange", "dark-green", "dark-cyan"};
        std::vector<int> pt{5, 7, 9, 11, 13};
        for (size_t i = 0; i < T.size(); ++i)
        {
            *plot << "'-' with points title '" << seq[i]->get_name() << " " << t_P[i].get_posix_time().date() << "'" << " lt rgb '" << lc[i % lc.size()] << "' pt " << pt[i % pt.size()] << " ps 2, ";
            *plot << "'-' with lines title '" << seq[i]->get_name() << " -> " << seq[i+1]->get_name() << " " << int(T[i] + 0.5) << " days'" << " lt rgb '" << lc[i % lc.size()] << "' lw 2, ";
        }
        *plot << "'-' with points title '" << seq.back()->get_name() << " " << t_P.back().get_posix_time().date() << "'" << " lt rgb '" << lc[T.size() % lc.size()] << "' pt " << pt[T.size() % pt.size()] << " ps 2" << std::endl;
        *plot << "0\t0\t0" << std::endl << "e" << std::endl;
        for (auto&& pl: sorted_planets)
        {
            plot_planet(*plot, *pl.second, departureEpoch);
        }
    }
    
    lambert_problem lp(r_P[0], r_P[1], T[0] * ASTRO_DAY2SEC, common_mu, 0, 0);
    if (plot) plot_lambert(*plot, lp, T.size() == 1);
    
    array3D v_beg = lp.get_v1()[0];
    array3D v_end = lp.get_v2()[0];
    array3D v_rel;
    
    diff(v_rel, v_beg, v_P[0]);
    double Vinf = norm(v_rel);
    
    VC.push_back( (norm(v_P[0]) / Vinf - 100.0) / 100.0 );
    
    DV.push_back( sqrt(2 * seq[0]->get_mu_self() / seq[0]->get_radius() + pow(Vinf, 2)) - sqrt(seq[0]->get_mu_self() / seq[0]->get_radius()) );
    
    if (os)
    {
		*os << std::fixed << std::setprecision(3);
        *os << "Voyage:    ";
        for (size_t i = 0; i < seq.size(); ++i)
        {
            *os << (i == 0 ? "" : " -> ") << seq[i]->get_name();
        }
        *os << std::endl;
        *os << std::endl;
        *os << "leg no. 1: " << seq[0]->get_name() << " to " << seq[1]->get_name() << std::endl;
        *os << "Departure: " << t_P[0].get_posix_time().date() << "  ( " << t_P[0].mjd() << " MJD )" << std::endl;
        *os << "Duration: " << std::setprecision(0) << T[0] << std::setprecision(3) << " days ( " << T[0] * ASTRO_DAY2YEAR << " years )" << std::endl;
        *os << "Vinf: " << Vinf / 1000 << " km/s" << std::endl;
        *os << "Eject DV: " << DV[0] / 1000 << " km/s" << std::endl;
        array3D vr_out;
        diff(vr_out, v_beg, v_P[0]);
        array3D ch_plane;
        cross(ch_plane, r_P[0], v_P[0]);
        array3D outward;
        cross(outward, v_P[0], ch_plane);
        *os << "TransX plan:" << std::endl;
		*os << std::setprecision(1);
        *os << "\tPrograde vel.  " << dot(vr_out, v_P[0]) / norm(v_P[0]) << " m/s" << std::endl;
        *os << "\tOutward vel.   " << dot(vr_out, outward) / norm(outward) << " m/s" << std::endl;
        *os << "\tCh. plane vel. " << -dot(vr_out, ch_plane) / norm(ch_plane) << " m/s" << std::endl;
		*os << std::setprecision(3);
    }
    
    for (size_t i = 1; i < T.size(); ++i)
    {
        array3D vr_in;
        diff(vr_in, v_end, v_P[i]);
        lambert_problem lp(r_P[i], r_P[i + 1], T[i] * ASTRO_DAY2SEC, common_mu, 0, 0);
        if (plot) plot_lambert(*plot, lp, i == T.size() - 1);
        v_beg = lp.get_v1()[0];
        v_end = lp.get_v2()[0];
        array3D vr_out;
        diff(vr_out, v_beg, v_P[i]);
        
        double curDV, curR;
        fb_calc(curDV, curR, vr_in, vr_out, seq[i]->get_mu_self());
        DV.push_back(curDV);
        RC.push_back( seq[i]->get_safe_radius() / curR - 1.0 );
        VC.push_back( (norm(v_P[i]) / norm(vr_out) - 100.0) / 100.0 );
        
        if (os)
        {
            *os << "\nleg no." << i+1 << ":      " << seq[i]->get_name() << " to " << seq[i+1]->get_name() << std::endl;
            *os << "Fly-by date:   " << t_P[i].get_posix_time().date() << "  ( " << t_P[i].mjd() << " MJD )" << std::endl;
            *os << "Duration:      " << std::setprecision(0) << T[i] << std::setprecision(3) << " days  (" << T[i] * ASTRO_DAY2YEAR << " years )" << std::endl;
            *os << "Fly-by radius: " << curR / seq[i]->get_radius() << " planetary radius ( " << curR / seq[i]->get_safe_radius() << " safe radius )" << std::endl;
			*os << std::setprecision(1);
            *os << "Fly-by DV:     " << DV[i] << " m/s" << std::endl;
			*os << std::setprecision(3);
            array3D ch_plane;
            cross(ch_plane, r_P[i], v_P[i]);
            array3D outward;
            cross(outward, v_P[i], ch_plane);
            array3D vr_proj;
            mul(vr_proj, ch_plane, dot(vr_out, ch_plane) / norm(ch_plane) / norm(ch_plane));
            array3D vr_pl;
            diff(vr_pl, vr_out, vr_proj);
            *os << "TransX plan:" << std::endl;
            array3D vrXvP;
            cross(vrXvP, vr_pl, v_P[i]);
            *os << "\tOutward angle " << acos(dot(vr_pl, v_P[i]) / norm(vr_pl) / norm(v_P[i])) * 180 / M_PI * (dot(vrXvP, ch_plane) > 0.0 ? 1.0 : -1.0) << std::endl;
            *os << "\tInc. angle    " << -asin(dot(vr_out, ch_plane) / norm(vr_out) / norm(ch_plane)) * 180 / M_PI << std::endl;
            *os << "\tVelocity      " << norm(vr_out) / 1000 << " km/s" << (DV[i] < 1.0 ? " (inherited)" : "") << std::endl;
        }
    }
    
    array3D Venc_rel;
    diff(Venc_rel, v_end, v_P.back());
    double Venc = norm(Venc_rel);
	double insertionRadius = insertR ? *insertR : seq.back()->get_safe_radius();
    double Vins = sqrt(2 * seq.back()->get_mu_self() / insertionRadius + Venc * Venc) - sqrt(2 * seq.back()->get_mu_self() / insertionRadius);
	double Vcirc = sqrt(2 * seq.back()->get_mu_self() / insertionRadius + Venc * Venc) - sqrt(seq.back()->get_mu_self() / insertionRadius);
	double Vsurf = sqrt(2 * seq.back()->get_mu_self() / seq.back()->get_radius() + Venc * Venc);
    
    if (os)
    {
        *os << "\nArrival at " << seq.back()->get_name() << std::endl;
        *os << "Arrival epoch:       " << t_P.back().get_posix_time().date() << "  ( " << t_P.back().mjd() << " MJD ) " << std::endl; 
        *os << "Arrival Vinf:        " << Venc / 1000 << " km/s" << std::endl;
        
        *os << "Arrival surface vel: " << Vsurf / 1000 << " km/s" << std::endl;
		*os << "Reference radius:    " << insertionRadius / 1000 << " km" << std::endl;
        *os << "Insertion DV:        " << Vins / 1000 << " km/s" << std::endl;
		*os << "Circularization DV:  " << Vcirc / 1000 << " km/s" << std::endl;
    }
    
    double DVcost = std::accumulate(DV.begin(), DV.end(), 0.0);    
	double ToF = std::accumulate(T.begin(), T.end(), 0.0) * ASTRO_DAY2YEAR;
	
    if (arrivalC == ArrivalConstraint::INSERT_LO)
    {
        DVcost += Vins;
    }
	else if (arrivalC == ArrivalConstraint::AEROBRAKE)
	{
		DVcost += std::max(Vsurf - ab_maxV, 0.0);
	}
	
    
    if (os)
    {
        *os << std::endl;
        *os << "Total mission time:  " << std::setprecision(0) << ToF / ASTRO_DAY2YEAR << std::setprecision(3) << " days ( " << ToF << " years )" << std::endl;
        *os << "Total DV cost: " << DVcost / 1000 << " km/s" << std::endl;
    }
    
	auto constrF = [](double s, double a){
		return a > 0.0 ? s + a : s;
	};
	
    DVcost *= exp( std::accumulate(VC.begin(), VC.end(), 0.0, constrF) );
	
    if (f)
    {
		f->assign(get_f_dimension(), 0.0);
		
        (*f)[0] = DVcost;
		if (optimizeToF) {
			(*f)[1] = ToF;
		}
    }
	
	if (c)
	{
    
		c->assign(get_c_dimension(), 0.0);

		if (zeroFB) (*c)[0] = std::accumulate(DV.begin() + 1, DV.end(), 0.0, constrF);
		
		(*c)[1] = std::accumulate(RC.begin(), RC.end(), 0.0, constrF);
		
		if (maxDV) (*c)[2] = DVcost / *maxDV - 1.0;
		
		(*c)[3] = ToF / maxToF - 1.0;
	}
}
