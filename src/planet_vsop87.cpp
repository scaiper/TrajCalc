#include "math_defines_hack.h"

#include "planet_vsop87.h"
#include "exceptions.h"

#include "core_functions/ic2par.h"
#include "core_functions/convert_anomalies.h"

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/lexical_cast.hpp>

#include <vector>
#include <map>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <stdexcept>

planet_vsop87::planet_vsop87(const VSOP87_Data& vsop87_data, const std::string& name, double radius, double safe_radius, double mu_self, size_t cache_size):
	vsop87_data(vsop87_data),
	cacheSize(cache_size)
{
	kep_toolbox::array6D orbital_elements;
	kep_toolbox::array3D r0, v0;
	kep_toolbox::epoch ref_epoch(2451545.0, kep_toolbox::epoch::JD);
	get_eph(ref_epoch, r0, v0);
	kep_toolbox::ic2par(r0,v0, ASTRO_MU_SUN, orbital_elements);
	orbital_elements[5] = kep_toolbox::e2m(orbital_elements[5], orbital_elements[1]);
	
	build_planet(ref_epoch, orbital_elements, ASTRO_MU_SUN, mu_self, radius, safe_radius, name);
}

planet_vsop87::planet_vsop87(const planet_vsop87& other) :
	planet(other),
	vsop87_data(other.vsop87_data),
	cacheSize(other.cacheSize)
{
}

void planet_vsop87::get_eph(const kep_toolbox::epoch& when, kep_toolbox::array3D &r, kep_toolbox::array3D &v) const
{
	double t = (when.jd() - 2451545.0) / 365250.0;
	
	if (!cache.get()) cache.reset(new cache_type(cacheSize));
	
	if (auto cachedPtr = cache->find(t))
	{
		r[0] = (*cachedPtr)[0];
		r[1] = (*cachedPtr)[1];
		r[2] = (*cachedPtr)[2];
		v[0] = (*cachedPtr)[3];
		v[1] = (*cachedPtr)[4];
		v[2] = (*cachedPtr)[5];
		return;
	}
	
	r[0] = vsop87_data.calc_Var(0, t);
	r[1] = vsop87_data.calc_Var(1, t);
	r[2] = vsop87_data.calc_Var(2, t);
	v[0] = vsop87_data.calc_dVar(0, t);
	v[1] = vsop87_data.calc_dVar(1, t);
	v[2] = vsop87_data.calc_dVar(2, t);
	
	cache->insert(t, {r[0], r[1], r[2], v[0], v[1], v[2]});
}

kep_toolbox::planet_ptr planet_vsop87::clone() const
{
	return kep_toolbox::planet_ptr(new planet_vsop87(*this));
}
