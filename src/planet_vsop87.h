#pragma once

#include "cache.h"
#include "vsop87.h"

#include <planet.h>

#include <boost/thread/tss.hpp>

#include <array>
#include <memory>

class planet_vsop87 : public kep_toolbox::planet
{
public:
	planet_vsop87(const VSOP87_Data& vsop_data, const std::string& name, double radius, double safe_radius, double mu_self, size_t cache_size = 4096);
	
	kep_toolbox::planet_ptr clone() const;
	
	planet_vsop87(const planet_vsop87& other);
	
	planet_vsop87(planet_vsop87&& other) = default;
	
	/// Computes the planet/system position and velocity w.r.t the Sun
	/**
		* \param[in] when Epoch in which ephemerides are required
		* \param[out] r Planet position at epoch (SI units)
		* \param[out] v Planet velocity at epoch (SI units)
		*/
	void get_eph(const kep_toolbox::epoch& when, kep_toolbox::array3D &r, kep_toolbox::array3D &v) const;

private:
	typedef Cache< double, std::array< double, 6 > > cache_type;

	VSOP87_Data vsop87_data;
	size_t cacheSize;
	mutable boost::thread_specific_ptr< cache_type > cache;
};
