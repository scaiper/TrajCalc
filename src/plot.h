#pragma once

#include "math_defines_hack.h"

#include <planet.h>
#include <astro_constants.h>
#include <lambert_problem.h>

#include <iostream>

void plot_planet(std::ostream& os, const kep_toolbox::planet& pl, const kep_toolbox::epoch& t0, unsigned N = 120, double units = ASTRO_AU);

void plot_lambert(std::ostream& os, const kep_toolbox::lambert_problem& lp, bool draw_dst = false, size_t sol = 0, unsigned N = 120, double units = ASTRO_AU);
