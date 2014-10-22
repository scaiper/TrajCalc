#pragma once

#include "mga.h"

#include <problem/base_meta.h>

#include <string>

/// Constrainted to Multi-Objective meta-problem
/**
 * Implements a meta-problem class that wraps some other constrained problems,
 * resulting in multi-objective problem.
 *
 * For a problem with c_n constraints, nobj + c_n objective functions, the first objectives functions are the original objective functions.
 *
 * Note: This constraints handling technique can only be used for <b>MINIMIZATION</b> problems.
 */

class con2mo : public pagmo::problem::base_meta
{
	public:
		explicit con2mo(const MGA& problem);

		pagmo::problem::base_ptr clone() const;
		std::string get_name() const;

	protected:
		std::string human_readable_extra() const;
		void objfun_impl(pagmo::fitness_vector &, const pagmo::decision_vector &) const;
};
