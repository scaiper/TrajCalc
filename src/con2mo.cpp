#include "con2mo.h"

#include <exceptions.h>

#include <cmath>
#include <algorithm>

static int __mo_dimension__(const pagmo::problem::base &original_problem)
{
	return original_problem.get_f_dimension() + original_problem.get_c_dimension();
}

con2mo::con2mo(const MGA &problem):
	pagmo::problem::base_meta(
		problem,
		problem.get_dimension(),
		problem.get_i_dimension(),
		problem.get_f_dimension() + problem.get_c_dimension(),
		0,
		0)
{
}

pagmo::problem::base_ptr con2mo::clone() const
{
	return pagmo::problem::base_ptr(new con2mo(*this));
}

void con2mo::objfun_impl(pagmo::fitness_vector &f, const pagmo::decision_vector &x) const
{
	auto prob = dynamic_cast<const MGA*>(m_original_problem.get());
	if (!prob) pagmo_throw(value_error, "Error: Wrong problem, should be MGA!"); 
	
	pagmo::constraint_vector c;
	pagmo::decision_vector original_f;
	prob->calc_objfun_constr(original_f, c, x);

	c_size_type number_of_constraints = c.size();
	c_size_type number_of_eq_constraints = number_of_constraints - prob->get_ic_dimension();
	
	const std::vector<double> &c_tol = prob->get_c_tol();
	
	// modify equality constraints to behave as inequality constraints:
	for(c_size_type i=0; i < number_of_constraints; i++) {
		if(i<number_of_eq_constraints){
			c[i] = std::abs(c[i]) - c_tol[i];
		}
		else{
			c[i] = c[i] - c_tol[i];
		}
	}

	f.assign(get_f_dimension(), 0.0);

	auto f_it = std::copy(original_f.begin(), original_f.end(), f.begin());
	std::transform(c.begin(), c.end(), f_it, [](double v){ return std::max(v, 0.0); });
}

std::string con2mo::human_readable_extra() const
{
	std::ostringstream oss;
	oss << m_original_problem->human_readable_extra() << std::endl;
	return oss.str();
}

std::string con2mo::get_name() const
{
	return m_original_problem->get_name() + " [con2mo]";
}
