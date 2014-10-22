#pragma once

#include "math_defines_hack.h"

#define PAGMO_NO_STD_VECTOR_STREAM_OVERLOADS
#include <planet.h>
#include <epoch.h>
#include <problem/base.h>

#include <string>
#include <vector>
#include <iostream>

class MGA : public pagmo::problem::base
{
private:
    enum ArrivalConstraint
    {
        IGNORE,
        INSERT_LO,
        AEROBRAKE
    };
    
public:
    MGA(std::vector<kep_toolbox::planet_ptr> seq,
        const kep_toolbox::epoch& t_begin, const kep_toolbox::epoch& t_end,
        double tof_min, double tof_max,
        bool zeroFB,
        bool optimizeToF);
        
    void setArrivalIgnore();
    void setArrivalInsert(boost::optional<double> R = boost::none);
    void setArrivalAerobrake(double max_Venc);
	
	void setMaxDV(boost::optional<double> maxDV);
    
    virtual void objfun_impl(pagmo::fitness_vector& f, const pagmo::decision_vector& x) const;
	virtual void compute_constraints_impl(pagmo::constraint_vector& c, const pagmo::decision_vector& x) const;
	
	void calc_objfun_constr(pagmo::fitness_vector& f, pagmo::constraint_vector& c, const pagmo::decision_vector& x) const;
	
    virtual pagmo::problem::base_ptr clone() const;
    virtual std::string human_readable_extra() const;
    
    void pretty(std::ostream& os, const pagmo::decision_vector& x) const;
    void plot(std::ostream& os, const pagmo::decision_vector& x) const;
    
private:
    void do_calc(const pagmo::decision_vector& x, pagmo::fitness_vector* f, std::ostream* os = nullptr, std::ostream* plot = nullptr, pagmo::constraint_vector* c = nullptr) const;

    std::vector<kep_toolbox::planet_ptr> seq;
    bool zeroFB;
    bool optimizeToF;
    ArrivalConstraint arrivalC;
    double ab_maxV;
    double common_mu;
	boost::optional<double> maxDV;
	double maxToF;
	boost::optional<double> insertR;
};
