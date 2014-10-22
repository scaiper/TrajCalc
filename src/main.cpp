#include "mga.h"
#include "planet_vsop87.h"
#include "con2mo.h"

#include <archipelago.h>
#include <algorithm/jde.h>
#include <pagmo.h>
#include <algorithm/cstrs_self_adaptive.h>
#include <topology/fully_connected.h>

#include <boost/format.hpp>

#include <iostream>
#include <fstream>

const double G = 6.67384e-11;

typedef std::vector< std::pair< pagmo::fitness_vector, pagmo::decision_vector > > front_type;

front_type calc_front(const pagmo::population& pop)
{
	front_type front;
	auto&& fronts = pop.compute_pareto_fronts();
	if (!pop.problem().feasibility_c( pop.get_individual(fronts.front().front()).cur_c )) {
		return front;
	}
	for (auto indIdx: fronts.front())
	{
		auto&& ind = pop.get_individual(indIdx);
		front.emplace_back(ind.cur_f, ind.cur_x);
	}
	
	std::sort(front.begin(), front.end());
	return front;
}

double cmp_fronts(const front_type& front1, const front_type& front2)
{
	double s = 0.0;
	for (auto&& p1: front1)
	{
		for (auto&& p2: front2)
		{
			if (std::mismatch(p1.first.begin(), p1.first.end(), p2.first.begin(), [](double a, double b){ return b <= a; }).first == p1.first.end())
			{
				s = std::inner_product(p1.first.begin(), p1.first.end(), p2.first.begin(), s, [](double a, double b){ return std::max(a, b); }, [](double a, double b){ return (a - b) / b; });
			}
		}
	}
	return s;
}

front_type merge_fronts(const front_type& front1, const front_type& front2)
{
	front_type result;
	
	for (auto&& p1: front1)
	{
		bool dominated = std::any_of(front2.begin(), front2.end(), [&](auto&& p2){
			return std::mismatch(p1.first.begin(), p1.first.end(), p2.first.begin(), [](double a, double b){ return b <= a; }).first == p1.first.end() && p1.first != p2.first;
		});
		if (!dominated) result.push_back(p1);
	}
	
	for (auto&& p1: front2)
	{
		bool dominated = std::any_of(front1.begin(), front1.end(), [&](auto&& p2){
			return std::mismatch(p1.first.begin(), p1.first.end(), p2.first.begin(), [](double a, double b){ return b <= a; }).first == p1.first.end();
		});
		if (!dominated) result.push_back(p1);
	}
	
	std::sort(result.begin(), result.end());
	result.erase(std::unique(result.begin(), result.end(), [](auto&& a, auto&& b){ return a.first == b.first; }), result.end());
	return result;
}

pagmo::population merge_pop(const pagmo::archipelago& archi, const pagmo::problem::base& prob)
{
	pagmo::population pop(prob);
    for (pagmo::archipelago::size_type i = 0; i < archi.get_size(); ++i)
    {
		for (auto&& ind: archi.get_island(i)->get_population()) {
			pop.push_back(ind.cur_x);
		}
    }
	
	return pop;
}

void plot_front(std::ostream& os, const front_type& front)
{
	os << "set grid" << std::endl;
	os << "plot '-'" << std::endl;
	
	os << std::setprecision(std::numeric_limits<double>::digits10 + 2);
	
	for (auto&& p: front)
	{
		bool first = true;
		for (auto&& v: p.first) {
			os << (first ? "" : "\t") << v;
			first = false;
		}
		os << std::endl;
	}
	
	os << "e" << std::endl;
}

struct PlanetData
{
	boost::optional< std::string > filename;
	boost::optional< std::string > name;
	boost::optional< double > radius;
	boost::optional< double > safe_radius;
	boost::optional< double > mu_self;
	
	bool isInitialized() const { return filename && name && radius && safe_radius && mu_self; }
};

std::vector<double> read_vector(std::istream& in)
{
	std::stringstream buf;
	while (in)
	{
		char c;
		in >> c;
		if (!in) break;
		if (c == '[' || c == ']' || c == ',') c = ' ';
		buf << c;
	}
	
	std::vector<double> result;
	while (buf)
	{
		double value;
		buf >> value;
		if (!buf) break;
		result.push_back(value);
	}
	return result;
}

void print_solution(const MGA& prob, const pagmo::decision_vector& x)
{
	prob.pretty(std::cout, x);
	{
		std::ofstream os("out.txt");
		prob.pretty(os, x);
	}
}

int main()
{
    try {
		std::map< std::string, planet_vsop87 > planets;
		{
			PlanetData tmpPlanet;
			
			std::ifstream in("planets.cfg");
            while (in.good())
            {
                std::string line;
                std::getline(in, line);
				
				std::string key, value;
				std::istringstream buf(line);
				buf >> key >> value;
				if ((in.fail() || key == "[planet]") && tmpPlanet.isInitialized())
				{
					VSOP87_Data vsop87_data(std::ifstream(*tmpPlanet.filename));
					planets.emplace(std::piecewise_construct,
					                std::forward_as_tuple(*tmpPlanet.name),
					                std::forward_as_tuple(vsop87_data, *tmpPlanet.name, *tmpPlanet.radius, *tmpPlanet.safe_radius, *tmpPlanet.mu_self));
					tmpPlanet = PlanetData{};
					
				}
                if (in.fail()) break;
				
				if (key == "file") {
					tmpPlanet.filename = std::move(value);
				} else if (key == "name") {
					tmpPlanet.name = std::move(value);
				} else if (key == "radius") {
					tmpPlanet.radius = boost::lexical_cast<double>(value);
				} else if (key == "safe_radius") {
					tmpPlanet.safe_radius = boost::lexical_cast<double>(value);
				} else if (key == "mass") {
					tmpPlanet.mu_self = boost::lexical_cast<double>(value) * G;
				}
			}
		}
		
        std::vector<kep_toolbox::planet_ptr> seq;
        double lowest_start_year = 1990.0;
        double highest_start_year = 2020.0;
        double min_tof = 0.0;
        double max_tof = 10.0;
        bool zeroFB = false;
        bool optimizeToF = false;
        double abMaxV;
		boost::optional<double> maxDV;
        int arrival = 0;
		boost::optional<double> insertR;
        
        unsigned gen = 100;
        unsigned threads = 4;
        unsigned population_per_thread = 80;
		
		std::vector<double> solution;
        
        {
            std::ifstream in("config.txt");
            while (in.good())
            {
                std::string line;
                std::getline(in, line);
                if (!in.fail())
                {
                    std::string key;
                    std::istringstream buf(line);
                    buf >> key;
                    if (key == "planets")
                    {
                        while (buf.good())
                        {
                            std::string value;
                            buf >> value;
                            if (!buf.fail())
                            {
                                seq.push_back(planets.at(value).clone());
							}
                        }
                    }
                    else if (key == "lowest_start_year")
                    {
                        buf >> lowest_start_year;
                    }
                    else if (key == "highest_start_year")
                    {
                        buf >> highest_start_year;
                    }
                    else if (key == "min_tof")
                    {
                        buf >> min_tof;
                    }
                    else if (key == "max_tof")
                    {
                        buf >> max_tof;
                    }
                    else if (key == "zero_fly_by")
                    {
                        std::string value;
                        buf >> value;
                        zeroFB = (value == "y" || value == "yes" || value == "1" || value == "true" || value == "on");
                    }
                    else if (key == "optimizeToF")
                    {
                        std::string value;
                        buf >> value;
                        optimizeToF = (value == "y" || value == "yes" || value == "1" || value == "true" || value == "on");
                    }
                    else if (key == "max_dv")
                    {
						double value;
                        buf >> value;
						if (buf) maxDV = value;
                    }
                    else if (key == "arrival")
                    {
                        std::string value;
                        buf >> value;
                        if (value == "ignore")
                        {
                            arrival = 0;
                        }
                        else if (value == "insert")
                        {
							double r;
                            arrival = 1;
							buf >> r;
							if (buf) insertR = r;
                        }
                        else if (value == "aerobrake")
                        {
                            arrival = 2;
                            buf >> abMaxV;
                        }
                    }
                    else if (key == "gen")
                    {
                        buf >> gen;
                    }
                    else if (key == "threads")
                    {
                        buf >> threads;
                    }
                    else if (key == "population_per_thread")
                    {
                        buf >> population_per_thread;
                    }
					else if (key == "solution")
					{
						solution = read_vector(buf);
					}
                }
            }
        }
        
        MGA prob(std::move(seq), (lowest_start_year - 2000.0) / ASTRO_DAY2YEAR, (highest_start_year - 2000.0) / ASTRO_DAY2YEAR, min_tof, max_tof, zeroFB, optimizeToF);
		prob.setMaxDV(maxDV);
		
        switch (arrival)
        {
        case 0: {
            prob.setArrivalIgnore();
        }
        break;
        case 1: {
			prob.setArrivalInsert(insertR);
        }
        break;
        case 2: {
            prob.setArrivalAerobrake(abMaxV);
        }
        break;
        }
        
		if (!solution.empty())
		{
			print_solution(prob, solution);
			return 0;
		}
		
		pagmo::algorithm::nsga2 algo(gen);
		con2mo q(prob);
        
        pagmo::topology::fully_connected topo;
        
		pagmo::archipelago archi(pagmo::algorithm::ms(algo, 30), q, threads, population_per_thread);
		std::cout << "Starting..." << std::endl;
		archi.evolve();
		archi.join();
		std::cout << "Calculating..." << std::endl;
		
		for (size_t i = 0; i < threads; ++i) archi.set_algorithm(i, algo);
		archi.set_topology(topo);
		
        bool prev_good = false;
        int evolve_count = 0;
        double champf_prev;
		bool feas = false;
		double delta = HUGE_VAL;
		front_type prev_front = calc_front(merge_pop(archi, prob));
        do
        {
			prev_good = feas && delta < 1.0e-6;
			
			champf_prev = delta;
            archi.evolve(5);
            archi.join();
            evolve_count += 5;
			
			front_type cur_front = calc_front(merge_pop(archi, prob));
			delta = cmp_fronts(prev_front, cur_front);
			feas = !cur_front.empty();
			cur_front = merge_fronts(prev_front, cur_front);
			
            std::cout << "After " << evolve_count << ": delta = " << delta << ", feasibility: " << feas << ", min DV: " << (cur_front.empty() ? 0.0 : cur_front.front().first[0])
			          << ", front size: " << cur_front.size()  << std::endl;
			
			if (cur_front.size() > prev_front.size() && cur_front.size() < 20) delta += cur_front.size() - prev_front.size();
			prev_front = std::move(cur_front);
        } while (!(feas && delta < 1.0e-6 && prev_good && (prev_front.size() >= 10 || !optimizeToF)) && (evolve_count < 500 || prev_front.empty()));
		
		if (optimizeToF)
		{
			{
				std::ofstream os("dv_ToF.plt");
				plot_front(os, prev_front);
			}
			{
				std::ofstream os("dv_ToF.txt");
				for (auto&& p: prev_front)
				{
					os << p.first << "\t" << p.second << std::endl;
				}
			}
		}

		if (!prev_front.empty())
		{
			auto&& min_x = prev_front[0].second;

			print_solution(prob, min_x);
			{
				std::ofstream os("out.plt");
				prob.plot(os, min_x);
			}
		}
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << std::endl;
    }
    
    return 0;
}
