#include "vsop87.h"

#include "math_defines_hack.h"

#include <astro_constants.h>

VSOP87_Data::VSOP87_Data(std::istream& in):
	abc(new std::vector< VarData >())
{
	while (in)
	{
		std::string str;
		std::getline(in, str);
		if (!in) break;
		
		if (str.length() >= 131 && str[1] == '1')
		{
			unsigned var = str[3] - '0' - 1;
			unsigned alpha = str[4] - '0';
			const char *b, *e;
			
			for (b = str.data() + 79; *b == ' '; ++b);
			e = str.data() + 96;
			double A = boost::lexical_cast<double>(b, e - b + 1);
			
			for (b = str.data() + 97; *b == ' '; ++b);
			e = str.data() + 110;
			double B = boost::lexical_cast<double>(b, e - b + 1);
			
			for (b = str.data() + 111; *b == ' '; ++b);
			e = str.data() + 130;
			double C = boost::lexical_cast<double>(b, e - b + 1);
			
			if (abc->size() <= var) abc->resize(var + 1);
			
			if ((*abc)[var].size() <= alpha) (*abc)[var].resize(alpha + 1);
			
			(*abc)[var][alpha].push_back({A, B, C});
		}
	}
}

VSOP87_Data::VSOP87_Data(std::istream&& in):
	VSOP87_Data(in)
{
}

double VSOP87_Data::calc_Var(unsigned var, double t) const
{
	double result = 0.0;
	double tn = 1.0;
	for (auto&& curABC: (*abc)[var])
	{
		double s = 0.0;
		for (auto&& d: curABC)
		{
			s += d.A * cos(d.B + d.C * t);
		}
		result += s * tn;
		tn *= t;
	}
	
	return result * ASTRO_AU;
}

double VSOP87_Data::calc_dVar(unsigned var, double t) const
{
	if ((*abc)[var].empty()) return 0.0;
	
	double result = 0.0;

	{
		for (auto&& d: (*abc)[var].front())
		{
			result -= d.A * d.C * sin(d.B + d.C * t);
		}
	}
	
	double tn_1 = 1.0;
	for (size_t alpha = 1; alpha < (*abc)[var].size(); ++alpha)
	{
		double s = 0.0;
		for (auto&& d: (*abc)[var][alpha])
		{
			s += d.A * (alpha * cos(d.B + d.C * t) - t * d.C * sin(d.B + d.C * t));
		}
		result += s * tn_1;
		tn_1 *= t;
	}
	
	return result * ASTRO_AU / 365250.0 / ASTRO_DAY2SEC;
}
