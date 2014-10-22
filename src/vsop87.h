#pragma once

#include <iostream>
#include <memory>
#include <string>
#include <vector>

class VSOP87_Data
{
private:
	struct ABC
	{
		double A;
		double B;
		double C;
		
		bool operator< (const ABC& other) const { return A < other.A; }
	};

	typedef std::vector< ABC > AlphaData;
	typedef std::vector< AlphaData > VarData;

public:
	VSOP87_Data(std::istream& in);
	VSOP87_Data(std::istream&& in);

	double calc_Var(unsigned var, double t) const;
	double calc_dVar(unsigned var, double t) const;
	
private:
	std::shared_ptr< std::vector< VarData > > abc;
};
