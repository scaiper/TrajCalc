#pragma once

#include <cmath>
#include <vector>
#include <iostream>

class F
{
public:
    F(double mu, double Vin2, double Vout2, double alpha): k1(Vin2 / mu), k2(Vout2 / mu), alpha(alpha) {}
    
    double operator() (double R) const
    {
        return asin(1.0 / (1.0 + R * k1)) + asin(1.0 / (1.0 + R * k2)) - alpha;
    }
    
    double der(double R) const
    {
        double z1 = 1.0 / (1.0 + R * k1);
        double z2 = 1.0 / (1.0 + R * k2);
        return -1.0 / sqrt(1.0 - z1*z1) * z1*z1 * k1 - 1.0 / sqrt(1.0 - z2*z2) * z2*z2 * k2;
    }
    
private:
    double k1;
    double k2;
    double alpha;
};

double solveNewton(const F& f, double x)
{
    double fx = f(x);
    while (fx > 1.0e-5)
    {
        double x_prev = x;
        x -= fx / f.der(x);
        if (x_prev == x) break;
        fx = f(x);
    }

    return x;
}

double calcR(double alpha, double mu, double Vin2, double Vout2)
{
    F f(mu, Vin2, Vout2, alpha);
    
    double r0 = (1.0 / sin(0.5 * alpha) - 1.0) / std::max(Vin2, Vout2) * mu;
    return solveNewton(f, r0);
}

template<class vettore3D>
inline void fb_calc(double& DV, double& R, const vettore3D& v_rel_in, const vettore3D& v_rel_out, double mu_self)
{
    double Vin2  = v_rel_in[0]*v_rel_in[0]+v_rel_in[1]*v_rel_in[1]+v_rel_in[2]*v_rel_in[2];
    double Vout2 = v_rel_out[0]*v_rel_out[0]+v_rel_out[1]*v_rel_out[1]+v_rel_out[2]*v_rel_out[2];
    
    double alpha = acos( (v_rel_in[0]*v_rel_out[0] + v_rel_in[1]*v_rel_out[1] + v_rel_in[2]*v_rel_out[2]) / sqrt(Vin2 * Vout2) );
    
    R = calcR(alpha, mu_self, Vin2, Vout2);
    DV = std::abs(sqrt(2 * mu_self / R + Vin2) - sqrt(2 * mu_self / R + Vout2));
    if (std::isnan(DV)) {
        DV = 0.0;
    }
}
