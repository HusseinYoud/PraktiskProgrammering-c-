#include "Spliner.h"
#include<cassert>
#include<cmath>

namespace pp{
    int binsearch(const vec& x, double z){
        assert (x.size() >= 2);
        assert (z >= x[0]);
        assert (z <= x[x.size()-1]);

        int i = 0;
        int j = x.size()-1;

        while (j - i > 1){
            int mid = (i+j)/2;

            if (z > x[mid]){
                i = mid;
            } else {
                j = mid;
            }
        }
        return i;
    }

    double linterp(const vec& x, const vec& y, double z){
        assert (x.size() == y.size());
        assert (x.size() >= 2);

        int i = binsearch(x, z);

        double dx = x[i+1] - x[i];
        double dy = y[i+1] - y[i];

        assert(dx > 0.0);

        return y[i] + dy/dx * (z - x[i]);
    }
    double linterp_intg(const vec& x, const vec& y, double z){
        assert (x.size() == y.size());
        assert (x.size() >= 2);
        
        int interval = binsearch(x, z);
        double integral = 0.0;

        for (int i = 0; i < interval; ++i){
            double dx = x[i+1] - x[i];
            assert (dx > 0.0);
            integral += 0.5 * (y[i] + y[i+1]) * dx;
        }

        double dx = x[interval +1] - x[interval];
        double dy = y[interval +1] - y[interval];
        
        assert(dx > 0.0);

        double slope = dy / dx;
        double dz = z-x[interval];
        
        integral += y[interval] * dz + 0.5 * slope * dz * dz;

        return integral;
        }
    qspline::qspline(const vec& xs, const vec& ys)
    : x(xs),
      y(ys),
      b(xs.size() > 1 ? xs.size() - 1 : 0),
      c(xs.size() > 1 ? xs.size() - 1 : 0)
{
    assert(xs.size() == ys.size());
    assert(xs.size() >= 2);

    const int n = x.size();

    vec h(n - 1);
    vec p(n - 1);

    for (int i = 0; i < n - 1; ++i) {
        h[i] = x[i + 1] - x[i];

        assert(h[i] > 0.0);

        p[i] = (y[i + 1] - y[i]) / h[i];
    }

    c[0] = 0.0;

    for (int i = 0; i < n - 2; ++i) {
        c[i + 1] =
            (p[i + 1] - p[i] - c[i] * h[i])
            / h[i + 1];
    }


    c[n - 2] /= 2.0;


    for (int i = n - 3; i >= 0; --i) {
        c[i] =
            (p[i + 1] - p[i] - c[i + 1] * h[i + 1])
            / h[i];
    }


    for (int i = 0; i < n - 1; ++i) {
        b[i] = p[i] - c[i] * h[i];
    }
}


int qspline::binsearch(double z) const
{
    return pp::binsearch(x, z);
}


double qspline::eval(double z) const
{
    int i = binsearch(z);

    double h = z - x[i];

    return y[i]
         + b[i] * h
         + c[i] * h * h;
}


double qspline::deriv(double z) const
{
    int i = binsearch(z);

    double h = z - x[i];

    return b[i] + 2.0 * c[i] * h;
}


double qspline::intg(double z) const
{
    int interval = binsearch(z);

    double integral = 0.0;


    for (int i = 0; i < interval; ++i) {
        double h = x[i + 1] - x[i];

        integral += y[i] * h
                  + 0.5 * b[i] * h * h
                  + c[i] * h * h * h / 3.0;
    }


    double h = z - x[interval];

    integral += y[interval] * h
              + 0.5 * b[interval] * h * h
              + c[interval] * h * h * h / 3.0;

    return integral;
}
        


}