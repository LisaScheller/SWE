//
// Created by lisa on 14.11.17.
//

#ifndef SWE1D_ADVECTIONGALERKIN_H
#define SWE1D_ADVECTIONGALERKIN_H
#define _USE_MATH_DEFINES
#include "../types.h"
#include <algorithm>
#include <cmath>
namespace scenarios
{

    class AdvectionGalerkin
    {

    private:
        /** Number of cells */
        const unsigned int m_size;

    public:
        AdvectionGalerkin(unsigned int size)
                : m_size(size)
        {
        }

        /**
         * @return Initial water height at pos
         */
        u getHeight(unsigned int i)
        {
            T delta_x = getCellSize();
            //Compute x-value from interval number
            T pos = delta_x*i;
            T posr = delta_x *(i+1);
            //Gaussian distribution with standard deviation 50/3 and expected value at m_size/2
            //Maximum (height 1) at m_size/2
            u res;
            //Standard deviation
            double sigma = 10.0;//3.0;
            //Get relative position in relation to m_size/2
            int x0 = 50;
            T xl = pos - x0;
            T xr = (posr) - x0;
            T *h = new T[2];
            h[0] = 1.0/(sigma*(sqrt(2.0*M_PI)))*exp((-1.0/2.0)*(xl/sigma)*(xl/sigma));
            h[1] = 1.0/(sigma*(sqrt(2.0*M_PI)))*exp((-1.0/2.0)*(xr/sigma)*(xr/sigma));

            res.u0 = 10*h[0]+0.1f;
            res.u1 = 10*h[1]+0.1f;
            return res;
        }

        /**
         * @return Cell size of one cell (= domain size/number of cells)
         */
        T getCellSize()
        {
            return 100.f / m_size;
        }
    };

}

#endif //SWE1D_ADVECTIONGALERKIN_H
