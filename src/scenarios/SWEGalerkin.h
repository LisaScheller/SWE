//
// Created by lisa on 17.11.17.
//

#ifndef SWE1D_SWEGALERKIN_H
#define SWE1D_SWEGALERKIN_H
#define _USE_MATH_DEFINES
#include "../types.h"
#include <algorithm>
#include <cmath>
namespace scenarios
{

    class SWEGalerkin
    {

    private:
        /** Number of cells */
        const unsigned int m_size;

    public:
        SWEGalerkin(unsigned int size)
                : m_size(size)
        {
        }

        /**
         * @return Initial water height at pos
         */
        q getHeight(unsigned int i)
        {
            q res;
            T delta_x = getCellSize();
            //Compute x-value from interval number
            T pos = delta_x*i;
            T factor = m_size/100.0;
            if (pos <= 50.0){
                res.h.u0 = 3.0;
                res.h.u1 = 3.0;
            }
            else {
                res.h.u0 = 1.0;
                res.h.u1 = 1.0;
            }
            res.hu.u0 = 0.0;
            res.hu.u1 = 0.0;
            return res;
        }
        q getHeightGauss(unsigned int i){
            q res;
            T delta_x = getCellSize();
            //Compute x-value from interval number
            T pos = delta_x*i;
            T posr = delta_x *(i+1);
            //Standard deviation
            double sigma = 10.0;//3.0;
            //Get relative position in relation to m_size/2
            int x0 = 50;
            T xl = pos - x0;
            T xr = (posr) - x0;
            vect h(2);
            h.at(0) = 1.0/(sigma*(sqrt(2.0*M_PI)))*exp((-1.0/2.0)*(xl/sigma)*(xl/sigma));
            h.at(1) = 1.0/(sigma*(sqrt(2.0*M_PI)))*exp((-1.0/2.0)*(xr/sigma)*(xr/sigma));
            res.h.u0 = 10.0*h.at(0)+0.1;
            res.h.u1 = 10.0*h.at(1)+0.1;
            res.hu.u0 = 0.0;
            res.hu.u1 = 0.0;
            return res;
        }

        /**
         * @return Cell size of one cell (= domain size/number of cells)
         */
        T getCellSize()
        {
            return 100.0 / m_size;
        }
    };

}
#endif //SWE1D_SWEGALERKIN_H
