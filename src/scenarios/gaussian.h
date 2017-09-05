//
// Created by lisa on 05.09.17.
//

#ifndef SWE1D_GAUSSIAN_H
#define SWE1D_GAUSSIAN_H
#define _USE_MATH_DEFINES
#include "../types.h"
#include <algorithm>
#include <cmath>

namespace scenarios
{

    class Gaussian
    {

    private:
        /** Number of cells */
        const unsigned int m_size;

    public:
        Gaussian(unsigned int size)
                : m_size(size)
        {
        }

        /**
         * @return Initial water height at pos
         */
        double getHeight(unsigned int pos)
        {
            //Gaussian distribution with standard deviation 50/3 and expected value at m_size/2
            //Maximum (height 1) at m_size/2

            //Standard deviation
            double sigma = 50.0/3.0;
            //Get relative position in relation to m_size/2
            int x0 = m_size/2;
            int x = pos - x0;
            double h = 1.0/(sigma*(sqrt(2.0*M_PI)))*exp((-1.0/2.0)*(x/sigma)*(x/sigma));
            //TODO Find correct factor for multiplication with h so that there are neither absurd oscillations nor no effect is visible
            return 20.0f*h;
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



#endif //SWE1D_GAUSSIAN_H
