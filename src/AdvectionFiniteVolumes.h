//
// Created by lisa on 26.11.17.
//

#ifndef SWE1D_ADVECTIONFINITEVOLUMES_H
#define SWE1D_ADVECTIONFINITEVOLUMES_H

#include "types.h"

class AdvectionFiniteVolumes{
private:

    //Left and right fluxes
    vect m_uNetUpdatesLeft;
    vect m_uNetUpdatesRight;
    //u
    vect m_u;

    int m_size;

    //Time variables
    T t0 = 0.0;
    //Advection velocity
    T a;

    T m_cellSize;
public:
    AdvectionFiniteVolumes(float a_def, vect h, unsigned int size, T cellsize):
            a(a_def),
            m_size(size),
            m_u(h.size()),
            m_uNetUpdatesRight(size+2,0.0),
            m_uNetUpdatesLeft(size+2,0.0),
            m_cellSize(cellsize)

    {
        for (unsigned int i = 1; i<size+1; i++){
            m_u.at(i) = h[i];

        }
        m_u[0] = m_u.at(1);
        m_u[size+1] = m_u.at(size);


    }

    ~AdvectionFiniteVolumes(){

    }
    //Methods
    void setBoundaryConditions();
    T computeLocalLaxFriedrichsFlux(T t);
    void updateUnknowns(T t);

    vect setU();

    T computeError(vect vector);

    vect updateUnknownsLocalLaxFriedrichs(T dt);
};


#endif //SWE1D_ADVECTIONFINITEVOLUMES_H
