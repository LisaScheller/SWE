//
// Created by lisa on 05.09.17.
//

#ifndef SWE1D_LAXFRIEDRICHSSOLVER_H
#define SWE1D_LAXFRIEDRICHSSOLVER_H
#include "types.h"

class LaxFriedrichsSolver{
private:
    vect m_h;
    vect m_hu;

    vect m_hNetUpdatesLeft;
    vect m_hNetUpdatesRight;


    vect m_huNetUpdatesLeft;
    vect m_huNetUpdatesRight;

    unsigned int m_size;

    T m_cellSize;

    static constexpr double g = 9.80665;
public:
    vect a_h;
    vect a_hu;
    LaxFriedrichsSolver(vect h, vect hu, vect ah, vect ahu, unsigned int size, T cellSize)
    : m_h(h),
    m_hu(hu),
    a_h(ah),
    a_hu(ahu),
    m_size(size),
    m_cellSize(cellSize),
    m_hNetUpdatesRight(size+2,0.0),
    m_hNetUpdatesLeft(size+2,0.0),
    m_huNetUpdatesRight(size+2,0.0),
    m_huNetUpdatesLeft(size+2,0.0)
    {



    }

    ~LaxFriedrichsSolver()
    {


    }
    //Standard computation of fluxes and timestep
    T computeLaxFriedrichsFlux();
    //Alternative computation of fluxes for use with updateUnknownsLaxFriedrichs2, needs
    //time of current timestep to compute provisoric dt
    T computeLaxFriedrichsFlux2(T t);

    void updateUnknownsLaxFriedrichs(T dt);
    vecu updateUnknownsLaxFriedrichs2(T dt);
    void updateUnknownsLaxFriedrichsDirect(T dt);

    T computeLocalLaxFriedrichsFlux(T t);

    vecu updateUnknownsLocalLaxFriedrichs(T dt);

    void solveAnalytically(T dt);


    vecu getAnalyticalSolution(T dt);
};

#endif //SWE1D_LAXFRIEDRICHSSOLVER_H
