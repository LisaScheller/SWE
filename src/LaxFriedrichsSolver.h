//
// Created by lisa on 05.09.17.
//

#ifndef SWE1D_LAXFRIEDRICHSSOLVER_H
#define SWE1D_LAXFRIEDRICHSSOLVER_H
#include "types.h"

class LaxFriedrichsSolver{
private:
    T *m_h;
    T *m_hu;

    T *m_hNetUpdatesLeft;
    T *m_hNetUpdatesRight;


    T *m_huNetUpdatesLeft;
    T *m_huNetUpdatesRight;

    unsigned int m_size;

    T m_cellSize;

    static constexpr double g = 9.80665;
public:
    T *a_h;
    T *a_hu;
    LaxFriedrichsSolver(T *h, T *hu, T *ah, T *ahu, unsigned int size, T cellSize)
    : m_h(h),
    m_hu(hu),
    a_h(ah),
    a_hu(ahu),
    m_size(size),
    m_cellSize(cellSize)
    {
        // Allocate net updates
        m_hNetUpdatesLeft = new T[size+1];
        m_hNetUpdatesRight = new T[size+1];
        m_huNetUpdatesLeft = new T[size+1];
        m_huNetUpdatesRight = new T[size+1];


    }

    ~LaxFriedrichsSolver()
    {
        // Free allocated memory
        delete [] m_hNetUpdatesLeft;
        delete [] m_hNetUpdatesRight;
        delete [] m_huNetUpdatesLeft;
        delete [] m_huNetUpdatesRight;

    }
    //Standard computation of fluxes and timestep
    T computeLaxFriedrichsFlux();
    //Alternative computation of fluxes for use with updateUnknownsLaxFriedrichs2, needs
    //time of current timestep to compute provisoric dt
    T computeLaxFriedrichsFlux2(T t);

    void updateUnknownsLaxFriedrichs(T dt);
    void updateUnknownsLaxFriedrichs2(T dt);
    void updateUnknownsLaxFriedrichsDirect(T dt);

    T computeLocalLaxFriedrichsFlux(T t);

    void updateUnknownsLocalLaxFriedrichs(T dt);

    void solveAnalytically(T dt);


};

#endif //SWE1D_LAXFRIEDRICHSSOLVER_H
