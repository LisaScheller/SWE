//
// Created by lisa on 05.09.17.
//
//TODO Find out why all these methods solve the dam break problem reasonably well but fail at the gaussian initial conditions
//TODO Assumption: Something in computing hu is still faulty
#include "types.h"
#include "LaxFriedrichsSolver.h"
#include "WavePropagation.h"
#include "cmath"
#include "algorithm"

T tPrev = 0;
T LaxFriedrichsSolver::computeLaxFriedrichsFlux() {
    float maxWaveSpeed = 0.f;
    //Loop over all edges
    for (unsigned int i = 1; i<m_size+2; i++){
        T maxEdgeSpeed;

        // Compute net updates (1-dim SWE, Leveque S. 255)
        m_hNetUpdatesRight[i] = m_hu[i];
        m_huNetUpdatesRight[i] = (m_hu[i]*m_hu[i])/m_h[i]+(0.5f*g*(m_h[i]*m_h[i]));
        m_hNetUpdatesLeft[i] = m_hu[i-1];
        m_huNetUpdatesLeft[i] = (m_hu[i-1]*m_hu[i-1])/m_h[i-1]+(0.5f*g*(m_h[i-1]*m_h[i-1]));

        //Compute edge speed

        if(m_hu[i-1] < 0.0){
            maxEdgeSpeed= -m_hu[i-1]+std::sqrt(g)*m_h[i-1];
        }else{
            maxEdgeSpeed=  m_hu[i-1]+std::sqrt(g)*m_h[i-1];
        }

        // Update maxWaveSpeed
        if (maxEdgeSpeed > maxWaveSpeed)
            maxWaveSpeed = maxEdgeSpeed;
    }
    // Compute CFL condition (with Courant number 0.4)

    T maxTimeStep = m_cellSize/maxWaveSpeed * .4f;

    return maxTimeStep;
}

T LaxFriedrichsSolver::computeLaxFriedrichsFlux2(T t) {
    float maxWaveSpeed = 0.f;
    T deltaT = t - tPrev;
    tPrev = t;
    //Loop over all edges
    for (unsigned int i = 1; i<m_size+2; i++){
        T maxEdgeSpeed;

        // Compute net updates (1-dim SWE, Leveque S. 255)
        m_hNetUpdatesLeft[i] = (0.5f*(m_hu[i-1]+m_hu[i]))-((m_cellSize/2*deltaT)*(m_h[i]-m_h[i-1]));
        m_hNetUpdatesRight[i] = (0.5f*(m_hu[i]+m_hu[i+1]))-((m_cellSize/2*deltaT)*(m_h[i+1]-m_h[i]));
        float fQLeft = ((m_hu[i-1]*m_hu[i-1])/m_h[i-1])+(0.5f*g*(m_h[i-1]*m_h[i-1]));
        float fQMiddle = ((m_hu[i]*m_hu[i])/m_h[i])+(0.5f*g*(m_h[i]*m_h[i]));
        float fQRight = ((m_hu[i+1]*m_hu[i+1])/m_h[i+1])+(0.5f*g*(m_h[i+1]*m_h[i+1]));
        m_huNetUpdatesLeft[i] = (0.5f*(fQLeft+fQMiddle))-((m_cellSize/2*deltaT)*(m_hu[i]-m_hu[i-1]));
        m_huNetUpdatesRight[i] = (0.5f*(fQMiddle+fQRight))-((m_cellSize/2*deltaT)*(m_hu[i+1]-m_hu[i]));


        //Compute edge speed

        if(m_hu[i-1] < 0.0){
            maxEdgeSpeed= -m_hu[i-1]+std::sqrt(g)*m_h[i-1];
        }else{
            maxEdgeSpeed=  m_hu[i-1]+std::sqrt(g)*m_h[i-1];
        }

        // Update maxWaveSpeed
        if (maxEdgeSpeed > maxWaveSpeed)
            maxWaveSpeed = maxEdgeSpeed;
    }
    // Compute CFL condition (with Courant number 0.4)

    T maxTimeStep = m_cellSize/maxWaveSpeed * .4f;

    return maxTimeStep;
}

void LaxFriedrichsSolver::updateUnknownsLaxFriedrichs(T dt){
    //Original method using the calculated Fluxes to compute h and hu at i
    //Loop over all inner cells
    //Leveque S. 71
    for (unsigned int i = 1; i < m_size+1; i++){
        m_h[i] = ((m_h[i-1]+m_h[i+1])/2) + dt/2*m_cellSize*(m_hNetUpdatesLeft[i]-m_hNetUpdatesRight[i]);
        m_hu[i] = ((m_hu[i-1]+m_hu[i+1])/2) + dt/2*m_cellSize*(m_huNetUpdatesLeft[i]-m_huNetUpdatesRight[i]);
}

}

void LaxFriedrichsSolver::updateUnknowsLaxFriedrichs2(T dt){
    //Alternative method using a combination of formula 4.21 (Leveque p. 71) for the fluxes
    //and formulas 4.6 and 4.7 (Leveque p. 66) to compute h and hu
    //Since the computation of the fluxes requires dt, the time of the current timestep
    //has to be given to the method calculating the fluxes
    //Loop over all inner cells
    for (unsigned int i = 1; i < m_size+1; i++){
        m_h[i] = ((m_h[i-1]+m_h[i+1])/2) + ((dt/m_cellSize)*(m_hNetUpdatesRight[i]-m_hNetUpdatesLeft[i]));
        m_hu[i] = ((m_hu[i-1]+m_hu[i+1])/2) + ((dt/m_cellSize)*(m_huNetUpdatesRight[i]-m_huNetUpdatesLeft[i]));
    }
}

void LaxFriedrichsSolver::updateUnknownsLaxFriedrichsDirect(T dt){
    //Alternative method that computes h and hu directly using formula 4.20 (Leveque p.71)
    //without needing an extra method for calculating the fluxes
    //Loop over all inner cells
    for (unsigned int i = 1; i < m_size+1; i++){
        //fQLeft = f(q(i-1)) fQRight = f(q(i+1))
        float fQLeft = ((m_hu[i-1]*m_hu[i-1])/m_h[i-1])+(0.5f*g*(m_h[i-1]*m_h[i-1]));
        float fQRight = ((m_hu[i+1]*m_hu[i+1])/m_h[i+1])+(0.5f*g*(m_h[i+1]*m_h[i+1]));
        m_h[i] = ((m_h[i-1]+m_h[i+1])/2) - ((dt/2.0f*m_cellSize)*(m_hu[i+1]-m_hu[i-1]));
        m_hu[i] = ((m_hu[i-1]+m_hu[i+1])/2) - ((dt/2.0f*m_cellSize)*(fQRight-fQLeft));
    }

}

