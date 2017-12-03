//
// Created by lisa on 26.11.17.
//

#include "AdvectionFiniteVolumes.h"
#include "types.h"
#include <cmath>
#include <algorithm>

void AdvectionFiniteVolumes::setBoundaryConditions(){
    m_u.at(0) = m_u.at(m_size); m_u.at(m_size+1) = m_u.at(1);
}

T AdvectionFiniteVolumes::computeLocalLaxFriedrichsFlux(T t){
    T maxWaveSpeed = 0.f;
    T deltaT = t - t0;
    t0 = t;
    //Loop over all intervalls
    for(unsigned int i = 1; i<m_size+1; i++){
        //Factor a for local Lax Friedrichs Method (aLeft = a_i-1/2, aRight = a_i+1/2
        T aLeft;
        T aRight;

        //Eigenvalue lambda = a
        aLeft  = std::abs(a);
        aRight = std::abs(a);

        //Compute fluxes for u

        m_uNetUpdatesRight.at(i) = 0.5f*(a*m_u.at(i)+a*m_u.at(i+1)-(aRight*(m_u.at(i+1)-m_u.at(i))));

        m_uNetUpdatesLeft.at(i) = 0.5f*(a*m_u.at(i-1)+a*m_u.at(i)-(aLeft*(m_u.at(i)-m_u.at(i-1))));

        // Update maxWaveSpeed
        maxWaveSpeed = std::max(std::max(maxWaveSpeed,aLeft),aRight);

    }
    // Compute CFL condition (delta_t <= c/a*(delta_x/2)) (DG book p. 97)
    T c = 0.4;
    T maxTimeStep = m_cellSize*(c/std::abs(a));

    return maxTimeStep;
}

vect AdvectionFiniteVolumes::updateUnknownsLocalLaxFriedrichs(T dt){
    vect res(m_size+2);
    for (unsigned int i = 0; i < m_size+2; i++){
        m_u.at(i) = m_u.at(i) - ((dt/m_cellSize)*(m_uNetUpdatesRight.at(i)-m_uNetUpdatesLeft.at(i)));
        res.at(i) = m_u.at(i);
    }
    return res;

}

vect AdvectionFiniteVolumes::setU() {
    vect result(m_u.size(), 0.0);
    for (unsigned int i = 0; i<m_u.size(); i++){
        result.at(i) = m_u.at(i);
    }
    return result;
}

T AdvectionFiniteVolumes::computeError(vect u0) {
    T res = 0.0;
    //T div = 0.0;
    for(unsigned int i = 1; i<m_size+1; i++){
        //res += std::abs((m_u.at(i) - u0.at(i))*(m_u.at(i) - u0.at(i)));
        //res += std::abs(m_u.at(i)- u0.at(i));
        res += std::abs((m_u.at(i) - u0.at(i))*(m_u.at(i) - u0.at(i)));
        //div += (u0.at(i)*u0.at(i));
    }
    //return std::sqrt(res);
    //return m_cellSize*res;
    return std::sqrt(m_cellSize*res);
    //return res/div;
}