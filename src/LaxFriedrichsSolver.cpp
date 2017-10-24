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
    T maxWaveSpeed = 0.f;
    //Loop over all edges
    for (unsigned int i = 1; i<m_size+2; i++){
        T maxEdgeSpeed;

        // Compute net updates (1-dim SWE, Leveque S. 255)
        m_hNetUpdatesRight[i] = m_hu[i];
        m_huNetUpdatesRight[i] = (m_hu[i]*m_hu[i])/m_h[i]+(0.5*g*(m_h[i]*m_h[i]));
        m_hNetUpdatesLeft[i] = m_hu[i-1];
        m_huNetUpdatesLeft[i] = (m_hu[i-1]*m_hu[i-1])/m_h[i-1]+(0.5*g*(m_h[i-1]*m_h[i-1]));

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

    T maxTimeStep = m_cellSize/maxWaveSpeed * .4;

    return maxTimeStep;
}

T LaxFriedrichsSolver::computeLaxFriedrichsFlux2(T t) {
    T maxWaveSpeed = 0.0;
    T deltaT = t - tPrev;
    tPrev = t;
    //Loop over all edges
    for (unsigned int i = 1; i<m_size+2; i++){
        T maxEdgeSpeed;

        // Compute net updates (1-dim SWE, Leveque S. 255)
        m_hNetUpdatesLeft[i] = (0.5*(m_hu[i-1]+m_hu[i]))-((m_cellSize/2*deltaT)*(m_h[i]-m_h[i-1]));
        m_hNetUpdatesRight[i] = (0.5*(m_hu[i]+m_hu[i+1]))-((m_cellSize/2*deltaT)*(m_h[i+1]-m_h[i]));
        T fQLeft = ((m_hu[i-1]*m_hu[i-1])/m_h[i-1])+(0.5*g*(m_h[i-1]*m_h[i-1]));
        T fQMiddle = ((m_hu[i]*m_hu[i])/m_h[i])+(0.5*g*(m_h[i]*m_h[i]));
        T fQRight = ((m_hu[i+1]*m_hu[i+1])/m_h[i+1])+(0.5*g*(m_h[i+1]*m_h[i+1]));
        m_huNetUpdatesLeft[i] = (0.5*(fQLeft+fQMiddle))-((m_cellSize/2*deltaT)*(m_hu[i]-m_hu[i-1]));
        m_huNetUpdatesRight[i] = (0.5*(fQMiddle+fQRight))-((m_cellSize/2*deltaT)*(m_hu[i+1]-m_hu[i]));


        //Compute edge speed
	T aLeft  = std::max(abs(m_hu[i-1]/m_h[i-1]) + sqrt(g*m_h[i-1]),abs(m_hu[i]/m_h[i]) + sqrt(g*m_h[i]));
	T aRight = std::max(abs(m_hu[i+1]/m_h[i+1]) + sqrt(g*m_h[i+1]),abs(m_hu[i]/m_h[i]) + sqrt(g*m_h[i]));

    maxWaveSpeed = std::fmax(std::fmax(maxWaveSpeed, aLeft),aRight);
	//maxWaveSpeed = std::max(std::max(maxWaveSpeed,aLeft),aRight);

    }
    // Compute CFL condition (with Courant number 0.4)

    T maxTimeStep = m_cellSize/maxWaveSpeed * .4;

    return maxTimeStep;
}

T LaxFriedrichsSolver::computeLocalLaxFriedrichsFlux(T t){
    T maxWaveSpeed = 0.f;
    T deltaT = t - tPrev;
    tPrev = t;
    //Loop over all edges
    for (unsigned int i = 1; i<m_size+2; i++){

        //Factor a for local Lax Friedrichs Method (aLeft = a_i-1/2, aRight = a_i+1/2
        T aLeft;
        T aRight;
        // //Compute Frobenius Norm of f'(Q_i-1) and f'(Q_i) and compare them to find the maximum
        // T frobLeft = std::sqrt(1+((-1*(m_hu[i-1]/m_h[i-1])*(m_hu[i-1]/m_h[i-1])+(g*m_h[i-1]))*(-1*(m_hu[i-1]/m_h[i-1])*(m_hu[i-1]/m_h[i-1])+(g*m_h[i-1])))+((2*m_hu[i-1]/m_h[i-1])*(2*m_hu[i-1]/m_h[i-1])));
        // T frobRight = std::sqrt(1+((-1*(m_hu[i]/m_h[i])*(m_hu[i]/m_h[i])+(g*m_h[i]))*(-1*(m_hu[i]/m_h[i])*(m_hu[i]/m_h[i])+(g*m_h[i])))+((2*m_hu[i]/m_h[i])*(2*m_hu[i]/m_h[i])));
        // if(frobLeft>frobRight){
        //     aLeft = frobLeft;
        // }
        // else aLeft = frobRight;
        // //Compute Frobenius Norm of f'(Q_i) and f'(Q_i+1) and compare them to find the maximum
        // T frobLeft2 = std::sqrt(1+((-1*(m_hu[i]/m_h[i])*(m_hu[i]/m_h[i])+(g*m_h[i]))*(-1*(m_hu[i]/m_h[i])*(m_hu[i]/m_h[i])+(g*m_h[i])))+((2*m_hu[i]/m_h[i])*(2*m_hu[i]/m_h[i])));
        // T frobRight2 = std::sqrt(1+((-1*(m_hu[i+1]/m_h[i+1])*(m_hu[i+1]/m_h[i+1])+(g*m_h[i+1]))*(-1*(m_hu[i+1]/m_h[i+1])*(m_hu[i+1]/m_h[i+1])+(g*m_h[i+1])))+((2*m_hu[i+1]/m_h[i+1])*(2*m_hu[i+1]/m_h[i+1])));
        // if(frobLeft2 > frobRight2){
        //     aRight = frobLeft2;
        // }
        // else aRight = frobRight2;

	//hu[i-1]
	aLeft  = std::max(abs(m_hu[i-1]/m_h[i-1]) + sqrt(g*m_h[i-1]),abs(m_hu[i]/m_h[i]) + sqrt(g*m_h[i]));
	aRight = std::max(abs(m_hu[i+1]/m_h[i+1]) + sqrt(g*m_h[i+1]),abs(m_hu[i]/m_h[i]) + sqrt(g*m_h[i])); 

        // Compute fluxes for h and hu
        m_hNetUpdatesLeft[i] = 0.5*(m_hu[i-1]+m_hu[i]-(aLeft*(m_h[i]-m_h[i-1])));
        m_hNetUpdatesRight[i] = 0.5*(m_hu[i]+m_hu[i+1]-(aRight*(m_h[i+1]-m_h[i])));
        T fQLeft = ((m_hu[i-1]*m_hu[i-1])/m_h[i-1])+(0.5*g*(m_h[i-1]*m_h[i-1]));
        T fQMiddle = ((m_hu[i]*m_hu[i])/m_h[i])+(0.5*g*(m_h[i]*m_h[i]));
        T fQRight = ((m_hu[i+1]*m_hu[i+1])/m_h[i+1])+(0.5*g*(m_h[i+1]*m_h[i+1]));
        m_huNetUpdatesLeft[i] = 0.5*(fQLeft+fQMiddle-(aLeft*(m_hu[i]-m_hu[i-1])));
        m_huNetUpdatesRight[i] = 0.5*(fQMiddle+fQRight-(aRight*(m_hu[i+1]-m_hu[i])));

        // Update maxWaveSpeed
	maxWaveSpeed = std::max(std::max(maxWaveSpeed,aLeft),aRight);
    }
    // Compute CFL condition (with Courant number 0.4)

    T maxTimeStep = m_cellSize/maxWaveSpeed * .4;

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

void LaxFriedrichsSolver::updateUnknownsLaxFriedrichs2(T dt){
    //Alternative method using a combination of formula 4.21 (Leveque p. 71) for the fluxes
    //and formulas 4.6 and 4.7 (Leveque p. 66) to compute h and hu
    //Since the computation of the fluxes requires dt, the time of the current timestep
    //has to be given to the method calculating the fluxes
    //Loop over all inner cells
    for (unsigned int i = 1; i < m_size+1; i++){
        /*m_h[i] = ((m_h[i-1]+m_h[i+1])/2) + ((dt/m_cellSize)*(m_hNetUpdatesRight[i]-m_hNetUpdatesLeft[i]));
        m_hu[i] = ((m_hu[i-1]+m_hu[i+1])/2) + ((dt/m_cellSize)*(m_huNetUpdatesRight[i]-m_huNetUpdatesLeft[i]));
        */
        m_h[i] = m_h[i] - ((dt/m_cellSize)*(m_hNetUpdatesRight[i]-m_hNetUpdatesLeft[i]));
        m_hu[i] = m_hu[i] - ((dt/m_cellSize)*(m_huNetUpdatesRight[i]-m_huNetUpdatesLeft[i]));
    }
}

void LaxFriedrichsSolver::updateUnknownsLocalLaxFriedrichs(T dt){
    for (unsigned int i = 1; i < m_size+1; i++){
        m_h[i] = m_h[i] - ((dt/m_cellSize)*(m_hNetUpdatesRight[i]-m_hNetUpdatesLeft[i]));
        m_hu[i] = m_hu[i] - ((dt/m_cellSize)*(m_huNetUpdatesRight[i]-m_huNetUpdatesLeft[i]));
    }
}

void LaxFriedrichsSolver::solveAnalytically(T dt){
    //Values for specific case h_l=11, h_r=10, u_l=u_r=0
    //Computed with maple
    T t = tPrev +dt;

    T h_l=11.0;
    T u_l=0.0;
    T h_r=10.0;
    T u_r=0.0;
    T s=10.26901816;
    T u_m=0.48339953;
    T h_m=10.49398975;
    T lambda_1_l=u_l-sqrt(g*h_l);
    T lambda_1_m=u_m-sqrt(g*h_m);
    for (unsigned int i =0; i<m_size; i++){
        //T x = i-50.0; //oder i+(i+1)/2?
        T x = (((i+(i+1))/2)/4)-50.0;
        if(x/t<=lambda_1_l){
            a_h[i]=h_l;
            a_hu[i]=a_h[i]*u_l;
        }
        else if(x/t>=s){
            a_h[i]=h_r;
            a_hu[i]=a_h[i]*u_r;
        }
        else if (x/t>lambda_1_l && x/t<lambda_1_m){
            a_h[i]=0.01133018014*(20.77239996-x/t)*(20.77239996-x/t);
            a_hu[i]=a_h[i]*(u_m-2.0*(sqrt(g*h_m)-sqrt(g*a_h[i])));
        }
        else {
            a_h[i]=h_m;
            a_hu[i]=a_h[i]*u_m;
        }
    }

}


void LaxFriedrichsSolver::updateUnknownsLaxFriedrichsDirect(T dt){
    //Alternative method that computes h and hu directly using formula 4.20 (Leveque p.71)
    //without needing an extra method for calculating the fluxes
    //Loop over all inner cells
    for (unsigned int i = 1; i < m_size+1; i++){
        //fQLeft = f(q(i-1)) fQRight = f(q(i+1))
        T fQLeft = ((m_hu[i-1]*m_hu[i-1])/m_h[i-1])+(0.5*g*(m_h[i-1]*m_h[i-1]));
        T fQRight = ((m_hu[i+1]*m_hu[i+1])/m_h[i+1])+(0.5*g*(m_h[i+1]*m_h[i+1]));
        m_h[i] = ((m_h[i-1]+m_h[i+1])/2) - ((dt/2.0*m_cellSize)*(m_hu[i+1]-m_hu[i-1]));
        m_hu[i] = ((m_hu[i-1]+m_hu[i+1])/2) - ((dt/2.0*m_cellSize)*(fQRight-fQLeft));
    }

}



