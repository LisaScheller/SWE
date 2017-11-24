//
// Created by lisa on 05.09.17.
//

#include <iostream>
#include "types.h"
#include "LaxFriedrichsSolver.h"
#include "WavePropagation.h"
#include "cmath"
#include "algorithm"

T tPrev = 0;
T LaxFriedrichsSolver::computeLaxFriedrichsFlux() {
    T maxWaveSpeed = 0.0;
    //Loop over all edges
    for (unsigned int i = 1; i<m_size+2; i++){
        T maxEdgeSpeed;

        // Compute net updates (1-dim SWE, Leveque S. 255)
        m_hNetUpdatesRight.at(i) = m_hu.at(i);
        m_huNetUpdatesRight.at(i) = (m_hu.at(i)*m_hu.at(i))/m_h.at(i)+(0.5*g*(m_h.at(i)*m_h.at(i)));
        m_hNetUpdatesLeft.at(i) = m_hu.at(i-1);
        m_huNetUpdatesLeft.at(i) = (m_hu.at(i-1)*m_hu.at(i-1))/m_h.at(i-1)+(0.5*g*(m_h.at(i-1)*m_h.at(i-1)));

        //Compute edge speed

        if(m_hu.at(i-1) < 0.0){
            maxEdgeSpeed= -m_hu.at(i-1)+std::sqrt(g)*m_h.at(i-1);
        }else{
            maxEdgeSpeed=  m_hu.at(i-1)+std::sqrt(g)*m_h.at(i-1);
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
    for (unsigned int i = 1; i<m_size+1; i++){
        T maxEdgeSpeed;

        // Compute net updates (1-dim SWE, Leveque S. 255)
        m_hNetUpdatesLeft.at(i) = (0.5*(m_hu.at(i-1)+m_hu.at(i)))-((m_cellSize/2*deltaT)*(m_h.at(i)-m_h.at(i-1)));
        m_hNetUpdatesRight.at(i) = (0.5*(m_hu.at(i)+m_hu.at(i+1)))-((m_cellSize/2*deltaT)*(m_h.at(i+1)-m_h.at(i)));
        T fQLeft = ((m_hu.at(i-1)*m_hu.at(i-1))/m_h.at(i-1))+(0.5*g*(m_h.at(i-1)*m_h.at(i-1)));
        T fQMiddle = ((m_hu.at(i)*m_hu.at(i))/m_h.at(i))+(0.5*g*(m_h.at(i)*m_h.at(i)));
        T fQRight = ((m_hu.at(i+1)*m_hu.at(i+1))/m_h.at(i+1))+(0.5*g*(m_h.at(i+1)*m_h.at(i+1)));
        m_huNetUpdatesLeft.at(i) = (0.5*(fQLeft+fQMiddle))-((m_cellSize/2*deltaT)*(m_hu.at(i)-m_hu.at(i-1)));
        m_huNetUpdatesRight.at(i) = (0.5*(fQMiddle+fQRight))-((m_cellSize/2*deltaT)*(m_hu.at(i+1)-m_hu.at(i)));


        //Compute edge speed
	T aLeft  = std::max(abs(m_hu.at(i-1)/m_h.at(i-1)) + sqrt(g*m_h.at(i-1)),abs(m_hu.at(i)/m_h.at(i)) + sqrt(g*m_h.at(i)));
	T aRight = std::max(abs(m_hu.at(i+1)/m_h.at(i+1)) + sqrt(g*m_h.at(i+1)),abs(m_hu.at(i)/m_h.at(i)) + sqrt(g*m_h.at(i)));

    maxWaveSpeed = std::max(std::max(maxWaveSpeed, aLeft),aRight);
	//maxWaveSpeed = std::max(std::max(maxWaveSpeed,aLeft),aRight);

    }
    // Compute CFL condition (with Courant number 0.4)

    T maxTimeStep = m_cellSize/maxWaveSpeed * .4;

    return maxTimeStep;
}

T LaxFriedrichsSolver::computeLocalLaxFriedrichsFlux(T t){
    T maxWaveSpeed = 0.0;
    T deltaT = t - tPrev;
    tPrev = t;
    //Loop over all edges
    for (unsigned int i = 1; i<m_size+1; i++){
        if (i == 49){
            int o = 3;
        }
        //Factor a for local Lax Friedrichs Method (aLeft = a_i-1/2, aRight = a_i+1/2
        T aLeft;
        T aRight;

	    aLeft  = std::max(abs(m_hu.at(i-1)/m_h.at(i-1)) + sqrt(g*m_h.at(i-1)),abs(m_hu.at(i)/m_h.at(i)) + sqrt(g*m_h.at(i)));
	    aRight = std::max(abs(m_hu.at(i+1)/m_h.at(i+1)) + sqrt(g*m_h.at(i+1)),abs(m_hu.at(i)/m_h.at(i)) + sqrt(g*m_h.at(i)));

        // Compute fluxes for h and hu
        m_hNetUpdatesLeft.at(i) = 0.5*(m_hu.at(i-1)+m_hu.at(i)-(aLeft*(m_h.at(i)-m_h.at(i-1))));
        m_hNetUpdatesRight.at(i) = 0.5*(m_hu.at(i)+m_hu.at(i+1)-(aRight*(m_h.at(i+1)-m_h.at(i))));
        T fQLeft = ((m_hu.at(i-1)*m_hu.at(i-1))/m_h.at(i-1))+(0.5*g*(m_h.at(i-1)*m_h.at(i-1)));
        T fQMiddle = ((m_hu.at(i)*m_hu.at(i))/m_h.at(i))+(0.5*g*(m_h.at(i)*m_h.at(i)));
        T fQRight = ((m_hu.at(i+1)*m_hu.at(i+1))/m_h.at(i+1))+(0.5*g*(m_h.at(i+1)*m_h.at(i+1)));
        m_huNetUpdatesLeft[i] = 0.5*(fQLeft+fQMiddle-(aLeft*(m_hu.at(i)-m_hu.at(i-1))));
        m_huNetUpdatesRight[i] = 0.5*(fQMiddle+fQRight-(aRight*(m_hu.at(i+1)-m_hu.at(i))));

        // Update maxWaveSpeed
	    maxWaveSpeed = std::max(std::max(maxWaveSpeed,aLeft),aRight);
    }
    // Compute CFL condition (with Courant number 0.04)
    T maxTimeStep = m_cellSize/maxWaveSpeed * 0.4;

    return maxTimeStep;
}


void LaxFriedrichsSolver::updateUnknownsLaxFriedrichs(T dt){
    //Original method using the calculated Fluxes to compute h and hu at i
    //Loop over all inner cells
    //Leveque S. 71
    for (unsigned int i = 1; i < m_size+1; i++){
        m_h.at(i) = ((m_h.at(i-1)+m_h.at(i+1))/2) + dt/2*m_cellSize*(m_hNetUpdatesLeft.at(i)-m_hNetUpdatesRight.at(i));
        m_hu.at(i) = ((m_hu.at(i-1)+m_hu.at(i+1))/2) + dt/2*m_cellSize*(m_huNetUpdatesLeft.at(i)-m_huNetUpdatesRight.at(i));
}

}

vecu LaxFriedrichsSolver::updateUnknownsLaxFriedrichs2(T dt){
    //Alternative method using a combination of formula 4.21 (Leveque p. 71) for the fluxes
    //and formulas 4.6 and 4.7 (Leveque p. 66) to compute h and hu
    //Since the computation of the fluxes requires dt, the time of the current timestep
    //has to be given to the method calculating the fluxes
    //Loop over all inner cells
    vecu res(m_size+2);
    for (unsigned int i = 1; i < m_size+1; i++){
        /*m_h[i] = ((m_h[i-1]+m_h[i+1])/2) + ((dt/m_cellSize)*(m_hNetUpdatesRight[i]-m_hNetUpdatesLeft[i]));
        m_hu[i] = ((m_hu[i-1]+m_hu[i+1])/2) + ((dt/m_cellSize)*(m_huNetUpdatesRight[i]-m_huNetUpdatesLeft[i]));
        */
        m_h.at(i) = m_h.at(i) - ((dt/m_cellSize)*(m_hNetUpdatesRight.at(i)-m_hNetUpdatesLeft.at(i)));
        m_hu.at(i) = m_hu.at(i) - ((dt/m_cellSize)*(m_huNetUpdatesRight.at(i)-m_huNetUpdatesLeft.at(i)));
        res.at(i).u0 = m_h.at(i);
        res.at(i).u1 = m_hu.at(i);
    }
    return res;
}

vecu LaxFriedrichsSolver::updateUnknownsLocalLaxFriedrichs(T dt){
    vecu res(m_size+2);
    for (unsigned int i = 1; i < m_size+1; i++){
        if (i==50){
            int y = 3;
        }
        m_h.at(i) = m_h.at(i) - ((dt/m_cellSize)*(m_hNetUpdatesRight.at(i)-m_hNetUpdatesLeft.at(i)));
        m_hu.at(i) = m_hu.at(i) - ((dt/m_cellSize)*(m_huNetUpdatesRight.at(i)-m_huNetUpdatesLeft.at(i)));
        res.at(i).u0 = m_h.at(i);
        res.at(i).u1 = m_hu.at(i);
    }
    return res;

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
    T factor = m_size/100;
    for (unsigned int i =0; i<m_size; i++){
        //T x = i-50.0; //oder i+(i+1)/2?
        int l = i/factor;
        int r = (i+1)/factor;
        T x = ((l+r)/2)-50.0;
        if(x/t<=lambda_1_l){
            a_h.at(i)=h_l;
            a_hu.at(i)=a_h.at(i)*u_l;
        }
        else if(x/t>=s){
            a_h.at(i)=h_r;
            a_hu.at(i)=a_h.at(i)*u_r;
        }
        else if (x/t>lambda_1_l && x/t<lambda_1_m){
            a_h.at(i)=0.01133018014*(20.77239996-x/t)*(20.77239996-x/t);
            a_hu.at(i)=a_h.at(i)*(u_m-2.0*(sqrt(g*h_m)-sqrt(g*a_h.at(i))));
        }
        else {
            a_h.at(i)=h_m;
            a_hu.at(i)=a_h.at(i)*u_m;
        }
    }

}

vecu LaxFriedrichsSolver::getAnalyticalSolution(T dt){
    //Values for specific case h_l=11, h_r=10, u_l=u_r=0
    //Computed with maple
    T t = tPrev +dt;

    T h_l=3.0;
    T u_l=0.0;
    T h_r=1.0;
    T u_r=0.0;
    /*T s=10.26901816;
    T u_m=0.48339953;
    T h_m=10.49398975;*/
    T u_m=2.332542824;
    T h_m=1.848576603;
    T c_m=std::sqrt(g*h_m);
    T lambda_1_l=u_l-sqrt(g*h_l);
    T lambda_1_m=u_m-sqrt(g*h_m);
    T factor = m_size/100.0;
    T c1 = std::sqrt(g*h_l);
    T a = -5.423991150;
    T b = -1.92517691;
    T s = 5.081313902;
    for (unsigned int i =0; i<m_size+2; i++){
        if (i == 47){
            int test = 38;
        }
        //T x = i-50.0; //oder i+(i+1)/2?
        //int l = i/factor;
        //int r = (i+1)/factor;
        //T x = ((l+r)/2)-50.0;
        T x = (i/factor) - 50.0;
        T sigma = x/t;
        if (sigma <= a/*-10.38615232*/){
            a_h.at(i) = h_l;
            a_hu.at(i) = a_h.at(i)*u_l;
        }
        else if (sigma > a/*-10.38615232*/ && sigma <= b/*-9.66105637*/){
            T test = c1-(sigma/2);
            T t2 = test*test;
            a_h.at(i)=(4/(9*g))*((c1-(sigma/2))*(c1-(sigma/2)));
            a_hu.at(i)=a_h.at(i)*((2/3)*(sigma+c1));

        }
        else if (sigma>b/*-9.66105637*/ && sigma <= s){
            a_h.at(i) = h_m;
            a_hu.at(i) = a_h.at(i)*(2*(c1-c_m));
        }
        else if (sigma>s){//1.012274368){
           a_h.at(i) = h_r;
           a_hu.at(i) = a_h.at(i)*u_r;
        }
        /*if((x/t)<=lambda_1_l){
            a_h.at(i)=h_l;
            //a_hu.at(i) = a_hu.at(i)/a_h.at(i);
            a_hu.at(i)=a_h.at(i)*u_l;
        }
        else if((x/t)>=s){
            a_h.at(i)=h_r;
            a_hu.at(i)=a_h.at(i)*u_r;
        }
        else if ((x/t)>lambda_1_l && (x/t)<lambda_1_m){
            T c3 = (1/9)*(2*c1-(x/t))*(2*c1-(x/t));
            a_h.at(i)= (c3*c3)/g;
            T u3 = (2/3)*(c1+(x/t));
            a_hu.at(i) = a_h.at(i)*u3;
            /*a_h.at(i)=(1/(9*g))*((2*std::sqrt(g*h_l))-(x/t))*((2*std::sqrt(g*h_l))-(x/t));
            //a_hu.at(i)=(u_m-(2*(std::sqrt(g*h_m)-std::sqrt(g*a_h.at(i)))))*a_h.at(i);
            a_hu.at(i)=(u_l-(2*(std::sqrt(g*10.51700522)-std::sqrt(g*a_h.at(i)))))*a_h.at(i);
            //a_hu.at(i)=u_m*a_h.at(i);
            //a_h.at(i)=0.01133018014*(20.77239996-(x/t))*(20.77239996-(x/t));
            //a_hu.at(i)=a_h.at(i)*(u_m-2.0*(sqrt(g*h_m)-sqrt(g*a_h.at(i))));
        }
        else {
            a_h.at(i)=h_m;
            a_hu.at(i)=a_h.at(i)*u_m;
        }*/
    }
    a_h.at(m_size+1) = a_h.at(m_size);
    a_hu.at(m_size+1) = a_hu.at(m_size);
    vecu res(m_size+2);
    for(unsigned int j = 0; j < m_size+2; j++){
        res.at(j).u0 = a_h.at(j);
        res.at(j).u1 = a_hu.at(j);
    }
    return res;
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



