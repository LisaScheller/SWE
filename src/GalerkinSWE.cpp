//
// Created by lisa on 17.11.17.
//

#include "GalerkinSWE.h"
#include <cmath>
#include <algorithm>
#include "types.h"
#include "NodalAdvection.h"
#include <vector>


q matMult(q a, q b){
    q res;
    res.h.u0 = a.h.u0*b.h.u0+a.h.u1*b.hu.u0;
    res.h.u1 = a.h.u0*b.h.u1+a.h.u1*b.hu.u1;
    res.hu.u0 = a.hu.u0*b.h.u0+a.hu.u1*b.hu.u0;
    res.hu.u1 = a.hu.u0*b.h.u1+a.hu.u1*b.hu.u1;
    return res;
}
u matVecMul(vect a,u b){
    //Mutliplication of 2x2 Matrix with 2x1 vector
    u res;
    res.u0 =0;
    res.u1 = 0;
    res.u0 = a[0]*b.u0+a[1]*b.u1;
    res.u1 = a[2]*b.u0+a[3]*b.u1;
    return res;
}
T GalerkinSWE::computeLocalLaxFriedrichsFluxes(T t) {
    T maxWaveSpeed = 0.0;
    T deltaT = t - t0;
    t0 = t;
    vecu hu (m_size+2);
    vecu h (m_size+2);
    for(int j = 0; j<m_size+2; j++){
        h.at(j).u0 = m_q.at(j).h.u0;
        h.at(j).u1 = m_q.at(j).h.u1;
        hu.at(j).u0 = m_q.at(j).hu.u0;
        hu.at(j).u1 = m_q.at(j).hu.u1;
    }
    //Loop over all intervalls
    for(unsigned int i = 1; i<m_size+1; i++){
        //Factor a for local Lax Friedrichs Method (aLeft = a_i-1/2, aRight = a_i+1/2
        u aLeft;
        u aRight;
        if (i == 200){
            int k = 0;
        }

        //Eigenvalues
        aLeft.u0  = std::max(abs(hu.at(i-1).u0/h.at(i-1).u0) + sqrt(g*h.at(i-1).u0),abs(hu.at(i).u0/h.at(i).u0) + sqrt(g*h.at(i).u0));
        aRight.u0 = std::max(abs(hu.at(i+1).u0/h.at(i+1).u0) + sqrt(g*h.at(i+1).u0),abs(hu.at(i).u0/h.at(i).u0) + sqrt(g*h.at(i).u0));
        aLeft.u1  = std::max(abs(hu.at(i-1).u1/h.at(i-1).u1) + sqrt(g*h.at(i-1).u1),abs(hu.at(i).u1/h.at(i).u1) + sqrt(g*h.at(i).u1));
        aRight.u1 = std::max(abs(hu.at(i+1).u1/h.at(i+1).u1) + sqrt(g*h.at(i+1).u1),abs(hu.at(i).u1/h.at(i).u1) + sqrt(g*h.at(i).u1));

        //T aL = std::max(aLeft.u0,aLeft.u1);
        //T aR = std::max(aRight.u0, aRight.u1);
        T aL = aLeft.u0;
        T aR = aRight.u1;
        //Compute fluxes for h and hu
        m_hNetUpdatesLeft.at(i) = 0.5*(hu.at(i-1).u1+hu.at(i).u0-(aL*(h.at(i).u0-h.at(i-1).u1)));
        m_hNetUpdatesRight.at(i) = 0.5*(hu.at(i).u1+hu.at(i+1).u0-(aR*(h.at(i+1).u0-h.at(i).u1)));
        T fQLeft = ((hu.at(i-1).u1*hu.at(i-1).u1)/h.at(i-1).u1)+(0.5*g*(h.at(i-1).u1*h.at(i-1).u1));
        T fQMiddleL = ((hu.at(i).u0*hu.at(i).u0)/h.at(i).u0)+(0.5*g*(h.at(i).u0*h.at(i).u0));
        T fQMiddleR = ((hu.at(i).u1*hu.at(i).u1)/h.at(i).u1)+(0.5*g*(h.at(i).u1*h.at(i).u1));
        T fQRight = ((hu.at(i+1).u0*hu.at(i+1).u0)/h.at(i+1).u0)+(0.5*g*(h.at(i+1).u0*h.at(i+1).u0));
        m_huNetUpdatesLeft.at(i) = 0.5*(fQLeft+fQMiddleL-(aL*(hu.at(i).u0-hu.at(i-1).u1)));
        m_huNetUpdatesRight.at(i) = 0.5*(fQMiddleR+fQRight-(aR*(hu.at(i+1).u0-hu.at(i).u1)));
        // Update maxWaveSpeed
        maxWaveSpeed = std::max(std::max(maxWaveSpeed,aL),aR);

    }
    // Compute CFL condition (delta_t <= c/a*(delta_x/2)) (DG book p. 97)
    T c = 0.04;
    T maxTimeStep = m_cellSize*(c/maxWaveSpeed);

    return maxTimeStep;
}

/*T GalerkinSWE::computeLocalLaxFriedrichsFluxes(T t){
    T maxWaveSpeed = 0.0;
    T deltaT = t - t0;
    t0 = t;
    //Loop over all intervalls
    for(unsigned int i = 1; i<m_size; i++){
        //Factor a for local Lax Friedrichs Method (aLeft = a_i-1/2, aRight = a_i+1/2
        u aLeft;
        u aRight;

        aLeft.u0  = std::max(abs(m_q.at(i-1).hu.u0/m_q.at(i-1).h.u0) + sqrt(g*m_q.at(i-1).h.u0),abs(m_q.at(i).hu.u0/m_q.at(i).h.u0) + sqrt(g*m_q.at(i).h.u0));
        aLeft.u1  = std::max(abs(m_q.at(i-1).hu.u1/m_q.at(i-1).h.u1) + sqrt(g*m_q.at(i-1).h.u1),abs(m_q.at(i).hu.u1/m_q.at(i).h.u1) + sqrt(g*m_q.at(i).h.u1));
        aRight.u0 = std::max(abs(m_q.at(i+1).hu.u0/m_q.at(i+1).h.u0) + sqrt(g*m_q.at(i+1).h.u0),abs(m_q.at(i).hu.u0/m_q.at(i).h.u0) + sqrt(g*m_q.at(i).h.u0));
        aRight.u1 = std::max(abs(m_q.at(i+1).hu.u1/m_q.at(i+1).h.u1) + sqrt(g*m_q.at(i+1).h.u1),abs(m_q.at(i).hu.u1/m_q.at(i).h.u1) + sqrt(g*m_q.at(i).h.u1));

        //Compute fluxes for h,hu
        //m_u[i][1]+m_u[i+1][0] ... m_u[i][1]  -m_u[i+1][0]
        m_hNetUpdatesLeft.at(i) = 0.5*(m_q.at(i-1).hu.u1+m_q.at(i).hu.u0-(aLeft.u0*(m_q.at(i).h.u0-m_q.at(i-1).h.u1)));
        m_hNetUpdatesRight.at(i) = 0.5*(m_q.at(i).hu.u1+m_q.at(i+1).hu.u0-(aRight.u1*(m_q.at(i+1).h.u0-m_q.at(i).h.u1)));

        T fQLeft = ((m_q.at(i-1).hu.u1*m_q.at(i-1).hu.u1)/m_q.at(i-1).h.u1)+(0.5*g*(m_q.at(i-1).h.u1*m_q.at(i-1).h.u1));
        T fQMiddleL = ((m_q.at(i).hu.u0*m_q.at(i).hu.u0)/m_q.at(i).h.u0)+(0.5*g*(m_q.at(i).h.u0*m_q.at(i).h.u0));
        T fQMiddleR = ((m_q.at(i).hu.u1*m_q.at(i).hu.u1)/m_q.at(i).h.u1)+(0.5*g*(m_q.at(i).h.u1*m_q.at(i).h.u1));
        T fQRight = ((m_q.at(i+1).hu.u0*m_q.at(i+1).hu.u0)/m_q.at(i+1).h.u0)+(0.5*g*(m_q.at(i+1).h.u0*m_q.at(i+1).h.u0));

        m_huNetUpdatesLeft.at(i) = 0.5*(fQLeft+fQMiddleL-(aLeft.u0*(m_q.at(i).hu.u0-m_q.at(i-1).hu.u1)));
        m_huNetUpdatesRight.at(i) = 0.5*(fQMiddleR+fQRight-(aRight.u1*(m_q.at(i+1).hu.u0-m_q.at(i).hu.u1)));

        // Update maxWaveSpeed
        maxWaveSpeed = std::max(std::max(maxWaveSpeed,aLeft.u0),aRight.u1);

    }
    // Compute CFL condition (delta_t = delta_x/(a*3))

    T maxTimeStep = m_cellSize/(std::abs(maxWaveSpeed)*3);

    return maxTimeStep;
}*/

void GalerkinSWE::computeTimeDerivative(){
    vecu m_qtH(m_size+2);
    vecu m_qtHu(m_size+2);
    for(unsigned int i = 1; i<m_size; i++) {
        //Compute S *(f(q))
        vect rightTermFirstH(2,0.0);
        vect rightTermFirstHu(2,0.0);

        for (unsigned int j = 0; j<2; j++){
            if (j==0){
                rightTermFirstH.at(0) += s.at(j)*m_q.at(i).hu.u0;
                rightTermFirstH.at(1) += s.at(2+j)*m_q.at(i).hu.u0;
                /*rightTermFirstHu.at(0) += s.at(j)*m_q.at(i).hu.u0*m_q.at(i).hu.u0;
                rightTermFirstHu.at(1) += s.at(2+j)*m_q.at(i).hu.u1*m_q.at(i).hu.u1;*/
                //rightTermFirstHu0.at(1) += s.at(2+j)*(((m_q.at(i).hu.u1*m_q.at(i).hu.u1)/m_q.at(i).h.u1)+(0.5*g*m_q.at(i).h.u1*m_q.at(i).h.u1))*m_q.at(i).hu.u1;
                //rightTermFirstHu0.at(0) += s.at(j)*(((m_q.at(i).hu.u0*m_q.at(i).hu.u0)/m_q.at(i).h.u0)+(0.5*g*m_q.at(i).h.u0*m_q.at(i).h.u0))*m_q.at(i).hu.u0;
                rightTermFirstHu.at(0) += s.at(j)*(((m_q.at(i).hu.u0*m_q.at(i).hu.u0)/m_q.at(i).h.u0)+(0.5*g*m_q.at(i).h.u0*m_q.at(i).h.u0));
                rightTermFirstHu.at(1) += s.at(2+j)*(((m_q.at(i).hu.u0*m_q.at(i).hu.u0)/m_q.at(i).h.u0)+(0.5*g*m_q.at(i).h.u0*m_q.at(i).h.u0));;

            }
            else {
                /*rightTermFirstH.at(1) += s.at(j+2)*(((m_q.at(i).hu.u0*m_q.at(i).hu.u0)/m_q.at(i).h.u0)+(0.5*g*m_q.at(i).h.u0*m_q.at(i).h.u0))*m_q.at(i).h.u0;
                rightTermFirstH.at(0) += s.at(j)*(((m_q.at(i).hu.u1*m_q.at(i).hu.u1)/m_q.at(i).h.u1)+(0.5*g*m_q.at(i).h.u1*m_q.at(i).h.u1))*m_q.at(i).h.u1;*/
                //rightTermFirstH.at(0) += s.at(j)*m_q.at(i).hu.u1*m_q.at(i).h.u1;
                rightTermFirstH.at(0) += s.at(j)*m_q.at(i).hu.u1;
                rightTermFirstH.at(1) += s.at(j+2)*m_q.at(i).hu.u1;
                rightTermFirstHu.at(0) += s.at(j)*(((m_q.at(i).hu.u1*m_q.at(i).hu.u1)/m_q.at(i).h.u1)+(0.5*g*m_q.at(i).h.u1*m_q.at(i).h.u1));
                rightTermFirstHu.at(1) += s.at(j+2)*(((m_q.at(i).hu.u1*m_q.at(i).hu.u1)/m_q.at(i).h.u1)+(0.5*g*m_q.at(i).h.u1*m_q.at(i).h.u1));
                //rightTermFirstH.at(1) += s.at(2+j)*m_q.at(i).hu.u0*m_q.at(i).h.u0;
                //rightTermFirstHu.at(1) += s.at(2+j)*(((m_q.at(i).hu.u0*m_q.at(i).hu.u0)/m_q.at(i).h.u0)+(0.5*g*m_q.at(i).h.u0*m_q.at(i).h.u0))*m_q.at(i).hu.u0;
                //rightTermFirstHu.at(0) += s.at(j)*(((m_q.at(i).hu.u1*m_q.at(i).hu.u1)/m_q.at(i).h.u1)+(0.5*g*m_q.at(i).h.u1*m_q.at(i).h.u1))*m_q.at(i).hu.u1;
            }

        }
        //s.o.
        //Compute n0*netupdatesRight*delta_x
        vect rightTermSecondH(2,0.0);
        //vect rightTermSecondH0(2,0.0);
        //vect rightTermSecondH1(2,0.0);
        vect rightTermSecondHu(2,0.0);

        for(unsigned int k = 0; k<2; k++){
            /*rightTermSecondH.at(0) += n0.at(k)*m_hNetUpdatesLeft.at(i);
            rightTermSecondH.at(1) += n0.at(2+k)*m_hNetUpdatesLeft.at(i);
            rightTermSecondHu.at(0) += n0.at(k)*m_huNetUpdatesLeft.at(i);
            rightTermSecondHu.at(1) += n0.at(2+k)*m_huNetUpdatesLeft.at(i);*/
            rightTermSecondH.at(0) = m_hNetUpdatesLeft.at(i);
            rightTermSecondH.at(1) = 0.0;
            rightTermSecondHu.at(0) = m_huNetUpdatesLeft.at(i);
            rightTermSecondHu.at(1) = 0.0;
        }

        if (i  == 49){
            int y = 0;
        }
        //Compute n1*netUpdatesLeft*delta_x
        vect rightTermThirdH(2,0.0);
        vect rightTermThirdHu(2,0.0);

        for(unsigned int l = 0; l<2; l++){
            /*rightTermThirdH.at(0) += n1.at(l)*m_hNetUpdatesRight.at(i);
            rightTermThirdH.at(1) += n1.at(2+l)*m_hNetUpdatesRight.at(i);
            rightTermThirdHu.at(0) += n1.at(l)*m_huNetUpdatesRight.at(i);
            rightTermThirdHu.at(1) += n1.at(2+l)*m_huNetUpdatesRight.at(i);*/
            rightTermThirdH.at(0) = 0.0;
            rightTermThirdH.at(1) = m_hNetUpdatesRight.at(i);
            rightTermThirdHu.at(0) = 0.0;
            rightTermThirdHu.at(1) = m_huNetUpdatesRight.at(i);
        }
        //Put right side first+second-third together

        u rightTermH;
        u rightTermHu;

        rightTermH.u0 = (rightTermFirstH.at(0)+rightTermSecondH.at(0)-rightTermThirdH.at(0))/m_cellSize;
        rightTermH.u1 = (rightTermFirstH.at(1)+rightTermSecondH.at(1)-rightTermThirdH.at(1))/m_cellSize;
        rightTermHu.u0 = (rightTermFirstHu.at(0)+rightTermSecondHu.at(0)-rightTermThirdHu.at(0))/m_cellSize;
        rightTermHu.u1 = (rightTermFirstHu.at(1)+rightTermSecondHu.at(1)-rightTermThirdHu.at(1))/m_cellSize;


        /*q m_inv;
        m_inv.h.u0 = inv_m.at(0);
        m_inv.h.u1 = inv_m.at(1);
        m_inv.hu.u0 = inv_m.at(2);
        m_inv.hu.u1 = inv_m.at(3);*/
        //m_qtH.at(i) = matVecMul(inv_m, rightTermH);
        //m_qtHu.at(i)=matVecMul(inv_m,rightTermHu);
        //m_qt.at(i) = matMult(m_inv, rightTerm);
        //m_qt.at(i).h = m_qtH.at(i);
        //m_qt.at(i).hu = m_qtHu.at(i);
        u ht;
        u hut;
        ht = matVecMul(inv_m, rightTermH);
        hut = matVecMul(inv_m, rightTermHu);
        m_qt.at(i).h.u0 = ht.u0;
        m_qt.at(i).h.u1 = ht.u1;
        m_qt.at(i).hu.u0 = hut.u0;
        m_qt.at(i).hu.u1 = hut.u1;

        /*m_qt.at(i).h.u0 = rightTermH.u0;
        m_qt.at(i).h.u1 = rightTermH.u1;
        m_qt.at(i).hu.u0 = rightTermHu.u0;
        m_qt.at(i).hu.u1 = rightTermHu.u1;*/

    }
    m_qt.at(0).h.u0 = m_qt.at(1).h.u0;
    m_qt.at(0).hu.u0 = m_qt.at(1).hu.u0;
    m_qt.at(0).h.u1 = m_qt.at(1).h.u1;
    m_qt.at(0).hu.u1 = m_qt.at(1).hu.u1;

    m_qt.at(m_size+1).h.u0 = m_qt.at(m_size).h.u0;
    m_qt.at(m_size+1).hu.u0 = m_qt.at(m_size).hu.u0;
    m_qt.at(m_size+1).h.u1 = m_qt.at(m_size).h.u1;
    m_qt.at(m_size+1).hu.u1 = m_qt.at(m_size).hu.u1;
}

void GalerkinSWE::computeEulerStep(T delta_t){
    q q0;
    q0.h.u0 = 0.0;
    q0.hu.u0 = 0.0;
    q0.h.u1 = 0.0;
    q0.hu.u1 = 0.0;
    vecq tmp(m_size+2,q0);
    for(unsigned int i = 0; i<m_size+2; i++){
        if (i==49){
            int o = 0;
        }
        tmp.at(i).h.u0 = m_q.at(i).h.u0 + (delta_t*m_qt.at(i).h.u0);
        tmp.at(i).h.u1 = m_q.at(i).h.u1 + (delta_t*m_qt.at(i).h.u1);
        tmp.at(i).hu.u0=m_q.at(i).hu.u0+(delta_t*m_qt.at(i).hu.u0);
        tmp.at(i).hu.u1=m_q.at(i).hu.u1+(delta_t*m_qt.at(i).hu.u1);
    }
    for(unsigned int j = 0; j<m_size+2; j++){
        m_q.at(j).h.u0 = tmp.at(j).h.u0;
        m_q.at(j).hu.u0 = tmp.at(j).hu.u0;
        m_q.at(j).h.u1 = tmp.at(j).h.u1;
        m_q.at(j).hu.u1 = tmp.at(j).hu.u1;
    }
    /*for (unsigned int j = 0; j<m_size+2; j++) {
        if (tmp.at(j).u0 > 0.1 && tmp.at(j).u1 > 0.1) {
            m_u.at(j).u0 = tmp.at(j).u0;
            m_u.at(j).u1 = tmp.at(j).u1;
        }
        else if (tmp.at(j).u0 <= 0.1 && tmp.at(j).u1 > 0.1){
            m_u.at(j).u1 = tmp.at(j).u1;
        }
        else if (tmp.at(j).u0 > 0.1 && tmp.at(j).u1 <= 0.1){
            m_u.at(j).u0 = tmp.at(j).u0;
        }
    }*/
}

void GalerkinSWE::setBoundaryConditions() {
    m_q.at(0).h = m_q.at(1).h; m_q.at(m_size+1).h = m_q.at(m_size).h;
    m_q.at(0).hu = m_q.at(1).hu; m_q.at(m_size+1).hu = m_q.at(m_size).hu;

}

vecq GalerkinSWE::setQ(){
    q q0;
    q0.h.u0 = 0.0;
    q0.hu.u0 = 0.0;
    q0.h.u1 = 0.0;
    q0.hu.u1 = 0.0;
    vecq result(m_q.size(), q0);
    for (unsigned int i = 0; i<m_q.size(); i++){
        result.at(i) = m_q.at(i);
    }
    return result;
}

