//
// Created by lisa on 17.11.17.
//

#ifndef SWE1D_GALERKINSWE_H
#define SWE1D_GALERKINSWE_H
#include "types.h"

class GalerkinSWE{
private:
    //Stiffness Matrix
    vect s;
    //Mass matrix
    vect m;
    //inverse mass matrix
    vect inv_m;
    //First N-Matrix
    vect n0;
    //Second N-Matrix
    vect n1;
    //Left and right fluxes
    vect m_hNetUpdatesLeft;
    vect m_hNetUpdatesRight;

    vect m_huNetUpdatesLeft;
    vect m_huNetUpdatesRight;
    //u
    vecq m_q;
    vecq aq;
    vecq m_q_bef;
    //Time derivative of u
    vecq m_qt;
    int m_size;

    T m_cellSize;

    T t0 = 0.0;
    double g = 9.80665;

public:
    GalerkinSWE(vecq q, unsigned int size, T cellsize):
            m_size(size),
            m_q(q.size()),
            m_q_bef(q.size()),
            aq(q.size()),
            m_qt(q.size()),
            m_hNetUpdatesRight(size+2,0.0),
            m_hNetUpdatesLeft(size+2,0.0),
            m_huNetUpdatesRight(size+2,0.0),
            m_huNetUpdatesLeft(size+2                                                                                                                                                                             ,0.0),
            m_cellSize(cellsize),
            s(4),
            m(4),
            inv_m(4),
            n0(4),
            n1(4)
    {   //Allocate matrices

        //vect m_uNetUpdatesLeft(size+1);
        //vect m_uNetUpdatesRight(size+1);

        /*for (unsigned int j = 0; j< size+1; j++){
            m_uNetUpdatesLeft[j]=0.0;
            m_uNetUpdatesRight[j]=0.0;
        }*/
        //vecu m_ut(size+2);
        //vecu m_u(size+2);
        for (unsigned int i = 0; i<size+2; i++){
            m_q.at(i).h.u0 = q.at(i).h.u0;
            m_q.at(i).h.u1 = q.at(i).h.u1;
            m_q.at(i).hu.u0 = q.at(i).hu.u0;
            m_q.at(i).hu.u1 = q.at(i).hu.u1;
            m_q_bef.at(i).h.u0 = q.at(i).h.u0;
            m_q_bef.at(i).h.u1 = q.at(i).h.u1;
            m_q_bef.at(i).hu.u0 = q.at(i).hu.u0;
            m_q_bef.at(i).hu.u1 = q.at(i).hu.u1;
            m_qt.at(i).h.u0 = 0.0;
            m_qt.at(i).h.u1 = 0.0;
            m_qt.at(i).hu.u0 = 0.0;
            m_qt.at(i).hu.u1 = 0.0;
        }


        //Define entries of matrices (0=00, 1 = 01, 2 = 10, 3 = 11)
        s[0] = -0.5;
        s[1] = -0.5;
        s[2] = 0.5;
        s[3] = 0.5;

        m[0] = 1.0f/3.0;
        m[1] = 1.0f/6.0;
        m[2] = 1.0f/6.0;
        m[3] = 1.0f/3.0;

        inv_m[0] = 4.0;
        inv_m[1] = -2.0;
        inv_m[2] = -2.0;
        inv_m[3] = 4.0;

        n0[0] = 1.0;
        n0[1] = 0.0;
        n0[2] = 0.0;
        n0[3] = 0.0;

        n1[0] = 0.0;
        n1[1] = 0.0;
        n1[2] = 0.0;
        n1[3] = 1.0;

    }

    ~GalerkinSWE(){
        /* //Free allocated memory
         delete [] s;
         delete [] m;
         delete [] inv_m;
         delete [] n0;
         delete [] n1;
         delete [] m_uNetUpdatesLeft;
         delete [] m_uNetUpdatesRight;
         delete [] m_ut;*/
    }
    //Methods
    void setBoundaryConditions();

    T computeLocalLaxFriedrichsFluxes(T t);
    void updateUnknowns(T t);

    void computeTimeDerivative();

    void computeEulerStep(T delta_t);

    vecq setQ();

    vect getExactSolution(T t);

    vecq getAnalyticalSolution(T dt);

    void computeHalfEulerStep(T delta_t);
};
#endif //SWE1D_GALERKINSWE_H
