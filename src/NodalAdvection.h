//
// Created by lisa on 18.10.17.
//

#ifndef SWE1D_NODALADVECTION_H
#define SWE1D_NODALADVECTION_H

#include "types.h"

class NodalAdvection{
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
    vect m_uNetUpdatesLeft;
    vect m_uNetUpdatesRight;
    //u
    vecu m_u;
    //Time derivative of u
    vecu m_ut;
    int m_size;

    //Time variables
    T t0 = 0.0;
    //Advection velocity
    T a;

    T m_cellSize;
public:
    NodalAdvection(float a_def, vecu h, unsigned int size, unsigned int cellsize):
            a(a_def),
            m_size(size),
            m_u(h.size()),
            m_ut(h.size()),
            m_uNetUpdatesRight(size,0.0),
            m_uNetUpdatesLeft(size,0.0),
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
                    m_u[i].u0 = h[i].u0;
                    m_u[i].u1 = h[i].u1;
                    m_ut.at(i).u0 = 0.0;
                    m_ut.at(i).u1 = 0.0;
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

    ~NodalAdvection(){
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
    T compute(T t);
    T computeLocalLaxFriedrichsFluxes(T t);
    void updateUnknowns(T t);

    void computeTimeDerivative();

    void computeEulerStep(T delta_t);

    vecu setH();

    vect getExactSolution(T t);
};

#endif //SWE1D_NODALADVECTION_H
