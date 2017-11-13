//
// Created by lisa on 18.10.17.
//

#ifndef SWE1D_NODALADVECTION_H
#define SWE1D_NODALADVECTION_H

#include "types.h"

class NodalAdvection{
private:
    //Stiffness Matrix
    T *s;
    //Mass matrix
    T *m;
    //inverse mass matrix
    T *inv_m;
    //First N-Matrix
    T *n0;
    //Second N-Matrix
    T *n1;
    // m_u[nElements + 2][2]
    T *m_u;
    T *m_uNetUpdatesLeft;
    T *m_uNetUpdatesRight;
    //Time derivative of u
    T **m_ut;
    T m_size;

    //Time variables
    T t0 = 0.0;
    //Advection velocity
    T a;

    T m_cellSize;
public:
    NodalAdvection(float a_def, T *h, unsigned int size, unsigned int cellsize):
            a(a_def),
            m_size(size),
            m_u(h),
            m_cellSize(cellsize)
            {   //Allocate matrices
                s = new T[4];
                m = new T[4];
                inv_m = new T[4];
                n0 = new T[4];
                n1 = new T[4];
                m_uNetUpdatesLeft = new T[size+1];
                m_uNetUpdatesRight = new T[size+1];
                T m_ut[size+2][4];

                //Define entries of matrices (0=00, 1 = 01, 2 = 10, 3 = 11)
                s[0] = -0.5f;
                s[1] = -0.5f;
                s[2] = 0.5f;
                s[3] = 0.5f;

                m[0] = 1.0f/3.0f;
                m[1] = 1.0f/6.0f;
                m[2] = 1.0f/6.0f;
                m[3] = 1.0f/3.0f;

                inv_m[0] = 4.0f;
                inv_m[1] = -2.0f;
                inv_m[2] = -2.0f;
                inv_m[3] = 4.0f;

                n0[0] = 1.0f;
                n0[1] = 0.0f;
                n0[2] = 0.0f;
                n0[3] = 0.0f;

                n1[0] = 0.0f;
                n1[1] = 0.0f;
                n1[2] = 0.0f;
                n1[3] = 1.0f;

            }

    ~NodalAdvection(){
        //Free allocated memory
        delete [] s;
        delete [] m;
        delete [] inv_m;
        delete [] n0;
        delete [] n1;
        delete [] m_uNetUpdatesLeft;
        delete [] m_uNetUpdatesRight;
        delete [] m_ut;
    }
    //Methods
    void setBoundaryConditions();
    T compute(T t);
    T computeLocalLaxFriedrichsFluxes(T t);
    void updateUnknowns(T t);

    void computeTimeDerivative();

    void computeEulerStep(T delta_t);
};

#endif //SWE1D_NODALADVECTION_H
