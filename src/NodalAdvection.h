//
// Created by lisa on 18.10.17.
//

#ifndef SWE1D_NODALADVECTION_H
#define SWE1D_NODALADVECTION_H

#include "types.h"

class NodalAdvection{
private:
    //Vertex coordinates (Size K+1)
    T *VX;
    //Element to Vertex (Kx2)
    T *EToV;
    //Legendre-Gauss-Lobatto quadrature points
    T *r;
    //Array containing physical coordinates of grid points (size N_pxK)
    T *x;
    //Face to Vertex matrix (2KxN)
    T *FToV;
    //Connectivity matrix
    T *FToF;
    //Element to element matrix
    T *EToE;
    //Element to face
    T *EtoF;
    //
    T *vMapM;
    T *vMapP;
    //Time variables
    T t0 = 0.0;
    T h;
public:
    NodalAdvection(unsigned int K, ):{};
    ~NodalAdvection(){}
    //Methods
    void setBoundaryConditions();
    T compute(T t);
    void updateUnknowns(T t);

//Order of approximation
static unsigned int N = 1;
};

#endif //SWE1D_NODALADVECTION_H
