//
// Created by lisa on 18.10.17.
//

#include <cmath>
#include <algorithm>
#include "types.h"
#include "NodalAdvection.h"


T* matmult(T* a, T* b, int n, int l, int m){
    //TODO implement algorithm for multiplication of sparse matrices
    //Mutliplication of Matrix A (nxm), Matrix b(mxl)
    T* res = new T[n*l];
    for(int i = 0; i<n; i++){
        for(int k = 0; k<l; k++){
            for(int j = 0; j<m; j++){
                res[i*n+k] += (a[i*m+j]*b[j*m+k]);
            }
        }
    }
    return res;
}

/*//Compute Jacobi Polynomial of order n at value x
T JacobiP(T x, T alpha, T beta, int n){
    T res = 0;
    if(n==0){
        res = 1;
    }
    else if (n==1){
        res = 0.5*(alpha-beta+(alpha+beta+2)*x);
    }
    else {
        T a_1n = 2*(n+1)*(n+alpha+beta+1)*(2*n+alpha+beta);
        T a_2n = (2*n+alpha+beta+1)*(alpha+alpha-beta*beta);
        T a_3n = (2*n+alpha+beta)*(2*n+alpha+beta+1)*(2*n+alpha+beta+2);
        T a_4n = 2*(n+alpha)*(n+beta)*(2*n+alpha+beta+2);
        res = (((a_2n+a_3n*x)*JacobiP(x, alpha, beta, n-1))-a_4n*JacobiP(x,alpha,beta,n-2))/a_1n;
    }
}

T* gradJacobiP(T* r, T alpha, T beta){
    //Evaluate derivative of Jacobi Polynomial
    T *res = new T[r.length];
    for (int i = 0; i<r.size; i++){
        res[i]=0;
    }
    if (NodalAdvection::N == 0){
        for (int i = 0; i<r.size; i++){
            res[i]=0;
        }
    }
    else {
        for(int i = 0; i<r.size(); i++){
            res[i]=std::sqrt(N*(N+alpha+beta+1))*JacobiP(r[i], alpha+1, beta+1, N-1);
        }
    }
    return res;
}*/

//Solve standard 1D advection equation u_t + a*u_x = 0
//Time derivative u_t = inv_m*(s*a*deltax*u+no*netUpdatesRight*deltax-n1*netUpdatesLeft*deltax)
//Euler step u_n+1 = u_n + deltat*u_t(t_n, u_n)
//CFL-Condition deltat = deltax/(a*3)

T NodalAdvection::computeLocalLaxFriedrichsFluxes(T t){
    T maxWaveSpeed = 0.f;
    T deltaT = t - t0;
    t0 = t;
    //Loop over all intervalls
    for(unsigned int i = 0; i<m_size+2; i++){
        //Factor a for local Lax Friedrichs Method (aLeft = a_i-1/2, aRight = a_i+1/2
        T aLeft;
        T aRight;

        //Eigenvalue lambda = a
        aLeft  = std::abs(a);
        aRight = std::abs(a);

        //Compute fluxes for u
        m_uNetUpdatesRight[i]=0.5f*(m_u[i]+m_u[i+1]-(aRight*(m_u[i+1]-m_u[i])));
        m_uNetUpdatesLeft[i]=0.5f*(m_u[i-1]+m_u[i]-(aLeft*(m_u[i]-m_u[i-1])));

        // Update maxWaveSpeed
        maxWaveSpeed = std::max(std::max(maxWaveSpeed,aLeft),aRight);

    }
    // Compute CFL condition (delta_t = delta_x/(a*3))

    T maxTimeStep = m_cellSize/(a*3);

    return maxTimeStep;
}

void NodalAdvection::computeTimeDerivative(){
    for(unsigned int i = 0; i<m_size+2; i++){
        //Compute S *(a*delta_x*u)
        T *rightTermFirst = new T[4];
        for (unsigned int j = 0; j<4; j++){
            rightTermFirst[j] = s[j]*a*m_u[i]*m_cellSize;
        }
        //Compute n0*netupdatesRight*delta_x
        T *rightTermSecond = new T[4];
        for(unsigned int k = 0; k<4; k++){
            rightTermSecond[k] = n0[k]*m_uNetUpdatesRight[i]*m_cellSize;
        }
        //Compute n1*netUpdatesLeft*delta_x
        T *rightTermThird = new T[4];
        for(unsigned int l = 0; l<4; l++){
            rightTermThird[l] = n1[l]*m_uNetUpdatesLeft[i]*m_cellSize;
        }
        //Put right side first+second-third together
        T *rightTerm = new T[4];
        for(unsigned int o = 0; o<4; o++){
            rightTerm[o] = rightTermFirst[o]+rightTermSecond[o]-rightTermThird[o];
        }
        m_ut[i]=matmult(inv_m,rightTerm,2,2,2);
    }
}

void NodalAdvection::computeEulerStep(T delta_t){
    for(unsigned int i = 0; i<m_size+2; i++){
        //Commentend out because of invalid operation scalar + matrix
        //m_u[i]=m_u[i]+(delta_t*m_ut[i]);
    }
}


