//
// Created by lisa on 18.10.17.
//

#include <cmath>
#include <algorithm>
#include "types.h"
#include "NodalAdvection.h"
#include <vector>


vect matmult(vect a, vect b, int n, int l, int m){

    //Mutliplication of Matrix A (nxm), Matrix b(mxl)
    vect res (n*l);
    //fill with zeros
    for (int a = 0; a<n*l; a++){
        res[a] = 0.0;
    }
    for(int i = 0; i<n; i++){
        for(int k = 0; k<l; k++){
            for(int j = 0; j<m; j++){
                res[i*n+k] += (a[i*m+j]*b[j*m+k]);
            }
        }
    }
    return res;
}

u matVecMult(vect a,u b){
    //Mutliplication of 2x2 Matrix with 2x1 vector
    u res;
    res.u0 =0;
    res.u1 = 0;
    res.u0 = a[0]*b.u0+a[1]*b.u1;
    res.u1 = a[2]*b.u0+a[3]*b.u1;
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
    for(unsigned int i = 1; i<m_size; i++){
        //Factor a for local Lax Friedrichs Method (aLeft = a_i-1/2, aRight = a_i+1/2
        T aLeft;
        T aRight;

        //Eigenvalue lambda = a
        aLeft  = std::abs(a);
        aRight = std::abs(a);

        //Compute fluxes for u
        //m_u[i][1]+m_u[i+1][0] ... m_u[i][1]  -m_u[i+1][0]

        m_uNetUpdatesRight.at(i) = 0.5f*(a*m_u.at(i).u1+a*m_u.at(i+1).u0-(aRight*(m_u.at(i+1).u0-m_u.at(i).u1)));
        //m_uNetUpdatesRight[i]=0.5f*(m_u[i]+m_u[i+1]-(aRight*(m_u[i+1]-m_u[i])));
        //m_u[i-1][1]+m_u[i][0] ... m_u[i-1][1]-m_u[i][0]
        m_uNetUpdatesLeft.at(i) = 0.5f*(a*m_u.at(i-1).u1+a*m_u.at(i).u0-(aLeft*(m_u.at(i).u0-m_u.at(i-1).u1)));
        //m_uNetUpdatesLeft[i]=0.5f*(m_u[i-1]+m_u[i]-(aLeft*(m_u[i]-m_u[i-1])));

        // Update maxWaveSpeed
        maxWaveSpeed = std::max(std::max(maxWaveSpeed,aLeft),aRight);

    }
    // Compute CFL condition (delta_t = delta_x/(a*3))

    T maxTimeStep = m_cellSize/(std::abs(a)*3);

    return maxTimeStep;
}

void NodalAdvection::computeTimeDerivative(){
    for(unsigned int i = 1; i<m_size; i++) {
        //Compute S *(a*delta_x*u)
        //only size 2
        vect rightTermFirst(2,0.0);

        for (unsigned int j = 0; j<2; j++){
        //for (unsigned int j = 0; j<4; j++){
            if (j==0){
                rightTermFirst.at(0) += s.at(j) *a * m_u.at(i).u0;
                rightTermFirst.at(1) += s.at(2+j) *a * m_u.at(i).u0;
            }
            else {
                rightTermFirst.at(0) += s.at(j) * a * m_u.at(i).u1 ;
                rightTermFirst.at(1) += s.at(2+j)* a * m_u.at(i).u1 ;
            }

        //rightTermFirst[j] = s[j]*a*m_u[i]*m_cellSize;

        }
        //s.o.
        //Compute n0*netupdatesRight*delta_x
        vect rightTermSecond(2,0.0);

        for(unsigned int k = 0; k<2; k++){
            //s.o.
            rightTermSecond.at(0) += n0.at(k)*m_uNetUpdatesLeft.at(i);
            rightTermSecond.at(1) += n0.at(2+k)*m_uNetUpdatesLeft.at(i);
            //rightTermSecond[k] = n0[k]*m_uNetUpdatesRight[i]*m_cellSize;
        }
        //s.o.
        //Compute n1*netUpdatesLeft*delta_x
        vect rightTermThird(2,0.0);

        //s.o.
        for(unsigned int l = 0; l<2; l++){
            //s.o.
            rightTermThird.at(0) += n1.at(l)*m_uNetUpdatesRight[i];
            rightTermThird.at(1) += n1.at(2+l)*m_uNetUpdatesRight[i];
            //rightTermThird[l] = n1[l]*m_uNetUpdatesLeft[i]*m_cellSize;
        }
        //Put right side first+second-third together
        //s.o.
        u rightTerm;
        rightTerm.u0 = (rightTermFirst.at(0)+rightTermSecond.at(0)-rightTermThird.at(0))/m_cellSize;
        rightTerm.u1 = (rightTermFirst.at(1)+rightTermSecond.at(1)-rightTermThird.at(1))/m_cellSize;

        /*u temp = matVecMult(inv_m,rightTerm);
        m_ut[i] = temp;*/

        m_ut.at(i) = matVecMult(inv_m,rightTerm);

        //m_ut[i]=matmult(inv_m,rightTerm,2,1,2);
    }
}

void NodalAdvection::computeEulerStep(T delta_t){
    u u0;
    u0.u0 = 0.0;
    u0.u1 = 0.0;
    vecu tmp(m_size+2,u0);
    for(unsigned int i = 0; i<m_size+2; i++){
        tmp.at(i).u0=m_u.at(i).u0+(delta_t*m_ut.at(i).u0);
        tmp.at(i).u1=m_u.at(i).u1+(delta_t*m_ut.at(i).u1);
    }
    /*for(unsigned int i = 0; i<m_size+2; i++){
        m_u.at(i).u0=m_u.at(i).u0+(delta_t*m_ut.at(i).u0);
        m_u.at(i).u1=m_u.at(i).u1+(delta_t*m_ut.at(i).u1);
    }*/
    for (unsigned int j = 0; j<m_size+2; j++) {
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
    }
}

void NodalAdvection::setBoundaryConditions() {
    m_u.at(0) = m_u.at(1); m_u.at(m_size+1) = m_u.at(m_size);
}

vecu NodalAdvection::setH(){
    u u0;
    u0.u0 = 0.0;
    u0.u1 = 0.0;
    vecu result(m_u.size(), u0);
    for (unsigned int i = 0; i<m_u.size(); i++){
        result.at(i) = m_u.at(i);
    }
    return result;
}


