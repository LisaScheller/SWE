//
// Created by lisa on 18.10.17.
//

#include <cmath>
#include "types.h"
#include "NodalAdvection.h"
int N = NodalAdvection::N;
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

//Compute Jacobi Polynomial of order n at value x
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
}

T* vandermonde1D(T*r){
    //Compute the Vandermonde Matrix, dimension lenght(r)*N+1
    T *V1D = new T[r.length*(N+1)];
    for(int j = 1; j<=N+1; j++){
        //TODO Find correct index of r
        V1D[2*j-2] = JacobiP(r[2*j-2],0,0,j-1);
        V1D[2*j-1] = JacobiP(r[2*j-1],0,0,j-1);
    }
    return V1D;

}

T* gradVandermonde1D(T* r){
    T *DVr = new T[r.length*(N+1)];
    for(int i = 0; i<N+1; i++){
        for(int j = 0; j<r.length; j++){
            DVr[i*(N+1)+j]=gradJacobiP(r[i*(N+1)+j],0,0,i);
        }
    }
    return DVr;
}
T* DMatrix(T*r, T* V){
    T* Vr = gradVandermonde1D(r);
    T* Dr = new T[];
    //TODO Compute with maple???
    for (int i = 0; i<Vr.size; i++){
        Dr[i] = Vr[i]/V[i];
    }
    return Dr;
}

T* Lift1D(){
    T* Emat = new T[Np*(Nfaces*Nfp)];
    Emat[1*(Nfaces*Nfp)+1]=1.0;
    Emat[Np*(Nfaces*Nfp)+2]=1.0;
    T* lift = new T[];
    lift = V*(invV*Emat);
}

T* computeX(){
    T* va = new T[K];
    T* vb = new T[K];
    //va first column of EToV, vb second column
    for(int i = 0; i<K; i++){
        va[i] = EToV[i*2+1];
        vb[i] = EToV[i*2+2];
    }
    T* ones = new T[Np];
    for (int i = 0; i<Np; i++){
        ones[i] = 1;
    }
    return ones*Vx(va)+0.5*(r+1)*(VX(vb)-VX(va));
}

T* geometricFactors1D(T* x, T* Dr){
    T *xr = Dr*x;
    J = xr;
    rx = 1/J;
}

T* normals1D(){
    T* nx = new T[(Nfp*Nfaces)*K];
    for(int i = 0; i<2; i++){
        for (int j = 0; j<K; j++){
            if(i==0){
                nx[i*k+j]=-1.0;
            }
            else nx[i*K+j]= 1.0;
        }
    }
    return nx;
}

//TODO Compute FToF and FToV

T* connect1D(){
    int NFaces = 2;
    int K = sizeof(EToV,1);
    int totalFaces = NFaces*K;
    int NV = K+1;
    T* vn = new T[1*2];
    
}
