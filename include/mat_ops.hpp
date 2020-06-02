
/**
 * @file mat_ops.hpp
 * @author Lee J. O'Riordan (lee.oriordan@ichec.ie)
 * @brief Templated methods to manipulate small matrices. Adapted from QNLP.
 * @version 0.2
 * @date 2020-06-01
 */

#ifndef MAT_OPS
#define MAT_OPS

#include <complex>
#include <vector>
#include <iostream>

/**
 * @brief Calculates the unitary matrix square root (U == VV, where V is returned)
 * 
 * @tparam Type ComplexDP or ComplexSP
 * @param U Unitary matrix to be rooted
 * @return Matrix V such that VV == U
 */
template <class Mat2x2Type>            
const Mat2x2Type matrixSqrt(const Mat2x2Type& U){
    Mat2x2Type V(U);
    std::complex<double> delta = U(0,0)*U(1,1) - U(0,1)*U(1,0);
    std::complex<double> tau = U(0,0) + U(1,1);
    std::complex<double> s = sqrt(delta);
    std::complex<double> t = sqrt(tau + 2.0*s);

    //must be a way to vectorise these; TinyMatrix have a scale/shift option?
    V(0,0) += s;
    V(1,1) += s;
    std::complex<double> scale_factor(1.,0.);
    scale_factor /= t;
    V(0,0) *= scale_factor; //(std::complex<double>(1.,0.)/t);
    V(0,1) *= scale_factor; //(1/t);
    V(1,0) *= scale_factor; //(1/t);
    V(1,1) *= scale_factor; //(1/t);

    return V;
}

/**
 * @brief Function to calculate the adjoint of an input matrix
 * 
 * @tparam Type ComplexDP or ComplexSP
 * @param U Unitary matrix to be adjointed
 * @return qhipster::TinyMatrix<Type, 2, 2, 32> U^{\dagger}
 */
template <class Mat2x2Type>            
Mat2x2Type adjointMatrix(const Mat2x2Type& U){
    Mat2x2Type Uadjoint(U);
    std::complex<double> tmp;
    tmp = Uadjoint(0,1);
    Uadjoint(0,1) = Uadjoint(1,0);
    Uadjoint(1,0) = tmp;
    Uadjoint(0,0) = std::conj(Uadjoint(0,0));
    Uadjoint(0,1) = std::conj(Uadjoint(0,1));
    Uadjoint(1,0) = std::conj(Uadjoint(1,0));
    Uadjoint(1,1) = std::conj(Uadjoint(1,1));
    return Uadjoint;
}
#endif