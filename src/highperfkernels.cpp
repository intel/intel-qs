#include "../include/highperfkernels.hpp"

namespace iqs {

#define Specialization(MainLoop_) \
  frac_of_state_accessed = 1.; \
  /* Handling many special cases: should be easy to add  \
     more other cases; overhead is small  */ \
  if (m01 == Type(0., 0.) && m10 == Type(0., 0.))  \
  {  \
    if (m00 == Type(1., 0.)) {  \
      frac_of_state_accessed = 0.5;  \
      if (m11 == Type(1., 0.)) {  \
        /* identify  */\
        /*  \
             |1  0|  \
         I = |    |  \
             |0  1|  \
        */  \
        frac_of_state_accessed = 0.0;  \
        label += "_Id";  \
      } else if (m11 == Type(-1., 0.)) {  \
        /*  \
             |1  0|  \
         Z = |    |  \
             |0 -1|  \
        */   \
        label += "_Z";  \
        MainLoop_(KeyLoop_000T,  \
                  DC,      DC,  \
                  DC,      -in1);  \
      } else if (m11 == Type(0., 1.)) {  \
        /*  \
             |1  0|  \
         S = |    |  \
             |0  i|  \
        */   \
        label += "_S";  \
        MainLoop_(KeyLoop_000T,  \
                  DC,           DC,  \
                  DC,      Type(0., 1.)*in1);  \
      } else {  \
        /*  \
             |1  0|  \
         U = |    |  \
             |0  c|  \
        */   \
        label += "_100c";  \
        MainLoop_(KeyLoop_000T,  \
                  DC,           DC,  \
                  DC,      m11*in1);  \
      }  \
    } else {  \
      /*  \
           |c1 0 |  \
       U = |     |  \
           |0  c2|  \
      */   \
      MainLoop_(KeyLoop_T00T,  \
                m00*in0,       DC,  \
                DC,      m11*in1);  \
    }  \
  } else if (m00 == Type(0., 0.) && m11 == Type(0., 0.)) {  \
      if (m01.real() == 0. && m10.real() == 0) {  \
        if (m01.imag() == -1. && m10.imag() == 1) {  \
          /*  \
               |0    -i |  \
           Y = |        |  \
               |i    0  |  \
          */   \
          label += "_Y";  \
          MainLoop_(KeyLoop_0TT0,  \
                     DC,  Type(0., -1.) * in1,  \
                     Type(0., 1.) * in0, DC);  \
        } else {  \
          /*  \
               |0       c1.im()|  \
           U = |               |  \
               |c2.im() 0      |  \
          */   \
           MainLoop_(KeyLoop_0TT0,  \
                     DC,  Type(0., m01.imag()) * in1,  \
                     Type(0., m10.imag()) * in0, DC);         \
        }  \
      } else if (m01.real() == 1. && m10.real() == 1. && m01.imag() == 0. && m10.imag() == 0.) {  \
       /*  \
            |0  1|  \
        X = |    |  \
            |1  0|  \
       */   \
        label += "_X";  \
        MainLoop_(KeyLoop_0TT0,  \
                  DC, in1,  \
                  in0,      DC);         \
      } else {  \
       /*  \
            |0   c1|  \
        U = |      |  \
            |c2  0 |  \
       */   \
        MainLoop_(KeyLoop_0TT0,  \
                  DC,  m01 * in1,  \
                  m10 * in0, DC);         \
      }  \
  } else {  \
      if (m00.imag() == 0. && m01.imag() == 0. && m10.imag() == 0. && m11.imag() == 0.) {  \
        /*  \
             |c00.re  c01.re|  \
         U = |              |  \
             |c10.re  c11.re|  \
        */   \
        label += "_H";  \
        MainLoop_(KeyLoop_TTTT,  \
                  Type(m00.real(), 0.) * in0, Type(m01.real(), 0.) * in1,  \
                  Type(m10.real(), 0.) * in0, Type(m11.real(), 0.) * in1);  \
      } else if (m00.imag() == 0. && m01.real() == 0. && m10.real() == 0. && m11.imag() == 0.) {  \
        /*  \
             |c00.re  c01.im|  \
         U = |              |  \
             |c10.im  c11.re|  \
        */   \
        MainLoop_(KeyLoop_TTTT,  \
                  Type(m00.real(), 0.) * in0,  Type(0., m01.imag()) * in1,  \
                  Type(0., m10.imag()) * in0,  Type(m11.real(), 0.) * in1);  \
      } else {  \
        /*  \
             |c00  c01|  \
         U = |              |  \
             |c10  c11|  \
        */   \
        MainLoop_(KeyLoop_TTTT,  \
                  m00 * in0, m01 * in1,  \
                  m10 * in0, m11 * in1);  \
      }  \
  }
  

// specialized loops with autovectorization
// TODO: Not using ICC vectorization: need to rewrite in intrincics
// #define SIMD _Pragma("simd vectorlength(2) assert")
#define SIMD 
#define KeyLoop_TTTT(simdpragma, from, to, indsht0, indsht1, state0, state1, t00, t01, t10, t11) \
  simdpragma                                            \
  _Pragma("vector aligned")                             \
  for (std::size_t ind0 = from; ind0 < to; ind0++) {    \
    std::size_t i0 = ind0 + indsht0;                    \
    std::size_t i1 = ind0 + indsht1;                    \
    Type in0 = state0[i0], in1 = state1[i1];     \
    state0[i0] = (t00) + (t01);                         \
    state1[i1] = (t10) + (t11);                         \
  }
  //__asm call LOOP_END
#define KeyLoop_T00T(simdpragma, from, to, indsht0, indsht1, state0, state1, t00, t01, t10, t11) \
  simdpragma                                            \
  _Pragma("vector aligned")                             \
  for (std::size_t ind0 = from; ind0 < to; ind0++) {    \
    std::size_t i0 = ind0 + indsht0;                    \
    std::size_t i1 = ind0 + indsht1;                    \
    Type in0 = state0[i0], in1 = state1[i1];     \
    state0[i0] = t00;                                   \
    state1[i1] = t11;                                   \
  }
  // __asm call LOOP_END
#define KeyLoop_0TT0(simdpragma, from, to, indsht0, indsht1, state0, state1, t00, t01, t10, t11) \
  simdpragma                                            \
  _Pragma("vector aligned")                             \
  for (std::size_t ind0 = from; ind0 < to; ind0++) {    \
    std::size_t i0 = ind0 + indsht0;                    \
    std::size_t i1 = ind0 + indsht1;                    \
    Type in0 = state0[i0], in1 = state1[i1];     \
    state0[i0] = t01;                                   \
    state1[i1] = t10;                                   \
  }
  // __asm call LOOP_END

#define KeyLoop_000T(simdpragma, from, to, indsht0, indsht1, state0, state1, t00, t01, t10, t11) \
  simdpragma                                            \
  _Pragma("vector aligned")                             \
  for (std::size_t ind0 = from; ind0 < to; ind0++) {    \
    std::size_t i0 = ind0 + indsht0;                    \
    std::size_t i1 = ind0 + indsht1;                    \
    Type in1 = state1[i1];                              \
    state1[i1] = t11;                                   \
  }


#define _Loop_SN_(KeyLoop, t00, t01, t10, t11) \
{\
  _Pragma("omp parallel for") \
  KeyLoop(SIMD, start, end, indsht0, indsht1, state0, state1, t00, t01, t10, t11) \
}
template <class Type>
__attribute__((noinline))
void Loop_SN(std::size_t start, std::size_t end, Type *state0, Type *state1,
             std::size_t indsht0, std::size_t indsht1, TM2x2<Type> const&m, 
             bool specialize, Timer *timer)
{
  Type m00 = m[0][0], 
       m01 = m[0][1], 
       m10 = m[1][0], 
       m11 = m[1][1];

  std::string label;
  double frac_of_state_accessed;
  double ttot = 0., tnov = 0., ttmp1, ttmp2;
  ttmp1 = sec();

  if(specialize == false)
  {
      frac_of_state_accessed = 1.;
      label = "general";
#pragma omp parallel for
      for (std::size_t i = start; i < end; i++)
      {
          std::size_t i0 = i + indsht0;
          std::size_t i1 = i + indsht1;
          Type in0 = state0[i0], in1 = state1[i1];
          state0[i0] = m00 * in0 + m01 * in1;
          state1[i1] = m10 * in0 + m11 * in1;
      }
  }
  else
  {
    Specialization(_Loop_SN_);
  }

  if (timer)
  {
      ttot = sec() - ttmp1;
      double datab = ((state0 == state1) ? 2.0 : 4.0) * 
                     frac_of_state_accessed * sizeof(state0[0]) * D(end - start);
      // printf("datab=%lf len=%lu time=%lf bw=%lf\n", datab, end-start, ttot, datab / ttot / 1e9);
      double flops = D(1L << 19) * 38.0;
      double gflops = flops / ttot / 1e9;
      // printf("label=%s ttot = %.4lfs bw = %.2lf GB/s\n",
      //        label.c_str(), ttot, datab / ttot / 1e9);

      timer->record_sn(ttot, datab / ttot);
  }
}
template 
__attribute__((noinline))
void Loop_SN(std::size_t start, std::size_t end, ComplexSP *state0, ComplexSP *state1,
             std::size_t indsht0, std::size_t indsht1, TM2x2<ComplexSP> const&m, 
             bool specialize, Timer *timer);
template 
__attribute__((noinline))
void Loop_SN(std::size_t start, std::size_t end, ComplexDP *state0, ComplexDP *state1,
             std::size_t indsht0, std::size_t indsht1, TM2x2<ComplexDP> const&m, 
             bool specialize, Timer *timer);


#define _Loop_DN_(KeyLoop, t00, t01, t10, t11) \
{ \
  /*label = std::string(xstr(t00)) + " " + std::string(xstr(t01)) + \
          " " + std::string(xstr(t10)) + " " + std::string(xstr(t11)); */\
  ttmp2 = sec();\
  if ((gend - gstart) / (UL(1) << (pos + UL(1))) >= nthreads) { \
    if (pos == 0) { \
      _Pragma("omp parallel for") \
      for (std::size_t group = gstart; group < gend; group += (UL(1) << (pos + UL(1)))) { \
        KeyLoop(, group, (group + (UL(1) << pos)), indsht0, indsht1, state0, state1, t00, t01, t10, t11) \
      } \
    } else { \
      _Pragma("omp parallel for") \
      for (std::size_t group = gstart; group < gend; group += (UL(1) << (pos + UL(1)))) { \
        KeyLoop(SIMD, group, (group + (UL(1) << pos)), indsht0, indsht1, state0, state1, t00, t01, t10, t11) \
      } \
    } \
  } else { \
    if (pos == 0) { \
      for (std::size_t group = gstart; group < gend; group += (UL(1) << (pos + UL(1)))) { \
        _Pragma("omp parallel for") \
        KeyLoop(, group, (group + (UL(1) << pos)), indsht0, indsht1, state0, state1, t00, t01, t10, t11) \
      } \
    } else { \
      for (std::size_t group = gstart; group < gend; group += (UL(1) << (pos + UL(1)))) { \
        _Pragma("omp parallel for") \
        KeyLoop(SIMD, group, (group + (UL(1) << pos)), indsht0, indsht1, state0, state1, t00, t01, t10, t11) \
      } \
    } \
  } \
  tnov = sec() - ttmp2; \
}


template <class Type>
__attribute__((noinline))
void Loop_DN(std::size_t gstart, std::size_t gend, std::size_t pos,
             Type *state0, Type *state1,
             std::size_t indsht0, std::size_t indsht1,
             TM2x2<Type> const&m, 
             bool specialize, Timer *timer)
{
  // TODO: Allow for case where state is not aligned: need SIMD ISA for un-aligned access.
  if (specialize) // FIXME: this condition on the assertions was added later. Need to be validated.
  {
      assert((UL(state0) % 256) == 0);
      assert((UL(state1) % 256) == 0);
  }
#if defined(__ICC) || defined(__INTEL_COMPILER)
  __assume_aligned(state0, 256);
  __assume_aligned(state1, 256);
#endif

  Type m00 = m[0][0],
       m01 = m[0][1],
       m10 = m[1][0],
       m11 = m[1][1];

  std::string label;
  double frac_of_state_accessed;
  double ttot = 0., tnov = 0., ttmp1, ttmp2;
  ttmp1 = sec();

  size_t nthreads = 1;
#ifdef _OPENMP
#pragma omp parallel 
  {
      nthreads = omp_get_num_threads();
  }
#endif
  // TODO: Add nthreads check to clamp to smaller number if too little work.
  // TODO: Generalize for AVX3 cases so we check for pos <=1 etc.

  if(specialize == false)
  {
      frac_of_state_accessed = 1.;
      label = "general";

      if((gend - gstart) / (1L << (pos+1L)) >= nthreads)
      {
#pragma omp parallel for 
          for(std::size_t group = gstart; group < gend; group += (1L << (pos + 1L)))
          {
              for(std::size_t ind0 = group; ind0 < group + (1L << pos); ind0++)
              {
                  std::size_t i0 = ind0 + indsht0;
                  std::size_t i1 = ind0 + indsht1;
                  Type in0 = state0[i0], in1 = state1[i1];
                  state0[i0] = m00 * in0 + m01 * in1;
                  state1[i1] = m10 * in0 + m11 * in1;
              }
            }
      }
      else
      {
          for(std::size_t group = gstart; group < gend; group += (1L << (pos + 1)))
          {
#pragma omp parallel for
              for(std::size_t ind0 = group; ind0 < group + (1L << pos); ind0++)
              {
                  std::size_t i0 = ind0 + indsht0;
                  std::size_t i1 = ind0 + indsht1;
                  Type in0 = state0[i0], in1 = state1[i1];
                  state0[i0] = m00 * in0 + m01 * in1;
                  state1[i1] = m10 * in0 + m11 * in1;
              }
          }
      }
  }
  else
  {
      Specialization(_Loop_DN_);
  }

  if(timer)
  {
      ttot = sec() - ttmp1;     
      double datab = 2.0 * frac_of_state_accessed * sizeof(state0[0]) * D(gend - gstart);
      double flops = D(1L << 19) * 38.0;
      double gflops = flops / ttot / 1e9;
      // printf("label=%s ttot = %.4lfs bw = %.2lf GB/s\n",
      //         label.c_str(), ttot, datab / ttot / 1e9);

      timer->record_dn(ttot, datab / ttot);
  }


}
template 
__attribute__((noinline))
void Loop_DN(std::size_t gstart, std::size_t gend, std::size_t pos,
             ComplexSP *state0, ComplexSP *state1,
             std::size_t indsht0, std::size_t indsht1,
             TM2x2<ComplexSP> const&m, 
             bool specialize, Timer *timer);
template 
__attribute__((noinline))
void Loop_DN(std::size_t gstart, std::size_t gend, std::size_t pos,
             ComplexDP *state0, ComplexDP *state1,
             std::size_t indsht0, std::size_t indsht1,
             TM2x2<ComplexDP> const&m, 
             bool specialize, Timer *timer);


template <class Type>
__attribute__((noinline))
void Loop_TN(Type *state,
             std::size_t c11, std::size_t c12, std::size_t c13,
             std::size_t c21, std::size_t c22, std::size_t c23,
             std::size_t c31, std::size_t c32, 
             std::size_t ind_shift, TM2x2<Type> const&m, bool specialize, Timer *timer)
{
  double ttmp1 = sec(), ttot = 0.;
  Type m00 = m[0][0],
       m01 = m[0][1],
       m10 = m[1][0],
       m11 = m[1][1];

  // std::cout << m00 << " " << m01 << " " << m10 << " " << m11 << std::endl;
  size_t nthreads = 1;
#ifdef _OPENMP
#pragma omp parallel 
  {
      nthreads = omp_get_num_threads();
  }
#endif

  if ((c12 - c11) / c13 >= nthreads)
  {
#pragma omp parallel for
      for (std::size_t l1 = c11; l1 < c12; l1 += c13)
      {
          for (std::size_t l2 = l1 + c21; l2 < l1 + c22; l2 += c23)
          {      
              for (std::size_t ind0 = l2 + c31; ind0 < l2 + c32; ind0++)
              {
                  std::size_t ind1 = ind0 + ind_shift;
                  Type in0 = state[ind0], in1 = state[ind1];
                  state[ind0] = m00 * in0 + m01 * in1;
                  state[ind1] = m10 * in0 + m11 * in1;
              }
          }
      }
  }
  else
  {
      for (std::size_t l1 = c11; l1 < c12; l1 += c13)
      {
          if ((l1 + c22 - l1 - c21) / c23 >= nthreads)
          {
#pragma omp parallel for
              for (std::size_t l2 = l1 + c21; l2 < l1 + c22; l2 += c23)
              {
                  for (std::size_t ind0 = l2 + c31; ind0 < l2 + c32; ind0++)
                  {
                      std::size_t ind1 = ind0 + ind_shift;
                      Type in0 = state[ind0], in1 = state[ind1];
                      state[ind0] = m00 * in0 + m01 * in1;
                      state[ind1] = m10 * in0 + m11 * in1;
                  }
              }
          }
          else
          {
              for (std::size_t l2 = l1 + c21; l2 < l1 + c22; l2 += c23)
              {
#pragma omp parallel for
                  for (std::size_t ind0 = l2 + c31; ind0 < l2 + c32; ind0++)
                  {
                      std::size_t ind1 = ind0 + ind_shift;
                      Type in0 = state[ind0], in1 = state[ind1];
                      state[ind0] = m00 * in0 + m01 * in1;
                      state[ind1] = m10 * in0 + m11 * in1;
                  }
              }
          }
      }
  }

  assert(((c12 - c11) % c13) == 0);
  assert(((c22 - c21) % c23) == 0);

  if (timer)
  {
    ttot = sec() - ttmp1;
    double datab =
      4.0 * sizeof(state[0]) * D((c12 - c11) / c13) * D((c22 - c21) / c23) * D(c32 - c31);
    double flops = D(1L << 19) * 38.0;
    double gflops = flops / ttot / 1e9;
    // printf("label=%s ttot = %.4lfs bw = %.2lf GB/s\n",
    //        "ScaleState", ttot, datab / ttot / 1e9);
    timer->record_tn(ttot, datab / ttot);
  }
}

template
__attribute__((noinline))
void Loop_TN(ComplexSP *state, std::size_t c11, std::size_t c12,
             std::size_t c13, std::size_t c21, std::size_t c22, 
             std::size_t c23, std::size_t c31, std::size_t c32, 
             std::size_t ind_shift, TM2x2<ComplexSP> const&m, bool specialize, Timer *timer);
template
__attribute__((noinline))
void Loop_TN(ComplexDP *state, std::size_t c11, std::size_t c12,
             std::size_t c13, std::size_t c21, std::size_t c22, 
             std::size_t c23, std::size_t c31, std::size_t c32, 
             std::size_t ind_shift, TM2x2<ComplexDP> const&m, bool specialize, Timer *timer);

template <typename Type>
void ScaleState(std::size_t start, std::size_t end, Type *state, 
                const Type &s, Timer *timer)
{
  double ttmp1 = sec(), ttot = 0.;
  if (s != Type(1., 0.)) {
#pragma omp parallel for    
    for (std::size_t i = start;  i < end; i++) state[i] *= s;
  }
  if (timer) {
    ttot = sec() - ttmp1;
    double datab = 2.0 * sizeof(state[0]) * D(end - start);
    double flops = D(1L << 19) * 38.0;
    double gflops = flops / ttot / 1e9;
    // printf("label=%s ttot = %.4lfs bw = %.2lf GB/s\n",
    //        "ScaleState", ttot, datab / ttot / 1e9);
    timer->record_sn(ttot, datab / ttot);
  }

}
template void ScaleState<ComplexSP>(std::size_t start, std::size_t end, 
                                    ComplexSP *state, const ComplexSP &s, Timer *timer);
template void ScaleState<ComplexDP>(std::size_t start, std::size_t end, 
                                    ComplexDP *state, const ComplexDP &s, Timer *timer);

} // close namespace iqs
