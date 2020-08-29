#include "../include/spec_kernels.hpp"

// Declare loops
#define PARALLEL_FOR_1D                                   \
  _Pragma("omp parallel for")                             \
  for (std::size_t ind0 = gstart; ind0 < gend; ind0++)

  #define SERIAL_FOR_1D                                  \
  for (std::size_t ind0 = gstart; ind0 < gend; ind0++)

#define PARALLEL_FOR_2D                                   \
  _Pragma("omp parallel for collapse(2)")                 \
  for(std::size_t group = gstart; group < gend;           \
      group += (1L << pos + 1L))                          \
    for(std::size_t ind0 = 0; ind0 < (1L << pos); ind0++)

#define SERIAL_FOR_2D                                     \
  for(std::size_t group = gstart; group < gend;           \
      group += (1L << pos + 1L))                          \
    for(std::size_t ind0 = 0; ind0 < (1L << pos); ind0++)

#define PARALLEL_FOR_3D                             \
  _Pragma("omp parallel for collapse(3)")           \
  for (std::size_t l1 = c11; l1 < c12; l1 += c13)   \
    for (std::size_t l2 = c21; l2 < c22; l2 += c23) \
      for (std::size_t l3 = c31; l3 < c32; l3++)		

#define SERIAL_FOR_3D                               \
  for (std::size_t l1 = c11; l1 < c12; l1 += c13)   \
    for (std::size_t l2 = c21; l2 < c22; l2 += c23) \
      for (std::size_t l3 = c31; l3 < c32; l3++)

// Declare gate bodies
#define HADAMARD_BODY_2D {                    \
    std::size_t i0 = ind0 + indsht0 + group;  \
    std::size_t i1 = ind0 + indsht1 + group;  \
    Type in0 = state0[i0], in1 = state1[i1];  \
    state0[i0] = (in0 + in1) * isqrt2;        \
    state1[i1] = (in0 - in1) * isqrt2;        \
}

#define RX_BODY_2D {                                                          \
    std::size_t i0 = ind0 + indsht0 + group;                                  \
    std::size_t i1 = ind0 + indsht1 + group;                                  \
    Type in0 = state0[i0], in1 = state1[i1];                                  \
    state0[i0] = cos_2 * in0 + Type(sin_2 * in1.imag(), msin_2 * in1.real()); \
    state1[i1] = Type(sin_2 * in0.imag(), msin_2 * in0.real()) + cos_2 * in1; \
}

#define RY_BODY_2D {                                \
    std::size_t i0 = ind0 + indsht0 + group;        \
    std::size_t i1 = ind0 + indsht1 + group;        \
    Type in0 = state0[i0], in1 = state1[i1];        \
    state0[i0] = cos_2 * in0 + msin_2 * in1;        \
    state1[i1] = sin_2 * in0 + cos_2 * in1;         \
}

#define RZ_BODY_2D {                                \
    std::size_t i0 = ind0 + indsht0 + group;        \
    std::size_t i1 = ind0 + indsht1 + group;        \
    Type in0 = state0[i0], in1 = state1[i1];        \
    state0[i0] = Type(cos_2, msin_2) * in0;         \
    state1[i1] = Type(cos_2, sin_2) * in1;          \
}

#define PX_BODY_2D {                          \
    std::size_t i0 = ind0 + indsht0 + group;  \
    std::size_t i1 = ind0 + indsht1 + group;  \
    Type in0 = state0[i0], in1 = state1[i1];  \
    state0[i0] = in1;                         \
    state1[i1] = in0;                         \
}

// 1q pauli y gate body
#define PY_BODY_2D {                            \
    std::size_t i0 = ind0 + indsht0 + group;    \
    std::size_t i1 = ind0 + indsht1 + group;    \
    Type in0 = state0[i0], in1 = state1[i1];    \
    state0[i0] = Type(in1.imag(), -in1.real()); \
    state1[i1] = Type(-in0.imag(), in0.real()); \
}

// 1q pauli z gate body
#define PZ_BODY_2D {                          \
    std::size_t i1 = ind0 + indsht1 + group;  \
    state1[i1] = -state1[i1];                 \
}

#define T_BODY_2D {                           \
    std::size_t i1 = ind0 + indsht1 + group;  \
    state1[i1] *= texp;                       \
}

// Controlled Gates
#define RX_BODY_3D {                                                           \
    std::size_t ind0 = l1 + l2 + l3;                                           \
    std::size_t ind1 = ind0 + index_shift;                                     \
    Type in0 = state[ind0], in1 = state[ind1];                                 \
    state[ind0] = cos_2 * in0 + Type(sin_2 * in1.imag(), msin_2 * in1.real()); \
    state[ind1] = Type(sin_2 * in0.imag(), msin_2 * in0.real()) + cos_2 * in1; \
}

#define RY_BODY_3D {                               \
    std::size_t ind0 = l1 + l2 + l3;               \
    std::size_t ind1 = ind0 + index_shift;         \
    Type in0 = state[ind0], in1 = state[ind1];     \
    state[ind0] = cos_2 * in0 + msin_2 * in1;      \
    state[ind1] = sin_2 * in0 + cos_2 * in1;       \
}

#define RZ_BODY_3D {                              \
    std::size_t ind0 = l1 + l2 + l3;              \
    std::size_t ind1 = ind0 + index_shift;        \
    Type in0 = state[ind0], in1 = state[ind1];    \
    state[ind0] = Type(cos_2, msin_2) * in0;     \
    state[ind1] = Type(cos_2, sin_2) * in1;      \
}

#define PX_BODY_3D {                        \
    std::size_t ind0 = l1 + l2 + l3;        \
    std::size_t ind1 = ind0 + index_shift;  \
    std::swap(state[ind0], state[ind1]);    \
}

#define PY_BODY_3D {                                 \
    std::size_t ind0 = l1 + l2 + l3;                 \
    std::size_t ind1 = ind0 + index_shift;           \
    Type in0 = state[ind0], in1 = state[ind1];       \
    state[ind0] = Type(in1.imag(), -in1.real());     \
    state[ind1] = Type(-in0.imag(), in0.real());     \
}

#define PZ_BODY_3D {                        \
    std::size_t ind0 = l1 + l2 + l3;        \
    std::size_t ind1 = ind0 + index_shift;  \
    state[ind1] = -state[ind1];             \
}

#define HADAMARD_BODY_3D {                      \
    std::size_t ind0 = l1 + l2 + l3;            \
    std::size_t ind1 = ind0 + index_shift;      \
    Type in0 = state[ind0], in1 = state[ind1];  \
    state[ind0] = (in0 + in1) * isqrt2;         \
    state[ind1] = (in0 - in1) * isqrt2;         \
}

#define CP_BODY_3D {                        \
    std::size_t ind0 = l1 + l2 + l3;        \
    std::size_t ind1 = ind0 + index_shift;  \
    state[ind1] *= pexp;                    \
}

// Kernels

template< typename Type >
__attribute__((noinline))
void Loop_SN(std::size_t gstart, std::size_t gend,
             Type *state0, Type *state1,
             std::size_t indsht0, std::size_t indsht1,
	           GateSpec1Q spec, Timer *timer, double angle)
{
  assert((UL(state0) % 256) == 0);
  assert((UL(state1) % 256) == 0);
#if defined(__ICC) || defined(__INTEL_COMPILER)
  __assume_aligned(state0, 256);
  __assume_aligned(state1, 256);
#endif

  // Declare constants
  const auto theta = static_cast<decltype(state0[0].real())>(angle);
  const decltype(theta) isqrt2 = 1 / std::sqrt(2);
  const decltype(theta) cos_2 = std::cos(theta / 2);
  const decltype(theta) sin_2 = std::sin(theta / 2);
  const decltype(theta) msin_2 = -sin_2;
  const Type texp = Type(std::cos(M_PI / 4), std::sin(M_PI / 4)); 

  constexpr size_t group = 0;

  size_t nthreads = 1;
#ifdef _OPENMP
#pragma omp parallel 
  nthreads = omp_get_num_threads();
#endif
  bool par = nthreads > 1;
  // printf("Executing here with %d threads...\n", nthreads);
  switch(spec) {

   case GateSpec1Q::Hadamard:
     if (par) { PARALLEL_FOR_1D HADAMARD_BODY_2D; }
     else { SERIAL_FOR_1D HADAMARD_BODY_2D; }
     break;

  case GateSpec1Q::RotationX:
    if (par) { PARALLEL_FOR_1D RX_BODY_2D; }
    else { SERIAL_FOR_1D RX_BODY_2D; }
    break;

  case GateSpec1Q::RotationY:
    if (par) { PARALLEL_FOR_1D RY_BODY_2D; }
    else { SERIAL_FOR_1D RY_BODY_2D; }
    break;

  case GateSpec1Q::RotationZ:
    if (par) { PARALLEL_FOR_1D RZ_BODY_2D; }
    else { SERIAL_FOR_1D RZ_BODY_2D; }
    break;

  case GateSpec1Q::PauliX:
    if (par) { PARALLEL_FOR_1D PX_BODY_2D; }
    else { SERIAL_FOR_1D PX_BODY_2D; }
    break;

  case GateSpec1Q::PauliY:
    if (par) { PARALLEL_FOR_1D PY_BODY_2D; }
    else { SERIAL_FOR_1D PY_BODY_2D; }
    break;
  
  case GateSpec1Q::PauliZ:
    if (par) { PARALLEL_FOR_1D PZ_BODY_2D; }
    else { SERIAL_FOR_1D PZ_BODY_2D; }
    break;

  case GateSpec1Q::T:
    if (par) { PARALLEL_FOR_1D T_BODY_2D; }
    else { SERIAL_FOR_1D T_BODY_2D; }
    break;

   default:
     break;
 }
}


template < typename Type >
__attribute__((noinline))
void Loop_DN(std::size_t gstart, std::size_t gend, std::size_t pos,
             Type *state0, Type *state1,
             std::size_t indsht0, std::size_t indsht1,
	           GateSpec1Q spec, Timer *timer, double angle)
{
  double ttmp1 = sec(), ttot = 0.;
  assert((UL(state0) % 256) == 0);
  assert((UL(state1) % 256) == 0);
#if defined(__ICC) || defined(__INTEL_COMPILER)
  __assume_aligned(state0, 256);
  __assume_aligned(state1, 256);
#endif

  // Declare constants
  const auto theta = static_cast<decltype(state0[0].real())>(angle);
  const decltype(theta) isqrt2 = 1 / std::sqrt(2);
  const decltype(theta) cos_2 = std::cos(theta / 2);
  const decltype(theta) sin_2 = std::sin(theta / 2);
  const decltype(theta) msin_2 = -sin_2;
  const Type texp = Type(std::cos(M_PI / 4), std::sin(M_PI / 4)); 

  size_t nthreads = 1;
#ifdef _OPENMP
#pragma omp parallel 
  nthreads = omp_get_num_threads();
#endif
  bool par = nthreads > 1;

  switch(spec) {

   case GateSpec1Q::Hadamard:
     if (par) { PARALLEL_FOR_2D HADAMARD_BODY_2D; }
     else { SERIAL_FOR_2D HADAMARD_BODY_2D; }
     break;

  case GateSpec1Q::RotationX:
    if (par) { PARALLEL_FOR_2D RX_BODY_2D; }
    else { SERIAL_FOR_2D RX_BODY_2D; }
    break;

  case GateSpec1Q::RotationY:
    if (par) { PARALLEL_FOR_2D RY_BODY_2D; }
    else { SERIAL_FOR_2D RY_BODY_2D; }
    break;

  case GateSpec1Q::RotationZ:
    if (par) { PARALLEL_FOR_2D RZ_BODY_2D; }
    else { SERIAL_FOR_2D RZ_BODY_2D; }
    break;

  case GateSpec1Q::PauliX:
    if (par) { PARALLEL_FOR_2D PX_BODY_2D; }
    else { SERIAL_FOR_2D PX_BODY_2D; }
    break;

  case GateSpec1Q::PauliY:
    if (par) { PARALLEL_FOR_2D PY_BODY_2D; }
    else { SERIAL_FOR_2D PY_BODY_2D; }
    break;
  
  case GateSpec1Q::PauliZ:
    if (par) { PARALLEL_FOR_2D PZ_BODY_2D; }
    else { SERIAL_FOR_2D PZ_BODY_2D; }
    break;

  case GateSpec1Q::T:
    if (par) { PARALLEL_FOR_2D T_BODY_2D; }
    else { SERIAL_FOR_2D T_BODY_2D; }
    break;

   default:
    // This should not happen!
    throw std::runtime_error("InvalidArgument: Loop_DN SpecializeV2 is called with GateSpec1Q::None!");
 }

 if(timer)
  {
      ttot = sec() - ttmp1;     
      double datab = 2.0 * sizeof(state0[0]) * D(gend - gstart);
      timer->record_dn(ttot, datab / ttot);
  }
}

template <class Type>
__attribute__((noinline))
void Loop_TN(Type *state,
             std::size_t c11, std::size_t c12, std::size_t c13,
             std::size_t c21, std::size_t c22, std::size_t c23,
             std::size_t c31, std::size_t c32, 
             std::size_t index_shift, GateSpec2Q spec,
             Timer *timer, double angle)
{
  double ttmp1 = sec(), ttot = 0.;
  assert((UL(state) % 256) == 0);
#if defined(__ICC) || defined(__INTEL_COMPILER)
  __assume_aligned(state, 256);
#endif
  size_t nthreads = 1;
#ifdef _OPENMP
#pragma omp parallel 
  nthreads = omp_get_num_threads();
#endif
  bool par = nthreads > 1;

   // Declare constants
  const auto theta = static_cast<decltype(state[0].real())>(angle);
  const decltype(theta) isqrt2 = 1 / std::sqrt(2);
  const decltype(theta) cos_2 = std::cos(theta / 2);
  const decltype(theta) sin_2 = std::sin(theta / 2);
  const decltype(theta) msin_2 = -sin_2;
  const Type pexp = Type(std::cos(theta), std::sin(theta));

  switch(spec) {

    case GateSpec2Q::CHadamard:
      if (par) { PARALLEL_FOR_3D HADAMARD_BODY_3D; }
      else { SERIAL_FOR_3D HADAMARD_BODY_3D; }
      break;

    case GateSpec2Q::CRotationX:
      if (par) { PARALLEL_FOR_3D RX_BODY_3D; }
      else { SERIAL_FOR_3D RX_BODY_3D; }
      break;

    case GateSpec2Q::CRotationY:
      if (par) { PARALLEL_FOR_3D RY_BODY_3D; }
      else { SERIAL_FOR_3D RY_BODY_3D; }
      break;

    case GateSpec2Q::CRotationZ:
      if (par) { PARALLEL_FOR_3D RZ_BODY_3D; }
      else { SERIAL_FOR_3D RZ_BODY_3D; }
      break;

    case GateSpec2Q::CPauliX:
      if (par) { PARALLEL_FOR_3D PX_BODY_3D; }
      else { SERIAL_FOR_3D PX_BODY_3D; }
      break;

    case GateSpec2Q::CPauliY:
      if (par) { PARALLEL_FOR_3D PY_BODY_3D; }
      else { SERIAL_FOR_3D PY_BODY_3D; }
      break;

    case GateSpec2Q::CPauliZ:
      if (par) { PARALLEL_FOR_3D PZ_BODY_3D; }
      else { SERIAL_FOR_3D PZ_BODY_3D; }
      break;
    
    case GateSpec2Q::CPhase:
      if (par) { PARALLEL_FOR_3D CP_BODY_3D; }
      else { SERIAL_FOR_3D CP_BODY_3D; }
      break;

    default: 
      // This should not happen!
      throw std::runtime_error("InvalidArgument: Loop_TN SpecializeV2 is called with GateSpec2Q::None!");
  }

  if (timer)
  {
    ttot = sec() - ttmp1;
    double datab =
      4.0 * sizeof(state[0]) * D((c12 - c11) / c13) * D((c22 - c21) / c23) * D(c32 - c31);
    timer->record_tn(ttot, datab / ttot);
  }
}

// Declarations
template void Loop_SN<ComplexSP>(std::size_t gstart, std::size_t gend,
                                 ComplexSP *state0, ComplexSP *state1,
                                 std::size_t indsht0, std::size_t indsht1,
                                 GateSpec1Q spec, Timer *timer, double angle);

template void Loop_SN<ComplexDP>(std::size_t gstart, std::size_t gend,
                                 ComplexDP *state0, ComplexDP *state1,
                                 std::size_t indsht0, std::size_t indsht1,
                                 GateSpec1Q spec, Timer *timer, double angle);

template void Loop_DN<ComplexSP>(std::size_t gstart, std::size_t gend, std::size_t pos,
                                 ComplexSP *state0, ComplexSP *state1,
                                 std::size_t indsht0, std::size_t indsht1,
                                 GateSpec1Q spec, Timer *timer, double angle);

template void Loop_DN<ComplexDP>(std::size_t gstart, std::size_t gend, std::size_t pos,
                                 ComplexDP *state0, ComplexDP *state1,
                                 std::size_t indsht0, std::size_t indsht1,
                                 GateSpec1Q spec, Timer *timer, double angle);

template void Loop_TN<ComplexSP>(ComplexSP *state,
                                 std::size_t c11, std::size_t c12, std::size_t c13,
                                 std::size_t c21, std::size_t c22, std::size_t c23,
                                 std::size_t c31, std::size_t c32, 
                                 std::size_t index_shift, GateSpec2Q spec, 
                                 Timer *timer, double angle);

template void Loop_TN<ComplexDP>(ComplexDP *state,
                                 std::size_t c11, std::size_t c12, std::size_t c13,
                                 std::size_t c21, std::size_t c22, std::size_t c23,
                                 std::size_t c31, std::size_t c32, 
                                 std::size_t index_shift, GateSpec2Q spec, 
                                 Timer *timer, double angle);
