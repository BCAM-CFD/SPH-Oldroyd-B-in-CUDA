/******************************************************
This code has been developed by Adolfo Vazquez-Quesada,
from the Department of Fundamental Physics at UNED, in
Madrid, Spain.
email: a.vazquez-quesada@fisfun.uned.es
********************************************************/

// ===========================================================
// Universal implementation universal of atomicAdd for double precission
// For old GPUs: compute capability < 6.0
// ===========================================================

#ifndef ATOMIC_ADD_DOUBLE_CUH
#define ATOMIC_ADD_DOUBLE_CUH

// Universal emulation based on atomicCAS
__device__ double atomicAdd_double_emulated(double* address, double val)
{
  unsigned long long int* addr_as_ull =
    (unsigned long long int*)address;

  unsigned long long int old = *addr_as_ull, assumed;

  do {
    assumed = old;
    old = atomicCAS(
		    addr_as_ull,
		    assumed,
		    __double_as_longlong(
					 __longlong_as_double(assumed) + val
					 )
		    );
  } while (assumed != old);

  return __longlong_as_double(old);
}

#endif

