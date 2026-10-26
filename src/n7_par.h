#ifndef NUMERICALS7_N7_PAR_H
#define NUMERICALS7_N7_PAR_H

#include <cstddef>
#include <fenv.h>
#include <RcppParallel.h>

// The parallel driver of the elementwise kernels, in the shape the toolkit's
// other two drivers have (distributions7's d7_par.h carries the reasoning in
// full). The decomposition is over the elements of the OUTPUT: element i is
// computed and written by one thread, so no reduction is ever split and the
// result is bit-identical at any thread count. The worker's loop is noinline
// and the sequential branch runs THROUGH the worker, so both branches execute
// one compiled copy; and the worker installs the calling thread's
// floating-point environment, without which some of the platform's own math
// routines return per-thread last bits.
//
// A body here must not touch the R API, and that excludes more of Rmath than
// allocation does: anything that can raise a warning kills the process from a
// worker thread, which puts the p/q family, lchoose and the Bessel functions
// out of reach. It is why log_bessel_k_cpp() is NOT threaded -- its hybrid
// branch calls R::bessel_k -- while log_bessel_i_cpp(), whose every branch is
// this file's own arithmetic, is.
namespace n7 {

#if defined(__GNUC__) || defined(__clang__)
#define N7_NOINLINE __attribute__((noinline))
#elif defined(_MSC_VER)
#define N7_NOINLINE __declspec(noinline)
#else
#define N7_NOINLINE
#endif

template <typename Body>
struct BodyWorker : public RcppParallel::Worker {
  const Body& body;
  fenv_t env;
  explicit BodyWorker(const Body& b) : body(b) { fegetenv(&env); }
  N7_NOINLINE void operator()(std::size_t begin, std::size_t end) {
    fesetenv(&env);
    for (std::size_t i = begin; i < end; ++i) body(i);
  }
};

// The count is passed to parallelFor rather than left to the process-level
// setting: resolveValue() prefers an explicit positive value, so a fit that
// sized the pool through local_threads() keeps that size and a caller that
// did not gets the count it asked for.
template <typename Body>
inline void par_for(std::size_t n, int threads, std::size_t threshold,
                    const Body& body) {
  BodyWorker<Body> w(body);
  if (threads > 1 && n >= threshold) {
    RcppParallel::parallelFor(0, n, w, 1, threads);
  } else {
    w(0, n);
  }
}

// The threshold of the kernels with a full transcendental body per element,
// the one distributions7 uses for its own. It was set for log I when its
// series cost several microseconds an element; since 0.24.0 log I costs
// about 0.2 and pays from about 500 elements (log_bessel.cpp).
constexpr std::size_t kMinCostly = 128;

} // namespace n7

#endif
