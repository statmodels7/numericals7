#ifndef NUMERICALS7_N7_DD_H
#define NUMERICALS7_N7_DD_H

#include <cmath>

// Double-double arithmetic (Dekker, Knuth; the algorithms of the QD
// library): a value is the unevaluated sum hi + lo with |lo| <= ulp(hi)/2,
// about 32 significant digits. Only the operations the kernels of this
// package use are defined.
namespace n7 {

struct dd {
  double hi, lo;
  dd(double h = 0.0, double l = 0.0) : hi(h), lo(l) {}
};

static inline dd dd_two_sum(double a, double b) {
  const double s = a + b, bb = s - a;
  return dd(s, (a - (s - bb)) + (b - bb));
}
static inline dd dd_quick(double a, double b) {
  const double s = a + b;
  return dd(s, b - (s - a));
}
static inline dd operator+(const dd& a, const dd& b) {
  dd s = dd_two_sum(a.hi, b.hi);
  const dd t = dd_two_sum(a.lo, b.lo);
  s.lo += t.hi;
  s = dd_quick(s.hi, s.lo);
  s.lo += t.lo;
  return dd_quick(s.hi, s.lo);
}
static inline dd operator-(const dd& a) { return dd(-a.hi, -a.lo); }
static inline dd operator-(const dd& a, const dd& b) { return a + (-b); }
// The exact product by Dekker's splitting rather than std::fma: R builds
// with -msse2 and no -mfma, so on x86 std::fma is a software routine, and
// it made the product the dominant cost (1.3 us against 0.2 per
// derivative of the Bessel ratio). The split overflows past 1e300, far
// beyond the magnitudes the kernels form.
static inline dd dd_two_prod(double a, double b) {
  const double p = a * b;
  const double ta = 134217729.0 * a, ah = ta - (ta - a), al = a - ah;
  const double tb = 134217729.0 * b, bh = tb - (tb - b), bl = b - bh;
  return dd(p, ((ah * bh - p) + ah * bl + al * bh) + al * bl);
}
static inline dd operator*(const dd& a, const dd& b) {
  dd p = dd_two_prod(a.hi, b.hi);
  p.lo += a.hi * b.lo + a.lo * b.hi;
  return dd_quick(p.hi, p.lo);
}
static inline dd operator/(const dd& a, const dd& b) {
  const double q1 = a.hi / b.hi;
  dd r = a - b * dd(q1);
  const double q2 = r.hi / b.hi;
  r = r - b * dd(q2);
  const double q3 = r.hi / b.hi;
  return dd_quick(q1, q2) + dd(q3);
}
static inline dd operator+(const dd& a, double b) { return a + dd(b); }
static inline dd operator-(double a, const dd& b) { return dd(a) - b; }
static inline dd operator*(double a, const dd& b) { return dd(a) * b; }

}  // namespace n7

#endif
