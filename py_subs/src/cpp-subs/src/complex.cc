#include "trm/complex.h"

Subs::Complex Subs::operator*(const Complex& c1, const Complex& c2) {
  return Complex(c1.real_*c2.real_ - c1.imag_*c2.imag_, c1.real_*c2.imag_ + c1.imag_*c2.real_);
}

