#include <quantum/impl.hpp>

namespace ChronusQ {

  template class Quantum<double>;
  template class Quantum<dcomplex>;

  // Instantiate converting constructors
  template Quantum<dcomplex>::Quantum(const Quantum<double> &, int);

  // Instantiate converting constructors
  template Quantum<dcomplex>::Quantum( Quantum<double> &&, int);

}; // namespace ChronusQ
