#include <quantum/impl.hpp>

namespace ChronusQ {

  template class Quantum<double>;
  template class Quantum<dcomplex>;

  // Instantiate copy constructors
  template Quantum<double>::Quantum(  const Quantum<double> &);
  template Quantum<dcomplex>::Quantum(const Quantum<dcomplex> &);
  template Quantum<dcomplex>::Quantum(const Quantum<double> &);
  // Instantiate copy ructors
  template Quantum<double>::Quantum(   Quantum<double> &&);
  template Quantum<dcomplex>::Quantum( Quantum<dcomplex> &&);
  template Quantum<dcomplex>::Quantum( Quantum<double> &&);

}; // namespace ChronusQ
