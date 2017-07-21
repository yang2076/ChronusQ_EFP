#include <wavefunction/impl.hpp>

namespace ChronusQ {

  template class WaveFunction<double>;
  template class WaveFunction<dcomplex>;

  // Instantiate copy constructors
  template WaveFunction<dcomplex>::WaveFunction(const WaveFunction<double> &,
    int);
  // Instantiate copy ructors
  template WaveFunction<dcomplex>::WaveFunction( WaveFunction<double> &&,int);

}; // namespace ChronusQ
