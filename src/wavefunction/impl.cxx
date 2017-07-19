#include <wavefunction/impl.hpp>

namespace ChronusQ {

  template class WaveFunction<double>;
  template class WaveFunction<dcomplex>;

  // Instantiate copy constructors
  template WaveFunction<double>::WaveFunction(  const WaveFunction<double> &);
  template WaveFunction<dcomplex>::WaveFunction(const WaveFunction<dcomplex> &);
  template WaveFunction<dcomplex>::WaveFunction(const WaveFunction<double> &);
  // Instantiate copy ructors
  template WaveFunction<double>::WaveFunction(   WaveFunction<double> &&);
  template WaveFunction<dcomplex>::WaveFunction( WaveFunction<dcomplex> &&);
  template WaveFunction<dcomplex>::WaveFunction( WaveFunction<double> &&);

}; // namespace ChronusQ
