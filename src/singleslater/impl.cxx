#include <singleslater/impl.hpp>

namespace ChronusQ {

  template class SingleSlater<double>;
  template class SingleSlater<dcomplex>;

  // Instantiate copy constructors
  template SingleSlater<double>::SingleSlater(  const SingleSlater<double> &);
  template SingleSlater<dcomplex>::SingleSlater(const SingleSlater<dcomplex> &);
  template SingleSlater<dcomplex>::SingleSlater(const SingleSlater<double> &);
  // Instantiate copy ructors
  template SingleSlater<double>::SingleSlater(   SingleSlater<double> &&);
  template SingleSlater<dcomplex>::SingleSlater( SingleSlater<dcomplex> &&);
  template SingleSlater<dcomplex>::SingleSlater( SingleSlater<double> &&);

}; // namespace ChronusQ
