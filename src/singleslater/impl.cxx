#include <singleslater/impl.hpp>

namespace ChronusQ {

  template class SingleSlater<double>;
  template class SingleSlater<dcomplex>;

  // Instantiate copy constructors
  template SingleSlater<dcomplex>::SingleSlater(const SingleSlater<double> &,
    int);

  // Instantiate copy ructors
  template SingleSlater<dcomplex>::SingleSlater( SingleSlater<double> &&, int);

}; // namespace ChronusQ
