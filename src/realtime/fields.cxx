/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2017 Li Research Group (University of Washington)
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  
 *  Contact the Developers:
 *    E-Mail: xsli@uw.edu
 *  
 */
#include <realtime/fields.hpp>

// This file instantiates the various field envelop definitions
namespace ChronusQ {

  /**
   *  \brief Definition of a step function field envelope
   *  
   *  \f[
   *    F(t) = \begin{cases}
   *              1 &  \mathrm{tOn} <= t <= \mathrm{tOff} \\
   *              0 &  \mathrm{else}
   *           \end{cases}
   *  \f]
   *
   *  Also handle the case where t is close to the time
   *  boundaries due to numerical precision.
   *
   *  \param [in] t Time point
   *  \returns      F(t)
   *
   */ 
  double StepField::getAmp(double t) {
    if( (t >= this->tOn) and (t <= this->tOff) or 
        (std::abs(t - this->tOn)  < 1e-10) or 
        (std::abs(t - this->tOff) < 1e-10) ) return 1.;
    else return 0.;
  };


}; // namespace ChronusQ