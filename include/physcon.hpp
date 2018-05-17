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
#ifndef __INCLUDE_PHYSCON_HPP
#define __INCLUDE_PHYSCON_HPP

namespace ChronusQ {

  // Physical Constants 
  // XXX: Could use some citations / more accurate values here
  constexpr double AngPerBohr    = 0.52917721092;
  constexpr double EBohrPerDebye = 0.393430307;
  constexpr double EVPerHartree  = 27.211396132;
  constexpr double NMPerHartree  = 45.56335;
  constexpr double SpeedOfLight  = 137.035999074;
  constexpr double FSPerAUTime   = 2.41884326505e-2;

};

#endif
