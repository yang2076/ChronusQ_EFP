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
#ifndef __INCLUDED_OPTIONS_HPP__
#define __INCLUDED_OPTIONS_HPP__

#include <cxxapi/input.hpp>

#include <molecule.hpp>
#include <basisset.hpp>

// Preprocessor directive to aid the digestion of optional input arguements
#define OPTOPT(x) try{ x; } catch(...) { ; }

namespace ChronusQ {

  // Function definitions ofr option parsing. See src/cxxapi/input/*opts.cxx
  // for documentation

  // Parse the options relating to the Molecule object
  Molecule CQMoleculeOptions(std::ostream &, CQInputFile &);

  // Parse the options relating to the BasisSet
  BasisSet CQBasisSetOptions(std::ostream &, CQInputFile &, Molecule &);
};

#endif
