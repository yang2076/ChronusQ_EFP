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
#ifndef __INCLUDED_BASISSET_DEF_HPP__
#define __INCLUDED_BASISSET_DEF_HPP__

#include <chronusq_sys.hpp>
#include <util/typedefs.hpp>
#include <memmanager.hpp>
#include <molecule.hpp>

#include <libint2/shell.h>

namespace ChronusQ {

  /**
   *  \brief The BasisSet struct. Contains information pertinant
   *  for a gaussian type orbital (GTO) basis set. 
   *
   *  Precomputes useful intermediates relating to the basis set
   *  as well as stores useful constants (nBasis, nPrimitive) and 
   *  maps (maps shell # -> basis function #, etc). Utilizes
   *  a storage format consistent with the Libint2  integral
   *  package (libint2::Shell)
   */ 
  struct BasisSet {

    std::string basisName;

    bool forceCart;

    size_t nBasis;      ///< Number of CGTO basis functions
    size_t nPrimitive;  ///< Number of primitive GTO functions
    size_t nShell;      ///< Number of CGTO basis shells
    size_t maxPrim;     ///< Max primitive dimension of basis
    size_t maxL;        ///< Max angular momentum of basis

    cartvec_t centers; ///< A list of centers that comprise the BasisSet

    std::vector<libint2::Shell> shells; ///< Basis shells

    std::vector<std::vector<double>> unNormCont;
      ///< Unnormalized basis coefficients
        
    std::vector<size_t> mapSh2Bf;  ///< Map Shell # -> BF #
    std::vector<size_t> mapSh2Cen; ///< Map Shell # -> Cen #
    std::vector<size_t> mapCen2BfSt; ///< Map Cen # -> Starting BF #

    // Disable default constructor
    BasisSet() = delete;


    // Path / Molecule constructor.
    // See src/basisset/basisset.cxx for documentation
    BasisSet(std::string _basisName, const Molecule &mol, 
      bool _forceCart = false, bool doPrint = true);
       
    /**
     *  Copy constructor.
     *
     *  Copies the contents of one BasisSet object to another
     */ 
    BasisSet(const BasisSet &) = default;


    /**
     *  Move constructor.
     *
     *  Moves the contents of one BasisSet object to another
     */ 
    BasisSet(BasisSet &&) = default;


    /**
     *  Copy assignment operator
     *
     *  Assigns one BasisSet object to another through a copy
     */ 
    BasisSet& operator=(const BasisSet &) = default;


    /**
     *  Move assignment operator
     *
     *  Assigns one BasisSet object to another through a move (rvalue reference)
     */ 
    BasisSet& operator=(BasisSet &&) = default;


    // Ouputs BasisSet info, see src/basisset/basisset.cxx for 
    // documentation
    friend std::ostream& operator<<(std::ostream &, const BasisSet&);

    // Misc functions
      
    std::vector<libint2::Shell> uncontractShells();
    void makeMapPrim2Cont(double *, double *, CQMemManager&);


    private:

      // Update BasisSet object member data (nBasis, etc).
      // See src/basisset/basisset.cxx for documentation
      void update();

  }; // BasisSet struct

// SS start  
  extern std::vector<std::vector<std::array<int,3>>> cart_ang_list; //list of cartesian angular momentum.
  extern std::vector<std::vector<double>> car2sph_matrix; //SS: external variable,
          //vector of transform matrix from cartesian to spherical
  
  void pop_cart_ang_list(); //SS:populate angular momentum list, can only be called once.   

  void pop_car2sph_matrix(); //populate transform matrix up to L=6

  std::complex<double> cart2sphCoeff(int, int, int, int, int); 
          //SS: calculate the matrix elements of cartesian to spherical transform

  double factorial(int);
  
  double doubleFact(int);

  double polyCoeff(int, int);

  void cart2sph_transform( int, int, std::vector<double>&, std::vector<double>& ); 

// SS end
 
}; // namespace ChronusQ

#endif
