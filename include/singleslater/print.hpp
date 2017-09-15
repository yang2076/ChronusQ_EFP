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
#ifndef __INCLUDED_SINGLESLATER_PRINT_HPP__
#define __INCLUDED_SINGLESLATER_PRINT_HPP__

#include <singleslater.hpp>
#include <util/matout.hpp>


namespace ChronusQ {

  template <typename T>
  void SingleSlater<T>::printFock(std::ostream &out) {

    size_t NB = aoints.basisSet().nBasis;

    prettyPrintSmart(out,"Fock (AO) Scalar",fock[SCALAR],NB,NB,NB);

    if( fock.size() > 1 )
      prettyPrintSmart(out,"Fock (AO) MZ",fock[MZ],NB,NB,NB);

    if( fock.size() > 2 ) {
      prettyPrintSmart(out,"Fock (AO) MY",fock[MY],NB,NB,NB);
      prettyPrintSmart(out,"Fock (AO) MX",fock[MX],NB,NB,NB);
    }


  }; // SingleSlater<T>::printFock

  template <typename T>
  void SingleSlater<T>::print1PDMOrtho(std::ostream &out) {

    size_t NB = aoints.basisSet().nBasis;

    prettyPrintSmart(out,"1PDM (Ortho) Scalar",onePDMOrtho[SCALAR],NB,NB,NB);

    if( onePDMOrtho.size() > 1 )
      prettyPrintSmart(out,"1PDM (Ortho) MZ",onePDMOrtho[MZ],NB,NB,NB);

    if( onePDMOrtho.size() > 2 ) {
      prettyPrintSmart(out,"1PDM (Ortho) MY",onePDMOrtho[MY],NB,NB,NB);
      prettyPrintSmart(out,"1PDM (Ortho) MX",onePDMOrtho[MX],NB,NB,NB);
    }


  }; // SingleSlater<T>::print1PDMOrtho

  template <typename T>
  void SingleSlater<T>::printGD(std::ostream &out) {

    size_t NB = aoints.basisSet().nBasis;

    prettyPrintSmart(out,"GD (AO) Scalar",GD[SCALAR],NB,NB,NB);

    if( GD.size() > 1 )
      prettyPrintSmart(out,"GD (AO) MZ",GD[MZ],NB,NB,NB);

    if( GD.size() > 2 ) {
      prettyPrintSmart(out,"GD (AO) MY",GD[MY],NB,NB,NB);
      prettyPrintSmart(out,"GD (AO) MX",GD[MX],NB,NB,NB);
    }


  }; // SingleSlater<T>::printGD


  template <typename T>
  void SingleSlater<T>::printJ(std::ostream &out) {

    size_t NB = aoints.basisSet().nBasis;

    prettyPrintSmart(out,"J (AO) Scalar",JScalar,NB,NB,NB);


  }; // SingleSlater<T>::printJ


  template <typename T>
  void SingleSlater<T>::printK(std::ostream &out) {

    size_t NB = aoints.basisSet().nBasis;

    prettyPrintSmart(out,"K (AO) Scalar",K[SCALAR],NB,NB,NB);

    if( K.size() > 1 )
      prettyPrintSmart(out,"K (AO) MZ",K[MZ],NB,NB,NB);

    if( K.size() > 2 ) {
      prettyPrintSmart(out,"K (AO) MY",K[MY],NB,NB,NB);
      prettyPrintSmart(out,"K (AO) MX",K[MX],NB,NB,NB);
    }


  }; // SingleSlater<T>::printK


  template <typename T>
  void SingleSlater<T>::printMiscProperties(std::ostream &out) {

    out << "\nCharge Analysis:\n" << bannerTop << "\n\n";

    out << std::setw(20) << std::left << "  Atom";
    out << std::setw(20) << std::right << "Mulliken Charges";
    out << std::setw(20) << std::right << "Lowdin Charges" << std::endl;

    out << std::right << bannerMid << std::endl;

    Molecule & mol = aoints.molecule();
    for(auto iAtm = 0; iAtm < mol.nAtoms; iAtm++) {

      // Get symbol
      std::map<std::string,Atom>::const_iterator it = 
      std::find_if(atomicReference.begin(),atomicReference.end(),
        [&](const std::pair<std::string,Atom> &st){ 
          return (st.second.atomicNumber == mol.atoms[iAtm].atomicNumber) and
                 (st.second.massNumber == mol.atoms[iAtm].massNumber);}
         );

      out << "  " << std::setw(18) << std::left <<
        (it == atomicReference.end() ? "X" : it->first);

      out << std::setprecision(5) << std::right;

      out << std::setw(20) << mullikenCharges[iAtm];
      out << std::setw(20) << lowdinCharges[iAtm];

      out << std::endl;
    }

    out << std::endl << bannerEnd << std::endl;

  }; // SingleSlater<T>::printMiscProperties

  template <typename T>
  void SingleSlater<T>::printMOInfo(std::ostream &out) {

    out << std::scientific << std::setprecision(4);

    out << "\n\n" << "SCF Results:\n" << BannerTop << "\n\n";

    out << "Orbital Eigenenergies " << (this->nC == 1 ? "(Alpha) " : "" )
        << "/ Eh\n" << bannerTop << "\n";

    size_t NO = (this->nC == 1 ? this->nOA : this->nO);
    for(auto i = 0ul; i < this->nC * aoints.basisSet().nBasis; i++) {

      if( i == 0 )
        out << "Occupied:\n";
      else if( i == NO )
        out << "\n\nVirtual:\n";

      out << std::setw(13) << this->eps1[i];

      if( i < NO and (i + 1) % 5 == 0 )  out << "\n";
      else if( i >= NO and ((i - NO) + 1) % 5 == 0 ) out << "\n";
    }

     
    out << "\n" << bannerEnd << "\n";

    if( this->nC == 1 and not this->iCS ) {
      out << "\n\nOrbital Eigenenergies (Beta) / Eh\n" << bannerTop << "\n";

      for(auto i = 0ul; i < aoints.basisSet().nBasis; i++) {

        if( i == 0 )
          out << "Occupied:\n";
        else if( i == this->nOB )
          out << "\n\nVirtual:\n";

        out << std::setw(13) << this->eps2[i];

        if( i < this->nOB and (i + 1) % 5 == 0 )  out << "\n";
        else if( i >= this->nOB and ((i - this->nOB) + 1) % 5 == 0 ) out << "\n";
      }

      out << "\n" << bannerEnd << "\n";
    }

    out << "\n" << BannerEnd << "\n\n";
  }; // SingleSlater<T>::printMOInfo

}; // namespace ChronusQ

#endif
