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
#include <cxxapi/options.hpp>
#include <cerr.hpp>

namespace ChronusQ {

  void CQSCFOptions(std::ostream &out, CQInputFile &input,
    SingleSlaterBase &ss) {

    // SCF section not required
    if( not input.containsSection("SCF") ) return;

    // Optionally parse guess

    // Energy convergence tolerance
    OPTOPT( ss.scfControls.eneConvTol = 
              input.getData<double>("SCF.ENETOL"); )

    // Energy convergence tolerance
    OPTOPT( ss.scfControls.denConvTol = 
              input.getData<double>("SCF.DENTOL"); )

    // Maximum SCF iterations
    OPTOPT( ss.scfControls.maxSCFIter = 
              input.getData<size_t>("SCF.MAXITER"); )


    // Toggle extrapolation in its entireity
    OPTOPT(
      ss.scfControls.doExtrap = 
        input.getData<bool>("SCF.EXTRAP");
    )

    // Handle DIIS options
    OPTOPT(
      bool doDIIS = input.getData<bool>("SCF.DIIS");
      if( not doDIIS )
        ss.scfControls.diisAlg = NONE;
    );

    // Number of terms for keep for DIIS
    OPTOPT( ss.scfControls.nKeep = 
              input.getData<size_t>("SCF.NKEEP"); )

    // Parse Damping options
      
    OPTOPT(
      ss.scfControls.doDamp = input.getData<bool>("SCF.DAMP");
    );

    OPTOPT(
      ss.scfControls.dampStartParam = 
        input.getData<double>("SCF.DAMPPARAM");
    );

    OPTOPT(
      ss.scfControls.dampError = 
        input.getData<double>("SCF.DAMPERROR");
    );


    // Handling eqivalences in the input options


    // Setting the damp param to 0. is equivalent to
    // turning damping off
    if( ss.scfControls.dampStartParam == 0. )
      ss.scfControls.doDamp = false;

    // Turning off both damping and DIIS is equivalent
    // to turning off extrapolation entirely
    if( not ss.scfControls.doDamp and 
        ss.scfControls.diisAlg == NONE )
      ss.scfControls.doExtrap = false;


  }; // CQSCFOptions

}; // namespace ChronusQ
