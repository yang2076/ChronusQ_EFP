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
#ifndef __INCLUDED_SINGLESLATERBASE_SCF_HPP__
#define __INCLUDED_SINGLESLATERBASE_SCF_HPP__

#include <singleslater.hpp>
#include <util/matout.hpp>

namespace ChronusQ {

  /**
   *  \brief Performs the self--consistant field procedure given a set of 
   *  orbitals.
   *
   *  \warning SCF procedure assumes that the 1PDM and orbital (mo1/2) storage
   *  has been populated in some way.
   */ 
  void SingleSlaterBase::SCF() {

    SCFInit();

    // Initialize type independent parameters
    bool isConverged = false;
    scfControls.dampParam = scfControls.dampStartParam;
    scfControls.doIncFock = scfControls.doIncFock and (aoints.cAlg == DIRECT);

    printSCFHeader(std::cout);

    for( scfConv.nSCFIter = 0; scfConv.nSCFIter < scfControls.maxSCFIter; 
         scfConv.nSCFIter++) {

      // Save current state of the wave function (method specific)
      saveCurrentState();

      // Exit loop on convergence
      if(isConverged) break;

      // Get new orbtials and densities from current state: 
      //   C/D(k) -> C/D(k + 1)
      getNewOrbitals();

      // Evaluate convergence
      isConverged = evalConver();

      // Print out iteration information
      printSCFProg(std::cout);

    }; // Iteration loop

    // Save current state of the wave function (method specific)
    saveCurrentState();

    SCFFin();


    //printSCFFooter(isConverged);
    if(not isConverged)
      CErr(std::string("SCF Failed to converged within ") + 
        std::to_string(scfControls.maxSCFIter) + 
        std::string(" iterations"));
    else {
      std::cout << std::endl << "SCF Completed: E("
                << refShortName_ << ") = " << std::fixed
                << std::setprecision(10) << this->totalEnergy
                << " Eh after " << scfConv.nSCFIter
                << " SCF Iterations" << std::endl;
    }
    std::cout << BannerEnd << std::endl;
    
  }; // SingleSlater<T>::SCF()

  
  void SingleSlaterBase::printSCFHeader(std::ostream &out) {
    out << BannerTop << std::endl;
    out << "Self Consistent Field (SCF) Settings:" << std::endl << std::endl;

    out << std::setw(38) << std::left << "  Reference:" << refLongName_ 
           << std::endl;

    out << std::setprecision(6) << std::scientific;
    out << std::setw(38) << std::left << "  Density Convergence Tolerence:" 
           <<  scfControls.denConvTol << std::endl;

    out << std::setw(38) << std::left << "  Energy Convergence Tolerence:" 
           <<  scfControls.eneConvTol << std::endl;

    out << std::setw(38) << std::left << "  Maximum Number of SCF Cycles:" 
           << scfControls.maxSCFIter << std::endl;



    if (scfControls.doExtrap) {
      if (scfControls.doDamp) {
        out << std::setw(38)   << std::left << "  Static Damping Factor:" 
               <<  scfControls.dampParam << std::endl;
        out << std::setw(38)   << std::left << "  Damping Error:" 
               <<  scfControls.dampError << std::endl;
      }

      if (scfControls.diisAlg != NONE) {
        out << std::setw(38) << std::left << "  DIIS Extrapolation Algorithm:";
        if (scfControls.diisAlg == CDIIS) out << "CDIIS";
        out << std::endl;

        out << std::left << "    * DIIS will track up to " 
            << scfControls.nKeep << " previous iterations" << std::endl;
      }
 
    } else {
        out << std::setw(38)   << std::left << "  SCF Algorithm:"
               <<  "Standard Roothaan-Hall" << std::endl;
    }


    if( scfControls.doIncFock ) {
      out << "\n  * Will Perform Incremental Fock Build -- Restarting Every "
          << scfControls.nIncFock << " SCF Steps\n";
         
    }




    out << std::endl << BannerMid << std::endl << std::endl;
    out << std::setw(16) << "SCF Iteration";
    out << std::setw(18) << "Energy (Eh)";
    out << std::setw(18) << "\u0394E (Eh)";
    out << std::setw(18) << " |\u0394P(S)|";
    if(not iCS or nC > 1)
      out << std::setw(18) << "  |\u0394P(M)|";
     
    out << std::endl;
    out << std::setw(16) << "-------------";
    out << std::setw(18) << "-----------";
    out << std::setw(18) << "-------";
    out << std::setw(18) << "-------";
    if(not iCS or nC > 1)
      out << std::setw(18) << "-------";
    out << std::endl;

  }; // SingleSlater<T>::printSCFHeader


  /**
   *  \brief Print the current convergence information of the SCF
   *  procedure
   */ 
  void SingleSlaterBase::printSCFProg(std::ostream &out) {

    // SCF Iteration
    out << "  SCFIt: " <<std::setw(6) << std::left << scfConv.nSCFIter + 1;

    // Current Total Energy
    out << std::setw(18) << std::fixed << std::setprecision(10)
                         << std::left << this->totalEnergy;

    out << std::scientific << std::setprecision(7);
    // Current Change in Energy
    out << std::setw(14) << std::right << scfConv.deltaEnergy;
    out << "   ";
    out << std::setw(13) << std::right << scfConv.RMSDenScalar;
    if(not iCS or nC > 1) {
      out << "   ";
      out << std::setw(13) << std::right << scfConv.RMSDenMag;
    }
  
    out << std::endl;
  }; // SingleSlater<T>::printSCFProg

}; // namespace ChronusQ

#endif
