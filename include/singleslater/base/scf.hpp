/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2018 Li Research Group (University of Washington)
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
#include <chronusqefp.hpp>
#define _DEBUG_SCF
namespace ChronusQ {

  /**
   *  \brief Performs the self--consistant field procedure given a set of 
   *  orbitals.
   *
   *  \warning SCF procedure assumes that the 1PDM and orbital (mo1/2) storage
   *  has been populated in some way.
   */ 
  void SingleSlaterBase::SCF(EMPerturbation &pert, EFPBase* EFP_1, bool EFP_bool) {
    
    SCFInit();
    // Initialize type independent parameters
    bool isConverged = false;
    scfControls.dampParam = scfControls.dampStartParam;

//XSLIC
//    scfControls.doIncFock = scfControls.doIncFock and (this->aoints.cAlg == DIRECT);
    scfControls.doIncFock = false;
                 
    if( scfControls.scfAlg == _NEWTON_RAPHSON_SCF )
      scfControls.doExtrap = false;
    
    
    // Compute initial properties
    double efp_wf_denp = 0.0 ;

    this->EFPEnergy = 0.0; 
    this->computeProperties(pert,efp_wf_denp);

    if( printLevel > 0 and MPIRank(comm) == 0 ) {
      printSCFHeader(std::cout,pert);
      printSCFProg(std::cout,false);
    }

    // EFP permanent multipole energy
    EFP_multipole_contri(EFP_1,EFP_bool);

    for( scfConv.nSCFIter = 0; scfConv.nSCFIter < scfControls.maxSCFIter; 
         scfConv.nSCFIter++) {

      // Save current state of the wave function (method specific)
      saveCurrentState();

      // Exit loop on convergence
      if(isConverged) break;
 

      // Get new orbtials and densities from current state: 
      //   C/D(k) -> C/D(k + 1)
      getNewOrbitals(pert, EFP_1, EFP_bool);
      
      // EFP induction energy
      double efp_wf_denp = 0.0;
      EFP_wf_dependent(EFP_1, EFP_bool, &efp_wf_denp);
      
      // Evaluate convergence
      isConverged = evalConver(pert,&efp_wf_denp);

      // Print out iteration information
      if( printLevel > 0 and (MPIRank(comm) == 0)) printSCFProg(std::cout);

    }; // Iteration loop

    // Save current state of the wave function (method specific)
    saveCurrentState();

    SCFFin();

    efp_wf_denp = 0.0;
    EFP_wf_dependent(EFP_1, EFP_bool, &efp_wf_denp);
    
    // Compute initial properties
    this->computeProperties(pert,efp_wf_denp);

#ifdef _DEBUG_SCF
    isConverged = true;
#endif
    //printSCFFooter(isConverged);
    if(not isConverged)
      CErr(std::string("SCF Failed to converged within ") + 
        std::to_string(scfControls.maxSCFIter) + 
        std::string(" iterations"));
    else if( printLevel > 0 ) {
      std::cout << std::endl << "SCF Completed : E("
                << refShortName_ << ") = " << std::fixed
                << std::setprecision(10) << this->totalEnergy
                << " Eh after " << scfConv.nSCFIter
                << " SCF Iterations" << std::endl;
      std::cout << "Energy details: " << std::endl
                << "One electron energy: " 
                << std::setprecision(10) << this->OBEnergy << std::endl
                << "Two electron energy: " 
                << std::setprecision(10) << this->MBEnergy << std::endl
                << "Nuclear repulsion energy: " 
                << std::setprecision(10) << this->NREnergy << std::endl
                << "EFP impact: " 
                << std::setprecision(10) << this->EFPEnergy << std::endl;
    }

    if( printLevel > 0 ) std::cout << BannerEnd << std::endl;
    if( printLevel > 0 ) {
      this->printEFPContribution(std::cout,EFP_1,EFP_bool);
      this->printMOInfo(std::cout);
      this->printMultipoles(std::cout);
      this->printSpin(std::cout);
      this->printMiscProperties(std::cout);
    }
    
  }; // SingleSlaterBase::SCF()


  
  void SingleSlaterBase::printSCFHeader(std::ostream &out, 
    EMPerturbation &pert) {

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


    out << std::setw(38)   << std::left << "  SCF Algorithm:";
    if( scfControls.scfAlg == _CONVENTIONAL_SCF )
      out << "Conventional SCF";
    else if( scfControls.scfAlg == _NEWTON_RAPHSON_SCF )
      out << "Newton-Raphson (2nd Order)";
    out << "\n";


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
 
    }


    if( scfControls.doIncFock ) {
      out << "\n  * Will Perform Incremental Fock Build -- Restarting Every "
          << scfControls.nIncFock << " SCF Steps\n";
         
    }

    // Field print
    if( pert.fields.size() != 0 ) {

      out << "\n\n  * SCF will be performed in the presence of an EM "
          << "perturbation:\n\n";

      for(auto &field : pert.fields) {

        auto amp = field->getAmp();

        out << "     * ";
        if( field->emFieldTyp == Electric ) out << "Electric";
        else                                out << "Magnetic";
        
        out << " ";
      
        if( field->size == 3 )        out << "Dipole";
        else if ( field->size == 6 )  out << "Quadrupole";
        else if ( field->size == 10 ) out << "Octupole";
        
        out << " Field: ";
        out << "{ ";
        for(auto i = 0; i < amp.size(); i++) {
          out << amp[i]; if(i != amp.size() - 1) out << ", ";
        }
        out << " }\n";

      }


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
  void SingleSlaterBase::printSCFProg(std::ostream &out,
    bool printDiff) {

    // SCF Iteration
    out << "  SCFIt: " <<std::setw(6) << std::left;

    if( printDiff ) out << scfConv.nSCFIter + 1;
    else            out << 0;

    // Current Total Energy
    out << std::setw(18) << std::fixed << std::setprecision(10)
                         << std::left << this->totalEnergy;

    if( printDiff ) {
      out << std::scientific << std::setprecision(7);
      // Current Change in Energy
      out << std::setw(14) << std::right << scfConv.deltaEnergy;
      out << "   ";
      out << std::setw(13) << std::right << scfConv.RMSDenScalar;
      if(not iCS or nC > 1) {
        out << "   ";
        out << std::setw(13) << std::right << scfConv.RMSDenMag;
      }
    }
  
    out << std::endl;
  }; // SingleSlater<T>::printSCFProg

}; // namespace ChronusQ

#endif
