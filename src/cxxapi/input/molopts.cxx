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
#include <physcon.hpp>
#include <cxxapi/options.hpp>
#include <cerr.hpp>

namespace ChronusQ {

  /**
   *  Construct a Molecule object using the input file.
   *
   *  \param [in] out   Output device for data output
   *  \param [in] input Input file datastructure
   *
   *  \returns Appropriate Molecule object for the input parameters
   */
  Molecule CQMoleculeOptions(std::ostream &out, CQInputFile &input) {

    Molecule mol;
    
    // Obtain Charge
    try { mol.charge = input.getData<int>("MOLECULE.CHARGE"); }
    catch (...) {
      CErr("Unable to set Molecular Charge!", out);
    }

    // Obtain Multiplicity
    try { mol.multip = input.getData<size_t>("MOLECULE.MULT"); }
    catch (...) {
      //CErr("Unable to set Molecular Spin Multiplicity!", out);
    }

    // Parse Geometry
    std::string geomStr;
    try { geomStr = input.getData<std::string>("MOLECULE.GEOM"); }
    catch (...) {
      CErr("Unable to find Molecular Geometry!", out);
    }



    std::istringstream geomStream; geomStream.str(geomStr);
    std::vector<std::string> tokens;
    std::vector<Atom> atoms;
    std::locale loc;

    // Loop over lines of geometry specification
    for(std::string line; std::getline(geomStream, line); ){
      split(tokens,line," \t");

      if( tokens.size() == 0 ) continue;

      if( std::any_of(tokens[0].begin(),tokens[0].end(),
        [&](char a) {return std::isdigit(a,loc); }) ){



        auto it = 
        std::find_if(atomicReference.begin(),atomicReference.end(),
          [&](std::pair<std::string,Atom> st){ 
            return st.second.atomicNumber == std::stoi(tokens[0]);}
           );

        atoms.emplace_back((it == atomicReference.end() ? "X" : it->first));

        
      } else
        atoms.emplace_back(tokens[0]);
      
      // Convert to Bohr
      atoms.back().coord[0] = std::stod(tokens[1]) / AngPerBohr;
      atoms.back().coord[1] = std::stod(tokens[2]) / AngPerBohr;
      atoms.back().coord[2] = std::stod(tokens[3]) / AngPerBohr;
      
    }

    // Set the Atoms vector in Molecule (calls Molecule::update())
    mol.setAtoms(atoms);

    // Output Molecule data
    out << mol << std::endl;

    return mol; // Return Molecule object (no intermediates)

  }; // CQMoleculeOptions

}; // namespace ChronusQ

