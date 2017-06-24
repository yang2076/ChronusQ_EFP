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

#include <basisset/basisset_def.hpp>
#include <basisset/reference.hpp>
#include <cxxapi/output.hpp>

namespace ChronusQ {

  // Basis set keyword map
  // TODO: Keywords for the Sapporo sets
  std::unordered_map<std::string,std::string> basisKeyword = {
    {  "STO-3G"         , "sto3g.gbs"                     },
    {  "STO-6G"         , "sto6g.gbs"                     },
    {  "3-21G"          , "3-21g.gbs"                     },
    {  "4-31G"          , "4-31g.gbs"                     },
    {  "6-31G"          , "6-31g.gbs"                     },
    {  "6-31G(D)"       , "6-31g*.gbs"                    },
    {  "6-31++G(D)"     , "6-31++g*.gbs"                  },
    {  "6-311G"         , "6-311g.gbs"                    },
    {  "6-311+G(D)"     , "6-311+g*.gbs"                  },
    {  "6-311+G(D,P)"   , "6-311+g**.gbs"                 },
    {  "6-311+G(2D,P)"  , "6-311+g_2d_p.gbs"              },
    {  "CC-PVDZ"        , "cc-pvdz.gbs"                   },
    {  "CC-PVTZ"        , "cc-pvtz.gbs"                   },
    {  "CC-PVQZ"        , "cc-pvqz.gbs"                   },
    {  "CC-PV5Z"        , "cc-pv5z.gbs"                   },
    {  "CC-PV6Z"        , "cc-pv6z.gbs"                   },
    {  "CC-PVDZ-RI"     , "cc-pvdz_ri.gbs"                },
    {  "DEF2-SVP"       , "def2-svp.gbs"                  },
    {  "DEF2-SVPD"      , "def2-svpd.gbs"                 },
    {  "DEF2-TZVP"      , "def2-tzvp.gbs"                 },
    {  ""  , "sapporo-dkh3-dzp-2012_all.gbs" },
    {  ""  , "sapporo-dkh3-dzp-2012_no.gbs"  },
    {  ""  , "sapporo-dkh3-dzp-2012_sp.gbs"  },
    {  ""  , "sapporo-dkh3-dzp_all.gbs"      },
    {  ""  , "sapporo-dkh3-dzp_no.gbs"       },
    {  ""  , "sapporo-dkh3-dzp_sp.gbs"       },
    {  ""  , "sapporo-dkh3-qzp-2012_all.gbs" },
    {  ""  , "sapporo-dkh3-qzp-2012_no.gbs"  },
    {  ""  , "sapporo-dkh3-qzp-2012_sp.gbs"  },
    {  ""  , "sapporo-dkh3-qzp_all.gbs"      },
    {  ""  , "sapporo-dkh3-qzp_no.gbs"       },
    {  ""  , "sapporo-dkh3-qzp_sp.gbs"       },
    {  ""  , "sapporo-dkh3-tzp-2012_all.gbs" },
    {  ""  , "sapporo-dkh3-tzp-2012_no.gbs"  },
    {  ""  , "sapporo-dkh3-tzp-2012_sp.gbs"  },
    {  ""  , "sapporo-dkh3-tzp_all.gbs"      },
    {  ""  , "sapporo-dkh3-tzp_no.gbs"       },
    {  ""  , "sapporo-dkh3-tzp_sp.gbs"       },
    {  ""  , "sapporo-dzp-2012_all.gbs"      },
    {  ""  , "sapporo-dzp-2012_no.gbs"       },
    {  ""  , "sapporo-dzp-2012_sp.gbs"       },
    {  ""  , "sapporo-dzp_all.gbs"           },
    {  ""  , "sapporo-dzp_no.gbs"            },
    {  ""  , "sapporo-dzp_sp.gbs"            },
    {  ""  , "sapporo-qzp-2012_all.gbs"      },
    {  ""  , "sapporo-qzp-2012_no.gbs"       },
    {  ""  , "sapporo-qzp-2012_sp.gbs"       },
    {  ""  , "sapporo-qzp_all.gbs"           },
    {  ""  , "sapporo-qzp_no.gbs"            },
    {  ""  , "sapporo-qzp_sp.gbs"            },
    {  ""  , "sapporo-tzp-2012_all.gbs"      },
    {  ""  , "sapporo-tzp-2012_no.gbs"       },
    {  ""  , "sapporo-tzp-2012_sp.gbs"       },
    {  ""  , "sapporo-tzp_all.gbs"           },
    {  ""  , "sapporo-tzp_no.gbs"            },
    {  ""  , "sapporo-tzp_sp.gbs"            }
  };


  /**
   *  Path / Molecule BasisSet constructor 
   *
   *  Constructs a BasisSet object and populates member data given
   *  a basis set name and a Molecule object
   *
   *  \param [in] basisName Either a basis file in G94 format
   *                        or a known keyword which maps to a basis
   *                        file
   *  \param [in] mol       Molecule for which to construct the BasisSet
   */
  BasisSet::BasisSet(std::string basisName, const Molecule &mol,
    bool _forceCart) {

    forceCart = _forceCart;

    std::string uppercase(basisName);

    // Convert the passed string to uppercase
    std::transform(basisName.begin(),basisName.end(),uppercase.begin(),
      [](unsigned char c){ return std::toupper(c); });  

    // Possibly find appropriate basis file for keyword
    if( basisKeyword.find(uppercase) != basisKeyword.end() )
      basisName = basisKeyword[basisName];

    // Generate the reference basis set of that keyword
    ReferenceBasisSet ref(basisName,_forceCart);

    // Update appropriate shell set and coefficients for the Molecule
    // object
    std::tie(shells,unNormCont) = std::move(ref.generateShellSet(mol));

    // Obtain a copy of the basis centers
    std::for_each(mol.atoms.begin(),mol.atoms.end(),
      [&]( const Atom &at ){ centers.emplace_back(at.coord); }
    );

    // Update the BasisSt member data
    update();
    
  }; // BasisSet::BasisSet 


  /**
   *  Updates the member data of a BasisSet object.
   *
   *  Computes nShell, nBasis and nPrimitive. Determines  maxL and maxPrim.
   *  Constructs the shell mappings to basis function number and center number
   */ 
  void BasisSet::update() {

    // Compute the number of shells
    nShell = shells.size();

    // Compute the number of Basis functions and primitives
    std::tie(nBasis,nPrimitive) = std::accumulate(shells.begin(),shells.end(),
      std::pair<size_t,size_t>{0,0},
      [](std::pair<size_t,size_t> init, libint2::Shell &sh) ->
         std::pair<size_t,size_t> {
        return { init.first + sh.size(), 
                 init.second + sh.size()*sh.alpha.size() };
      }
    );

    // Determine the maximum angular momentum
    maxL = std::max_element(shells.begin(), shells.end(),
      [](libint2::Shell &sh1, libint2::Shell &sh2){
        return sh1.contr[0].l < sh2.contr[0].l;
      }
    )->contr[0].l;

    // Determine the maximum contraction depth
    maxPrim = std::max_element(shells.begin(), shells.end(),
      [](libint2::Shell &sh1, libint2::Shell &sh2){
        return sh1.alpha.size() < sh2.alpha.size();
      }
    )->alpha.size();

    // Create basis maps
    // Maps Sh # -> BF #
    // Maps Sh # -> Center index
    size_t run(0);
    std::for_each(shells.begin(),shells.end(),
      [&](libint2::Shell &sh) {
        auto cenIt = std::find_if(centers.begin(),centers.end(),
          [&](cart_t &cen){ return sh.O == cen; }
        );
        mapSh2Cen.emplace_back(std::distance(centers.begin(),cenIt));

        mapSh2Bf.emplace_back(run);
        run += sh.size(); 
        
      }
    );
    
  }; // BasisSet::update



  /**
   *  Outputs relevant information for the BasisSet object
   *  to a specified output.
   *
   *  \param [in/out] out Ouput device
   *  \paral [in]     mol BasisSet object to output.
   */ 
  std::ostream& operator<<(std::ostream &out, const BasisSet &basis) {

    out << std::endl;
    out << "Basis Set Information:" << std::endl << BannerTop
        << std::endl << std::endl;
    
    out << std::left;

    out << "  " << std::setw(25) << "NBasis" << basis.nBasis 
        << std::endl;
    out << "  " << std::setw(25) << "NPrimitive" << basis.nPrimitive 
        << std::endl;
    out << "  " << std::setw(25) << "NShell" << basis.nShell 
        << std::endl;
    out << "  " << std::setw(25) << "Max Primitive" << basis.maxPrim 
        << std::endl;
    out << "  " << std::setw(25) << "Max L" << basis.maxL << std::endl;

    out << std::endl << std::endl;

    out << "  " << "Shell Information:" << std::endl << bannerMid << std::endl;

    out << "  " << "  " << std::left;
    out << std::setw(5) << "#" ;
    out << std::setw(5) << "L" ;
    out << std::setw(15) << std::right << "Exponents";
    out << std::setw(30) << std::right << "Normalized Contraction";
  //out << std::setw(30) << std::right << "Unnormalized Contraction";
    out << std::endl << std::endl;
    
    for(auto iShell = 0; iShell < basis.nShell; iShell++){
      out << "  " << "  " << std::left << std::setprecision(4) 
          << std::scientific;
      out << std::setw(5) << iShell ;
      out << std::setw(5) << basis.shells[iShell].contr[0].l << std::right;
      out << std::setw(15) << basis.shells[iShell].alpha[0];
      out << std::setw(30) << basis.shells[iShell].contr[0].coeff[0];
    //out << std::setw(30) << basis.unNormCont[iShell][0];

      out << std::endl;
      for(auto iPrim = 1; iPrim < basis.shells[iShell].alpha.size(); iPrim++){
        out << "  " << "  " << std::setw(5+5) << " ";
        out << std::setw(15) << basis.shells[iShell].alpha[iPrim];
        out << std::setw(30) << basis.shells[iShell].contr[0].coeff[iPrim];
      //out << std::setw(30) << basis.unNormCont[iShell][iPrim];

        out << std::endl;
      }
      out << std::endl;
    }
    out << std::endl << BannerEnd << std::endl;

    return out; // return std::ostream reference

  }; // BasisSet::operator<<

}; // namespace ChronusQ

