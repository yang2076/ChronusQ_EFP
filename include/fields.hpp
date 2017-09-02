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
#ifndef __INCLUDED_FIELDS_HPP__
#define __INCLUDED_FIELDS_HPP__

#include <chronusq_sys.hpp>
#include <util/typedefs.hpp>

namespace ChronusQ {


  enum FieldGauge {
    Length
  };

  enum EMFieldTyp {
    Electric,
    Magnetic
  };


  struct EMFieldBase {

    EMFieldTyp emFieldTyp;
    FieldGauge fieldGauge;

    EMFieldBase(EMFieldTyp em = Electric, FieldGauge fg = Length):
      fieldGauge(fg), emFieldTyp(em){ };

    virtual std::valarray<double> getAmp() = 0;

  };

  /**
   *  \brief Tensorial specification of an EM field perturbation.
   *
   *  Handles the specification of various descriptions of the EM field
   *  perturbation, such as the multipole expansion.
   *
   *  Expectes field amplidudes of the following:
   *  1. Dipole operator     - [x, y, z]
   *  2. Quadrupole operator - [xx, xy, xz, yy, yz, zz]
   *  ect...
   */ 
  template <typename _VecTyp>
  struct EMField : public EMFieldBase {

    _VecTyp ampVec;

    EMField()                 = delete;
    EMField( const EMField& ) = default;
    EMField( EMField&& )      = default;

    EMField(EMFieldTyp em, FieldGauge fg,
      const _VecTyp &x): ampVec(x), EMFieldBase{em,fg}{ };

    EMField(EMFieldTyp em, const _VecTyp &x): 
      ampVec(x), EMFieldBase{em,Length}{ };

    EMField(const _VecTyp &x): EMField(Electric,Length,x){ };

    std::valarray<double> getAmp() {
      return std::valarray<double>(ampVec.begin(),ampVec.size());
    }

  };

  using DipoleField     = EMField<std::array<double,3>>;
  using QuadrupoleField = EMField<std::array<double,6>>;
  using OctupoleField   = EMField<std::array<double,10>>;

  /**
   *  \brief Cast a templated EMField shared_ptr to
   *  Base class.
   *
   *  dynamic_pointer_cast EMField<> -> EMFieldBase
   *
   *  \param [in] x Shared pointer for a EMField object
   *  \returns      Shared pointer for a EMFieldBase object
   */ 
  template < typename _UnitVector >
  std::shared_ptr<EMFieldBase> 
    cast(std::shared_ptr<EMField<_UnitVector>> x) {
    return std::dynamic_pointer_cast<
             EMFieldBase,EMField<_UnitVector>>(x);
  };
  
  /**
   *  \brief Full specification of a static EM perturbation
   *
   *  General to multiple fields.
   */ 
  struct EMPerturbation {

    std::vector<
      std::shared_ptr<EMFieldBase>
        > fields; ///< Fields for the EM perturbation


    /**
     *  \brief Add a field to the perturbation.
     *
     *  \param [in] f    Tensorial field specification
     */ 
    inline void addField(std::shared_ptr<EMFieldBase> f) {
      fields.push_back(f);
    };

    template <typename _Field>
    inline void addField(std::shared_ptr<_Field> f){ addField(cast(f)); };

    template <typename _Field>
    inline void addField(_Field f){ addField(std::make_shared<_Field>(f)); }

    template <typename _UV>
    inline void addField(EMFieldTyp em, FieldGauge fg, const _UV &uv) {
      addField(EMField<_UV>(em,fg,uv));    
    }

    template <typename _UV>
    inline void addField(EMFieldTyp em, const _UV &uv) {
      addField(EMField<_UV>(em,uv));    
    }



    /**
     *  \brief Obtain the tensorial field amplitude 
     *  as the sum of all the field contributions.
     *
     *  FIXME: This function only works if all fields in perturbation
     *  are of the same type (i.e. electric dipole, etc).
     *
     *  \param [in]
     */ 
    std::valarray<double> getAmp(){

      assert(fields.size() != 0);

      std::valarray<double> amp = fields[0]->getAmp();
      for(auto i = 1; i < fields.size(); i++)
        amp += fields[i]->getAmp();

      return amp;

    }
    
  }; // struct TDEMPerturbation

};

#endif
