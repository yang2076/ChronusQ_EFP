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
#include <aointegrals.hpp>


// Template for a collective operation on the members of an 
// AOIntegrals object
#define AOIntegrals_COLLECTIVE_OP(OP_MEMBER, OP_OP, OP_VEC_OP) \
    OP_MEMBER(threshSchwartz_); \
    OP_MEMBER(cAlg_); \
    OP_MEMBER(orthoType_); \
    \
    /* Copy over meta  */ \
    OP_OP(double,schwartz); \
    OP_OP(double,ortho1); \
    OP_OP(double,ortho2); \
    \
    /* 1-e Integrals */ \
    OP_OP(double,overlap); \
    OP_OP(double,kinetic); \
    OP_OP(double,potential); \
    OP_VEC_OP(double,lenElecDipole); \
    OP_VEC_OP(double,lenElecQuadrupole); \
    OP_VEC_OP(double,lenElecOctupole); \
    OP_VEC_OP(double,velElecDipole); \
    OP_VEC_OP(double,velElecQuadrupole); \
    OP_VEC_OP(double,velElecOctupole); \
    OP_VEC_OP(double,coreH); \
    \
    /* 2-e Integrals */ \
    OP_OP(double,ERI)



// Dummy function
#define DUMMY(X) 
#define DUMMY2(X,Y) 


// Deallocation functions
#define DEALLOC_OP(typ,PTR) if(PTR != nullptr) memManager_.free(PTR);

#define DEALLOC_VEC_OP(typ, VEC_PTR) \
  for(auto i = 0; i < VEC_PTR.size(); i++) \
    if(VEC_PTR[i] != nullptr) DEALLOC_OP(typ,VEC_PTR[i]); \
  VEC_PTR.clear();




// Copy functions
#define COPY_OTHER_MEMBER(X) X = other.X;

#define COPY_OTHER_OP(typ,PTR) \
  if(other.PTR != nullptr) { \
    size_t OPSZ = memManager_.getSize<typ>(other.PTR); \
    PTR = memManager_.malloc<typ>(OPSZ); \
    std::copy_n(other.PTR, OPSZ, PTR); \
  } 

#define COPY_OTHER_VEC_OP(typ, VEC_PTR) \
  for(auto i = 0; i < other.VEC_PTR.size(); i++) \
    if(other.VEC_PTR[i] != nullptr) { \
      size_t OPSZ = memManager_.getSize<typ>(other.VEC_PTR[i]); \
      VEC_PTR.emplace_back(memManager_.malloc<typ>(OPSZ)); \
      std::copy_n(other.VEC_PTR[i], OPSZ, VEC_PTR.back()); \
    } 



// Move functions
#define MOVE_OTHER_OP(typ,PTR) \
  if(other.PTR != nullptr) { \
    size_t OPSZ = memManager_.getSize<typ>(other.PTR); \
    PTR = memManager_.malloc<typ>(OPSZ); \
    std::copy_n(other.PTR, OPSZ, PTR); \
    DEALLOC_OP(typ,other.PTR); \
  } 

#define MOVE_OTHER_VEC_OP(typ, VEC_PTR) \
  for(auto i = 0; i < other.VEC_PTR.size(); i++) \
    if(other.VEC_PTR[i] != nullptr) { \
      size_t OPSZ = memManager_.getSize<typ>(other.VEC_PTR[i]); \
      VEC_PTR.emplace_back(memManager_.malloc<typ>(OPSZ)); \
      std::copy_n(other.VEC_PTR[i], OPSZ, VEC_PTR.back()); \
      DEALLOC_OP(typ,other.VEC_PTR[i]); \
    } \
  other.VEC_PTR.clear();

    

namespace ChronusQ {

  /**
   *  Constructs an AOIntegrals object to another by copy.
   *
   *  \param [in] other AOIntegrals object to copy
   */ 
  AOIntegrals::AOIntegrals(const AOIntegrals &other) :
    AOIntegrals(other.memManager_, other.molecule_, other.basisSet_){

    AOIntegrals_COLLECTIVE_OP(COPY_OTHER_MEMBER,COPY_OTHER_OP,
      COPY_OTHER_VEC_OP);

  }; // AOIntegrals::AOIntegrals(const AOIntegrals &other)


  /**
   *  Constructs an AOIntegrals object to another by move.
   *
   *  \warning Deallocates the passed AOIntegrals object
   *
   *  \param [in] other AOIntegrals object to move
   */ 
  AOIntegrals::AOIntegrals(AOIntegrals &&other) :
    AOIntegrals(other.memManager_, other.molecule_, other.basisSet_){

    AOIntegrals_COLLECTIVE_OP(COPY_OTHER_MEMBER,MOVE_OTHER_OP,
      MOVE_OTHER_VEC_OP);


  }; // AOIntegrals::AOIntegrals(AOIntegrals &&other)



  /**
   *  Deallocates the internal memory an AOIntegrals object
   */ 
  void AOIntegrals::dealloc() {

    AOIntegrals_COLLECTIVE_OP(DUMMY,DEALLOC_OP,DEALLOC_VEC_OP);

  }; // AOIntegrals::dealloc()




}; // namespace ChronusQ
