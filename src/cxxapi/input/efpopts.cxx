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
#include <physcon.hpp>
#include <chronusq_sys.hpp>
#include <cxxapi/options.hpp>
#include <cxxapi/procedural.hpp>
#include <cerr.hpp>
#include <memmanager.hpp>
#include "efp.h"

namespace ChronusQ {
  
  /**
   *
   *  Check valid keywords in the section.
   *
  */
  void CQEFP_VALID( std::ostream &out, CQInputFile &input ) {
  // Allowed keywords
    std::vector<std::string> allowedKeywords = {
      "TERMS",
      "COORD_TYPE",
      "FRAGLIB_PATH",
      "EFP_PARAM_FILE",
      "DISP_DAMP",
      "ELEC_DAMP",
      "POL_DAMP",
      "POL_DRIVER",
      "ENABLE_PBC",
      "ENABLE_CUTOFF",
      "SWF_CUTOFF",
      "FRAG_NUM",
      "FRAGMENTS"
    };
  
  // Specified keywords
    std::vector<std::string> efpKeywords = input.getDataInSection("EFP");

  // Make sure all of the efpKeywords in allowedKeywords
    for( auto &keyword : efpKeywords ){
      auto ipos = std::find(allowedKeywords.begin(),allowedKeywords.end(),keyword);
      if( ipos == allowedKeywords.end() )
        CErr("Keyword EFP." + keyword + " is not recognized",std::cout);//Error
    }
  }
  

   /**
   *  Construct a Molecule object using the input file.
   *
   *  \param [in] out   Output device for data output
   *  \param [in] input Input file datastructure
   *
   *  \returns Appropriate Molecule object for the input parameters
   */
  ::efp_opts CQEFPOptions(std::ostream &out, CQInputFile &input){
    ::efp_opts efp_options;

    // Obtain terms
    try { auto t1 = input.getData<std::string>("EFP.TERMS");
          efp_options.terms = (unsigned int) (std::atoi(t1.c_str()));
        }
    catch (...) {
      CErr("Unable to set EFP terms!", out);
    }

    // Obtain disp_damp
    try { auto t1 = input.getData<std::string>("EFP.DISP_DAMP");
          if(t1 == "OVERLAP" ){
            efp_options.disp_damp = EFP_DISP_DAMP_OVERLAP;
          }
          if(t1 == "TT" )
            efp_options.disp_damp = EFP_DISP_DAMP_TT;
          if(t1 == "OFF" )
            efp_options.disp_damp = EFP_DISP_DAMP_OFF;
        }
    catch (...) {
      //CErr("Default dispersion damping!", out);
    }


    try { auto t1 = input.getData<std::string>("EFP.ELEC_DAMP"); 
          if(t1 == "SCREEN" )
            efp_options.elec_damp = EFP_ELEC_DAMP_SCREEN;
          if(t1 == "OVERLAP" )
            efp_options.elec_damp = EFP_ELEC_DAMP_OVERLAP;
          if(t1 == "OFF" )
            efp_options.elec_damp = EFP_ELEC_DAMP_OFF;
        }
    catch (...) {
      //CErr("Default electrostatic damping!", out);
    }

    
    try { auto t1 = input.getData<std::string>("EFP.POL_DAMP");
          if(t1 == "TT" )
            efp_options.pol_damp = EFP_POL_DAMP_TT;
          if(t1 == "OFF" )
            efp_options.pol_damp = EFP_POL_DAMP_OFF;
        }
    catch (...) {
      //CErr("Default polarization damping!", out);
    }

    
    try { auto t1 = input.getData<std::string>("EFP.POL_DRIVER");
          if(t1 == "ITERATIVE" )
            efp_options.pol_driver = EFP_POL_DRIVER_ITERATIVE;
          if(t1 == "DIRECT" )
            efp_options.pol_driver = EFP_POL_DRIVER_DIRECT;
        }
    catch (...) {
      //CErr("Unable to find polarization driver!", out);
    }

    
    try { efp_options.enable_pbc = input.getData<int>("EFP.ENABLE_PBC"); }
    catch (...) {
      //CErr("Unable to create PBC!", out);
    }

      
    try { efp_options.enable_cutoff = input.getData<int>("EFP.ENABLE_CUTOFF"); }
    catch (...) {
      //CErr("Unable to create cutoff!", out);
    }
    
    
    try { efp_options.swf_cutoff = input.getData<double>("EFP.SWF_CUTOFF"); }
    catch (...) {
      //CErr("Unable to find cutoff!", out);
    }
    return efp_options;
  };
  
  
  std::string CQEFPFilePath(std::ostream &out, CQInputFile &input){
    std::string efp_file_path;

    
    try { efp_file_path = input.getData<std::string>("EFP.EFP_PARAM_FILE"); }
    catch (...) {
      //CErr("Unable to ascertain the typr of coordination!", out);
    }
    return efp_file_path;
  };

  void CQFragCoord(std::ostream &out, CQInputFile &input, CQMemManager& mem_, struct Fragment_ifm* frag_ifm){
  
    size_t t = 0;

    try { t = input.getData<int>("EFP.FRAG_NUM"); }
    catch (...) {
      CErr("Unable to confirm the number of fragments", out);
    }
    // the number of the fragments
     
    int m = 0;
    try { auto t1 = input.getData<std::string>("EFP.COORD_TYPE");
          if(t1 == "XYZABC" ){
            m = 6;  
            (frag_ifm->Coord_type) = EFP_COORD_TYPE_XYZABC;
          }
          else if(t1 == "POINTS" ){
            (frag_ifm->Coord_type) = EFP_COORD_TYPE_POINTS;
            m = 9;
          }
          else if(t1 == "ROTMAT" ){
            (frag_ifm->Coord_type) = EFP_COORD_TYPE_ROTMAT;
            m = 12;
          }
        }
    catch (...) {
      CErr("Unable to ascertain the type of coordination!", out);
    }   

    std::string geomStr;
    try { geomStr = input.getData<std::string>("EFP.FRAGMENTS"); }
    catch (...) {
      CErr("Unable to find Fragment Geometry!", out);
    }
    std::istringstream geomStream; geomStream.str(geomStr);
    std::vector<std::string> tokens;
    std::locale loc;
    int k = 0;
    auto frag_matrix = mem_.template malloc<double>(m*t); 
    // Loop over lines of geometry specification
    for(std::string line; std::getline(geomStream, line); ){
      split(tokens,line," \t");
      
      if( tokens.size() == 0 ) continue;

      std::string fragSymb = tokens[0];
      bool hasDig = std::any_of(fragSymb.begin(),fragSymb.end(),
        [&](char a) {return std::isdigit(a,loc); });
      bool hasAlpha = std::any_of(fragSymb.begin(),fragSymb.end(),
        [&](char a) {return std::isalpha(a,loc); });

      
      bool isAtNum   = hasDig   and not hasAlpha;
      bool isComIso  = hasAlpha and not hasDig  ;
      bool isSpecIso = hasDig   and     hasAlpha;
      if( isComIso ) {
        (frag_ifm->fragments).push_back(tokens[1]);
      }
      if( isAtNum ) {
        if( frag_ifm->Coord_type == EFP_COORD_TYPE_XYZABC ){
          for(int j = 0; j < 6; j++){
            auto k1 = std::stod(tokens[j]);
            *(frag_matrix+j+6*k) = k1;
          }
          k += 1;             
        }
        if( frag_ifm->Coord_type == EFP_COORD_TYPE_POINTS ){
          
          for(int j = 0; j < 3; j++){
            auto k1 = std::stod(tokens[j]);
            *(frag_matrix+j+3*k) = k1;
          }
          k += 1;
        }
        if( frag_ifm->Coord_type == EFP_COORD_TYPE_ROTMAT ){
          
          for(int j = 0; j < 3; j++){
            auto k1 = std::stod(tokens[j]);
            *(frag_matrix+j+3*k) = k1;
          }
          k += 1;
        }
      }
      frag_ifm->Frag_coord = const_cast<double*>(frag_matrix);
    }
  };
  
  std::shared_ptr<EFPBase> CQEFPControl(std::ostream &out,
    CQInputFile &input, std::shared_ptr<AOIntegralsBase> &aoints,
    std::shared_ptr<SingleSlaterBase> &ss, bool EFP_bool){
    
    std::shared_ptr<EFPBase> efp_base;
    if(EFP_bool == true){
      auto e1 = CQEFPOptions(out,input);
      auto e2 = CQEFPFilePath(out,input);
      Fragment_ifm e3;
    
      auto slater0 = dynamic_cast<SingleSlater<double,double>* >(ss.get());
      auto slater1 = dynamic_cast<SingleSlater<double,dcomplex>* >(ss.get());
      auto slater2 = dynamic_cast<SingleSlater<dcomplex,double>* >(ss.get());
      auto slater3 = dynamic_cast<SingleSlater<dcomplex,dcomplex>* >(ss.get());
      if(slater0 != NULL){
        auto aoint_cast = std::dynamic_pointer_cast<AOIntegrals<double> >(aoints);
        CQFragCoord(out,input,aoint_cast->memManager_,&e3);     
        auto efp_trans  = std::make_shared<EFP<double,double> >(*aoint_cast,ss.get(),aoint_cast->memManager_);
        (efp_trans.get())->Initialize(&e1,e2,&e3);
        efp_base = std::dynamic_pointer_cast<EFPBase>(efp_trans);
        efp_base->type_choice = 0;
      }
      if(slater1 != NULL){
        auto aoint_cast = std::dynamic_pointer_cast<AOIntegrals<dcomplex> >(aoints);
        CQFragCoord(out,input,aoint_cast->memManager_,&e3);     
        auto efp_trans  = std::make_shared<EFP<dcomplex,double> >(*aoint_cast,ss.get(),aoint_cast->memManager_);
        (efp_trans.get())->Initialize(&e1,e2,&e3);
        efp_base = std::dynamic_pointer_cast<EFPBase>(efp_trans);
        efp_base->type_choice = 1;
      }
      if(slater2 != NULL){
        auto aoint_cast = std::dynamic_pointer_cast<AOIntegrals<double> >(aoints);
        CQFragCoord(out,input,aoint_cast->memManager_,&e3);     
        auto efp_trans  = std::make_shared<EFP<double,dcomplex> >(*aoint_cast,ss.get(),aoint_cast->memManager_);
        (efp_trans.get())->Initialize(&e1,e2,&e3);
        efp_base = std::dynamic_pointer_cast<EFPBase>(efp_trans);
        efp_base->type_choice = 2;
      }
      if(slater3 != NULL){
        auto aoint_cast = std::dynamic_pointer_cast<AOIntegrals<dcomplex> >(aoints);
        CQFragCoord(out,input,aoint_cast->memManager_,&e3);     
        auto efp_trans  = std::make_shared<EFP<dcomplex,dcomplex> >(*aoint_cast,ss.get(),aoint_cast->memManager_);
        (efp_trans.get())->Initialize(&e1,e2,&e3);
        efp_base = std::dynamic_pointer_cast<EFPBase>(efp_trans);
        efp_base->type_choice = 3;
      }
    }
    else
      efp_base = NULL;
    return efp_base; 


  };
};//using ChronusQ namespace

