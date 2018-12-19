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
#ifndef __INCLUDED_UTIL_FILES_HPP__
#define __INCLUDED_UTIL_FILES_HPP__

#include <chronusq_sys.hpp>
#include <cxxapi/input.hpp>
#include <cerr.hpp>

#include "H5Cpp.h"

template <typename T>
inline H5::CompType H5PredType() {

  if( std::is_same<double,T>::value )
    return H5::CompType(
      H5::DataType(H5::PredType::NATIVE_DOUBLE).getId()
    );
  else if( std::is_same<uint64_t,T>::value ) 
    return H5::CompType(
      H5::DataType(H5::PredType::NATIVE_UINT64).getId()
    );
  else if( std::is_same<int32_t,T>::value ) 
    return H5::CompType(
      H5::DataType(H5::PredType::NATIVE_INT32).getId()
    );
  else if( std::is_same<dcomplex,T>::value ) {
    typedef struct {
      double re;
      double im;
    } complex_t;

    H5::CompType complexType(sizeof(complex_t));

    complexType.insertMember(
      "r",HOFFSET(complex_t,re),H5::PredType::NATIVE_DOUBLE
    );
    complexType.insertMember(
      "i",HOFFSET(complex_t,im),H5::PredType::NATIVE_DOUBLE
    );
    return complexType;
  }
  else ChronusQ::CErr();

};

#define OpenH5File(file,type) \
  assert(not fName_.empty()); \
  H5::H5File file(fName_,type); \
  exists_ = true;

#define CreateH5File(file) \
  OpenH5File(file,H5F_ACC_TRUNC)

#define CreateDataSet(file,obj,data,type,nDims,dims) \
  OpenH5File(file,H5F_ACC_RDWR); \
  H5::DataSpace space(nDims,dims); \
  H5::DataSet obj = file.createDataSet(data,type,space);

#define OpenDataSet(file,obj,data) \
  OpenH5File(file,H5F_ACC_RDWR); \
  H5::DataSet obj = file.openDataSet(data);


namespace ChronusQ {

  class SafeFile {
  
    std::string fName_;
    bool        exists_;
  
    public:

      // Defaulted ctors
      SafeFile(const SafeFile &) = default;

      // String ctor
      SafeFile( const std::string &fName = "", 
        bool exists = false ) :
        fName_(fName), exists_(exists) { }
  

      // Member functions

      inline bool exists() const { return exists_; }
      inline std::string fName() const{ return fName_; }
      inline void setFile(const std::string &name) { fName_ = name; }

      inline void createFile() {
        CreateH5File(file);
      }; 

      inline void createGroup(const std::string &group) {
        OpenH5File(file,H5F_ACC_RDWR)

        // Create nested groups if need be
        std::vector<std::string> tokens;
        split(tokens,group,"/");

        // Remove empty tokens
        std::remove_if(tokens.begin(),tokens.end(),
          [](std::string &x) -> bool{ return x.empty(); }
        );

        std::string tmpName;
        for(auto &str : tokens) {
          H5::Group tmpGroup = 
            file.openGroup(tmpName.empty() ? "/" : tmpName);
          tmpName += "/" + str;

          try { tmpGroup.openGroup(str); } 
          catch(...) { tmpGroup.createGroup(str); }
        }
      }

      template <typename T>
      inline void createDataSet(const std::string &data,
        const std::vector<hsize_t> &dims) {

        CreateDataSet(file,obj,data,H5PredType<T>(),
          dims.size(),&dims[0]);

      };


      template <typename T>
      void readData(const std::string &dataSet, T* data) {
        OpenDataSet(file,obj,dataSet);
        obj.read(data,H5PredType<T>());
      };

      size_t sizeOfData(const std::string &dataSet) {

        size_t datSize;
        OpenDataSet(file,obj,dataSet);

        H5::DataSpace space = obj.getSpace();
        datSize = space.getSelectNpoints();
  
        return datSize;
        
      }; //sizeOfData

      template <typename T>
      void writeData(const std::string &dataSet, T* data) {
        OpenDataSet(file,obj,dataSet);
        obj.write(data,H5PredType<T>());
      };


      template <typename T>
      void safeWriteData(const std::string &dataSet, T* data,
        const std::vector<hsize_t> &dims) {

        try {
          writeData(dataSet,data);

        // DataSet doesnt exist
        } catch(...) {

          try { this->template createDataSet<T>(dataSet,dims); }
          catch(...) {

            // Separate Group from DataSet
            auto sPos = dataSet.rfind("/");
            if( sPos == std::string::npos ) throw;
            else {
              createGroup(dataSet.substr(0,sPos));
              this->template createDataSet<T>(dataSet,dims);
            }

          }

          writeData(dataSet,data);

        }

      };

      std::vector<hsize_t> getDims(const std::string &dataSet) {

        std::vector<hsize_t> dims;
        try {
          OpenDataSet(file,obj,dataSet);

          H5::DataSpace space = obj.getSpace();
          int rank = space.getSimpleExtentNdims();
   
          dims.resize(rank);

          space.getSimpleExtentDims(&dims[0],NULL);
        } catch(...){ 
          // Return empty vector 
        }
  
        return dims;

      }


      H5T_class_t getTypeClass(const std::string &dataSet) {

        H5T_class_t type;

        try {

          OpenDataSet(file,obj,dataSet);

          type = obj.getTypeClass();

        } catch(...){ 
          // Return empty vector 
        }


        return type;

      }



  }; // class SafeFile

}; // namespace ChronusQ

#endif
