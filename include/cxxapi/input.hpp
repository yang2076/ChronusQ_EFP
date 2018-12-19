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
#ifndef _INCLUDED_INPUT_HPP_
#define _INCLUDED_INPUT_HPP_

#include <chronusq_sys.hpp>
#include "efp.h"

namespace ChronusQ {


  struct Fragment_ifm{
    std::vector<std::string> fragments;
    const double* Frag_coord;
    enum efp_coord_type Coord_type;
  };
  /**
   *  \brief A class to handle the parsing a data fetching from a 
   *  ChronusQ input file.
   */
  class CQInputFile {
  
    std::shared_ptr<std::ifstream> inFile_ = nullptr;  ///< Input file
  
    std::unordered_map<std::string,
      std::unordered_map<std::string,std::string>> dict_; 
    ///< Input data fields partitioned by section headings 
  
  
  
  
  
    // Parses the input file
    // (See src/cxxapi/input/parse/cxxapi.cxx for documentation)
    void parse();
  
    // Splits query string on "."
    // (See src/cxxapi/input/parse/cxxapi.cxx for documentation)
    static std::pair<std::string,std::string> 
      splitQuery(const std::string&);
  
    /**
     *  std::ofstream constructor.
     *
     *  Sets and parses input file from std::ofstream object
     *  \param [in] inFile  File object to parse
     */  
    CQInputFile(std::shared_ptr<std::ifstream> inFile) :
      inFile_(inFile){ parse(); }
  
  
  
  public:
  
    // Disable default, copy and move constructors and assignment operators
    CQInputFile()                               = delete;
    CQInputFile(const CQInputFile &)            = delete;
    CQInputFile(CQInputFile &&)                 = delete;
    CQInputFile& operator=(const CQInputFile &) = delete; 
    CQInputFile& operator=(CQInputFile &&)      = delete; 
  
    /**
     *  Filename constructor.
     *
     *  Sets and parses CQ input file given a file name
     *  \param [in] inFileName  Name of CQ input file
     */ 
    CQInputFile(std::string inFileName) :
      CQInputFile(std::make_shared<std::ifstream>(inFileName)){ }
  
  
  
  
  
    /**
     *  \brief Template function which returns the value of a data field
     *  from the input file in a specified datatype given a formatted 
     *  query string.
     *
     *  i.e.
     *
     *  Input entry:
     *    [SCF]
     *    DENTOL = 1E-6
     *    
     *  Query
     *    double tol = input.getData<double>("SCF.DENTOL");
     *
     *  This example returns the value of the string data field "SCF.DENTOL"
     *  as a double precision number. Various specializations of this function
     *  exist for various datatypes
     *
     *  \param [in] s Formatted query string to be parsed
     *  \return       Value of query data field as specified datatype
     */ 
    template <typename T> T getData(std::string) ; 
  
  
  
  
    /**
     *  Checks whether or not the parsed CQ input file contains
     *  a query section.
     *
     *  \paral  [in] str Query string of a section heading
     *  \return      True if input file contains that heading
     */ 
    inline bool containsSection(std::string str) const {
      return dict_.find(str) != dict_.end();
    }
  
    /**
     *  Checks whether or not the parsed CQ input file contains
     *  a query data field.
     *
     *  \paral  [in] str Query string of a data field (includes section heading)
     *  \return      True if input file contains that data field
     */ 
    inline bool containsData(std::string str) const {
      auto pr = splitQuery(str);
      if( not containsSection(pr.first) ) return false;
      return dict_.at(pr.first).find(pr.second) != dict_.at(pr.first).end();
    }
  




    inline std::vector<std::string> getDataInSection( std::string section )  {

      std::vector<std::string> datasets;

      if( containsSection(section) ) {

        for(auto & data : dict_[section])
          datasets.emplace_back(data.first);

      }

      return datasets;

    }

  }; // CQInputFile class
  
  
  // Misc string functions
  
  /**
   *  Trim a string of left trailing whitespace
   *
   *  \param [in/out] s std::string to be trimmed
   */
  static inline std::string& trim_left(std::string &s) {
      s.erase(s.begin(), std::find_if(s.begin(), s.end(),
              std::not1(std::ptr_fun<int, int>(std::isspace))));
      return s;
  }; // trim_left
  
  
  /**
   *  Trim a string of right trailing whitespace
   *
   *  \param [in/out] s std::string to be trimmed
   */
  static inline std::string& trim_right(std::string &s) {
      s.erase(std::find_if(s.rbegin(), s.rend(),
              std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
      return s;
  }; // trim_right
  
  
  /**
   *  Trim a string of trailing whitespace from both ends
   *
   *  \param [in/out] s std::string to be trimmed
   */
  static inline std::string &trim(std::string &s) {
      return trim_left(trim_right(s));
  }; // trim
  
  /**
   *  Splits a string into tokens  based on a demiliter
   *
   *  \param [out] tokens     std::vector of std::string objects which hold
   *                          the split tokens
   *  \param [in]  str        std::string to split
   *  \param [in]  delimiters Delimiters on which to split str
   */
  static inline void split(std::vector<std::string>& tokens, 
    const std::string& str, const std::string& delimiters = " ") {
  
      tokens.clear();
      // Skip delimiters at beginning.
      std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
      // Find first "non-delimiter".
      std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
  
      while (std::string::npos != pos || std::string::npos != lastPos) {
          // Found a token, add it to the vector.
          tokens.push_back(str.substr(lastPos, pos - lastPos));
          // Skip delimiters.  Note the "not_of"
          lastPos = str.find_first_not_of(delimiters, pos);
          // Find next "non-delimiter"
          pos = str.find_first_of(delimiters, lastPos);
      }
  }; // split

}; // namespace ChronusQ

#endif
