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




/**
 *
 *  This file acts as a generic template for ChronusQ UTs through
 *  the Boost.Test framework. When compiling, one must add a definition
 *  of BOOST_TEST_MODULE in the compiler invocation
 *
 *  CXX -DBOOST_TEST_MODULE=ModName 
 *
 */




#include <cxxapi/boilerplate.hpp>

#include <boost/test/impl/unit_test_main.ipp>
#include <boost/test/impl/framework.ipp>
#include <boost/test/unit_test.hpp>
using namespace boost::unit_test;


/**
 *  \brief A global fixture for ChronusQ UTs through the
 *  Boost.Test framework.
 *
 *  Ensures that ChronusQ::initialize and ChronusQ::finalize
 *  only get called once throughout the test module
 */
struct CQConfig {

  CQConfig() { ChronusQ::initialize(); }
  ~CQConfig() { ChronusQ::finalize(); }

};

BOOST_GLOBAL_FIXTURE( CQConfig );



