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

#include <aointegrals.hpp>
#include <cqlinalg.hpp>
#include <cqlinalg/svd.hpp>
#include <cqlinalg/blasutil.hpp>
#include <physcon.hpp>
#include <util/matout.hpp>
#include <util/threads.hpp>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Core>

// Debug directives
//#define _DEBUGORTHO
//#define _DEBUGERI
//#define _DEBUGGIAOERI //SS
//#define _DEBUGGIAOONEE //SS 


namespace ChronusQ {

  typedef std::vector<libint2::Shell> shell_set; 

  // EFP INTEGRAL ENGINE
  template <typename IntsT>
  std::vector<IntsT*> AOIntegrals<IntsT>::OneEDriverEFP(
    libint2::Operator op, shell_set& shells, int deriv, 
	std::array<double,3> xyz) {

    // Determine the number of basis functions for the passed shell set
    size_t NB = std::accumulate(shells.begin(),shells.end(),0,
      [](size_t init, libint2::Shell &sh) -> size_t {
        return init + sh.size();
      }
    );

    size_t NBSQ = NB*NB;


    // Determine the maximum angular momentum of the passed shell set
    int maxL = std::max_element(shells.begin(), shells.end(),
      [](libint2::Shell &sh1, libint2::Shell &sh2){
        return sh1.contr[0].l < sh2.contr[0].l;
      }
    )->contr[0].l;

    // Determine the maximum contraction depth of the passed shell set
    int maxPrim = std::max_element(shells.begin(), shells.end(),
      [](libint2::Shell &sh1, libint2::Shell &sh2){
        return sh1.alpha.size() < sh2.alpha.size();
      }
    )->alpha.size();

    // Determine the number of OpenMP threads
    int nthreads = GetNumThreads();

    // Create a vector of libint2::Engines for possible threading
    std::vector<libint2::Engine> engines(nthreads);

    // Initialize the first engine for the integral evaluation
    engines[0] = libint2::Engine(op,maxPrim,maxL,deriv);
    engines[0].set_precision(0.0);


    // If engine is V, define nuclear charges
    std::vector<std::pair<double,std::array<double,3>>> q;
    q.push_back( { static_cast<double>(1), xyz } );
    engines[0].set_params(q);

    // Copy over the engines to other threads if need be
    for(size_t i = 1; i < nthreads; i++) engines[i] = engines[0];


    // Determine the number of operators
    std::vector<IntsT*> mats( engines[0].results().size() );

    std::vector<
      Eigen::Map<
        Eigen::Matrix<IntsT,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>
      > 
    > matMaps;
    for( auto i = 0; i < mats.size(); i++ ) {
      mats[i] = memManager_.template malloc<IntsT>(NBSQ);
      std::fill_n(mats[i],NBSQ,0.);
      matMaps.emplace_back(mats[i],NB,NB);
    }


    #pragma omp parallel
    {
      int thread_id = GetThreadID();

      const auto& buf_vec = engines[thread_id].results();
      size_t n1,n2;

      // Loop over unique shell pairs
      for(size_t s1(0), bf1_s(0), s12(0); s1 < shells.size(); bf1_s+=n1, s1++){ 
        n1 = shells[s1].size(); // Size of Shell 1
      for(size_t s2(0), bf2_s(0); s2 <= s1; bf2_s+=n2, s2++, s12++) {
        n2 = shells[s2].size(); // Size of Shell 2

        // Round Robbin work distribution
        #ifdef _OPENMP
        if( s12 % nthreads != thread_id ) continue;
        #endif

        // Compute the integrals       
        engines[thread_id].compute(shells[s1],shells[s2]);

        // If the integrals were screened, move on to the next batch
        if(buf_vec[0] == nullptr) continue;

        // Place integral blocks into their respective matricies
        // XXX: USES EIGEN
        for(auto iMat = 0; iMat < buf_vec.size(); iMat++){
          Eigen::Map<
            const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,
              Eigen::RowMajor>>
            bufMat(buf_vec[iMat],n1,n2);

          matMaps[iMat].block(bf1_s,bf2_s,n1,n2) = bufMat.template cast<IntsT>();
        }

      } // Loop over s2 <= s1
      } // Loop over s1

    } // end OpenMP context


    // Symmetrize the matricies 
    // XXX: USES EIGEN
    for(auto nMat = 0; nMat < matMaps.size(); nMat++) 
      matMaps[nMat] = matMaps[nMat].template selfadjointView<Eigen::Lower>();


    return mats;

  }; // AOIntegrals::OneEDriver
  template std::vector<double*> AOIntegrals<double>::OneEDriverEFP(
    libint2::Operator op, shell_set& shells, int deriv,
	std::array<double,3> xyz);
  template std::vector<dcomplex*> AOIntegrals<dcomplex>::OneEDriverEFP(
    libint2::Operator op, shell_set& shells, int deriv,
	std::array<double,3> xyz);

  template <typename IntsT>
  IntsT* AOIntegrals<IntsT>::computeEFPContributions_Coulomb(
	double* xyz_efp, double* mult, size_t* n_mult) {
        size_t NBS = basisSet_.nBasis;
        // Allocate an NBasis**2 size matrix to return
        auto resultMatrix = memManager_.template malloc<IntsT>(NBS*NBS);
        memset(resultMatrix, 0, NBS*NBS*sizeof(IntsT));

        int nthreads = GetNumThreads();

        if(envE0.size() == 0 and envE1.size() == 0 and envE2.size() == 0)
          computeEFP_integrals(xyz_efp, n_mult, 0);        
        std::vector<std::array<IntsT,3> > dipole(nthreads);
        std::vector<std::array<IntsT,3> > dipole_1(nthreads);
        std::vector<std::array<IntsT,6> > quadrupole(nthreads);
	std::vector<std::array<IntsT,6> > quadrupole_1(nthreads);

        #pragma omp parallel for
        for(int i = 0; i < *n_mult; i++){
          int thread_id = GetThreadID();
           
          // Get OneEDriverEFP for a given multipole site (deriv=0) [V] 
          IntsT monopole = *(mult+20*i);
          // LAPACK DSCAL call on resulting integral [qV]
          Scale(NBS*NBS, monopole, envE0[i], 1);
          // Dipole
          // Get OneEDriverEFP (deriv=1) [E]
          
          dipole[thread_id][0] =  *(mult+1+20*i);
          dipole[thread_id][1] =  *(mult+2+20*i);
          dipole[thread_id][2] =  *(mult+3+20*i);
          
          for(int j = 0; j < NBS*NBS; j++){
			dipole_1[thread_id][0] = *(envE1[i][0]+j);
			dipole_1[thread_id][1] = *(envE1[i][1]+j);
			dipole_1[thread_id][2] = *(envE1[i][2]+j);
              auto temp = InnerProd<IntsT,IntsT,IntsT>(3, &(dipole[thread_id][0]), 1, &(dipole_1[thread_id][0]), 1);

              #pragma omp critical
              *(resultMatrix+j) +=  *(envE0[i]+j) + temp;
          }

          // Quadrupole
          // Get OneEDriverEFP (deriv=2) [E']
          quadrupole[thread_id] = {    
              *(mult+4+20*i),   *(mult+7+20*i), *(mult+8+20*i), 
                                *(mult+5+20*i), *(mult+9+20*i),
                                                *(mult+6+20*i)
          };
          // LAPACK DTRMM [OE']
		  // Maybe allocate via memManager_
          for(int j = 0; j < NBS*NBS; j++){
			quadrupole_1[thread_id][0] = *(envE2[i][0]+j);
			quadrupole_1[thread_id][1] = *(envE2[i][1]+j);
			quadrupole_1[thread_id][2] = *(envE2[i][2]+j);
			quadrupole_1[thread_id][3] = *(envE2[i][3]+j);
			quadrupole_1[thread_id][4] = *(envE2[i][4]+j);
			quadrupole_1[thread_id][5] = *(envE2[i][5]+j);
            auto trans_ =(quadrupole_1[thread_id][0]*quadrupole[thread_id][0]
                         +quadrupole_1[thread_id][3]*quadrupole[thread_id][3]
                         +quadrupole_1[thread_id][5]*quadrupole[thread_id][5]
                         +quadrupole_1[thread_id][1]*quadrupole[thread_id][1]*2.
                         +quadrupole_1[thread_id][2]*quadrupole[thread_id][2]*2.
                         +quadrupole_1[thread_id][4]*quadrupole[thread_id][4]*2.)/3.0;
              #pragma omp critical
              *(resultMatrix+j) +=  trans_;
          } 
          Scale(NBS*NBS, IntsT(1.)/monopole, envE0[i], 1);

        }
        
        return resultMatrix;
        
  };
  template double* AOIntegrals<double>::computeEFPContributions_Coulomb(
	double* xyz_efp, double* mult, size_t* n_mult); 
  template dcomplex* AOIntegrals<dcomplex>::computeEFPContributions_Coulomb(
	double* xyz_efp, double* mult, size_t* n_mult);

 
  // Induced dipole contributions
  template <typename IntsT>
  IntsT* AOIntegrals<IntsT>::computeEFPContributions_Pol(
	double* xyz_pol, double* pol_aver, size_t* n_pol) {
        size_t NBS = basisSet_.nBasis;
        auto resultMatrix = memManager_.template malloc<IntsT>(NBS*NBS);
        memset(resultMatrix, 0, NBS*NBS*sizeof(IntsT));
        time_t currentClockTime;
        time(&currentClockTime);
        std::cout << ctime(&currentClockTime) << std::endl;

        int nthreads = GetNumThreads();

        std::vector<std::array<IntsT,3> > dipole(nthreads);
        std::vector<std::array<IntsT,3> > dipole_1(nthreads);

        if(indE1.size() == 0)
          computeEFP_integrals(xyz_pol, n_pol, 1); 
        #pragma omp parallel for
        for(int i = 0; i < *n_pol ; i++){
          
          int thread_id = GetThreadID();

          dipole[thread_id][0] = *(pol_aver+3*i);
          dipole[thread_id][1] = *(pol_aver+1+3*i);
          dipole[thread_id][2] = *(pol_aver+2+3*i);

          
          // LAPACK DDOT [uE]
		  // Maybe want to allocate via memManager_
          for(int j = 0; j < NBS*NBS; j++){
                dipole_1[thread_id][0] = *(indE1[i][0]+j);
                dipole_1[thread_id][1] = *(indE1[i][1]+j);
                dipole_1[thread_id][2] = *(indE1[i][2]+j);

              auto temp = InnerProd<IntsT,IntsT,IntsT>(3, &(dipole[thread_id][0]), 1, &(dipole_1[thread_id][0]), 1);
              #pragma omp critical
              *(resultMatrix+j) += temp;  
          } 
        
        }
        time(&currentClockTime);
        std::cout << ctime(&currentClockTime) << std::endl;
        return resultMatrix;
  };
  template double* AOIntegrals<double>::computeEFPContributions_Pol(
	double* xyz_pol, double* pol_aver, size_t* n_pol);
  template dcomplex* AOIntegrals<dcomplex>::computeEFPContributions_Pol(
	double* xyz_pol, double* pol_aver, size_t* n_pol);

  template <typename IntsT>
  void AOIntegrals<IntsT>::computeEFP_integrals(
	double* xyz, size_t* num, int option){
    std::array<double,3> coord_;
    size_t NBS = basisSet_.nBasis;
    
    if(option == 0){
      for(int i = 0; i < *num; i++){
        coord_[0] = *(xyz+3*i);
        coord_[1] = *(xyz+1+i*3);
        coord_[2] = *(xyz+2+i*3);
      
        auto E0V = OneEDriverEFP(libint2::Operator::nuclear, basisSet().shells, 0, coord_); 
        envE0.push_back(E0V[0]);
        auto E1V = OneEDriverEFP(libint2::Operator::nuclear, basisSet().shells, 1, coord_);
        envE1.push_back(std::vector<IntsT*>{E1V[6], E1V[7], E1V[8]});
        auto E2V = OneEDriverEFP(libint2::Operator::nuclear, basisSet().shells, 2, coord_);
        envE2.push_back(std::vector<IntsT*>{E2V[39], E2V[40], E2V[41], E2V[42], E2V[43], E2V[44]});
      }
    }
    else if(option == 1){
      for(int i = 0; i < *num; i++){
        coord_[0] = *(xyz+i*3);
        coord_[1] = *(xyz+1+i*3);
        coord_[2] = *(xyz+2+i*3);
        
        auto E1V = OneEDriverEFP(libint2::Operator::nuclear, basisSet().shells, 1, coord_);
        indE1.push_back(std::vector<IntsT*>{E1V[6],E1V[7],E1V[8]});
      }
    }
  }; 
  template void AOIntegrals<double>::computeEFP_integrals(
	double* xyz, size_t* num, int option);
  template void AOIntegrals<dcomplex>::computeEFP_integrals(
	double* xyz, size_t* num, int option);
  
}; // namespace ChronusQ



