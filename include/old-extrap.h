/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2016 Li Research Group (University of Washington)
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

template<typename F>
class DIIS {
  // Useful Typedefs
  typedef Eigen::Matrix<F,Dynamic,Dynamic,ColMajor> FMat;
  typedef Eigen::Map<FMat> FMap;

  std::size_t               nExtrap_;
  std::size_t               eDim_;
  std::vector<F>            coeffs_;
  std::vector<H5::DataSet*> eFiles_;
  const std::function<void(H5::DataSet*,std::size_t,F*)> &read_;
  const std::function<F*(std::size_t)> &allocF_;
  const std::function<void(F*,std::size_t)> &freeF_;
  std::size_t               iCurrent_;
  bool doPrune_;

  bool doNeg_;

public:

  DIIS() = delete;
  DIIS(std::size_t nExtrap, std::size_t eDim,
    const std::function<void(H5::DataSet*,std::size_t,F*)> &read,
    const std::function<F*(std::size_t)> &allocF,
    const std::function<void(F*,std::size_t)> &freeF,
    std::vector<H5::DataSet*> &eFiles, 
    std::size_t iCurrent, bool doPrune = false, bool doNeg = true):
    nExtrap_(nExtrap), eDim_(eDim), read_(read), 
    eFiles_(eFiles), allocF_(allocF),freeF_(freeF),
    iCurrent_(iCurrent), doPrune_(doPrune), doNeg_(doNeg){

    coeffs_.resize(nExtrap_+1);
  }

  bool extrapolate(std::ostream &out=cout);


  F* coeffs() { return coeffs_.data(); }
};

template<typename F>
bool DIIS<F>::extrapolate(std::ostream &out){

  int N(nExtrap_ + 1);
  int NRHS(1); 
  int INFO(1);
  int LWORK(5*N);
  int indRemove;
  char NORM = 'O';
  double RCOND(0);
  bool InvFail(false);

  std::vector<F> BMem(N*N);
  std::vector<F> tempcoeff(N);
  FMap B(BMem.data(),N,N);
  std::vector<int> iPiv(N), iWork(N); //TODO: do we need iPiv since we don't use it?
  std::vector<F> BDiag;
  Eigen::Matrix<F,Dynamic,Dynamic,ColMajor> ReducedB(N,N);
  Eigen::Matrix<F,Dynamic,Dynamic,ColMajor> OriginalB(N,N);

  F* SCRATCH1 = allocF_(eDim_);
  F* SCRATCH2 = allocF_(eDim_);

  FMap S1(SCRATCH1,eDim_,1);
  FMap S2(SCRATCH2,eDim_,1);

  B.setZero();
  for(auto iFile : eFiles_) {
    for(auto j = 0ul; j < nExtrap_; j++){
      read_(iFile,j,SCRATCH1);
      for(auto k = 0ul; k <= j; k++){
        read_(iFile,k,SCRATCH2);
        B(j,k) += S1.frobInner(S2.conjugate());
      }
    }
  }

  B = B.template selfadjointView<Lower>();

  double fact = -1.0;
  if(!doNeg_) fact = 1.0;

  for(auto l = 0; l < nExtrap_; l++){
    B(nExtrap_,l) = fact;
    B(l,nExtrap_) = fact;
  }
  B(nExtrap_,nExtrap_) = 0.0;

  // Save B so it can be used later if the extrapolation space is being
  // pruned
  ReducedB = B;
  OriginalB = B;

  // Loop over solution of linear equation checking that the calculated
  // coefficients meet certain criteria and removing rows/cols from the 
  // matrix and resolving if necessary
  int NSaved = N;
  double coeffLimit = 10.0;
  bool goodDIIS = false;
  while (!goodDIIS) {

    std::fill(tempcoeff.begin(),tempcoeff.end(),0.0);
    tempcoeff[NSaved-1] = fact;

    double ANORM = B.template lpNorm<1>();
    std::vector<F> WORK(LWORK);

    if(std::is_same<F,dcomplex>::value) {
      std::vector<double> RWORK(LWORK);
      // Linear Solve
      zgesv_(&NSaved,&NRHS,reinterpret_cast<dcomplex*>(B.data()),&NSaved,
          iPiv.data(),reinterpret_cast<dcomplex*>(tempcoeff.data()),&NSaved,&INFO);

      InvFail = (INFO != 0);

      // Obtain condition number from LU given from ZGESV
      zgecon_(&NORM,&NSaved,reinterpret_cast<dcomplex*>(B.data()),&NSaved,
          &ANORM,&RCOND,reinterpret_cast<dcomplex*>(WORK.data()),RWORK.data(),
          &INFO);

    } else if(std::is_same<F,double>::value) {
      // Linear Solve
      dgesv_(&NSaved,&NRHS,reinterpret_cast<double*>(B.data()),&NSaved,
          iPiv.data(),reinterpret_cast<double*>(tempcoeff.data()),&NSaved,&INFO);

      InvFail = (INFO != 0);

      // Obtain condition number from LU given from DGESV
      dgecon_(&NORM,&NSaved,reinterpret_cast<double*>(B.data()),&NSaved,
          &ANORM,&RCOND,reinterpret_cast<double*>(WORK.data()),iWork.data(),
          &INFO);
    }

    // Check that the coefficients are below the max allowed value. If they're
    // too high remove the row/column from B corresponding to the Fock matrix
    // with the highest energy. B cannot have a dimension smaller than 3x3.
    auto max_coeff = std::max_element(tempcoeff.begin(), tempcoeff.end()-1,
      [](const F& a,const F& b) {
        return std::abs(a) < std::abs(b); 
      });
    if (std::abs(*max_coeff) > coeffLimit and NSaved > 2 and doPrune_) {
//    out << "WARNING COEFF IS TOO HIGH = " << *max_coeff << endl;

      // Find which matrix to remove and update order list to be used again if the
      // redimensioned B matrix still has a problem
      // TODO: add check not to remove current fock matrix
      // TODO: check for iCurrent doesn't account for changing extrapolation length
      indRemove = 0; 
      for (auto i = 1; i < NSaved-1; i++){
        if ( std::abs(ReducedB(i,i)) > std::abs(ReducedB(indRemove,indRemove)) and 
             i != iCurrent_) {
          indRemove = i; } }
      BDiag.push_back(ReducedB(indRemove,indRemove));

      // Remove a row and column from Reduced B and put it back in B
      ReducedB.block(  indRemove  ,0,NSaved-indRemove-1,NSaved)   
      = ReducedB.block(indRemove+1,0,NSaved-indRemove-1,NSaved);
      ReducedB.conservativeResize(NSaved-1,NSaved);
      ReducedB.block(  0,indRemove  ,NSaved-1,NSaved-indRemove-1) 
      = ReducedB.block(0,indRemove+1,NSaved-1,NSaved-indRemove-1);
      ReducedB.conservativeResize(NSaved-1,NSaved-1);
      new (&B) FMap (BMem.data(),NSaved-1,NSaved-1);
      B = ReducedB;

      NSaved -= 1;

    } else {
      goodDIIS=true; 
    }

  } // end while

  freeF_(SCRATCH1,eDim_);
  freeF_(SCRATCH2,eDim_);

  // Build the final coefficient vector placing zeros in place for the 
  // matrices removed from the extrapolation
  if (doPrune_) {
    int j = 0;
    for (auto i = 0; i < nExtrap_; i++){
      auto result1 = std::find(BDiag.begin(),BDiag.end(),OriginalB(i,i));
      if (result1 != BDiag.end()) {
        coeffs_[i] = 0.; 
      } else {
        coeffs_[i] = tempcoeff[j];
        j++;
      }
    }
/*
    for (auto i = 0; i < BDiag.size(); i++){
      out << "BDiag = " << BDiag[i] << endl; }
    for (auto i = 0; i < nExtrap_; i++){
      out << "coeff = " << coeffs_[i] << endl; }
*/
  } else {
    for (auto i = 0; i < nExtrap_; i++){
      coeffs_[i] = tempcoeff[i];
    }
  }

  return !InvFail;
};

