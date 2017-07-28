#include <basisset/basisset_util.hpp>

namespace ChronusQ {

  /**
   *  Level 1 Basis Set Evaluation Function
   *  Evaluates a shell set over a specified number of cartesian points.
   *  All distances of the points from each shell origin are computed and used
   *  in the Level 2 Basis Set Evaluation Function
   *  \param [in] memManager CQ memory manager (for allocate distances inside)
   *  \param [in] typ        Type of evaluation to perform (gradient, etc)
   *  \param [in] shells     Shell set for evaluation(vector of libint2::Shell).
   *  \param [in] pts        Raw storage of the cartesian points (dimension 3 * npts)
   *  \param [in] npts       Number of cartesian points to evaluate.
   *  \param [in/out] fEval  eval Storage for the shell set evaluation, f(ixyz,iSh,ipt). 
   *                         Variable dimensions(npts*NBasisEff*1 or 4) allocated outside.
   *                         This storage will have all values of the functions in the shell, 
   *                         for each shell, for each point. If requested there will be appended 
   *                         f/dx(ixyz,iSh,ipt), f/dy(ixyz,iSh,ipt), f/dz(ixyz,iSh,ipt) values of 
   *                         the functions in the shell, for each shell(nShSize), for each point(npts)
   */ 
  void evalShellSet(CQMemManager &memManager,SHELL_EVAL_TYPE typ, std::vector<libint2::Shell> &shells, 
    double *pts, size_t npts, double *fEval) {
    size_t NBasisEff = 0;
    size_t nShSize = shells.size();
    double * r   = memManager.malloc<double>(3*npts*nShSize);
    double * rSq = memManager.malloc<double>(npts*nShSize);
    //figure the size (number of basis) of the all shells inputed that need to be evaluated
    // and store in NBasisEff. It will be used for defining the pointers later on.
    for (auto iSh = 0; iSh < nShSize; iSh++){
    NBasisEff += shells[iSh].size();
    }
    for (auto ipts = 0; ipts < npts; ipts++){
      double xp = *(pts + ipts*3);
      double yp = *(pts + 1 + ipts*3);
      double zp = *(pts + 2 + ipts*3);
      //we need this counter to keep track of the nbasis evaluated for all shells (at a given point)
      for (auto iSh = 0; iSh < nShSize; iSh++){
        r[   iSh*3 + ipts*3*nShSize] =  xp - shells[iSh].O[0];
        r[1+ iSh*3 + ipts*3*nShSize] =  yp - shells[iSh].O[1];
        r[2+ iSh*3 + ipts*3*nShSize] =  zp - shells[iSh].O[2];
        rSq[ iSh + ipts*nShSize] = r[    iSh*3 + ipts*3*nShSize]*r[    iSh*3 + ipts*3*nShSize] + 
                                       r[1 + iSh*3 + ipts*3*nShSize]*r[1 + iSh*3 + ipts*3*nShSize] + 
                                       r[2 + iSh*3 + ipts*3*nShSize]*r[2 + iSh*3 + ipts*3*nShSize]; 
        /*
        //Debug Printing
        std::cout << ipts << " " << iSh << " xdiff1 " << r[    iSh*3 + ipts*3*nShSize]  << 
                                           " ydiff1 " << r[1 + iSh*3 + ipts*3*nShSize] << 
                                           " zdiff1 " << r[2 + iSh*3 + ipts*3*nShSize]  <<
                                           " rSq1"    << rSq[  iSh   + ipts*nShSize] << std::endl;
        */
      }
    }
    //Call to Level 2 Basis Set Evaluation
    evalShellSet(typ,shells,rSq,r,npts,NBasisEff,fEval); 
  }; // evalShellSet Level 3

  /**
   *  Level 2 Basis Set Evaluation Function
   *  Evaluates a shell set over a specified number of cartesian points. This function requires a precomputed
   *  set of the distances and their x,y,z component for each point from each shell origin in the shells vector.
   *  \param [in] typ        Type of evaluation to perform (gradient, etc)
   *  \param [in] shells     Shell set for evaluation(vector of libint2::Shell).
   *  \param [in] rSq        Raw storage of overall distances between each point and the shell origin (precomputed outside)
   *                         Dimension (nshells * npts)
   *  \param [in] r          Raw storage of overall x,y,z of the vector between each point and the shell origin 
   *                         (precomputed outside). Dimension (3*nshells*npts) 
   *  \param [in] npts       Number of cartesian points to evaluate.
   *  \param [in] NBasisEff  Number of basis function to be evaluated given all the shells in input.
   *  \param [in/out] fEval  eval Storage for the shell set evaluation, f(ixyz,iSh,ipt). 
   *                         Variable dimensions(npts*NBasisEff*1 or 4) allocated outside.
   *                         This storage will have all values of the functions in the shell, 
   *                         for each shell, for each point. If requested there will be appended 
   *                         f/dx(ixyz,iSh,ipt), f/dy(ixyz,iSh,ipt), f/dz(ixyz,iSh,ipt) values of 
   *                         the functions in the shell, for each shell(nShSize), for each point(npts)
   */ 
  void evalShellSet(SHELL_EVAL_TYPE typ, std::vector<libint2::Shell> &shells, double* rSq,double *r, size_t npts, 
    size_t NBasisEff, double *fEval) {

    size_t nShSize = shells.size();
    size_t IOff =  npts*NBasisEff;
    for (auto ipts = 0; ipts < npts; ipts++){
      auto Ic = 0;
      for (auto iSh = 0; iSh < nShSize; iSh++){
        double * fStart    = fEval+Ic+ipts*NBasisEff;
        double  xdiff     = r[iSh*3 + ipts*3*nShSize];
        double  ydiff     = r[1 + iSh*3 + ipts*3*nShSize];
        double  zdiff     = r[2 + iSh*3 + ipts*3*nShSize];
        double  rSqVal    = rSq[iSh + ipts*nShSize];
        std::array<double,3> rVal ({xdiff,ydiff,zdiff});
        /*
        //Debug Printing
        std::cout <<"Distance " << ipts <<" " <<rVal[0] << " "<< rVal[1] <<" "<<rVal[2]<<" "
        <<std::sqrt(rSqVal) << std::endl;
        */
        evalShellSet(typ,shells[iSh],rSqVal,rVal,fStart,IOff); 
        Ic += shells[iSh].size();
      } // loop over shells
    } // loop over points
  }; // evalShellSet Level 2


  /**
   *   Level 3 Basis Set Evaluation Function
   *  
   *   Evaluates a single shell over a single cartesian point. This function requires a precomputed
   *   the distance and its x,y,z components for the point from the shell origin. An offset
   *   to properly store the results can be used.
   *   \param [in] typ        Type of evaluation to perform (gradient, etc)
   *   \param [in] shell      Shell for evaluation(libint2::Shell).
   *   \param [in] rSq        Raw storage of square distance between the point and the shell origin (precomputed outside)
   *   \param [in] r          Raw storage of x,y,z of the vector between each point and the shell origin 
   *                          (precomputed outside).  
   *   \param [in/out] fEval  eval Storage for the shell set evaluation, f(ixyz,iSh,ipt). 
   *                          Variable dimensions(npts*NBasisEff*1 or 4) allocated outside.
   *                          This storage will have all values of the functions in the shell, 
   *                          for each shell, for each point. If requested there will be appended 
   *                          f/dx(ixyz,iSh,ipt), f/dy(ixyz,iSh,ipt), f/dz(ixyz,iSh,ipt) values of 
   *                          the functions in the shell, for each shell(nShSize), for each point(npts)
   *   \param [in] IOff       OffSet to properly store the basis set. 
   */ 
  void evalShellSet(SHELL_EVAL_TYPE typ, const libint2::Shell &shell,double rSq, const std::array<double,3> &xyz, 
    double *fEval, size_t IOff) {
    auto L         = shell.contr[0].l;
    auto shSize    = ((L+1)*(L+2))/2; 
    double * f     = fEval ;
    double * dx = f   + IOff;
    double * dy = dx  + IOff;
    double * dz = dy  + IOff;
    auto contDepth = shell.alpha.size(); 
    double alpha(0.0);
    double expFactor(0.0);
    double expArg(0);
    double tmpcoef,tmpalpha;
    int lx,ly,lz, ixyz;
    double tmpxyz;
    double tmpdx;
    double tmpdy;
    double tmpdz;
    // Generating the expArgument, expFactotr and the
    // alpha (for derivatives later on) and store them
    // in temp variables
    for(auto k = 0; k < contDepth; k++){
      tmpcoef = shell.contr[0].coeff[k];
      tmpalpha = shell.alpha[k];
      expArg = std::exp(-tmpalpha*rSq);
      expFactor += tmpcoef * expArg;
      if (typ == GRADIENT) { 
        // quantities for derivatives
        tmpcoef *= tmpalpha;
        alpha += tmpcoef * expArg;
      }
    } 
    if (typ == GRADIENT)  { alpha *= 2;}
    for(auto i = 0u, I = 0u; i <= L; i++) {
      lx = L - i;
      for( auto j = 0u; j <= i; j++, I++) {
        ly = i - j;
        lz = L - lx - ly;
        tmpxyz= 1.0;
        tmpdx = 0.0;
        tmpdy = 0.0;
        tmpdz = 0.0;
        for(ixyz = 0; ixyz < lx-1; ixyz++) tmpxyz *= xyz[0];
        for(ixyz = 0; ixyz < ly-1; ixyz++) tmpxyz *= xyz[1];
        for(ixyz = 0; ixyz < lz-1; ixyz++) tmpxyz *= xyz[2];
        f[I]  =  tmpxyz;
        if (typ == GRADIENT) {
        // Derivatives
          if(lx> 0) {tmpdx = -expFactor * lx;}
          if(ly> 0) {tmpdy = -expFactor * ly;}
          if(lz> 0) {tmpdz = -expFactor * lz;}
           
          dx[I] = tmpxyz*tmpdx;
          dy[I] = tmpxyz*tmpdy;
          dz[I] = tmpxyz*tmpdz;
    
          // finishing up        
          if(lx> 0) {f[I]  *= xyz[0]; dy[I] *=xyz[0];dz[I] *=xyz[0];}
          if(ly> 0) {f[I]  *= xyz[1]; dx[I] *=xyz[1];dz[I] *=xyz[1];}
          if(lz> 0) {f[I]  *= xyz[2]; dx[I] *=xyz[2];dy[I] *=xyz[2];}
    
          dx[I] += f[I] * xyz[0] * alpha;
          dy[I] += f[I] * xyz[1] * alpha;
          dz[I] += f[I] * xyz[2] * alpha;
          f[I]  *= expFactor;
    
        } else{
        // Only basis (not GGA)
          if(lx> 0) {f[I]  *= xyz[0];}
          if(ly> 0) {f[I]  *= xyz[1];}
          if(lz> 0) {f[I]  *= xyz[2];}
          f[I]  *= expFactor;
        }
        /*
          // Debug Printing
          std::cout << I <<" "<< lx << " " 
            << ly << " "<<lz <<"  f(pt) "<< f[I] <<" " <<std::endl;
          std::cout << I<<" "<< lx << " "  
            << ly << " "<<lz <<" dx(pt) "<< dx[I] << std::endl;
          std::cout << I<<" "<< lx << " " 
            << ly << " "<<lz <<" dy(pt) "<< dy[I] << std::endl;
          std::cout << I<<" "<< lx << " " 
            << ly << " "<<lz <<" dz(pt) "<< dz[I] << std::endl;
        */
      } //loop overj, j[0,i]
    } //loop over i, i[0,L] this to loop required to build the lx,ly,lz combination given L
  }; // evalShellSet Level3

  
  /**
   *  \brief only for debug 
   *  it tests 5 pts and the evaluation of the f,dx,dy,dz at those points for a vector of shells
   *
   */ 
/* DBWY: use a template for UT later
  void testEval(CQMemManager &memManager,double *SCR, std::vector<libint2::Shell> &vshells){
    std::vector<std::array<double,3>>  testpts;
    testpts.push_back({0,0,0});
    testpts.push_back({0.1,0,0});
    testpts.push_back({0,0.1,0});
    testpts.push_back({0,0.,0.1});
    testpts.push_back({1.,0.5,0.1});
    size_t npts = testpts.size();
    std::cout <<"inside testEval" <<std::endl;
    evalShellSet(memManager,GRADIENT,vshells,&testpts[0][0],npts,SCR);
  }; // testEval
*/

}; // namespace ChronusQ
