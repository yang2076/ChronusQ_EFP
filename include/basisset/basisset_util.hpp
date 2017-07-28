
#include <libint2/shell.h>
#include <memmanager.hpp>

namespace ChronusQ {


  /**
   *  \brief Specification for shellset evaluation.
   *
   *  Currently supports GRADIENT and NOGRAD (simple eval) evaluation
   *  types.
   */ 
  enum SHELL_EVAL_TYPE {
    GRADIENT,
    NOGRAD
  };
 
  /**
   *  \brief Level 1 Basis Set Evaluation Function
   *  \brief Evaluates a shell set over a specified number of cartesian points.
   */ 
  void evalShellSet(CQMemManager &,SHELL_EVAL_TYPE, std::vector<libint2::Shell> &, double *, size_t, double *);

  inline void evalShellSet(CQMemManager &memManager,SHELL_EVAL_TYPE typ, std::vector<libint2::Shell> &shells, 
    std::vector<std::array<double,3>> &pts, double *eval){

    evalShellSet(memManager,typ,shells,&(pts[0][0]),pts.size(),eval);

  }; // evalShellSet (over vector of arrays)

  /**
   *  \brief Level 2 Basis Set Evaluation Function.
   *  \brief Evaluates a shell set over a specified number of cartesian points. This 
   *  \brief function requires a precomputed set of distances and their x,y,z components 
   *  \brief for each point from each shell origin in the shells vector..
   */ 
  void evalShellSet(SHELL_EVAL_TYPE, std::vector<libint2::Shell> &, double *, double *, size_t,size_t, double*);

  /**
   *  \brief Level 3 Basis Set Evaluation Function
   *  \brief Evaluates a single shell over a single cartesian point. This function requires a precomputed
   *  \brief distance and its x,y,z components for the point from the shell origin. An offset
   *  \brief to properly store the results can be used..
   */ 
  void evalShellSet(SHELL_EVAL_TYPE,const libint2::Shell&,double,const std::array<double,3>&, double *, size_t);

}; // namespace ChronusQ
