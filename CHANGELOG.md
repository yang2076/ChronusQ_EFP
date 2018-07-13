                     CHRONUSQ RELEASE HISTORY


  This file simply provides a brief overview of the release history of the
  ChronusQ program as released by the Li Research Group at the University
  of Washington. For the most up-to-date record of functionality, etc,
  please refer to the [ChronusQ Wiki](https://github.com/liresearchgroup/chronusq_public/wiki).


  FORMAT: YYYY-MM-DD

  - 2018-07-13 0.2.0 (BETA)
    - Full integration of GIAO basis set into SCF and RT modules
    - Added RESPONSE module
      - Supports the PolarizationPropagator (TDDFT/TDHF) and ParticleParticlePropagator (pp-RPA/pp-TDA/hh-TDA)
      - RESIDUE -> eigen decomposition
      - (D)FDR -> (damped) frequency depenent response
      - GPLHR for partial diagonalization
        - Supports arbitrary energy domain for diagonalization
    - Full integration of MPI functionality throughout (CQ_ENABLE_MPI)
      - Using MXX for C++11 MPI bindings
      - Using CXXBLACS for C++ interface to BLACS type functionality
      - Integration of ScaLAPACK into RESPONSE module
    - Bump Libint -> 2.4.2
    - Bump Libxc  -> 4.0.4
    - Added support for coverage checks (CQ_ENABLE_COVERAGE)
    - Various logic checks / bug fixes
  

  - 2017-09-01: 0.1.0 (BETA)
    - Complete overhaul of ChronusQ development stream (new repo)
    - Currently tested functionality:
      - Full support for Hartree-Fock and Kohn-Sham references
      - SCF general to Restricted (R), Unrestricted (U), and Generalized (G) references
      - Real-time propagation of R/U/G references
      - X2C Relativistic references (both HF and KS)
      - Full integration of OpenMP in performance critical code
      - Support for both INCORE (full ERI) and DIRECT integral contraction schemes
