#ifndef _SUPERLU_SUPERLUOPTIONS_HPP_
#define _SUPERLU_SUPERLUOPTIONS_HPP_

#include  "environment.hpp"
#include <string>

namespace PEXSI{
  /// @struct SuperLUOptions
  /// @brief A thin interface for passing parameters to set the SuperLU
  /// options.  
  ///
  struct SuperLUOptions{
    /// @brief Number of processors for parallel symbolic factorization.
    ///
    /// numProcSymbFact should be a power of 2, and is only useful when
    ///
    /// ColPerm = "PARMETIS"
    Int              numProcSymbFact;
    /// @brief The maximum pipeline depth. 
    ///
    /// @todo
    /// This option should not be here and should be moved into PMatrix.
    Int              maxPipelineDepth; 

    /// @brief Option of matrix permutation strategy.
    ///
    /// The following options of permutation strategy is available (case
    /// sensitive):
    ///
    /// - "MMD_AT_PLUS_A": Multiple minimum degree ordering. This is the
    /// default option.
    /// - "METIS_AT_PLUS_A": Sequential ordering using METIS. This
    /// requires the usage of METIS package.
    /// - "PARMETIS": Parallel ordering. This requires the usage of
    /// ParMETIS/PT-SCOTCH package.
    ///
    std::string      ColPerm;

    // Member functions to setup the default value
    SuperLUOptions(): numProcSymbFact(0), maxPipelineDepth(-1), ColPerm("MMD_AT_PLUS_A") {}
  };

}

#endif
