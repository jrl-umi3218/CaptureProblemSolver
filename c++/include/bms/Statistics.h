#pragma once

#include <vector>

#ifdef USE_STATS
# define STATISTICS(x) x
#else
# define STATISTICS(x) (void)0
#endif

namespace cps
{
  namespace stats
  {
    struct LSStats
    {
      void reset();

      int iter;
      int activation;
      int deactivation;
      int rankLoss;
      int activeConstraints;
    };

    struct SQPStats
    {
      void reset();

      std::vector<int> lineSearchSteps; //number of steps for each iter
      std::vector<LSStats> lsStats;
    };

    inline void LSStats::reset()
    {
      iter = 0;
      activation = 0;
      deactivation = 0;
      rankLoss = 0;
      activeConstraints = 0;
    }

    inline void SQPStats::reset()
    {
      lineSearchSteps.clear();
      lsStats.clear();
    }
  }
}