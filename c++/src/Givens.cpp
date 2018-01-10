#include <bms/Givens.h>

namespace bms
{
  Givens::Givens()
    : Givens(0,0,1,0)
  {
  }

  Givens::Givens(Index i, Index j, double c, double s)
    : i_(i), j_(j), Jt_(c,-s)
  {
  }

  Givens::Givens(Index i, double c, double s)
    : Givens(i,i+1,c,s)
  {
  }

  void Givens::extend(Index incr)
  {
    i_ += incr;
    j_ += incr;
  }
}
