#include "Givens.h"

namespace bms
{
  Givens::Givens()
    : Givens(0,0,1,0)
  {
  }

  Givens::Givens(int i, int j, double c, double s)
    : i_(i), j_(j), Jt_(c,-s)
  {
  }

  Givens::Givens(int i, double c, double s)
    : Givens(i,i+1,c,s)
  {
  }

  void Givens::extend(int incr)
  {
    i_ += incr;
    j_ += incr;
  }
}