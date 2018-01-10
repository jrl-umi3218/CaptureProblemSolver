#include <Eigen/Dense>
#include <iosfwd>

class toMatlab
{
public:
  template<typename Derived>
  toMatlab(const Eigen::DenseBase<Derived>& M)
    : mat(M)
  {}

private:
  Eigen::MatrixXd mat;

  friend std::ostream& operator<< (std::ostream&, const toMatlab&);
};

inline std::ostream& operator<< (std::ostream& o, const toMatlab& tom)
{
  if (tom.mat.cols() == 1)
  {
    Eigen::IOFormat fmt(Eigen::StreamPrecision, 0, ", ", ";", "", "", "[", "]");
    o << tom.mat.format(fmt);
  }
  else
  {
    Eigen::IOFormat fmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
    o << tom.mat.format(fmt);
  }
  return o;
}