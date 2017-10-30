#include "Problem.h"

#include <fstream>
#include <map>
#include <sstream>
#include <vector>

namespace
{
  // trim from left
  inline std::string& ltrim(std::string& s, const char* t = " \t\n\r\f\v")
  {
    s.erase(0, s.find_first_not_of(t));
    return s;
  }

  // trim from right
  inline std::string& rtrim(std::string& s, const char* t = " \t\n\r\f\v")
  {
    s.erase(s.find_last_not_of(t) + 1);
    return s;
  }

  // trim from left & right
  inline std::string& trim(std::string& s, const char* t = " \t\n\r\f\v")
  {
    return ltrim(rtrim(s, t), t);
  }

  double parseDouble_(const std::string& s)
  {
    std::stringstream ss(s);
    double d;
    ss >> d;
    if (ss.fail())
      throw std::runtime_error("Failed to read double value");

    return d;
  }

  Eigen::VectorXd parseVector_(const std::string& s)
  {
    std::stringstream ss(s);
    double d;
    char c;
    std::vector<double> v;
    
    ss >> c;
    while (c != ']')
    {
      ss >> d;
      ss >> c;
      if (ss.fail())
        throw std::runtime_error("Failed to read vector value");

      v.push_back(d);
    }

    return Eigen::Map<Eigen::VectorXd>(&v[0], static_cast<Eigen::DenseIndex>(v.size()));
  }

  double parseDouble(const std::map<std::string, std::string>& table, const std::string& key, bool optional = false, double defaultValue = 0)
  {
    if (table.count(key))
      return parseDouble_(table.at(key));
    else
    {
      if (optional)
        return defaultValue;
      else
        throw std::runtime_error("No element " + key + " found in the file");
    }
  }

  Eigen::VectorXd parseVector(const std::map<std::string, std::string>& table, const std::string& key, bool optional = false)
  {
    if (table.count(key))
      return parseVector_(table.at(key));
    else
    {
      if (optional)
        return Eigen::VectorXd();
      else
        throw std::runtime_error("No element " + key + " found in the file");
    }
  }
}

namespace bms
{
  void RawProblem::read(const std::string& filepath)
  {
    std::ifstream aif(filepath);
    if (aif.is_open())
    {
      std::map<std::string, std::string> table;
      for (std::string line; std::getline(aif, line); )
      {
        auto i = line.find("=");
        if (i != std::string::npos && i>0)
        {
          std::string substr = line.substr(0, i);
          std::string token = trim(substr);
          auto j = line.find(";", i + 1);
          if (j != std::string::npos)
          {
            substr = line.substr(i + 1, j - i - 1);
            std::string value = trim(substr);
            table[token] = value;
          }
          else
            throw std::runtime_error("error in reading line\n" + line);
        }
      }

      g = parseDouble(table, "g");
      lmin = parseDouble(table, "lambda_min");
      lmax = parseDouble(table, "lambda_max");
      delta = parseVector(table, "Delta");
      wi_min = parseDouble(table, "omega_i_min");
      wi_max = parseDouble(table, "omega_i_max");
      zi = parseDouble(table, "z_bar");
      dzi = parseDouble(table, "zd_bar");
      zf = parseDouble(table, "z_f");
      Phi_ = parseVector(table, "Phi", true);
    }
    else
    {
      throw std::runtime_error("Unable to open" + filepath);
    }
  }

  Problem::Problem(const RawProblem& pb)
    : lso_(pb.delta)
    , lc_(pb.lmin*pb.delta, pb.lmax*pb.delta, pb.wi_min*pb.wi_min, pb.wi_max*pb.wi_max)
    , bc_(pb.delta, pb.zi / pb.g, pb.dzi / pb.g)
    , raw_(pb)
  {
    double d = pb.delta[0] * pb.g / pb.zf;
    lc_.changeBounds(0, d, d);
  }

  const LeastSquareObjective& Problem::objective() const
  {
    return lso_;
  }

  LeastSquareObjective& Problem::objective()
  {
    return lso_;
  }

  const BoundenessConstraint& Problem::nonLinearConstraint() const
  {
    return bc_;
  }

  BoundenessConstraint& Problem::nonLinearConstraint()
  {
    return bc_;
  }

  const LinearConstraints& Problem::linearConstraints() const
  {
    return lc_;
  }

  LinearConstraints& Problem::linearConstraints()
  {
    return lc_;
  }

  void Problem::set_zf(double zf)
  {
    raw_.zf = zf;
    computeAndSetBounds0();
  }

  void Problem::set_zi(double zi)
  {
    raw_.zi = zi;
    computeAndSetAlpha();
  }

  void Problem::set_dzi(double dzi)
  {
    raw_.dzi = dzi;
    computeAndSetb();
  }

  void Problem::set_lambda_min(double lmin)
  {
    raw_.lmin = lmin;
    computeAndSetZonotopeBounds();
  }

  void Problem::set_lambda_max(double lmax)
  {
    raw_.lmax = lmax;
    computeAndSetZonotopeBounds();
  }

  void Problem::set_lambdas(double lmin, double lmax)
  {
    raw_.lmin = lmin;
    raw_.lmax = lmax;
    computeAndSetZonotopeBounds();
  }

  void Problem::set_wi_min(double wi_min)
  {
    raw_.wi_min = wi_min;
    computeAndSetBoundsN();
  }

  void Problem::set_wi_max(double wi_max)
  {
    raw_.wi_max = wi_max;
    computeAndSetBoundsN();
  }

  void Problem::set_wi(double wi_min, double wi_max)
  {
    raw_.wi_min = wi_min;
    raw_.wi_max = wi_max;
    computeAndSetBoundsN();
  }

  void Problem::computeAndSetBounds0()
  {
    double d = raw_.delta[0] * raw_.g / raw_.zf;
    lc_.changeBounds(0, d, d);
  }

  void Problem::computeAndSetZonotopeBounds()
  {
    lc_.changeBounds(raw_.lmin*raw_.delta, raw_.lmax*raw_.delta);
    computeAndSetBounds0();
  }

  void Problem::computeAndSetBoundsN()
  {
    lc_.changeBounds(raw_.delta.size(), raw_.wi_min*raw_.wi_min, raw_.wi_max*raw_.wi_max);
  }

  void Problem::computeAndSetAlpha()
  {
    bc_.setAlpha(raw_.zi / raw_.g);
  }

  void Problem::computeAndSetb()
  {
    bc_.setb(raw_.dzi / raw_.g);
  }

}
