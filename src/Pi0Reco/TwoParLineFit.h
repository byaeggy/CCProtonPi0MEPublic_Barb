#ifndef TwoParLineFit_h_seen
#define TwoParLineFit_h_seen

#include "OneParLineFit.h"

// Fitting function for line fit in two dimensions
class TwoParLineFit : public ROOT::Minuit2::FCNBase {
  public:
  TwoParLineFit(const std::vector<LineFit::Point>& data)
      : fData(data), useError_(false)  {}

    void UseError(bool v) { useError_ = v;}
    
    double operator() (const std::vector<double>& p) const {
        double k = p[0];
        double b = p[1];
        
        double chi2 = 0.0;
        double cov = 1.0;
        double total_weight = 0.0;
        for (std::vector<LineFit::Point>::const_iterator i = fData.begin();
             i != fData.end(); ++i) {
            const double    x_i = i->GetX();
            const double    y_i = i->GetY();
            const double yerr_i = i->GetYerr();
            const double w      = i->GetWeight();
            total_weight += w;

            double delta = y_i - f(k,b,x_i);

            if (useError_) cov = 1./(yerr_i*yerr_i);
            chi2 += delta*delta*cov*w;
        }

        chi2 /= total_weight;
        
        return chi2;

    }
        /* See FCNBase documentation */
    double Up() const { return 1.;}
    
private:
    double f(double k, double b, double x) const {
        return k*x + b;
    }

    std::vector<LineFit::Point> fData;
    bool useError_;
};


#endif
