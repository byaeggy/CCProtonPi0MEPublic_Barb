#ifndef OneParLineFit_h
#define OneParLineFit_h

#include <vector>

#include <TFitterMinuit.h>

namespace LineFit {

    class Point {
      public:
        Point() {}
        Point(double x, double y, double yerr, double w) : fX(x), fY(y), fYerr(yerr), fWeight(w) {}
        
        double GetX() const { return fX; }
        double GetY() const { return fY; }
        double GetYerr() const { return fYerr; }
        double GetWeight() const { return fWeight; }
        
      private:
        double fX;
        double fY;
        double fYerr;
        double fWeight;
    };

}


/// Force the fitted line to go through a fixed point
class OneParLineFit : public ROOT::Minuit2::FCNBase {
public:
  OneParLineFit(const std::vector<LineFit::Point>& data)
      : fData(data), useError_(false)  {}

    double GetFCN() const { return fval_; }

    void UseError(bool v) { useError_ = v;}
    
    void SetFixedPoint(double x, double y) {
        x0_ = x;
        y0_ = y;
    }

    double operator() (const std::vector<double>& p) const {
        double k = p[0];

        double chi2 = 0.0;
        double cov = 1.0;
        for (std::vector<LineFit::Point>::const_iterator i = fData.begin();
             i != fData.end(); ++i) {
            const double    x_i = i->GetX();
            const double    y_i = i->GetY();
            const double yerr_i = i->GetYerr();
            const double w      = i->GetWeight(); 

            double delta = y_i - f(k,y0_,x_i-x0_);

            if (useError_) cov = 1./(yerr_i*yerr_i);
            chi2 += delta*delta*cov*w;
        }

        fval_ = chi2;
        
        return chi2;
    }
        /* See FCNBase documentation */
    double Up() const { return 1.;}
    
private:
    double f(double k, double b, double x) const {
        return k*x + b;
    }

    std::vector<LineFit::Point> fData;
    double x0_;
    double y0_;
    mutable double fval_;
    bool useError_;
};




#endif
