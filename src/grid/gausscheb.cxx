#include <grid/quadrature.hpp>

namespace ChronusQ {


  std::pair<double,double> unitBoundTransform(double lowBound, double upBound, double pt, double wgt) {

    // Both upper and lower bounds are finite: Map (-1,1) -> (a,b)
    if( upBound  != std::numeric_limits<double>::infinity() and 
        lowBound != -std::numeric_limits<double>::infinity() ) {
        
      // Factor jacobian into the weights
      // dx' = ((b-a)/2) * dx
      wgt *= (upBound - lowBound) / 2;

      // Perform the coordinate shift
      // x' = ((b-a) / 2)*x + (a+b)/2
      pt = (upBound - lowBound) * pt / 2 + (upBound + lowBound) / 2;

    // Lower bound is finite, upper is infinite: Map (-1,1) -> (a,\inf)
    } else if ( lowBound != -std::numeric_limits<double>::infinity() ) {

      // Factor jacobian into the weights
      // dx' = \frac{2}{(1-x)^2} dx
      wgt *= 2.0  / ( (1 - pt) * (1 - pt) );

      // Perform the coordinate shift
      // x' = a + \frac{(1+x)}{(1-x)}
      pt   = lowBound  + (1 + pt) / (1 - pt);

    // Upper bound is finite, low is infinite: Map (-1,1) -> (-\inf,b)
    } else if ( upBound != std::numeric_limits<double>::infinity() ) {
  
      // Factor jacobian into the weights
      // dx' = - \frac{2}{(1+x)^2} dx
      wgt *= 2.0  / ( (1 + pt) * (1 + pt) );

      // Perform the coordinate shift
      // x' = b - \frac{(1-x)}{(1+x)}
      pt   = upBound  - (1 - pt) / (1 + pt);
 
    // Both infinite: Map (-1,1) -> (-\inf,\inf)
    } else { 

      // Factor jacobian into the weights
      // dx' = \frac{1+x^2}{(x^2 - 1)^2} dx
      wgt *=  (1 + pt*pt) / ((pt*pt - 1)*(pt*pt - 1));

      // Perform the coordinate shift
      // x' = x / (1+x) / (1-x)
      pt = pt / (1 + pt) / (1 - pt);

    }

    return std::make_pair(pt,wgt);

  };

  void GaussChebFst::generateQuadrature() {

    Quadrature::generateQuadrature();

    std::cerr << "GaussChebFst" << std::endl;

    assert(upBound > lowBound);

    for(size_t i = 0; i < this->nPts; i++) {

      // Generate raw points and weights on (-1,1)
      double pt  = std::cos( (2.0*(i+1)-1.0) / (2*this->nPts) * M_PI );
      double wgt = (M_PI / this->nPts);

      // As we're integrating f(x) not \frac{f(x)}{\sqrt{1 - x^2}}
      // we must factor the \sqrt{1-x^2} into the weights

      wgt *= std::sqrt(1 - pt * pt);



      // Transform points and populate arrays
      std::tie(pts[i],weights[i]) = unitBoundTransform(lowBound,upBound,pt,wgt);

    }; // Loop over points

  };



  void GaussChebSnd::generateQuadrature() {

    Quadrature::generateQuadrature();

    assert(upBound > lowBound);

    for(size_t i = 0; i < this->nPts; i++) {

      // Generate raw points and weights on (-1,1)
      double pt = std::cos(M_PI * (i+1) / (this->nPts +1));
      double wgt = M_PI / (this->nPts + 1) * 
        std::pow(
            std::sin(M_PI * (i+1) / (this->nPts +1)),
            2.0
        );

      // As we're integrating f(x) not f(x)\sqrt{1 - x^2}
      // we must factor the \sqrt{1-x^2} into the weights

      wgt /= std::sqrt(1 - pt * pt);



      // Transform points and populate arrays
      std::tie(pts[i],weights[i]) = unitBoundTransform(lowBound,upBound,pt,wgt);

    }; // Loop over points

  };

};
