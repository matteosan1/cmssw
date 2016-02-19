#include "PhysicsTools/TagAndProbe/interface/RooCBExGaussShape.h" 

ClassImp(RooCBExGaussShape) 

RooCBExGaussShape::RooCBExGaussShape(const char *name, const char *title, 
				     RooAbsReal& _m,
				     RooAbsReal& _m0,
				     RooAbsReal& _sigma,
				     RooAbsReal& _alpha,
				     RooAbsReal& _n,
				     RooAbsReal& _sigma_2
				     ) :
RooAbsPdf(name,title), 
  m("m","m",this,_m),
  m0(" m0"," m0",this,_m0),
  sigma(" sigma"," sigma",this,_sigma),
  alpha(" alpha"," alpha",this,_alpha),
  n(" n"," n",this,_n),
  sigma_2(" sigma_2"," sigma_2",this,_sigma_2)
{}

RooCBExGaussShape::RooCBExGaussShape(const RooCBExGaussShape& other, const char* name):
  RooAbsPdf(other,name), 
  m("m",this,other.m),
  m0(" m0",this,other. m0),
  sigma(" sigma",this,other. sigma),
  alpha(" alpha",this,other. alpha),
  n(" n",this,other. n),
  sigma_2(" sigma_2",this,other. sigma_2)
{}


Double_t RooCBExGaussShape::evaluate() const 
{ 
  Double_t rval=0;

  Double_t t = (m-m0)/sigma;
  Double_t t0 = (m-m0)/sigma_2;
  if (alpha < 0){ 
    t = -t;
    t0 = -t0;
  }

  Double_t absAlpha = fabs((Double_t)alpha);

  if (t>0) {
    rval= exp(-0.5*t0*t0);
  }
  else if (t > -absAlpha) {
    rval= exp(-0.5*t*t);
  }
  else {
    Double_t a =  TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
    Double_t b= n/absAlpha - absAlpha; 
    rval= a/TMath::Power(b - t, n);
  }

  return rval;
} 
