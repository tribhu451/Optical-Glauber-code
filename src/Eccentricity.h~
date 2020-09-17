#include <iostream>
#include <TF1.h>
#include <TF2.h>
#include <TMath.h>

using namespace std;

class Eccentricity{

public:
      void Give_Impact_Parameter(double anum);
      void Give_Norm_const(double aconst);
      void GiveTA(double aThetaA);
      void GiveTB(double aThetaB);
      void GivePA(double aPhiA)  ;
      void GivePB(double aPhiB)  ;
      void eccentricity();
      void Give_n (int an);
      double Eccentricityp;
      double PartPlane;
      double Get_Eccentricity();
      double Get_Participant_Plane();
      void Set_Inputs(double ,double ,double , double , double , double , double , double , double );

private:
      double pi=TMath::Pi();
      int n;
      double mThetaA;
      double mPhiA;
      double mThetaB;
      double mPhiB;
      double b;
      double Norm_const;

      double normalised_function(double* x,double* p);
      double Integrand1(double* x,double* p);
      double Integrand2(double* x,double* p);
      double Integrand3(double* x,double* p);

        double beta2;
        double beta4;
        double rsq;
        double a;
        double R;
        double sigma;
        double A;
        double B;
        double n_pp;
        double X_hard;
      


};
