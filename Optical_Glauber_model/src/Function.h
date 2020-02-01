# include<iostream>
#include<TF1.h>
#include<TF2.h>
#include<TMath.h>
#include "Math/Functor.h"
#include "Math/Integrator.h"
#include "Math/IntegratorMultiDim.h"

using namespace std;

class Glauber{

        private:

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
        
        double mThetaA;
        double mPhiA;
        double mThetaB;
        double mPhiB;
        double b;
        double Norm_const;

        double Modf_WoodSaxon(const double *x);
        double Norm_Function(double* x,double* p);
        double TA_TB(double* x,double* p);
        double TA_TB2(double* x,double* p);
        double TA,TB;

        public:
       
        double pi=TMath::Pi();
        double NColl();
        double NPart();
        double Multiplicity(double npart,double ncoll);
        void GiveTA(double aThetaA);
        void GiveTB(double aThetaA);
        void GivePA(double aThetaA);
        void GivePB(double aThetaA);
        void Giveb(double b);
        void Give_Norm_const(double aNorm);
        double Get_Norm_Constant();
        void Set_Inputs(double ,double ,double , double , double , double , double , double , double );
};

