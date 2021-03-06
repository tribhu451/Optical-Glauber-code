#include  "Function.h"

using namespace std;


double Glauber::Modf_WoodSaxon(const double* x){

        double theta=0;
        double phi=0;

        double A1 = (beta2/4.0)*(TMath::Sqrt(5.0/pi));
        double A2 = ((3.0*beta4)/16.0)*(TMath::Sqrt(1.0/pi));
        double ZP = (- ( (TMath::Sin(theta))*x[0] ) -  ((TMath::Cos(theta))*(TMath::Sin(phi))*x[1])  +  
                                           ((TMath::Cos(theta))*(TMath::Cos(phi))*x[2]));            
        double rsq = ((x[0]*x[0])+(x[1]*x[1])+(x[2]*x[2])) ;
        double RP = R* (        1+
                        (A1*(((3*ZP*ZP)/rsq)-1))+
                        A2*(    ((35*TMath::Power(ZP,4))/(rsq*rsq)) - ((30*TMath::Power(ZP,2))/rsq) + 3 ) 
               );


        double r = TMath::Sqrt(rsq);
        return 1/(TMath::Exp((r - RP)/a)+1);
}



double Glauber::Get_Norm_Constant(){


        ROOT::Math::Functor wf(this, &Glauber::Modf_WoodSaxon,3);
        double min[3] = {-3*R,-3*R,-3*R};
        double max[3] = {3*R,3*R,4*R};
        ROOT::Math::IntegratorMultiDim pd;
        pd.SetFunction(wf);
        double val = pd.Integral(min,max);
        return 1.0/val;

}




double Glauber::Norm_Function(double* x,double* p){

        double A1 = (beta2/4.0)*(TMath::Sqrt(5.0/pi));
        double A2 = ((3.0*beta4)/16.0)*(TMath::Sqrt(1.0/pi));
        double ZP = (- ( (TMath::Sin(p[2]))*(TMath::Cos(p[3]))*p[0] ) -  ((TMath::Sin(p[3]))*p[1])  +  
                                            ((TMath::Cos(p[2]))*(TMath::Cos(p[3]))*x[0]));            
        double rsq = ((x[0]*x[0])+p[0]*p[0]+p[1]*p[1]);
        double RP = R* ( 1+
                        (A1*(((3*ZP*ZP)/rsq)-1))+
                        A2*(    ((35*TMath::Power(ZP,4))/(rsq*rsq)) - ((30*TMath::Power(ZP,2))/rsq) + 3 ) 

                            );
        double r = TMath::Sqrt(rsq);
        return (Norm_const)/(TMath::Exp((r - RP)/a)+1);

}


double Glauber::TA_TB(double* x,double* p){


        TF1* f1;
        f1=new TF1 ("TA",this,&Glauber::Norm_Function,-3*R,3*R,4);
        f1->SetParameter(0,x[0]+(b/2.0)); 
        f1->SetParameter(1,x[1]); 
        f1->SetParameter(2,mThetaA); 
        f1->SetParameter(3,mPhiA); 
        double TA = f1->Integral(-3*R,3*R);
        f1->SetParameter(0,x[0]-(b/2.0)); 
        f1->SetParameter(2,mThetaB); 
        f1->SetParameter(3,mPhiB); 
        double TB = f1->Integral(-3*R,3*R);
        delete f1;
        return  A*B*sigma*TA*TB;

}

double Glauber::TA_TB2(double* x,double* p){


        TF1* f1;
        f1=new TF1 ("TA",this,&Glauber::Norm_Function,-3*R,3*R,4);
        f1->SetParameter(0,x[0]+(b/2.0)); 
        f1->SetParameter(1,x[1]); 
        f1->SetParameter(2,mThetaA); 
        f1->SetParameter(3,mPhiA); 
        double TA = f1->Integral(-3*R,3*R);

        f1->SetParameter(0,x[0]-(b/2.0)); 
        f1->SetParameter(2,mThetaB); 
        f1->SetParameter(3,mPhiB); 
        double TB = f1->Integral(-3*R,3*R);
        double va= ((A*TA*(1-(TMath::Power(1-TB*sigma,B)))) + (B*TB*(1-(TMath::Power(1-TA*sigma,A)))));
        delete f1;
        return  va;

}




double Glauber::NColl(){

        TF2* f2;
        f2= new TF2("TAB",this,&Glauber::TA_TB,-3*R,3*R,-3*R,3*R,0);
        double T_AB = f2->Integral(-3*R,3*R,-3*R,3*R);
        double N_coll= (T_AB);
        delete f2;
        return N_coll;
}

double Glauber::NPart(){
        TF2 *f3;
        f3=new TF2("T_AB",this,&Glauber::TA_TB2,-3*R,3*R,-3*R,3*R,0);
        double TAAB=f3->Integral(-3*R,3*R,-3*R,3*R);
        double N_Part=TAAB;
        delete f3;
        return N_Part;
}

double Glauber::Multiplicity(double npart,double ncoll){

        double N_ch=( (1-X_hard) * ( npart / 2.0 ) ) + ( X_hard * ncoll );
        return N_ch;
}

void Glauber::GiveTA(double aThetaA){  mThetaA=aThetaA;}
void Glauber::GiveTB(double aThetaB){  mThetaB=aThetaB;}
void Glauber::GivePA(double aPhiA)  {  mPhiA=aPhiA;    }
void Glauber::GivePB(double aPhiB)  {  mPhiB=aPhiB;    }
void Glauber::Giveb(double ab)       {  b=ab;           }
void Glauber::Give_Norm_const(double aNorm){Norm_const=aNorm;}
void Glauber::Set_Inputs(double abeta2,double abeta4,double aa, double aR, double asigma, double aA, double aB, double an_pp, double aX_hard){beta2= abeta2; beta4=abeta4;a=aa;R=aR;sigma=asigma;A=aA;B=aB;n_pp=an_pp;X_hard=aX_hard;}



