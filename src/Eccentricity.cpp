#include "Eccentricity.h"


using namespace std;

double Eccentricity::normalised_function(double* x,double* p){
           
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
        return (Norm_const)/(TMath::Exp((r-RP)/a)+1);

}

double Eccentricity::Integrand1(double* x,double* p){

        TF1* f1;
        f1=new TF1 ("TA",this,&Eccentricity::normalised_function,-3*R,3*R,4);
        f1->SetParameter(0,x[0]+(b/2.0)-xavg); 
        f1->SetParameter(1,x[1]-yavg); 
        f1->SetParameter(2,mThetaA); 
        f1->SetParameter(3,mPhiA); 
        double T_A = f1->Integral(-3*R,2.5*R,1.0E-10);
        f1->SetParameter(0,x[0]-(b/2.0)-xavg); 
        f1->SetParameter(2,mThetaB); 
        f1->SetParameter(3,mPhiB); 
        double T_B = f1->Integral(-3*R,2.5*R,1.0E-10);
        double N_coll=A*B*sigma*T_A*T_B;
        double N_Part= ((A*T_A*(1-(TMath::Power(1-(T_B*sigma),B)))) + (B*T_B*(1-(TMath::Power(1-(T_A*sigma),A)))));
        double N_ch=((1-X_hard)*n_pp*(N_Part/2))+(X_hard*n_pp*N_coll);
        delete f1;
        return TMath::Power( x[0]*x[0]+x[1]*x[1] ,(n)/2.0) * N_ch;


}

double Eccentricity::Integrand2(double* x,double* p){
        TF1* f2;
        f2=new TF1 ("TA",this,&Eccentricity::normalised_function,-3*R,3*R,4);
        f2->SetParameter(0,x[0]+(b/2.0)-xavg); 
        f2->SetParameter(1,x[1]-yavg); 
        f2->SetParameter(2,mThetaA); 
        f2->SetParameter(3,mPhiA); 
        double T_A = f2->Integral(-3*R,2.5*R,1.0E-10);
        f2->SetParameter(0,x[0]-(b/2.0)-xavg); 
        f2->SetParameter(2,mThetaB); 
        f2->SetParameter(3,mPhiB); 
        double T_B = f2->Integral(-3*R,2.5*R,1.0E-10);
        double N_coll=A*B*sigma*T_A*T_B;         
        double N_Part= ((A*T_A*(1-(TMath::Power(1-(T_B*sigma),B)))) + (B*T_B*(1-(TMath::Power(1-(T_A*sigma),A)))));
        double N_ch=((1-X_hard)*n_pp*(N_Part/2))+(X_hard*n_pp*N_coll);
        delete f2;
        return TMath::Power( x[0]*x[0]+x[1]*x[1] ,(n)/2.0) * TMath::Cos(n*(TMath::ATan2(x[1],x[0]) ) ) * N_ch;

}

double Eccentricity::Integrand3(double* x,double* p){
        TF1* f2;
        f2=new TF1 ("TA",this,&Eccentricity::normalised_function,-3*R,3*R,4);
        f2->SetParameter(0,x[0]+(b/2.0)-xavg); 
        f2->SetParameter(1,x[1]-yavg); 
        f2->SetParameter(2,mThetaA); 
        f2->SetParameter(3,mPhiA); 
        double T_A = f2->Integral(-3*R,2.5*R,1.0E-10);
        f2->SetParameter(0,x[0]-(b/2.0)-xavg); 
        f2->SetParameter(2,mThetaB); 
        f2->SetParameter(3,mPhiB); 
        double T_B = f2->Integral(-3*R,2.5*R,1.0E-10);
        double N_coll=A*B*sigma*T_A*T_B;         
        double N_Part= ((A*T_A*(1-(TMath::Power(1-(T_B*sigma),B)))) + (B*T_B*(1-(TMath::Power(1-(T_A*sigma),A)))));
        double N_ch=((1-X_hard)*n_pp*(N_Part/2))+(X_hard*n_pp*N_coll);
        delete f2;
        return TMath::Power( x[0]*x[0]+x[1]*x[1] ,(n)/2.0) * TMath::Sin(n*(TMath::ATan2(x[1],x[0]) ) ) * N_ch;

}

void Eccentricity::eccentricity(){
        TF2 *f3;
        f3=new TF2("1st",this,&Eccentricity::Integrand2,-4*R,4*R,-4*R,4*R,0);
        double Numerator1=f3->Integral(-3*R,3.1*R,-3*R,3.1*R,1.0E-03);

        TF2 *f5;
        f5=new TF2("13rd",this,&Eccentricity::Integrand3,-4*R,4*R,-4*R,4*R,0);
        double Numerator2=f5->Integral(-3*R,3.1*R,-3*R,3.1*R,1.0E-03);



        TF2 *f4;
        f4=new TF2("2nd",this,&Eccentricity::Integrand1,-4*R,3*R,-4*R,4*R,0);
        double Denominator=f4->Integral(-3*R,3.1*R,-3*R,3.1*R,1.0E-03);

        double Ans1=-(Numerator1/Denominator);
        double Ans2=-(Numerator2/Denominator);


        Eccentricityp = TMath::Sqrt((Ans2*Ans2)+(Ans1*Ans1));
        PartPlane = (1.0/n)*TMath::ATan2(Ans2,Ans1);

          
        delete f3;
        delete f4;
        delete f5;


}

void Eccentricity::Give_Impact_Parameter(double anum) {b=anum;}
void Eccentricity::Give_Norm_const(double aconst) {Norm_const=aconst;}
void Eccentricity::GiveTA(double aThetaA){  mThetaA=aThetaA;}
void Eccentricity::GiveTB(double aThetaB){  mThetaB=aThetaB;}
void Eccentricity::GivePA(double aPhiA)  {  mPhiA=aPhiA;    }
void Eccentricity::GivePB(double aPhiB)  {  mPhiB=aPhiB;    }
void Eccentricity::Give_n(int an)  {  n=an;    }
double Eccentricity::Get_Eccentricity(){return Eccentricityp;} 
double Eccentricity::Get_Participant_Plane(){return PartPlane;} 
void Eccentricity::Set_Inputs(double abeta2,double abeta4,double aa, double aR, double asigma, double aA, double aB, double an_pp, double aX_hard){beta2= abeta2; beta4=abeta4;a=aa;R=aR;sigma=asigma;A=aA;B=aB;n_pp=an_pp;X_hard=aX_hard;}
void Eccentricity::Set_CM(double axavg, double ayavg) {xavg=axavg;yavg=ayavg;}

