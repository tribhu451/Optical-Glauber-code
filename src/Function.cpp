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
                        A2*(    ((35*TMath::Power(ZP,4))/(rsq*rsq)) - ((30*TMath::Power(ZP,2))/rsq) + 3    ) 
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
        return 1/val;

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
        double TA = f1->Integral(-3*R,3*R,1.0E-09);
        f1->SetParameter(0,x[0]-(b/2.0)); 
        f1->SetParameter(2,mThetaB); 
        f1->SetParameter(3,mPhiB); 
        double TB = f1->Integral(-3*R,3*R,1.0E-09);
        delete f1;
        return  TA*TB;

}

double Glauber::TA_TB2(double* x,double* p){


        TF1* f1;
        f1=new TF1 ("TA",this,&Glauber::Norm_Function,-3*R,3*R,4);
        f1->SetParameter(0,x[0]+(b/2.0)); 
        f1->SetParameter(1,x[1]); 
        f1->SetParameter(2,mThetaA); 
        f1->SetParameter(3,mPhiA); 
        double TA = f1->Integral(-3*R,3*R,1.0E-9);

        f1->SetParameter(0,x[0]-(b/2.0)); 
        f1->SetParameter(2,mThetaB); 
        f1->SetParameter(3,mPhiB); 
        double TB = f1->Integral(-3*R,3*R,1.0E-9);
        double va= ((A*TA*(1-(TMath::Power(1-TB*sigma,B)))) + (B*TB*(1-(TMath::Power(1-TA*sigma,A)))));
        delete f1;
        return  va;

}


double Glauber::energyxy(double* x,double* p){


        TF1* f1;
        f1=new TF1 ("TA",this,&Glauber::Norm_Function,-3*R,3*R,4);
        f1->SetParameter(0,x[0]+(b/2.0)); 
        f1->SetParameter(1,x[1]); 
        f1->SetParameter(2,mThetaA); 
        f1->SetParameter(3,mPhiA); 
        double TA = f1->Integral(-3*R,3*R,1.0E-9);

        f1->SetParameter(0,x[0]-(b/2.0)); 
        f1->SetParameter(2,mThetaB); 
        f1->SetParameter(3,mPhiB); 
        double TB = f1->Integral(-3*R,3*R,1.0E-9);
	
        double va=(((1-X_hard)*(n_pp/2.0))*(((A*TA*(1-(TMath::Power(1-TB*sigma,B)))) + (B*TB*(1-(TMath::Power(1-TA*sigma,A)))))))+(X_hard*n_pp*A*B*sigma*TA*TB);
        delete f1;
        return  va;

}

double Glauber::xavgfn(double* x,double* p){


        TF1* f1;
        f1=new TF1 ("TA",this,&Glauber::Norm_Function,-3*R,3*R,4);
        f1->SetParameter(0,x[0]+(b/2.0)); 
        f1->SetParameter(1,x[1]); 
        f1->SetParameter(2,mThetaA); 
        f1->SetParameter(3,mPhiA); 
        double TA = f1->Integral(-3*R,3*R,1.0E-9);

        f1->SetParameter(0,x[0]-(b/2.0)); 
        f1->SetParameter(2,mThetaB); 
        f1->SetParameter(3,mPhiB); 
        double TB = f1->Integral(-3*R,3*R,1.0E-9);
	
        double va=(x[0])*((((1-X_hard)*(n_pp/2.0))*(((A*TA*(1-(TMath::Power(1-TB*sigma,B)))) + (B*TB*(1-(TMath::Power(1-TA*sigma,A)))))))+(X_hard*n_pp*A*B*sigma*TA*TB));
        delete f1;
        return  va;

}

double Glauber::yavgfn(double* x,double* p){


        TF1* f1;
        f1=new TF1 ("TA",this,&Glauber::Norm_Function,-3*R,3*R,4);
        f1->SetParameter(0,x[0]+(b/2.0)); 
        f1->SetParameter(1,x[1]); 
        f1->SetParameter(2,mThetaA); 
        f1->SetParameter(3,mPhiA); 
        double TA = f1->Integral(-3*R,3*R,1.0E-9);

        f1->SetParameter(0,x[0]-(b/2.0)); 
        f1->SetParameter(2,mThetaB); 
        f1->SetParameter(3,mPhiB); 
        double TB = f1->Integral(-3*R,3*R,1.0E-9);
	
        double va=x[1]*((((1-X_hard)*(n_pp/2.0))*(((A*TA*(1-(TMath::Power(1-TB*sigma,B)))) + (B*TB*(1-(TMath::Power(1-TA*sigma,A)))))))+(X_hard*n_pp*A*B*sigma*TA*TB));
        delete f1;
        return  va;

}



double Glauber::xavg(){

        TF2* f2;
        f2= new TF2("TAB",this,&Glauber::xavgfn,-4*R,4*R,-4*R,4*R,0);
        double one = f2->Integral(-3*R,3*R,-3*R,3*R,1.0E-06);
	TF2* f3;
        f3= new TF2("Txc",this,&Glauber::energyxy,-4*R,4*R,-4*R,4*R,0);
        double two = f3->Integral(-3*R,3*R,-3*R,3*R,1.0E-06);
	double Xaverage=one/two;
        Xavg=Xaverage;
        delete f2;
	delete f3;
        return Xaverage;
}

double Glauber::yavg(){

        TF2* f2;
        f2= new TF2("TAB",this,&Glauber::yavgfn,-4*R,4*R,-4*R,4*R,0);
        double one = f2->Integral(-3*R,3*R,-3*R,3*R,1.0E-06);
	TF2* f3;
        f3= new TF2("Txc",this,&Glauber::energyxy,-4*R,4*R,-4*R,4*R,0);
        double two = f3->Integral(-3*R,3*R,-3*R,3*R,1.0E-06);
	double Yaverage=one/two;
        Yavg=Yaverage;
        delete f2;
	delete f3;
        return Yaverage;
}




double Glauber::NColl(){

        TF2* f2;
        f2= new TF2("TAB",this,&Glauber::TA_TB,-4*R,4*R,-4*R,4*R,0);
        double T_AB = f2->Integral(-3*R,3*R,-3*R,3*R,1.0E-06);
        double N_coll= A*B*sigma*(T_AB);
        delete f2;
        return N_coll;
}

double Glauber::NPart(){
        TF2 *f3;
        f3=new TF2("T_AB",this,&Glauber::TA_TB2,-4*R,4*R,-4*R,4*R,0);
        double TAAB=f3->Integral(-3*R,3*R,-3*R,3*R,1.0E-06);
        double N_Part=TAAB;
        delete f3;
        return N_Part;
}

double Glauber::Multiplicity(double npart,double ncoll){

        double N_ch=((1-X_hard)*n_pp*(npart/2))+(X_hard*n_pp*ncoll);
        return N_ch;
}

void Glauber::GiveTA(double aThetaA){  mThetaA=aThetaA;}
void Glauber::GiveTB(double aThetaB){  mThetaB=aThetaB;}
void Glauber::GivePA(double aPhiA)  {  mPhiA=aPhiA;    }
void Glauber::GivePB(double aPhiB)  {  mPhiB=aPhiB;    }
void Glauber::Giveb(double ab)       {  b=ab;           }
void Glauber::Give_Norm_const(double aNorm){Norm_const=aNorm;}
void Glauber::Set_Inputs(double abeta2,double abeta4,double aa, double aR, double asigma, double aA, double aB, double an_pp, double aX_hard){beta2= abeta2; beta4=abeta4;a=aa;R=aR;sigma=asigma;A=aA;B=aB;n_pp=an_pp;X_hard=aX_hard;}


//***************************************************





// cross checking the CM




double Glauber::energyxy_(double* x,double* p){


        TF1* f1;
        f1=new TF1 ("TA",this,&Glauber::Norm_Function,-3*R,3*R,4);
        f1->SetParameter(0,x[0]+(b/2.0)); 
        f1->SetParameter(1,x[1]); 
        f1->SetParameter(2,mThetaA); 
        f1->SetParameter(3,mPhiA); 
        double TA = f1->Integral(-3*R,3*R,1.0E-9);

        f1->SetParameter(0,x[0]-(b/2.0)); 
        f1->SetParameter(2,mThetaB); 
        f1->SetParameter(3,mPhiB); 
        double TB = f1->Integral(-3*R,3*R,1.0E-9);
	
        double va=(((1-X_hard)*(n_pp/2.0))*(((A*TA*(1-(TMath::Power(1-TB*sigma,B)))) + (B*TB*(1-(TMath::Power(1-TA*sigma,A)))))))+(X_hard*n_pp*A*B*sigma*TA*TB);
        delete f1;
        return  va;

}

double Glauber::xavgfn_(double* x,double* p){


        TF1* f1;
        f1=new TF1 ("TA",this,&Glauber::Norm_Function,-3*R,3*R,4);
        f1->SetParameter(0,x[0]+(b/2.0)); 
        f1->SetParameter(1,x[1]); 
        f1->SetParameter(2,mThetaA); 
        f1->SetParameter(3,mPhiA); 
        double TA = f1->Integral(-3*R,3*R,1.0E-9);

        f1->SetParameter(0,x[0]-(b/2.0)); 
        f1->SetParameter(2,mThetaB); 
        f1->SetParameter(3,mPhiB); 
        double TB = f1->Integral(-3*R,3*R,1.0E-9);
	
        double va=(x[0]-Xavg)*((((1-X_hard)*(n_pp/2.0))*(((A*TA*(1-(TMath::Power(1-TB*sigma,B)))) + (B*TB*(1-(TMath::Power(1-TA*sigma,A)))))))+(X_hard*n_pp*A*B*sigma*TA*TB));
        delete f1;
        return  va;

}

double Glauber::yavgfn_(double* x,double* p){


        TF1* f1;
        f1=new TF1 ("TA",this,&Glauber::Norm_Function,-3*R,3*R,4);
        f1->SetParameter(0,x[0]+(b/2.0)); 
        f1->SetParameter(1,x[1]); 
        f1->SetParameter(2,mThetaA); 
        f1->SetParameter(3,mPhiA); 
        double TA = f1->Integral(-3*R,3*R,1.0E-9);

        f1->SetParameter(0,x[0]-(b/2.0)); 
        f1->SetParameter(2,mThetaB); 
        f1->SetParameter(3,mPhiB); 
        double TB = f1->Integral(-3*R,3*R,1.0E-9);
	
        double va=(x[1]-Yavg)*((((1-X_hard)*(n_pp/2.0))*(((A*TA*(1-(TMath::Power(1-TB*sigma,B)))) + (B*TB*(1-(TMath::Power(1-TA*sigma,A)))))))+(X_hard*n_pp*A*B*sigma*TA*TB));
        delete f1;
        return  va;

}



double Glauber::xavg_(){

        TF2* f2;
        f2= new TF2("TAB",this,&Glauber::xavgfn_,-4*R,4*R,-4*R,4*R,0);
        double one = f2->Integral(-3*R,3*R,-3*R,3*R,1.0E-06);
	TF2* f3;
        f3= new TF2("Txc",this,&Glauber::energyxy_,-4*R,4*R,-4*R,4*R,0);
        double two = f3->Integral(-3*R,3*R,-3*R,3*R,1.0E-06);
	double Xaverage=one/two;
        delete f2;
	delete f3;
        return Xaverage;
}

double Glauber::yavg_(){

        TF2* f2;
        f2= new TF2("TAB",this,&Glauber::yavgfn_,-4*R,4*R,-4*R,4*R,0);
        double one = f2->Integral(-3*R,3*R,-3*R,3*R,1.0E-06);
	TF2* f3;
        f3= new TF2("Txc",this,&Glauber::energyxy_,-4*R,4*R,-4*R,4*R,0);
        double two = f3->Integral(-3*R,3*R,-3*R,3*R,1.0E-06);
	double Yaverage=one/two;
        delete f2;
	delete f3;
        return Yaverage;
}















