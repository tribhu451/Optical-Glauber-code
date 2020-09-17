#include<iostream>
#include<TMath.h>
#include<TRandom3.h>
#include<TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <fstream>
#include "Function.h"
#include "Eccentricity.h"

int main(){

//------------------------------INPUTS----------------------------------------
      int Event_No=10;
      double beta2=0.28;
      double beta4=0.093;
      double a = 0.55;
      double R = 6.81;
      double sigma=4.2;
      double A=238.0;
      double B=238.0;
      double n_pp=2.49;
      double X_hard=0.13;
//----------------------------------------------------------------------------





//--------------------------Random Seed Set-----------------------------------
      TRandom3* t1= new TRandom3();
      t1->SetSeed(0);
      TRandom3* t2= new TRandom3();
      t2->SetSeed(0);

      TRandom3* rsduse=new TRandom3();
      rsduse->SetSeed(0);
      long kft=rsduse->GetSeed();
      gRandom->SetSeed(kft);
      TF1* f1= new TF1("f1","TMath::Abs(TMath::Sin(x))",0.0,TMath::Pi());


      TRandom3* rsduse1=new TRandom3();
      rsduse1->SetSeed(0);
      long kft1=rsduse1->GetSeed();
      gRandom->SetSeed(kft1);     
      TF1* f2= new TF1("f2","x",0.0,22.0);    
//---------------------------------------------------------------------------

//-----------------------------FILE setup------------------------------------
       ofstream myfile;
       myfile.open("optical_Uranium.dat");
       TFile* File2= new TFile ("Optical_Uranium.root","recreate");
//---------------------------------------------------------------------------



//-----------------------------TREE setup------------------------------------
       double Impact;double E2; double Psi2; double E3; double Psi3;
       double E4; double Psi4; double E5; double Psi5;  double E6; double Psi6;
       double ota;double otb;double opa; double opb;double Nch;int Npart;int Ncoll;
       TTree* tree1= new TTree("Tree_opt_U","Tree for U+U at 200 Gev");
       tree1-> Branch("N_Part",&Npart,"NPart/I");
       tree1->Branch("N_Coll",&Ncoll,"NColl/I");
       tree1->Branch("NCh",&Nch,"Nch/D");
       tree1->Branch("Impact_parameter",&Impact,"b/D");
       tree1->Branch("Ori_TA",&ota,"ota/D");
       tree1->Branch("Ori_PA",&opa,"opa/D");
       tree1->Branch("Ori_TB",&otb,"otb/D");
       tree1->Branch("Ori_PB",&opb,"opb/D");
       tree1->Branch("E2",&E2,"E2/D");
       tree1->Branch("Psi_2",&Psi2,"Psi2/D");
       tree1->Branch("E3",&E3,"E3/D");
       tree1->Branch("Psi_3",&Psi3,"Psi3/D");
       tree1->Branch("E4",&E4,"E4/D");
       tree1->Branch("Psi_4",&Psi4,"Psi4/D");
       tree1->Branch("E5",&E5,"E5/D");
       tree1->Branch("Psi_5",&Psi5,"Psi5/D");
       tree1->Branch("E6",&E6,"E6/D");
       tree1->Branch("Psi_6",&Psi6,"Psi6/D");
//-----------------------------------------------------------------------------


      double Eccen[8],PPlane[8];

      Glauber h;
      Eccentricity E;

      h.Set_Inputs(beta2,beta4,a,R,sigma,A,B,n_pp,X_hard);
      E.Set_Inputs(beta2,beta4,a,R,sigma,A,B,n_pp,X_hard);

      int event=0;
      do {
      double ThetaA=f1->GetRandom(0.0,TMath::Pi());
      double ThetaB=f1->GetRandom(0.0,TMath::Pi());
      double PhiA=2*TMath::Pi()*t1->Rndm();
      double PhiB=2*TMath::Pi()*t2->Rndm();
      double b=f2->GetRandom(0.0,20.0);   
      h.GiveTA(ThetaA);
      h.GiveTB(ThetaB);
      h.GivePA(PhiA);
      h.GivePB(PhiB);
      h.Giveb(b);

      E.GiveTA(ThetaA);
      E.GiveTB(ThetaB);
      E.GivePA(PhiA);  
      E.GivePB(PhiB);  
      E.Give_Impact_Parameter(b);

      double Norm_constant=h.Get_Norm_Constant();
      h.Give_Norm_const(Norm_constant);
      E.Give_Norm_const(Norm_constant);

      double N_Coll=h.NColl();
      double N_Part=h.NPart();
      double Multiplicity= h.Multiplicity(N_Part,N_Coll);

      double XAVG=h.xavg();
      double YAVG=h.yavg();
      //cout<<"X average= "<<XAVG<<"   Y Average= "<<YAVG<<endl;
      E.Set_CM(XAVG,YAVG);
      //XAVG=h.xavg_();
      //YAVG=h.yavg_();
      //cout<<"X average= "<<XAVG<<"   Y Average= "<<YAVG<<endl;
      for (int hrmo=2; hrmo<7; hrmo=hrmo+1){
      E.Give_n(hrmo); 
      E.eccentricity();
      Eccen[hrmo]=E.Get_Eccentricity();
      PPlane[hrmo]= E.Get_Participant_Plane();
                                      }

      if(N_Part > 0.1 ){ 
      myfile<<b<<"    "<< ThetaA<<"   "<< PhiA <<"    "<<ThetaB<<"    "<<PhiB<<"    "<<N_Coll<<"   "<<N_Part<<"    "<<Multiplicity<<"   "<<Eccen[2]<<"   "<<PPlane[2]<<"   "<<Eccen[3]<<"   "<<PPlane[3]<<"   "<<Eccen[4]<<"   "<<PPlane[4]<<"   "<<Eccen[5]<<"   "<<PPlane[5]<<"   "<<Eccen[6]<<"   "<<PPlane[6]<<"\n";

 //myfile<<b<<"    "<< ThetaA<<"   "<< PhiA <<"    "<<ThetaB<<"    "<<PhiB<<"    "<<N_Coll<<"   "<<N_Part<<"    "<<Multiplicity<<"   "<<Eccen[2]<<"   "<<PPlane[2]<<"   "<<Eccen[4]<<"   "<<PPlane[4]<<"\n";

      Npart=N_Part;
      Ncoll=N_Coll;
      Nch=Multiplicity;
      Impact=b;
      ota=ThetaA;
      otb=ThetaB;
      opa=PhiA;
      opb=PhiB;
      E2=Eccen[2];
      Psi2=PPlane[2];
      E3=Eccen[3];
      Psi3=PPlane[3];
      E4=Eccen[4];
      Psi4=PPlane[4];
      E5=Eccen[5];
      Psi5=PPlane[5];
      E6=Eccen[6];
      Psi6=PPlane[6];
      tree1->Fill();
      event=event+1;
      }// IF loop end
      
      if((event/10)*10 == event){cout<<event<<"\n";}
      }while(event<Event_No);  // event loop      

      myfile.close();
      tree1->Write();
      File2->Close();                      
      return 0;

}
