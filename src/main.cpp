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
      int Event_No=2000;
      double beta2=0.0;
      double beta4=0.0;
      double a = 0.45;
      double R = 6.66;
      double sigma=6.4;
      double A=208.0;
      double B=208.0;
      double X_hard=0.14;
      double n_pp = 3.92; // no use
//----------------------------------------------------------------------------





//--------------------------Random Seed Set-----------------------------------
      TRandom3* t1= new TRandom3();
      t1->SetSeed(0);
      TRandom3* t2= new TRandom3();
      t2->SetSeed(0);


//-----------------------------FILE setup------------------------------------
       ofstream myfile;
       myfile.open("opt_glau_pb_pb.dat");
       TFile* File2= new TFile ("opt_glau_pb_pb.root","recreate");
//---------------------------------------------------------------------------



//-----------------------------TREE setup------------------------------------
       double Impact;
       double Nch;double Npart;double Ncoll;
       TTree* tree1= new TTree("tree_opt_pb","Tree for pb+pb at 2760 GeV");
       tree1-> Branch("N_Part",&Npart,"NPart/D");
       tree1->Branch("N_Coll",&Ncoll,"NColl/D");
       tree1->Branch("NCh",&Nch,"Nch/D");
       tree1->Branch("Impact_parameter",&Impact,"b/D");
 
       Glauber h;


       h.Set_Inputs(beta2,beta4,a,R,sigma,A,B,n_pp,X_hard);


      int event=0;
      do {

      double test = t1->Rndm();
      double b;
      b = t2->Rndm();
      if(test < b ) 
       { b *= 17.0 + 0.001; }
      else 
       {continue ; }

 
      h.GiveTA(0.0);
      h.GiveTB(0.0);
      h.GivePA(0.0);
      h.GivePB(0.0);
      h.Giveb(b);


      double Norm_constant=h.Get_Norm_Constant();
      h.Give_Norm_const(Norm_constant);


      double N_Coll=h.NColl();
      double N_Part=h.NPart();
      double Multiplicity= h.Multiplicity(N_Part,N_Coll);


      if(N_Part > 2.0 ){ 
      myfile << b << "\t" << N_Coll << "\t" << N_Part << "\t" << Multiplicity << endl;
      Npart=N_Part;
      Ncoll=N_Coll;
      Nch=Multiplicity;
      Impact=b;
      tree1->Fill();
      event=event+1;
      }// IF loop end
      
      cout<<event<<"\n";
      }while(event<Event_No);  // event loop      

      myfile.close();
      tree1->Write();
      File2->Close();                      
      return 0;

      delete t1;
      delete t2;
}






