// root -b -q "centrality_determin.C(\"au1.root\")"

#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TObject.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TROOT.h"

using namespace std;

void centrality_determin()
{
  
  
  // inputs [
  int bins=10000;
  double min_bin_val = 0.;
  double max_bin_val = 500.;
  string out_file_name = "pb_pb_min_bias_cent_info.txt" ;    
  
  double nch; double imp_b; double npart; double ncoll ;
  double max_b; double max_npart; double max_ncoll;
  double npart_should; double ncoll_should; double nch_should;
  
  double xcut[20]={0.0};
  
  const int pr_f = 15;
  double percent[pr_f] = {0, 5, 10 , 20, 30, 40, 50, 60, 70,80};
  //double percent[pr_f] = {0,  6,  15, 25 , 35 , 45 , 55 };
  
  TFile* f1=new TFile("Pb_Pb_2760_xhard_0.14_min_bias.root","read");
  TH1F* h1 = new TH1F("h1","dNchdEta_count",bins,min_bin_val,max_bin_val);
  TTree* TreeAu =new TTree();
  f1->GetObject("tree_opt_pb",TreeAu);
  TreeAu->SetBranchAddress("NCh",&nch);
  TreeAu->SetBranchAddress("Impact_parameter",&imp_b);
  TreeAu->SetBranchAddress("N_Part",&npart);
  TreeAu->SetBranchAddress("N_Coll",&ncoll);
   
  double maxch = 0 ; 
  int Entries=TreeAu->GetEntries();        
  for(int j=0;j<Entries;j++)
    {
       TreeAu->GetEntry(j);
       if(nch > max_bin_val) 
	 {cout<<"nch > maximum bin range : "<<nch<<endl; exit(1);}
       h1->Fill(nch);
       if (nch > maxch){ maxch = nch ; max_b = imp_b ; max_npart=npart; max_ncoll=ncoll;} 
       if(imp_b < 0.1 ){npart_should = npart; ncoll_should = ncoll; nch_should = nch;}
    }
  
  int n1=0;
  double Sum=0.0;
  int per_count = 0 ;
  xcut[n1] = maxch ;
 
  for(int j=bins; j>=1; j--)
    {
      double aSum = h1->GetBinContent(j);
      //if(aSum == 0){cout<<"Bin contetnt zero at "<<j<<"th bin"<<endl; }
      double perc=percent[per_count+1] - percent[per_count];   // perc=10 means 0-10% , 10-20%, 20-30%....                                     
      if(Sum>(Entries*perc*0.01))
	{
	  n1=n1+1; 
	  double p0=h1->GetBinCenter(j);                                       
	  xcut[n1]=p0;
	  Sum=0.0;
	  per_count += 1;
	  if(bins == 1){break;}  
	}
      else
	{
	  if(bins == 1){break;}  
	  Sum=Sum+aSum;
	}
    }
  
  
  for(int i=0 ; i <= n1 ; i++)
    {
      cout<<percent[i]<<"-"<<percent[i+1]<<" %\t"<<"MULT :\t"<<xcut[i]<<"-"<<xcut[i+1]<<endl;
    }
  
 
  



  ///////////////////////////////////////////
  // average npart , average b calculation //
  ///////////////////////////////////////////
 
  double tevents[pr_f]={0.0};
  double avg_npart[pr_f]={0.0};
  double avg_b[pr_f] = {0.0};
  double avg_nch[pr_f] = {0.0};
  
  
  double avg_npart2[pr_f]={0.0};
  double avg_b2[pr_f] = {0.0};
  double avg_nch2[pr_f] = {0.0};
  
  
  for(int j=0;j<Entries;j++)
    {
      TreeAu->GetEntry(j);
      for(int i=0 ; i <= n1 ; i++)
	{
	  if(nch > xcut[i+1] && nch < xcut[i]) 
	    {
	      avg_npart[i] += npart ; 
	      avg_npart2[i] += npart*npart ;
	      
	      avg_b[i] += imp_b; 
	      avg_b2[i] += imp_b*imp_b ;
	      
	      avg_nch[i] += nch ; 
	      avg_nch2[i] += nch*nch ; 
	      
	      tevents[i] += 1.0 ; 
	    } 
	}
      
    } // entries loop
  



  
  ////////////////////////////////////////
  // print in terminal for confirmation //
  ////////////////////////////////////////

  cout << "\n\n" <<endl;
  
  for(int i=0 ; i < n1 ; i++)
    {
      cout<<"avg_npart : "<< avg_npart[i]/tevents[i] << " +/- " << sqrt( avg_npart2[i]/tevents[i] - pow( avg_npart[i]/tevents[i] , 2 ) ) << endl;
    }  
  
  cout << "\n\n" <<endl;
  
  for(int i=0 ; i < n1 ; i++)
    {
      cout<<"avg_b : "<< avg_b[i]/tevents[i] << " +/- " << sqrt( avg_b2[i]/tevents[i] - pow( avg_b[i]/tevents[i] , 2 ) ) << endl;
    }  
  
  cout << "\n\n" <<endl;
  
  for(int i=0 ; i < n1 ; i++)
    {
      cout<<"avg_nch : "<< avg_nch[i]/tevents[i] << " +/- " << sqrt( avg_nch2[i]/tevents[i] - pow( avg_nch[i]/tevents[i] , 2 ) ) << endl;
    }  
  
  



  /////////////////////////////
  // take the output in file //
  /////////////////////////////

  std::ofstream file_out;
  file_out.open(out_file_name.c_str());
  file_out<<"#\t"<<"centrality\t"<<"avg_b\t"<<"error_b\t"<<"avg_npart\t"<<"error_npart\t"
          <<"avg_nch\t"<<"error_nch\t"<<"avg_nc_by_npart_pair"<<endl;
  file_out<<std::setprecision(6)<<std::scientific;
  for(int i=0 ; i < n1 ; i++)
    {
      
      file_out<<(percent[i+1]+percent[i])/2.0 << "\t"
	      << avg_b[i]/tevents[i]          << "\t"
	      << sqrt( avg_b2[i]/tevents[i] - pow( avg_b[i]/tevents[i] , 2 ) ) << "\t"
	      << avg_npart[i]/tevents[i] << "\t"
	      << sqrt( avg_npart2[i]/tevents[i] - pow( avg_npart[i]/tevents[i] , 2 ) ) << "\t"
	      << avg_nch[i]/tevents[i] << "\t"
	      << sqrt( avg_nch2[i]/tevents[i] - pow( avg_nch[i]/tevents[i] , 2 ) ) << "\t"
	      << avg_nch[i]/tevents[i] /(0.5*(avg_npart[i]/tevents[i])) << endl;
      
    }
  
  file_out.close();
  
  
  
  
}




