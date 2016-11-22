#include <iostream>
#include <fstream>
#include <string>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TGaxis.h>

#include "multiplicity.h"

#ifdef LXPLUS
#define hadron_mult_pdf "/afs/cern.ch/user/n/nipierre/workspace/djangoh/hadron_mult.pdf"
#define pion_mult_pdf "/afs/cern.ch/user/n/nipierre/workspace/djangoh/pion_mult.pdf"
#define kaon_mult_pdf "/afs/cern.ch/user/n/nipierre/workspace/djangoh/kaon_mult.pdf"
#endif /*LXPLUS*/
#ifdef CCAGE
#define hadron_mult_pdf "/sps/compass/npierre/djangoh/out_mult/hadron_mult.pdf"
#define pion_mult_pdf "/sps/compass/npierre/djangoh/out_mult/pion_mult.pdf"
#define kaon_mult_pdf "/sps/compass/npierre/djangoh/out_mult/kaon_mult.pdf"
#endif /*CCAGE*/

using namespace std;

int main()
{
  ifstream in("luevents.dat");
  string dummys;
  double dummyd;
  int dummyi;
  vector<Events> Events_Tree;
  Events Events_storage;
#ifdef DEBUG
  int count = 1;
#endif /*DEBUG*/

#ifdef DEBUG
  cout << endl;
  cout << count << "\t";
#endif /*DEBUG*/

  while(in>>dummys)
  {

    Events Events_storage;
#ifdef DEBUG
    cout << dummys << "\t";
#endif /*DEBUG*/
    for(int i=0;i<2;i++)
    {
      in >> dummys;
#ifdef DEBUG
      cout << dummys << "\t";
#endif /*DEBUG*/
    }
#ifdef DEBUG
    cout << endl;
#endif /*DEBUG*/
    for(int i=0;i<10;i++)
    {
      in >> dummys;
#ifdef DEBUG
      cout << dummys << "\t";
#endif /*DEBUG*/
    }

    in >> dummys;
#ifdef DEBUG
    cout << endl;
    cout << dummys << "\t";
#endif /*DEBUG*/

    Events_storage.I.push_back(atoi(dummys.c_str()));

    while(dummys.compare("sum:"))
    {
      in >> dummys;
#ifdef DEBUG
      cout << dummys << "\t";
#endif /*DEBUG*/
      Events_storage.particle.push_back(dummys);

      if(!dummys.compare("u"))
      {
        in >> dummys;
#ifdef DEBUG
        cout << dummys << "\t";
#endif /*DEBUG*/
      }
      else if(!dummys.compare("diquark"))
      {
        in >> dummys;
#ifdef DEBUG
        cout << dummys << "\t";
#endif /*DEBUG*/
      }

      for(int i=0;i<3;i++)
      {
        in >> dummyi;
#ifdef DEBUG
        cout << dummyi << "\t";
#endif /*DEBUG*/
        if(i==0)
          Events_storage.KS.push_back(dummyi);
        else if(i==1)
          Events_storage.KF.push_back(dummyi);
        else
          Events_storage.orig.push_back(dummyi);
      }
      for(int i=0;i<5;i++)
      {
        in >> dummyd;
#ifdef DEBUG
        cout << dummyd << "\t";
#endif /*DEBUG*/
        if(i==0)
          Events_storage.p_x.push_back(dummyd);
        else if(i==1)
          Events_storage.p_y.push_back(dummyd);
        else if(i==2)
          Events_storage.p_z.push_back(dummyd);
        else if(i==3)
          Events_storage.E.push_back(dummyd);
        else
          Events_storage.m.push_back(dummyd);
      }

      in >> dummys;
#ifdef DEBUG
      cout << endl;
      cout << dummys << "\t";
#endif /*DEBUG*/
    }

    in >> dummyd;
#ifdef DEBUG
    cout << dummyd << "\t";
#endif /*DEBUG*/

    for(int i=0;i<5;i++)
    {
      in >> dummyd;
#ifdef DEBUG
      cout << dummyd << "\t";
#endif /*DEBUG*/
    }

    in >> dummyi;
#ifdef DEBUG
    cout << endl;
    cout << dummyi << "\t" << "xBj: ";
#endif /*DEBUG*/

    for(int i=0;i<3;i++)
    {
      in >> dummyd;
#ifdef DEBUG
      cout << dummyd << "\t";
#endif /*DEBUG*/
      if(i==0)
      {
        Events_storage.xBj = dummyd;
#ifdef DEBUG
        cout << "y: ";
#endif /*DEBUG*/
      }
      else if(i==1)
      {
        Events_storage.y = dummyd;
#ifdef DEBUG
        cout << "Q2: ";
#endif /*DEBUG*/
      }
      else
        Events_storage.Q2 = dummyd;
    }

#ifdef DEBUG
    cout << "\n\n" << endl;
#endif /*DEBUG*/

    Events_Tree.push_back(Events_storage);

  }

  int hp[9][5][12], hm[9][5][12], pip[9][5][12], pim[9][5][12],
      kp[9][5][12], km[9][5][12], pp[9][5][12], pm[9][5][12],
      DIS[9][5];
  int xbin, ybin, zbin;
  double nu, z;

  for(int i=0;i<9;i++)
  {
    for(int j=0;j<5;j++)
    {
      DIS[i][j] = 0;
      for(int k=0;k<12;k++)
      {
        hp[i][j][k] = 0;
        hm[i][j][k] = 0;
        pip[i][j][k] = 0;
        pim[i][j][k] = 0;
        kp[i][j][k] = 0;
        km[i][j][k] = 0;
        pp[i][j][k] = 0;
        pm[i][j][k] = 0;
      }
    }
  }

  for(int i=0;i<int(Events_Tree.size());i++)
  {
    if(!(0.004<Events_Tree[i].xBj && Events_Tree[i].xBj<0.4)) continue;
    if(!(0.1<Events_Tree[i].y && Events_Tree[i].y<0.7)) continue;

    nu = Events_Tree[i].y*160;

    if(0.004<Events_Tree[i].xBj && Events_Tree[i].xBj<0.01) xbin = 0;
    else if(0.01<=Events_Tree[i].xBj && Events_Tree[i].xBj<0.02) xbin = 1;
    else if(0.02<=Events_Tree[i].xBj && Events_Tree[i].xBj<0.03) xbin = 2;
    else if(0.03<=Events_Tree[i].xBj && Events_Tree[i].xBj<0.04) xbin = 3;
    else if(0.04<=Events_Tree[i].xBj && Events_Tree[i].xBj<0.06) xbin = 4;
    else if(0.06<=Events_Tree[i].xBj && Events_Tree[i].xBj<0.1) xbin = 5;
    else if(0.1<=Events_Tree[i].xBj && Events_Tree[i].xBj<0.14) xbin = 6;
    else if(0.14<=Events_Tree[i].xBj && Events_Tree[i].xBj<0.18) xbin = 7;
    else xbin = 8;

    if(0.1<Events_Tree[i].y && Events_Tree[i].y<0.15) ybin = 0;
    else if(0.15<=Events_Tree[i].y && Events_Tree[i].y<0.2) ybin = 1;
    else if(0.2<=Events_Tree[i].y && Events_Tree[i].y<0.3) ybin = 2;
    else if(0.3<=Events_Tree[i].y && Events_Tree[i].y<0.5) ybin = 3;
    else ybin = 4;

    DIS[xbin][ybin]++;

    for(int j=0;j<int(Events_Tree[i].KF.size());j++)
    {
#ifdef DEBUG
      cout << j << " " << Events_Tree[i].KF[j] << endl;
#endif /*DEBUG*/
      if(abs(Events_Tree[i].KF[j]) == 211)
      {

        if(nu>0) z = sqrt(pow(Events_Tree[i].p_x[j],2)+pow(Events_Tree[i].p_y[j],2)+pow(Events_Tree[i].p_z[j],2)+pow(fM_pi,2))/nu;
        else z = 0;
#ifdef DEBUG
        cout << z << endl;
#endif /*DEBUG*/
        if(!(0.2<z && z<0.85)) continue;

        if(0.2<z && z<0.25) zbin = 0;
        else if(0.25<=z && z<0.30) zbin = 1;
        else if(0.30<=z && z<0.35) zbin = 2;
        else if(0.35<=z && z<0.40) zbin = 3;
        else if(0.40<=z && z<0.45) zbin = 4;
        else if(0.45<=z && z<0.50) zbin = 5;
        else if(0.50<=z && z<0.55) zbin = 6;
        else if(0.55<=z && z<0.60) zbin = 7;
        else if(0.60<=z && z<0.65) zbin = 8;
        else if(0.65<=z && z<0.70) zbin = 9;
        else if(0.70<=z && z<0.75) zbin = 10;
        else zbin = 11;

        if(Events_Tree[i].KF[j]>0)
        {
          hp[xbin][ybin][zbin]++;
          pip[xbin][ybin][zbin]++;
        }
        else
        {
          hm[xbin][ybin][zbin]++;
          pim[xbin][ybin][zbin]++;
        }
      }
      else if(abs(Events_Tree[i].KF[j]) == 321)
      {

        if(nu>0) z = sqrt(pow(Events_Tree[i].p_x[j],2)+pow(Events_Tree[i].p_y[j],2)+pow(Events_Tree[i].p_z[j],2)+pow(fM_K,2))/nu;
        else z = 0;
#ifdef DEBUG
        cout << z << endl;
#endif /*DEBUG*/
        if(!(0.2<z && z<0.85)) continue;

        if(0.2<z && z<0.25) zbin = 0;
        else if(0.25<=z && z<0.30) zbin = 1;
        else if(0.30<=z && z<0.35) zbin = 2;
        else if(0.35<=z && z<0.40) zbin = 3;
        else if(0.40<=z && z<0.45) zbin = 4;
        else if(0.45<=z && z<0.50) zbin = 5;
        else if(0.50<=z && z<0.55) zbin = 6;
        else if(0.55<=z && z<0.60) zbin = 7;
        else if(0.60<=z && z<0.65) zbin = 8;
        else if(0.65<=z && z<0.70) zbin = 9;
        else if(0.70<=z && z<0.75) zbin = 10;
        else zbin = 11;

        if(Events_Tree[i].KF[j]>0)
        {
          hp[xbin][ybin][zbin]++;
          kp[xbin][ybin][zbin]++;
        }
        else
        {
          hm[xbin][ybin][zbin]++;
          km[xbin][ybin][zbin]++;
        }
      }
      else if(abs(Events_Tree[i].KF[j]) == 2212)
      {

        if(nu>0) z = sqrt(pow(Events_Tree[i].p_x[j],2)+pow(Events_Tree[i].p_y[j],2)+pow(Events_Tree[i].p_z[j],2)+pow(fM_p,2))/nu;
        else z = 0;
#ifdef DEBUG
        cout << z << endl;
#endif /*DEBUG*/
        if(!(0.2<z && z<0.85)) continue;

        if(0.2<z && z<0.25) zbin = 0;
        else if(0.25<=z && z<0.30) zbin = 1;
        else if(0.30<=z && z<0.35) zbin = 2;
        else if(0.35<=z && z<0.40) zbin = 3;
        else if(0.40<=z && z<0.45) zbin = 4;
        else if(0.45<=z && z<0.50) zbin = 5;
        else if(0.50<=z && z<0.55) zbin = 6;
        else if(0.55<=z && z<0.60) zbin = 7;
        else if(0.60<=z && z<0.65) zbin = 8;
        else if(0.65<=z && z<0.70) zbin = 9;
        else if(0.70<=z && z<0.75) zbin = 10;
        else zbin = 11;

        if(Events_Tree[i].KF[j]>0)
        {
          hp[xbin][ybin][zbin]++;
          pp[xbin][ybin][zbin]++;
        }
        else
        {
          hm[xbin][ybin][zbin]++;
          pm[xbin][ybin][zbin]++;
        }
      }
      else if(Events_Tree[i].KF[j] == 130)
      {

        if(nu>0) z = sqrt(pow(Events_Tree[i].p_x[j],2)+pow(Events_Tree[i].p_y[j],2)+pow(Events_Tree[i].p_z[j],2)+pow(fM_pi,2))/nu;
        else z = 0;
#ifdef DEBUG
        cout << z << " " << endl;
#endif /*DEBUG*/
        if(!(0.2<z && z<0.85)) continue;

        if(0.2<z && z<0.25) zbin = 0;
        else if(0.25<=z && z<0.30) zbin = 1;
        else if(0.30<=z && z<0.35) zbin = 2;
        else if(0.35<=z && z<0.40) zbin = 3;
        else if(0.40<=z && z<0.45) zbin = 4;
        else if(0.45<=z && z<0.50) zbin = 5;
        else if(0.50<=z && z<0.55) zbin = 6;
        else if(0.55<=z && z<0.60) zbin = 7;
        else if(0.60<=z && z<0.65) zbin = 8;
        else if(0.65<=z && z<0.70) zbin = 9;
        else if(0.70<=z && z<0.75) zbin = 10;
        else zbin = 11;

        if(Events_Tree[i].KF[j]>0)
        {
          hp[xbin][ybin][zbin]++;
        }
        else
        {
          hm[xbin][ybin][zbin]++;
        }
      }
    }
  }

  TCanvas c1("Hadron_Multiplicity","Hadron_Multiplicity",3200,1600);
  TCanvas c2("Pion_Multiplicity","Pion_Multiplicity",3200,1600);
  TCanvas c3("Kaon_Multiplicity","Kaon_Multiplicity",3200,1600);

  c1.Range(-10,-1,10,1);
  c1.SetLeftMargin(0.5);
  c1.SetRightMargin(0.15);
  c1.SetTopMargin(0.5);
  c1.SetBottomMargin(0.15);

  c1.SetFillColor(0);
  c2.SetFillColor(0);
  c3.SetFillColor(0);

  c1.Divide(9,5,0,0);
  c2.Divide(9,5,0,0);
  c3.Divide(9,5,0,0);

  TGraph* H[2][9][5];
  TGraph* Pi[2][9][5];
  TGraph* K[2][9][5];

  double z_range[12] = {.225,.275,.325,.375,.425,.475,.525,.575,.625,.675,.725,.8};
  double z_bin_width[12] = {.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.05,.1};

  for(int i=0;i<6;i++)
  {
    if(!i)
    {
      for(int j=0;j<10;j++)
      {
        c1.cd(j+1+i*10);
        gPad->SetFillStyle(0);
      }
    }
    else
    {
      c1.cd(1+i*10);
      gPad->SetFillStyle(0);
    }
  }

  for(int c=0;c<2;c++)
  {
    for(int i=0;i<9;i++)
    {
      for(int j=0;j<5;j++)
      {
        std::vector<double> p_a;
        std::vector<double> k_a;
        std::vector<double> h_a;

        if(c)
        {
          for(int k=0;k<12;k++)
          {
            p_a.push_back((DIS[i][j] ? double(pip[i][j][k]/(DIS[i][j]*z_bin_width[k])) : 0));
            k_a.push_back((DIS[i][j] ? double(kp[i][j][k]/(DIS[i][j]*z_bin_width[k])) : 0));
            h_a.push_back((DIS[i][j] ? double(hp[i][j][k]/(DIS[i][j]*z_bin_width[k])) : 0));
            #ifdef DEBUG
            cout << c << " " << i << " " << j << " " << k << endl;
            cout << "Pi: " << (DIS[i][j] ? double(pip[i][j][k]/(DIS[i][j]*z_bin_width[k])) : 0) << endl;
            cout << "K: " << (DIS[i][j] ? double(kp[i][j][k]/(DIS[i][j]*z_bin_width[k])) : 0) << endl;
            cout << "h: " << (DIS[i][j] ? double(hp[i][j][k]/(DIS[i][j]*z_bin_width[k])) : 0) << endl;
            #endif /*DEBUG*/
          }
        }
        else
        {
          for(int k=0;k<12;k++)
          {
            p_a.push_back((DIS[i][j] ? double(pim[i][j][k]/(DIS[i][j]*z_bin_width[k])) : 0));
            k_a.push_back((DIS[i][j] ? double(km[i][j][k]/(DIS[i][j]*z_bin_width[k])) : 0));
            h_a.push_back((DIS[i][j] ? double(hm[i][j][k]/(DIS[i][j]*z_bin_width[k])) : 0));
            #ifdef DEBUG
            cout << c << " " << i << " " << j << " " << k << endl;
            cout << "Pi: " << (DIS[i][j] ? double(pim[i][j][k]/(DIS[i][j]*z_bin_width[k])) : 0) << endl;
            cout << "K: " << (DIS[i][j] ? double(km[i][j][k]/(DIS[i][j]*z_bin_width[k])) : 0) << endl;
            cout << "h: " << (DIS[i][j] ? double(hm[i][j][k]/(DIS[i][j]*z_bin_width[k])) : 0) << endl;
            #endif /*DEBUG*/
          }
        }

        H[c][i][j] = new TGraph(int(h_a.size()),&(z_range[0]),&(h_a[0]));
        Pi[c][i][j] = new TGraph(int(p_a.size()),&(z_range[0]),&(p_a[0]));
        K[c][i][j] = new TGraph(int(k_a.size()),&(z_range[0]),&(k_a[0]));

        H[c][i][j]->SetMarkerColor(fMarkerColor[c]);
        Pi[c][i][j]->SetMarkerColor(fMarkerColor[c]);
        K[c][i][j]->SetMarkerColor(fMarkerColor[c]);

        H[c][i][j]->SetMarkerSize(2);
        Pi[c][i][j]->SetMarkerSize(2);
        K[c][i][j]->SetMarkerSize(2);

        H[c][i][j]->SetMarkerStyle(fMarkerStyle[c]);
        Pi[c][i][j]->SetMarkerStyle(fMarkerStyle[c]);
        K[c][i][j]->SetMarkerStyle(fMarkerStyle[c]);

        H[c][i][j]->GetYaxis()->SetTitle("");
        Pi[c][i][j]->GetYaxis()->SetTitle("");
        K[c][i][j]->GetYaxis()->SetTitle("");

        H[c][i][j]->GetXaxis()->SetTitle("");
        Pi[c][i][j]->GetXaxis()->SetTitle("");
        K[c][i][j]->GetXaxis()->SetTitle("");

        H[c][i][j]->SetTitle("");
        Pi[c][i][j]->SetTitle("");
        K[c][i][j]->SetTitle("");

        c1.cd(i+1+j*9);
        if(H[c][i][j])
        {
          if(!c)
          {
            H[c][i][j]->Draw("SAMEPA");
            H[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
            H[c][i][j]->SetMinimum(0.0);
            H[c][i][j]->SetMaximum(4.0);
          }
          else
          {
            H[c][i][j]->Draw("SAMEP");
            H[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
            H[c][i][j]->SetMinimum(0.0);
            H[c][i][j]->SetMaximum(4.0);
          }
        }
        c1.Update();

        c2.cd(i+1+j*9);
        if(Pi[c][i][j])
        {
          if(!c)
          {
            Pi[c][i][j]->Draw("SAMEPA");
            Pi[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
            Pi[c][i][j]->SetMinimum(0.0);
            Pi[c][i][j]->SetMaximum(4.0);
          }
          else
          {
            Pi[c][i][j]->Draw("SAMEP");
            Pi[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
            Pi[c][i][j]->SetMinimum(0.0);
            Pi[c][i][j]->SetMaximum(4.0);
          }
        }
        c2.Update();

        c3.cd(i+1+j*9);
        if(K[c][i][j])
        {
          if(!c)
          {
            K[c][i][j]->Draw("SAMEPA");
            K[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
            K[c][i][j]->SetMinimum(0.0);
            K[c][i][j]->SetMaximum(4.0);
          }
          else
          {
            K[c][i][j]->Draw("SAMEP");
            K[c][i][j]->GetXaxis()->SetLimits(0.1,0.9);
            K[c][i][j]->SetMinimum(0.0);
            K[c][i][j]->SetMaximum(4.0);
          }
        }
        c3.Update();
      }
    }
  }

  c1.cd(0);

  TGaxis *axis1 = new TGaxis(-9.0,0.8,9.0,0.8,0.004,0.4,510,"");
  axis1->SetName("axis1");
  axis1->Draw();
  TGaxis *axis2 = new TGaxis(-9.5,0.95,-9.5,-0.95,0,10,510,"");
  axis2->SetName("axis2");

  c1.Update();
  c2.Update();
  c3.Update();

  c1.Print(hadron_mult_pdf);
  c2.Print(pion_mult_pdf);
  c3.Print(kaon_mult_pdf);

  in.close();

  return 0;
}
