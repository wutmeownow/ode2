#include "RK.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include <iostream>
#include <cstdio>

using namespace std;

// differential equation to be solved
double fun1(double x, double y){
  (void)x;              // prevent unused variable warning
  return -2*y;          // f = y'(x,y) = -2 * y(x)  
}                       // solution: y(x) = 3 * exp(-2*x) ; with initial condition y(0)=3

double fun2(double x, double y){
  return -y/x-2/(x*x);  // f = y'(x,y) = -y(x)/x - 2/x^2 
}                       // -2*log(|x|)/x+2/x  ; with initial condition y(0)=2

// diff eqs from https://openstax.org/books/calculus-volume-2/pages/4-5-first-order-linear-equations
double prob245(double x, double y){
  return (1./(1.+TMath::Power(x,2)))*(y-1); // y'(x,y)= 1/(1+x^2) * (y-1)
}

double prob246(double x, double y){
  return (1./x)*(y+2.*x*TMath::Log(x)); // y'(x,y) = (1/x)*(y + 2xlnx)
}


double sol245(double* xPtr, double par[]){        
  double x = *xPtr; 
  return 1.-TMath::Exp(TMath::ATan(x));
}

double sol246(double* xPtr, double par[]){        
  double x = *xPtr; 
  return x*TMath::Power(TMath::Log(x),2) + 5*x;
}

// function to make error graphs from chatgpt
TGraph makeErrorGraph(const TGraph& g_method, const TF1& f_exact) {
    int n = g_method.GetN();
    // cout<<n<<endl;
    TGraph g_error;
    g_error.SetTitle("Percent Error; x;");

    for (int i = 0; i < n; ++i) {
        double x, y;
        g_method.GetPoint(i, x, y);
        double y_exact = f_exact.Eval(x);
        double err = 0;
        if (std::abs(y_exact) > 1e-12){
          err = std::abs(y - y_exact)/std::abs(y);
        } 
        g_error.SetPoint(i, x, err);
    }

    return g_error;
}

int main(int argc, char **argv){
  TApplication theApp("App", &argc, argv); // init ROOT App for displays

  // solve our DEQ using RK1 or RK2 methods!
  // Two examples are given.  Choose a function for testing
  // TGraph tg1=RK1Solve(fun1,3,30,0,3);                     // initial condition y(0)=3
  // TGraph tg2=RK2Solve(fun1,3,30,0,3);
  // TF1 fun_sol=TF1("fun_sol","3*exp(-2*x)",0,3);           // exact solution
  //TGraph tg1=RK1Solve(fun2,2,100,1,100);                // initial condition y(1)=2
  //TGraph tg2=RK2Solve(fun2,2,100,1,100);
  //TF1 fun_sol=TF1("fun_sol","-2*log(x)/x+2/x",1,100);   // exact solution

  // my odes and solutions
  const double xmax = 10;
  const double nsteps = xmax*10;
  TGraph tg1=RK1Solve(prob245,0,nsteps,0,xmax);                     // initial condition y(0)=0
  TGraph tg2=RK2Solve(prob245,0,nsteps,0,xmax);
  TGraph tg4=RK4Solve(prob245,0,nsteps,0,xmax);
  TF1 fun_sol=TF1("fun_sol",sol245,0,xmax);           // exact solution
  // TGraph tg1=RK1Solve(prob246,5,30,1,4);                     // initial condition y(1)=5
  // TGraph tg2=RK2Solve(prob246,5,30,1,4);
  // TF1 fun_sol=TF1("fun_sol",sol246,1,4);           // exact solution

  // ******************************************************************************
  // ** this block is useful for supporting both high and std resolution screens **
  // UInt_t dh = gClient->GetDisplayHeight()/2;   // fix plot to 1/2 screen height  
  //UInt_t dw = gClient->GetDisplayWidth();
  // UInt_t dw = 1.1*dh;
  // ******************************************************************************

  // percent error graphs
  TGraph g_err1 = makeErrorGraph(tg1, fun_sol);
  TGraph g_err2 = makeErrorGraph(tg2, fun_sol);
  TGraph g_err4 = makeErrorGraph(tg4, fun_sol);

  // std::cout << "Error graph points: " << g_err1.GetN() << ", "<< g_err2.GetN() << ", " << g_err4.GetN() << std::endl;


  // Style and draw
  g_err1.SetLineColor(kRed);
  g_err2.SetLineColor(kBlue);
  g_err4.SetLineColor(kGreen+2);
  g_err1.SetLineWidth(2);
  g_err2.SetLineWidth(2);
  g_err4.SetLineWidth(2);


  TCanvas *c1 = new TCanvas("c1","DEQ solutions",1600,800);
  c1->Divide(2,1);
  c1->cd(1);

  // tg1.SetMarkerSize(0.015*dh/8);  // size scale: 1 = 8 pixels, so here we choose the size to be 1.5% of the window height
  // tg2.SetMarkerSize(0.015*dh/8);
  tg1.SetMarkerStyle(kFullTriangleUp);
  tg2.SetMarkerStyle(kFullTriangleDown);
  tg4.SetMarkerStyle(kFullDotMedium);
  tg1.SetMarkerColor(kRed);
  tg2.SetMarkerColor(kGreen-2);
  tg4.SetMarkerColor(kBlue);
  fun_sol.SetLineColor(kBlack);
  fun_sol.SetLineStyle(2);
  
  // plot the results
  tg1.SetTitle("ODE demo;x;y");
  tg1.Draw("AP");
  tg2.Draw("P");
  tg4.Draw("P");
  fun_sol.Draw("same");
  
  TLegend *tl = new TLegend(0.6,0.7,0.9,0.9);
  tl->AddEntry(&tg1,"RK1 Solution","p");
  tl->AddEntry(&tg2,"RK2 Solution","p");
  tl->AddEntry(&tg4,"RK4 Solution","p");
  tl->AddEntry(&fun_sol,"Exact Solution","l");
  tl->Draw();

  // errors
  c1->cd(2);
  gPad->SetLogy();
  g_err1.Draw("AL");
  g_err2.Draw("L SAME");
  g_err4.Draw("L SAME");

  TLegend* leg = new TLegend(0.55, 0.75, 0.88, 0.88);
  leg->AddEntry(&g_err1, "RK1 (Euler)", "l");
  leg->AddEntry(&g_err2, "RK2", "l");
  leg->AddEntry(&g_err4, "RK4", "l");
  leg->Draw();


  c1->Draw();
  c1->Update();
  c1->Print("OED_cpp.png");

  // retreive the data from the graphs and write to a file
  FILE *fp=fopen("RKdemo.dat","w");
  double *x, *y1, *y2, *y4;
  x=tg1.GetX();
  y1=tg1.GetY();
  y2=tg2.GetY();
  y4=tg4.GetY();
  fprintf(fp,"#%8s %9s %9s %9s %9s\n","x","RK1","RK2","RK4","Exact");
  for (int i=0; i<tg1.GetN(); i++){
    fprintf(fp,"%9.4lf %9.4lf %9.4lf %9.4lf %9.4lf\n",x[i],y1[i],y2[i],y4[i],fun_sol.Eval(x[i]));
  }
  fclose(fp);
  
  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
  theApp.Run();
}

