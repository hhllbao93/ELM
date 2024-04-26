close all;
clear all;

 z =[  0.00710063521,...
  0.0279249996,...
  0.0622585751,...
  0.118865065,...
  0.2121934,...
  0.3660658,...
  0.619758487,...
  1.03802705,...
  1.72763526,...
  2.8646071,...
  4.73915672,...
  7.82976627,...
  12.9253206,...
  21.3264694,...
  35.1776199] ;
D=2.e-6;
u=1.e-7;


for kk = 1 : 600
    time=1.e-3+(kk-1).*86400;
    y1=trpoint_source_analytic_ade(z,u,D,time);
    plot(y1,-z);hold on;
    pause(0.1);
end