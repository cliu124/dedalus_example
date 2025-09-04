clear all;
close all;
clc;

Q=10^6;

%solve the k_c_max that minimize the critical Rayleigh number, and get the
%corresponding critical Rayleigh number
syms k;
eqn = 2*k^6+3*k^4*pi^2-pi^6-Q*pi^4;
k_c=solve(eqn,k);
real_ind=find(imag(k_c)==0);
k_c_max=max(double(k_c(real_ind)));
Ra_c=(pi^2+k_c_max^2)/k_c_max^2*((pi^2+k_c_max^2)^2+pi^2*Q);


%solve the k_c for given Ra and Q
Ra=1.5*10^7;
syms k;
eqn_Ra = -Ra*k^2+ (pi^2+k^2)*((pi^2+k^2)^2+pi^2*Q);
k_c_Ra=solve(eqn_Ra,k);
real_ind=find(imag(k_c_Ra)==0);
k_c_Ra_max=(double(k_c_Ra(real_ind)));



