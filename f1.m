function [f1]=f1(M)
%  ijk=-1; 
%  i=(-1)^0.5; j=(-1)^0.5; k=(-1)^0.5;
 j=(-1)^0.5;
 f1 = (1-exp(-real(M)))./(1+exp(-real(M)))+j./(1+exp(-imag(M)));
 %2*tanh(real(M))+2*j*tanh(imag(M));
 %f1 = 0.2*(1-exp(-real(M)))./(1+exp(-real(M)))+0.2*j./(1+exp(-imag(M)));