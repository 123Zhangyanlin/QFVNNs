function [f2]=f2(N)
  j=(-1)^0.5;
 % f2 = 2*tanh(real(N))+2*j*tanh(imag(N));
  f2 = (1-exp(-real(N)))./(1+exp(-real(N)))+j./(1+exp(-imag(N)));