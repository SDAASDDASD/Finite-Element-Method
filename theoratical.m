clc;clear all;
Hc=1/3;Hf=1/3;Ef=6.9e10;Ec=Ef*10^(-4);b=1;
A0=Ec*b*Hc^3/12+2*Ef*b*((Hc/2+Hf)^3-(Hc/2)^3)/3;
F=-100;L=50;
x1=linspace(0,L/2,200);
w1=-F.*x1.^3/(12*A0)+F*L^2.*x1/(16*A0);
x2=linspace(L/2,L,200);
w2=F.*x2.^3/(12*A0)-F*L.*x2.^2/(4*A0)+3*F*L^2.*x2/(16*A0)-F*L^3/(48*A0);
plot(x1,w1,'b',x2,w2,'b','linewidth',2)
xlabel('Distance / *10^{-2} m');
ylabel('Deflection / m');
title('Deflection-Distance')