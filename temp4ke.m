clc;clear all;
syms x y lx Hc l Ec G Hf Ef
lx=x/l;
N1=1-3*lx^2+2*lx^3;
N2=lx-2*lx^2+lx^3;
N3=3*lx^2-2*lx^3;
N4=lx^3-lx^2;
B1=y*[0 0 -1/l 0 0 1/l];
B2=[diff(N1,x),diff(N2,x),1-lx,diff(N3,x),diff(N4,x),lx];
B3=[y*diff(N1,x,2),y*diff(N2,x,2),Hc/(2*l),y*diff(N3,x,2),y*diff(N4,x,2),-Hc/(2*l)];
Ke1=B1'*B1*Ec;
Ke2=B2'*B2*G/2;
Ke3=B3'*B3*Ef*2;
Ke=int(int(Ke1,y,-Hc/2,Hc/2),x,0,l)+int(int(Ke2,y,-Hc/2,Hc/2),x,0,l)+int(int(Ke3,y,0,Hf),x,0,l);
%Ke=int(int(Ke2,y,-Hc/2,Hc/2),x,0,l);
save Ke4beam.mat Ke