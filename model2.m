clc;clear
load Ke4beam.mat
subplot(2,2,1)
Hc=1/3;Hf=1/3;Ef=6.9e10;Ec=Ef*10^(0);v=0.3;G=Ec/2*(1+v);b=1;
n=200;l=50/n;
Ke=subs(Ke);
K=zeros(3*(n+1),3*(n+1));
F=zeros(3*(n+1),1);
%%%%%%%%%%%%%%%集成总体刚度矩阵
for i=1:n
    K((3*i-2):(3*i+3),(3*i-2):(3*i+3))=Ke+K((3*i-2):(3*i+3),(3*i-2):(3*i+3));
end
%%%%%%%%%%%%%%%%引入位移边界条件
K(1,1)=1;
K(1,2:3*(n+1))=0;
K(2:3*(n+1),1)=0;
K(3*n+1,3*n+1)=1;
for i=1:3*(n+1)
    if i~=3*n+1
        K(3*n+1,i)=0;
        K(i,3*n+1)=0;
    end
end
%%%%%%%%%%%%%%%%%引入应力边界条件
F(3*(n/2+1)-2,1)=-100;
u=K\F;
x=linspace(0,50,n+1);
w=u(1:3:3*n+1);
plot(x,w,'b','linewidth',1.5)
hold on

A0=Ec*b*Hc^3/12+2*Ef*b*((Hc/2+Hf)^3-(Hc/2)^3)/3;
f=-100;L=50;
x1=linspace(0,L/2,200);
w1=-f.*x1.^3/(12*A0)+f*L^2.*x1/(16*A0);
x2=linspace(L/2,L,200);
w2=f.*x2.^3/(12*A0)-f*L.*x2.^2/(4*A0)+3*f*L^2.*x2/(16*A0)-f*L^3/(48*A0);
plot(x1,w1,'r',x2,w2,'r','linewidth',1.5)

xlabel('Distance / *10^{-2} m');
ylabel('Deflection / m');
title('E_c = E_f')
legend('FEM','CLT')

clc;clear
load Ke4beam.mat
subplot(2,2,2)
Hc=1/3;Hf=1/3;Ef=6.9e10;Ec=Ef*10^(-1);v=0.3;G=Ec/2*(1+v);b=1;
n=200;l=50/n;
Ke=subs(Ke);
K=zeros(3*(n+1),3*(n+1));
F=zeros(3*(n+1),1);
%%%%%%%%%%%%%%%集成总体刚度矩阵
for i=1:n
    K((3*i-2):(3*i+3),(3*i-2):(3*i+3))=Ke+K((3*i-2):(3*i+3),(3*i-2):(3*i+3));
end
%%%%%%%%%%%%%%%%引入位移边界条件
K(1,1)=1;
K(1,2:3*(n+1))=0;
K(2:3*(n+1),1)=0;
K(3*n+1,3*n+1)=1;
for i=1:3*(n+1)
    if i~=3*n+1
        K(3*n+1,i)=0;
        K(i,3*n+1)=0;
    end
end
%%%%%%%%%%%%%%%%%引入应力边界条件
F(3*(n/2+1)-2,1)=-100;
u=K\F;
x=linspace(0,50,n+1);
w=u(1:3:3*n+1);
plot(x,w,'b','linewidth',1.5)
hold on

A0=Ec*b*Hc^3/12+2*Ef*b*((Hc/2+Hf)^3-(Hc/2)^3)/3;
f=-100;L=50;
x1=linspace(0,L/2,200);
w1=-f.*x1.^3/(12*A0)+f*L^2.*x1/(16*A0);
x2=linspace(L/2,L,200);
w2=f.*x2.^3/(12*A0)-f*L.*x2.^2/(4*A0)+3*f*L^2.*x2/(16*A0)-f*L^3/(48*A0);
plot(x1,w1,'r',x2,w2,'r','linewidth',1.5)

xlabel('Distance / *10^{-2} m');
ylabel('Deflection / m');
title('E_c = E_f*10^{-1}')
legend('FEM','CLT')

clc;clear
load Ke4beam.mat
subplot(2,2,3)
Hc=1/3;Hf=1/3;Ef=6.9e10;Ec=Ef*10^(-2);v=0.3;G=Ec/2*(1+v);b=1;
n=200;l=50/n;
Ke=subs(Ke);
K=zeros(3*(n+1),3*(n+1));
F=zeros(3*(n+1),1);
%%%%%%%%%%%%%%%集成总体刚度矩阵
for i=1:n
    K((3*i-2):(3*i+3),(3*i-2):(3*i+3))=Ke+K((3*i-2):(3*i+3),(3*i-2):(3*i+3));
end
%%%%%%%%%%%%%%%%引入位移边界条件
K(1,1)=1;
K(1,2:3*(n+1))=0;
K(2:3*(n+1),1)=0;
K(3*n+1,3*n+1)=1;
for i=1:3*(n+1)
    if i~=3*n+1
        K(3*n+1,i)=0;
        K(i,3*n+1)=0;
    end
end
%%%%%%%%%%%%%%%%%引入应力边界条件
F(3*(n/2+1)-2,1)=-100;
u=K\F;
x=linspace(0,50,n+1);
w=u(1:3:3*n+1);
plot(x,w,'b','linewidth',1.5)
hold on

A0=Ec*b*Hc^3/12+2*Ef*b*((Hc/2+Hf)^3-(Hc/2)^3)/3;
f=-100;L=50;
x1=linspace(0,L/2,200);
w1=-f.*x1.^3/(12*A0)+f*L^2.*x1/(16*A0);
x2=linspace(L/2,L,200);
w2=f.*x2.^3/(12*A0)-f*L.*x2.^2/(4*A0)+3*f*L^2.*x2/(16*A0)-f*L^3/(48*A0);
plot(x1,w1,'r',x2,w2,'r','linewidth',1.5)

xlabel('Distance / *10^{-2} m');
ylabel('Deflection / m');
title('E_c = E_f*10^{-2}')
legend('FEM','CLT')

clc;clear
load Ke4beam.mat
subplot(2,2,4)
Hc=1/3;Hf=1/3;Ef=6.9e10;Ec=Ef*10^(-3);v=0.3;G=Ec/2*(1+v);b=1;
n=200;l=50/n;
Ke=subs(Ke);
K=zeros(3*(n+1),3*(n+1));
F=zeros(3*(n+1),1);
%%%%%%%%%%%%%%%集成总体刚度矩阵
for i=1:n
    K((3*i-2):(3*i+3),(3*i-2):(3*i+3))=Ke+K((3*i-2):(3*i+3),(3*i-2):(3*i+3));
end
%%%%%%%%%%%%%%%%引入位移边界条件
K(1,1)=1;
K(1,2:3*(n+1))=0;
K(2:3*(n+1),1)=0;
K(3*n+1,3*n+1)=1;
for i=1:3*(n+1)
    if i~=3*n+1
        K(3*n+1,i)=0;
        K(i,3*n+1)=0;
    end
end
%%%%%%%%%%%%%%%%%引入应力边界条件
F(3*(n/2+1)-2,1)=-100;
u=K\F;
x=linspace(0,50,n+1);
w=u(1:3:3*n+1);
plot(x,w,'b','linewidth',1.5)
hold on

A0=Ec*b*Hc^3/12+2*Ef*b*((Hc/2+Hf)^3-(Hc/2)^3)/3;
f=-100;L=50;
x1=linspace(0,L/2,200);
w1=-f.*x1.^3/(12*A0)+f*L^2.*x1/(16*A0);
x2=linspace(L/2,L,200);
w2=f.*x2.^3/(12*A0)-f*L.*x2.^2/(4*A0)+3*f*L^2.*x2/(16*A0)-f*L^3/(48*A0);
plot(x1,w1,'r',x2,w2,'r','linewidth',1.5)

xlabel('Distance / *10^{-2} m');
ylabel('Deflection / m');
title('E_c = E_f*10^{-3}')
legend('FEM','CLT')
%p=find(w==min(w));
%text(x(p),w(p),['(',num2str(w(p)),')'],'color','r')
%plot([x(p),x(p)],[0,w(p)],':b','linewidth',1.5);