%%%%%%求出单元刚度矩阵（2DQ8单元）
syms x y lx ly E v 
qx=2*x/lx;qy=2*y/ly;
N1b=(1-qx)*(1-qy)/4;
N2b=(1+qx)*(1-qy)/4;
N3b=(1+qx)*(1+qy)/4;
N4b=(1-qx)*(1+qy)/4;
N5=(1-qx^2)*(1-qy)/2;
N6=(1+qx)*(1-qy^2)/2;
N7=(1-qx^2)*(1+qy)/2;
N8=(1-qx)*(1-qy^2)/2;
N1=N1b-N5/2-N8/2;
N2=N1b-N5/2-N6/2;
N3=N1b-N6/2-N7/2;
N4=N1b-N7/2-N8/2;
N1m=[N1,0;0,N1];
N2m=[N2,0;0,N2];
N3m=[N3,0;0,N3];
N4m=[N4,0;0,N4];
N5m=[N5,0;0,N5];
N6m=[N6,0;0,N6];
N7m=[N7,0;0,N7];
N8m=[N8,0;0,N8];
N=[N1m,N2m,N3m,N4m,N5m,N6m,N7m,N8m];
N1n=[diff(N1,x),0;0,diff(N1,y);diff(N1,y),diff(N1,x)];
N2n=[diff(N2,x),0;0,diff(N2,y);diff(N2,y),diff(N2,x)];
N3n=[diff(N3,x),0;0,diff(N3,y);diff(N3,y),diff(N3,x)];
N4n=[diff(N4,x),0;0,diff(N4,y);diff(N4,y),diff(N4,x)];
N5n=[diff(N5,x),0;0,diff(N5,y);diff(N5,y),diff(N5,x)];
N6n=[diff(N6,x),0;0,diff(N6,y);diff(N6,y),diff(N6,x)];
N7n=[diff(N7,x),0;0,diff(N7,y);diff(N7,y),diff(N7,x)];
N8n=[diff(N8,x),0;0,diff(N8,y);diff(N8,y),diff(N8,x)];
B=[N1n,N2n,N3n,N4n,N5n,N6n,N7n,N8n];
D=E/(1-v^2)*[1,v,0;v,1,0;0,0,(1-v)/2];
Ke=int(int(B'*D*B,x,-lx/2,lx/2),y,-ly/2,ly/2);
save Ke.mat Ke