clc
clear all

a = 0.45 ;                     % crack length
L = 1 ;
D = 2 ;
xCr   = [0 D/2; a D/2];
xTip  = [a D/2];
seg   = xCr(2,:) - xCr(1,:);   % tip segment
alpha = atan2(seg(2),seg(1));  % inclination angle
QT    =[cos(alpha) sin(alpha); -sin(alpha) cos(alpha)];


fid = fopen('m0XY.prn','r');
node = fscanf(fid,'%d,%g,%g\n',[3,inf]);
node = node';
fclose(fid);

fid = fopen('m0Top.prn','r');
element = fscanf(fid,'%d,%d,%d,%d,%d\n',[5,inf]);
element = element';
fclose(fid);

numnode = size(node,1);
numelem = size(element,1);

x0  = xCr(1,1); y0 = xCr(1,2);
x1  = xCr(2,1); y1 = xCr(2,2);
t   = 1/norm(seg)*seg;
for i = 1 : numnode
    x = node(i,2);
    y = node(i,3);
    l   = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)) ;
    phi = (y0-y1)*x + (x1-x0)*y + (x0*y1-x1*y0);
    ls(i,1) = phi/l;            % normal LS
    ls(i,2) = ([x y]-xTip)*t';  % tangent LS
end

enrich_node = zeros(numnode,1);

count1 = 0;
count2 = 0;
for iel = 1 : numelem
    sctr = element(iel,2:5);
    phi  = ls(sctr,1);
    psi  = ls(sctr,2);
    if ( max(phi)*min(phi) < 0 )
        if max(psi) < 0
            count1 = count1 + 1 ; % ah, one split element
            split_elem(count1) = iel;
            enrich_node(sctr)   = 1;
        elseif max(psi)*min(psi) < 0
            count2 = count2 + 1 ; % ah, one tip element
            tip_elem(count2) = iel;
            enrich_node(sctr)   = 2;
        end
    end
end
split_nodes = find(enrich_node == 1);
tip_nodes   = find(enrich_node == 2);

TopType = element;
TopX = [];
TopNoX = [];
TopXTypeX = [];
for iel = 1 : numelem
    sctr = element(iel,2:5);
    TopType(iel,2:5) = enrich_node(sctr,1);
    if (any(enrich_node(sctr,1)))
    TopX = [TopX;element(iel,:)];
    TopXTypeX = [TopXTypeX;TopType(iel,:)];
    else
     TopNoX = [TopNoX;element(iel,:)];
    end
end

fid = fopen('TopX','w');
fprintf(fid,'%d,%d,%d,%d,%d\n',TopX');
fclose(fid);

fid = fopen('TopNoX','w');
fprintf(fid,'%d,%d,%d,%d,%d\n',TopNoX');
fclose(fid);

fid = fopen('TopXTypeX','w');
fprintf(fid,'%d,%d,%d,%d,%d\n',TopXTypeX');
fclose(fid);

ncrack = 1;
maxCP = 2;
nelemX = size(TopX,1);
nnodeX = size(split_nodes,1)+size(tip_nodes,1);
fid = fopen('GGinfoX','w');
fprintf(fid,'%d,%d,%d,%d',ncrack,maxCP,nelemX,nnodeX);
fclose(fid);

fid = fopen('GGXYC','w');
fprintf(fid,'%d\n',maxCP);
fprintf(fid,'%g,%g\n',x0,y0);
fprintf(fid,'%g,%g\n',x1,y1);
fclose(fid);


SETNodeX2dof = [];
SETNodeX4dof = [];
SETNodeX10dof = [];
for i =1:size(TopX,1)
    strc = TopX(i,2:5);
    for j =1:size(strc,2)
        if enrich_node(strc(j)) == 0
            SETNodeX2dof = union(SETNodeX2dof,strc(j));
        elseif enrich_node(strc(j)) == 1
            SETNodeX4dof = union(SETNodeX4dof,strc(j));
        elseif enrich_node(strc(j)) == 2
            SETNodeX10dof = union(SETNodeX10dof,strc(j));
        end
    end
end

fid = fopen('SETNodeX2dof','w');
fprintf(fid,'%d,%d,%d,%d,%d,%d,%d,%d\n',SETNodeX2dof');
fclose(fid);

fid = fopen('SETNodeX4dof','w');
fprintf(fid,'%d,%d,%d,%d,%d,%d,%d,%d\n',SETNodeX4dof');
fclose(fid);

fid = fopen('SETNodeX10dof','w');
fprintf(fid,'%d,%d,%d,%d,%d,%d,%d,%d\n',SETNodeX10dof');
fclose(fid);

EnodeX  = union(SETNodeX2dof,union(SETNodeX4dof,SETNodeX10dof));
GGnodeX = zeros(size(EnodeX,2),4);
GGnodeX(:,1) = EnodeX';
GGnodeX(:,2) = enrich_node(EnodeX,1);
GGnodeX(:,3:4) = ls(EnodeX,:);
fid = fopen('GGnodeX','w');
fprintf(fid,'%d,%d,%g,%g\n',GGnodeX');
fclose(fid);

