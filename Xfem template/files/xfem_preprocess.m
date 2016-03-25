type_elem = zeros(size(element,1),size(xCr,2));        % 存放单元类型
elem_crk = zeros(size(element,1),4);                          % 存放裂纹与单元边界交点坐标和裂尖点坐标                                                                                        
xTip = zeros(size(element,1),2);                                   % 存放单元内裂尖点坐标
xVertex = zeros(size(element,1),2);                              % 存放单元内裂纹顶点坐标，单元内最多只有一个裂纹顶点
xJertex = zeros(size(element,1),4);                               % 存放裂纹与单元边界交点坐标和交叉点坐标
enrich_node = zeros(size(node,1),size(xCr,2));           % 存放结点加强方式
for kk = 1:size(xCr,2)                                                       % 对所有裂纹循环
    for iel = 1:size(element,1)                                           % 对所有单元循环
        sctr = element(iel,2:5);                                                 % 当前单元连接信息
        vv = node(sctr,2:3);                                                      % 当前单元结点坐标
        crk_int = [];
        intes = 0; flag1 = 0; flag2 = 0;
        for kj = 1:size(xCr(kk).coor,1)-1                               % 对裂纹线段循环
            q1 = xCr(kk).coor(kj,:); q2 = xCr(kk).coor(kj+1,:);     % 裂纹线段的两个端点坐标
            sctrl = [sctr sctr(1,1)];
            for iedge = 1:size(sctr,2)                                            % 对单元边界循环
                nnode1 = sctrl(iedge); nnode2 = sctrl(iedge+1);  % 单元边界上的两个结点号
                p1 = node(nnode1,2:3); p2 = node(nnode2,2:3);          % 单元边界上的两个结点坐标
                intersect = segments_int_2d(p1,p2,q1,q2);          % 计算两个线段的交点
                % intersect(1) = 0, 两线段不相交；intersect(1) = 1，两线段相交
                % intersect(2:3) 为交点坐标
                intes = intes + intersect(1);
                if intersect(1) > 0
                    crk_int = [crk_int intersect(2) intersect(3)];
                    flag1 = inhull(xCr(kk).coor(kj,:),vv,[],-1e-8);           % 判断xCr(kk).coor(kj,:)是否在单元内，flag1 = 1，在单元内
                    flag2 = inhull(xCr(kk).coor(kj+1,:),vv,[],-1e-8);       % 判断xCr(kk).coor(kj+1,:)是否在单元内，flag2 = 1，在单元内
                    xCr_element(iel,:) = xCr(kk).coor(kj,:)*flag1 + xCr(kk).coor(kj+1,:)*flag2;
                end
            end
        end
        
        
        %判断单元类型
        if (((intes == 2)&&(flag1 == 0))&&(flag2 == 0))         %贯穿单元
            type_elem(iel,kk) = 2;
            elem_crk(iel,:) = crk_int;                                         %裂纹与单元边界的交点坐标
        end
        if (((flag1 ==1 )||(flag2 == 1))&&(intes == 2))            %含顶点的单元
            type_elem(iel,kk) = 2;
            elem_crk(iel,:) = crk_int;                                          %裂纹与单元边界的交点坐标
            xVertex(iel,:) = xCr_element(iel,:);                          %单元内裂纹顶点坐标
        end
        if (intes == 1)                                                                %裂尖单元或交叉单元
            type_elem(iel,kk) = 1;
            xTip(iel,:) = xCr_element(iel,:);                                %裂尖点坐标
            for kk1 =1:size(xCr,2)                                              %判断裂尖是否在其他裂纹上
                if (kk1 ~= kk)
                    for kj1 = 1:size(xCr(kk1).coor,1)-1
                        q1 = xCr(kk1).coor(kj1,:); q2 = xCr(kk1).coor(kj1+1,:);
                        [ival] = point_in_line(xTip(iel,:),q1,q2);        %判断xTip(iel,:)是否在q1-q2线上
                        if (ival == 1)                                                   %点在线上
                            type_elem(iel,kk) = 4;                              %交叉单元
                            xJertex(iel,:) = [crk_int xTip(iel,1) xTip(iel,2)];
                        end
                    end
                end
            end
           if (type_elem(iel,kk) ~= 4)                                         %裂尖单元
               type_elem(iel,kk) = 1;
               elem_crk(iel,:) = [crk_int xTip(iel,1) xTip(iel,2)];
           end
        end
    end
end


%判断结点加强方式
for kk = 1:size(xCr,2)
    for iel = 1:size(element,1)
        sctr = element(iel,2:5);
        if type_elem(iel,kk) == 1                       %裂尖单元
            enrich_node(sctr,kk) = 1;                  %采用裂尖分支函数加强
        elseif type_elem(iel,kk) ==2                 %贯穿单元
            for in = 1:length(sctr)                         %对单元结点数循环
                if enrich_node(sctr(in),kk) == 0     %没被加强的结点
                    [Aw,Awp] = support_area(sctr(in),iel,type_elem,elem_crk,xVertex,kk);
                    %计算结点支撑域面积Aw,以及位于裂纹一侧支撑域面积Awp
                    if (abs(Awp/Aw)>1e-4)&&(abs((Aw-Awp)/Aw)>1e-4)
                        enrich_node(strc(in),kk) = 2;        %阶跃函数加强
                    end
                end
            end
        elseif type_elem(iel,kk) == 3                        %含顶点单元
            for in = 1:length(sctr)
                if enrich_node(sctr(in),kk) == 0            %没被加强的结点
                    [Aw,Awp] = support_area(sctr(in),iel,type_elem,elem_crk,xVertex,kk);
                    if (abs(Awp/Aw)>1e-4)&&(abs((Aw-Awp)/Aw)>1e-4)
                        enrich_node(strc(in),kk) = 2;        %阶跃函数加强
                    end
                end
            end
        elseif type_elem(iel,kk) == 4                       %交叉单元
            for in = 1:length(sctr)
                enrich_node(strc(in),kk) = 3;               %连接函数加强
            end
        end
    end
end


%采用线增函数消除混合单元所增加的裂尖分支函数加强结点
for kk = 1:size(xCr,2)
    for iel = 1:size(element,1)
        sctr = element(iel,2:5);
        if (ismember(1,enrich_node(sctr,kk)) ~= 0)      %单元内含裂尖分支函数加强的结点
            for j = 1:length(sctr)
                if (enrich_node(sctr(j),kk) == 0)                 %没被加强的结点
                    enrich_node(sctr(j),kk) = 11;                  %裂尖分支函数加强
                elseif (enrich_node(sctr(j),kk) == 2)          %阶跃函数加强的结点
                    enrich_node(sctr(j),kk) = 22;                  %裂尖分支函数和阶跃函数共同加强
                end
            end
        end
    end
end