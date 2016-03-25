type_elem = zeros(size(element,1),size(xCr,2));        % ��ŵ�Ԫ����
elem_crk = zeros(size(element,1),4);                          % ��������뵥Ԫ�߽罻��������Ѽ������                                                                                        
xTip = zeros(size(element,1),2);                                   % ��ŵ�Ԫ���Ѽ������
xVertex = zeros(size(element,1),2);                              % ��ŵ�Ԫ�����ƶ������꣬��Ԫ�����ֻ��һ�����ƶ���
xJertex = zeros(size(element,1),4);                               % ��������뵥Ԫ�߽罻������ͽ��������
enrich_node = zeros(size(node,1),size(xCr,2));           % ��Ž���ǿ��ʽ
for kk = 1:size(xCr,2)                                                       % ����������ѭ��
    for iel = 1:size(element,1)                                           % �����е�Ԫѭ��
        sctr = element(iel,2:5);                                                 % ��ǰ��Ԫ������Ϣ
        vv = node(sctr,2:3);                                                      % ��ǰ��Ԫ�������
        crk_int = [];
        intes = 0; flag1 = 0; flag2 = 0;
        for kj = 1:size(xCr(kk).coor,1)-1                               % �������߶�ѭ��
            q1 = xCr(kk).coor(kj,:); q2 = xCr(kk).coor(kj+1,:);     % �����߶ε������˵�����
            sctrl = [sctr sctr(1,1)];
            for iedge = 1:size(sctr,2)                                            % �Ե�Ԫ�߽�ѭ��
                nnode1 = sctrl(iedge); nnode2 = sctrl(iedge+1);  % ��Ԫ�߽��ϵ���������
                p1 = node(nnode1,2:3); p2 = node(nnode2,2:3);          % ��Ԫ�߽��ϵ������������
                intersect = segments_int_2d(p1,p2,q1,q2);          % ���������߶εĽ���
                % intersect(1) = 0, ���߶β��ཻ��intersect(1) = 1�����߶��ཻ
                % intersect(2:3) Ϊ��������
                intes = intes + intersect(1);
                if intersect(1) > 0
                    crk_int = [crk_int intersect(2) intersect(3)];
                    flag1 = inhull(xCr(kk).coor(kj,:),vv,[],-1e-8);           % �ж�xCr(kk).coor(kj,:)�Ƿ��ڵ�Ԫ�ڣ�flag1 = 1���ڵ�Ԫ��
                    flag2 = inhull(xCr(kk).coor(kj+1,:),vv,[],-1e-8);       % �ж�xCr(kk).coor(kj+1,:)�Ƿ��ڵ�Ԫ�ڣ�flag2 = 1���ڵ�Ԫ��
                    xCr_element(iel,:) = xCr(kk).coor(kj,:)*flag1 + xCr(kk).coor(kj+1,:)*flag2;
                end
            end
        end
        
        
        %�жϵ�Ԫ����
        if (((intes == 2)&&(flag1 == 0))&&(flag2 == 0))         %�ᴩ��Ԫ
            type_elem(iel,kk) = 2;
            elem_crk(iel,:) = crk_int;                                         %�����뵥Ԫ�߽�Ľ�������
        end
        if (((flag1 ==1 )||(flag2 == 1))&&(intes == 2))            %������ĵ�Ԫ
            type_elem(iel,kk) = 2;
            elem_crk(iel,:) = crk_int;                                          %�����뵥Ԫ�߽�Ľ�������
            xVertex(iel,:) = xCr_element(iel,:);                          %��Ԫ�����ƶ�������
        end
        if (intes == 1)                                                                %�ѼⵥԪ�򽻲浥Ԫ
            type_elem(iel,kk) = 1;
            xTip(iel,:) = xCr_element(iel,:);                                %�Ѽ������
            for kk1 =1:size(xCr,2)                                              %�ж��Ѽ��Ƿ�������������
                if (kk1 ~= kk)
                    for kj1 = 1:size(xCr(kk1).coor,1)-1
                        q1 = xCr(kk1).coor(kj1,:); q2 = xCr(kk1).coor(kj1+1,:);
                        [ival] = point_in_line(xTip(iel,:),q1,q2);        %�ж�xTip(iel,:)�Ƿ���q1-q2����
                        if (ival == 1)                                                   %��������
                            type_elem(iel,kk) = 4;                              %���浥Ԫ
                            xJertex(iel,:) = [crk_int xTip(iel,1) xTip(iel,2)];
                        end
                    end
                end
            end
           if (type_elem(iel,kk) ~= 4)                                         %�ѼⵥԪ
               type_elem(iel,kk) = 1;
               elem_crk(iel,:) = [crk_int xTip(iel,1) xTip(iel,2)];
           end
        end
    end
end


%�жϽ���ǿ��ʽ
for kk = 1:size(xCr,2)
    for iel = 1:size(element,1)
        sctr = element(iel,2:5);
        if type_elem(iel,kk) == 1                       %�ѼⵥԪ
            enrich_node(sctr,kk) = 1;                  %�����Ѽ��֧������ǿ
        elseif type_elem(iel,kk) ==2                 %�ᴩ��Ԫ
            for in = 1:length(sctr)                         %�Ե�Ԫ�����ѭ��
                if enrich_node(sctr(in),kk) == 0     %û����ǿ�Ľ��
                    [Aw,Awp] = support_area(sctr(in),iel,type_elem,elem_crk,xVertex,kk);
                    %������֧�������Aw,�Լ�λ������һ��֧�������Awp
                    if (abs(Awp/Aw)>1e-4)&&(abs((Aw-Awp)/Aw)>1e-4)
                        enrich_node(strc(in),kk) = 2;        %��Ծ������ǿ
                    end
                end
            end
        elseif type_elem(iel,kk) == 3                        %�����㵥Ԫ
            for in = 1:length(sctr)
                if enrich_node(sctr(in),kk) == 0            %û����ǿ�Ľ��
                    [Aw,Awp] = support_area(sctr(in),iel,type_elem,elem_crk,xVertex,kk);
                    if (abs(Awp/Aw)>1e-4)&&(abs((Aw-Awp)/Aw)>1e-4)
                        enrich_node(strc(in),kk) = 2;        %��Ծ������ǿ
                    end
                end
            end
        elseif type_elem(iel,kk) == 4                       %���浥Ԫ
            for in = 1:length(sctr)
                enrich_node(strc(in),kk) = 3;               %���Ӻ�����ǿ
            end
        end
    end
end


%������������������ϵ�Ԫ�����ӵ��Ѽ��֧������ǿ���
for kk = 1:size(xCr,2)
    for iel = 1:size(element,1)
        sctr = element(iel,2:5);
        if (ismember(1,enrich_node(sctr,kk)) ~= 0)      %��Ԫ�ں��Ѽ��֧������ǿ�Ľ��
            for j = 1:length(sctr)
                if (enrich_node(sctr(j),kk) == 0)                 %û����ǿ�Ľ��
                    enrich_node(sctr(j),kk) = 11;                  %�Ѽ��֧������ǿ
                elseif (enrich_node(sctr(j),kk) == 2)          %��Ծ������ǿ�Ľ��
                    enrich_node(sctr(j),kk) = 22;                  %�Ѽ��֧�����ͽ�Ծ������ͬ��ǿ
                end
            end
        end
    end
end