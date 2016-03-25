function [ intersect ] = segments_int_2d( p1,p2,q1,q2 )
%二维线段求交点
%p1,p2,q1,q2为四个点坐标
% intersect(1) = 0, 两线段不相交；intersect(1) = 1，两线段相交
% intersect(2:3) 为交点坐标
intersect  = zeros(1,3);
flag1 = 0; flag2 = 0;
vector1 = p2 - p1;
vector2 = q2 - q1;
if (vector1(1)*vector2(2) == vector2(1)*vector1(2))                  %两线段平行
    intersect = [0,0,0];
else
    A = [vector1(2)      -vector1(1)
            vector2(2)      -vector2(1)];
    b = [p1(1)*p2(2)-p1(2)*p2(1)
            q1(1)*q2(2)-q1(2)*q2(1)];
    r = A\b;
    r = r';
    vector3 = r - p1;
    vector4 = r - p2;
    if abs(norm(vector3)+norm(vector4) - norm(vector1))<1e-4
        flag1 = 1;                                             %交点在单元边界上
        intersect(2:3) = [r(1),r(2)];                   %方便裂尖单元调用
    end
    vector5 = r - q1;
    vector6 = r - q2;
    if abs(norm(vector5)+norm(vector6) - norm(vector2))<1e-4
        flag2 = 1;
    end
    if (flag1 && flag2)
        intersect(1) = 1;
        intersect(2:3) = [r(1),r(2)];
    end
end

