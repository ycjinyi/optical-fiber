classdef FiberLmDistribution
%此类用于计算得出光纤端面理想出射光通量分布

    properties
        %光纤的半径(m)
        R;
        %最大出射角度
        U;
        %虚拟光源的高度
        H;  
    end

    methods
        function obj = FiberLmDistribution()
        end

        %根据径向距离x, 轴向距离h, 求解出射角度
        function theta = thetaCompute(obj, h, x)
            theta = atan(x / (obj.H + h));
        end

        %根据当前的轴向距离h(m), 径向距离x(m)计算光照度lm(lm)
        function lm = lmCompute(obj, h, x)
            theta = obj.thetaCompute(h, x);
            %超过最大出射角范围的照度为0
            if theta > obj.U
                lm = 0;
                return;
            end
            lm = power(cos(theta), 4) / (2 * power(h + obj.H, 2));
        end

        %该函数计算在光纤端面的不同轴向距离和径向距离下的光照度变化
        %输入参数R代表光纤的半径(m), U代表最大出射角(rad), S代表总的光通量大小(lm)
        %xPoints代表到纤芯的径向距离点(m), 1行n列
        %hpoints代表到纤芯的轴向距离点(m), 1行m列
        %返回值lmMatrix, m行n列, 代表在不同轴向距离下, 光照度随径向距离变化的情况
        function [lmMatrix] = illuminanceCompute(obj, R, U, S, ...
                xPoints, hPoints)
            obj.R = R;
            obj.U = U;
            %根据当前的参数计算虚拟光源的高度
            obj.H = R / tan(U);
            %准备返回值并进行计算
            xNumber = size(xPoints, 2);
            hNumber = size(hPoints, 2);
            lmMatrix = zeros(hNumber, xNumber);
            for i = 1: hNumber
                for j = 1: xNumber
                    lmMatrix(i, j) = obj.lmCompute(hPoints(1, i), xPoints(1, j));
                end 
            end
            %乘以归一化系数
            coff = S / (pi * power(sin(U), 2));
            lmMatrix = lmMatrix * coff;
        end

    end
end