classdef OptTool < handle
%此类作为基类用于提供光纤探测理论计算的所有工具函数
%注意,中间所有关于角度的计算都是弧度制

    properties
        %折射率
        %冰实部
        INR;
        %冰虚部
        INI;
        %水实部
        WNR;
        %水虚部
        WNI;
        %空气实部
        ANR;
    end

    methods
        function obj = OptTool()
            %在构造函数中,需要读取冰水折射率实部和虚部的数据
            load('BaseData.mat');
            obj.INR = iceNR;
            obj.INI = iceNI;
            obj.WNR = waterNR;
            obj.WNI = waterNI;
            obj.ANR = airNR;
        end

        %此函数将角度转换为弧度
        function rad = angle2Rad(~, angle)
            rad = mod(angle, 360) * pi / 180;
        end

        %此函数将弧度转换为角度
        function angle = ran2Angle(~, rad)
            angle = mod(rad, 2 * pi) * 180 / pi;
        end
       
        %此函数根据入射角和折射角计算反射比
        function refCoff = ref(~, inTheta, outTheta)
            %需要排除特殊情况
            if inTheta == 0 && outTheta == 0
                return 
            end
            %inTheta为入射角,outTheta为折射角,refCoff为反射比
            sum = inTheta + outTheta;
            diff = inTheta - outTheta;
            sDiff = power(sin(diff), 2);
            sSum = power(sin(sum), 2);
            tDiff = power(tan(diff), 2);
            tSum = power(tan(sum), 2);
            refCoff = (sDiff / sSum + tDiff / tSum) / 2;
        end

        %此函数根据斯涅尔定律求出射角
        function outTheta = snell(~, refIn, refOut, inTheta) 
            %refIn为入射介质的折射率,refOut为出射介质的折射率
            %inTheta为入射角度,outTheta为出射角
            outTheta = asin(refIn *  sin(inTheta) / refOut);
        end

        %此函数根据波长和折射率虚部计算吸收系数
        function absCoff = absCompute(~, NI, lambda)
            %NI为折射率虚部, lambda为波长单位为nm
            absCoff = 4 * pi * NI *  1e9 / lambda;
        end

        %该函数用于寻找目标波段的折射率
        function [N] = findN(~, lambda, NData)
            %NData第一列为波段,单位nm,第二列为折射率
            %需要转换为um
            lambda = lambda * 1e-3;
            N = 0;
            %考虑不在范围内的异常情况
            if(NData(1, 1) > lambda || NData(end, 1) < lambda) 
                return;
            end
            %找到第一个大于当前波段的下标
            l = 1;
            r = size(NData, 1);
            while(l < r)
                m = floor((l + r) / 2);
                if(NData(m, 1) >= lambda) 
                    r = m;
                else 
                    l = m + 1;
                end
            end
            if(NData(r, 1) == lambda) 
                N = NData(r, 2);
                return;
            end
            %此时需要进行插值
            px0 = NData(r - 1, 1);
            py0 = NData(r - 1, 2);
            px1 = lambda;
            px2 = NData(r, 1);
            py2 = NData(r, 2);
            N = py0 + (px1 - px0) * (py2 - py0) / (px2 - px0);
        end

        %该函数用于计算入射光纤的偏斜光能够在光纤中传输的临界角
        %x为光源到目标纤芯的距离,r为光纤的半径
        %phi为光源到光纤中心点和光源到目标点连线的夹角
        %n为入射光所在介质的折射率,NA为光纤子午平面内的数值孔径
        function [theta] = criticalAngle(~, x, r, phi, n, NA)
            coff = sqrt(1 - power(x * sin(phi) / r, 2));
            theta = asin(NA / (n * coff));
        end

        %该函数用于判断离散点是否存在于光纤接收范围内
        %x1, y1为当前光纤纤芯的坐标位置
        %r为光纤半径
        %nowX, nowPhi为当前离散点的位置
        function [contains] = acJudg(~, x1, y1, r, nowX, nowPhi)
            x2 = nowX * cos(nowPhi);
            y2 = nowX * sin(nowPhi);
            dis = power(x1 - x2, 2) + power(y1 - y2, 2);
            contains = dis <= power(r, 2);
        end

    end
end