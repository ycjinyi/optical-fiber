classdef SingleLayer < OptTool
%此类继承自OptCompute类,并进行单层情况下接收光通量的计算

    properties
        %介质的折射率实部和虚部表
        NR1;
        NI1;
        %外部介质的折射率实部表
        NR2;

        %-----计算时需要使用到的参数-----

        %光源的波长(nm)
        lambda;
        %特定波长下的折射率数据
        nr1;
        ni1;
        nr2;
        %接收光纤的半径(m)
        R;
        %dx(m),dphi(rad)代表网格的大小
        dx;
        dphi;
        %加快计算的集合
        rMap;

        %-----计算时的统计数据-----
        %包含的点数
        incluCount;
        %需要计算的点数
        ncaluCount;
        %实际计算的点数
        rcaluCount;
    end

    methods
        function obj = SingleLayer(NR1, NI1, NR2)
            obj = obj@OptTool();
            %在子类的构造函数中决定介质的属性
            obj.NR1 = NR1;
            obj.NI1 = NI1;
            obj.NR2 = NR2;
        end

        %根据距离x、厚度h, 求解出射角度
        function inTheta = thetaCompute(~, x, h)
            inTheta = atan(x / (2 * h));
        end

        %根据厚度h, 出射角度theta, 计算距离x
        function x = disCompute(~, theta, h)
            x = 2 * h * tan(theta);
        end

        %此函数计算dx = alpha * dtheta 中的alpha
        function alpha = alphaCompute(~, h, inTheta)
            alpha = 2 * h / power(cos(inTheta), 2);
        end

        %该函数根据光源波长lambda(nm)得到介质的折射率实部和虚部
        function [nr1, ni1, nr2] = NCoffCompute(obj, lambda)
            nr1 = obj.findN(lambda, obj.NR1);
            ni1 = obj.findN(lambda, obj.NI1);
            nr2 = obj.findN(lambda, obj.NR2);
        end

        %该函数计算单点的光通量
        %该函数根据距离x(m), 网格大小dphi(rad),dx(m), 波长lambda(nm)
        %介质的折射率实部和虚部nr1, ni1, 外部介质折射率实部nr2
        %厚度h(m), 最大出射角度U(rad), 光源的光通量S(lm)
        %求得在dxdphi下的光通量dFlux(lm)
        function dFlux = dFluxCompute(obj, inTheta, h)
            %计算光线进入外部介质中的折射角
            airTheta = obj.snell(obj.nr1, obj.nr2, inTheta);
            %根据上述角度计算界面的反射系数
            r = obj.ref(inTheta, airTheta);
            %计算介质中的吸收系数
            k = obj.absCompute(obj.ni1, obj.lambda);
            %计算位置x处,dxdphi下的光通量
            coff1 = cos(inTheta) * sin(inTheta) * r;
            coff2 = exp(-2 * k * h / cos(inTheta));
            dFlux = coff1 * coff2;
        end

        %该函数计算单根光纤之间的光通量
        %返回以dphi(rad), dx(m)为网格大小, 在厚为h(m)
        %最大出射角度为U(rad), 光源的光通量为S(lm), 光源波长为lambda(nm)
        %接收光纤的中心位置距离光源x(m), 位于光源相角为phi(rad) 
        %接收半径为R(m)的条件下的接收光通量
        function flux = fluxCompute(obj, x, phi, minU, maxU, h)
            %计算可能包含的网格点范围
            %当前光纤相对于光源的半角
            theta = asin(obj.R / x);
            %角度范围
            kphiLower = floor((phi - theta + 2 * pi) / obj.dphi);
            KphiUpper = ceil((phi + theta + 2 * pi) / obj.dphi);
            %计算当前光纤的坐标位置
            x1 = x * cos(phi);
            y1 = x * sin(phi);
            %计算包含的网格点, 并将网格点中的dFlux累加
            flux = 0;
            %计算边界条件
            minD = obj.disCompute(minU, h);
            maxD = obj.disCompute(maxU, h);
            kxLower = floor(max(x - obj.R, minD) / obj.dx);
            kxUpper = ceil(min(x + obj.R, maxD) / obj.dx);
            for i = kxLower: 1: kxUpper
                nowX = obj.dx * i;
                %计算当前距离下的入射角
                inTheta = obj.thetaCompute(nowX, h);
                %在x位置固定后,dphi的变化不会影响光通量的计算
                %但是可能处于临界角度之外
                count = 0;
                % for j = kphiLower: 1: KphiUpper
                %     nowPhi = obj.dphi * j;
                %     x2 = nowX * cos(nowPhi);
                %     y2 = nowX * sin(nowPhi);
                %     dis = power(x1 - x2, 2) + power(y1 - y2, 2);
                %     if dis > power(obj.R, 2)
                %         continue;
                %     end
                %     count = count + 1;
                % end

                %查找在接收范围内的下界
                l = kphiLower;
                r = KphiUpper;
                while l < r
                    j = floor((l + r) / 2);
                    if obj.acJudg(x1, y1, obj.R, nowX, obj.dphi * j)
                        r = j;
                    else 
                        l = j + 1;
                    end
                end
                k1 = l;
                %查找在接收范围内的上界
                l = kphiLower;
                r = KphiUpper;
                while l < r
                    j = ceil((l + r) / 2);
                    if obj.acJudg(x1, y1, obj.R, nowX, obj.dphi * j)
                        l = j;
                    else 
                        r = j - 1;
                    end
                end
                k2 = l;
                if k1 <= k2
                    count = k2 - k1 + 1;
                end
                
                if count > 0
                    obj.incluCount = obj.incluCount + count;
                    obj.ncaluCount = obj.ncaluCount + 1;
                    %只有该距离没有被计算过才需要计算
                    if ~isKey(obj.rMap, nowX)
                        obj.rMap(nowX) = obj.dFluxCompute(inTheta, h);
                        obj.rcaluCount = obj.rcaluCount + 1;
                    end
                    dflux = obj.rMap(nowX);
                    flux = flux + count * dflux;
                end
            end
        end

        %该函数计算任意的发光光源和接收光纤的组合下,每个接收光纤最终接收到的光通量
        %posMatrix为光源和接收光纤中心的距离矩阵, [a, b, c] 3个维度
        %a为1代表距离关系，a为2代表相角关系
        %b代表对应的发射光纤,c代表对应的接收光纤
        %注意所有的光源都是统一参数的
        %lambda为光源的波长(nm)
        %S为1行n列,代表n个发射光纤的光通量(lm)
        %U代表n个光纤的最小入射角、最大入射角(rad),分别是第1、2行
        %n代表光源所在介质的折射率
        %R代表接收光纤的半径(m)
        %dx(m),dphi(rad)代表网格的大小
        %hPoints代表厚度点,1行n列
        %返回值fluxMatrix的列代表不同接收光纤的光通量,行代表不同介质厚度下的情况
        %fluxMatrix的列数为接收光纤的个数,行数为上层厚度的情况数目
        %idx表示当前计算的波段序号, points代表所有需要计算的波段数目*厚度数目,用于输出计算进度
        function [fluxMatrix, ic, nc, rc] = fluxMatrixCompute(obj, posMatrix,...
                lambda, hPoints, S, U, n, R, dx, dphi, idx, points, ideal)
            obj.incluCount = 0;
            obj.ncaluCount = 0;
            obj.rcaluCount = 0;
            %首先根据传入的参数设置计算时需要的参数
            obj.lambda = lambda;
            obj.R = R;
            obj.dx = dx;
            obj.dphi = dphi;
            %计算光源的参数
            sNumber = size(posMatrix, 2);
            %先根据lambda求出介质的折射率实部和虚部
            [obj.nr1, obj.ni1, obj.nr2] = obj.NCoffCompute(obj.lambda);
            if ~ideal
                %将入射角度数据转换为出射角度数据
                for i = 1: size(U, 1)
                    for j = 1: size(U, 2)
                        U(i, j) = obj.snell(n, obj.nr1, U(i, j));
                    end
                end
            end
            %分配返回值
            hNumber = size(hPoints, 2);
            rNumber = size(posMatrix, 3);
            fluxMatrix = zeros(hNumber, rNumber);
            %当前已经计算了的点数
            count = (idx - 1) * hNumber;
            %分别计算不同厚度下的情况
            for H = 1: hNumber
                %计算数据只在相同波段下的相同厚度下生效，厚度变化后将失效
                %因此在每次计算前需要重置
                obj.rMap = containers.Map("KeyType", 'double', "ValueType", 'double');
                %临时保存结果的矩阵
                tempMatrix = zeros(sNumber, rNumber);
                %计算发射光纤和接收光纤两两之间的响应
                for i = 1: sNumber
                    if S(1, i) == 0
                        continue;
                    end
                    for j = 1: rNumber
                        tempMatrix(i, j) = obj.fluxCompute(posMatrix(1, i, j), ...
                            posMatrix(2, i, j), U(1, i), U(2, i), hPoints(1, H));
                    end
                    %根据不同光源的属性调整光通量
                    tempMatrix(i, :) = tempMatrix(i, :) * S(1, i) ...
                        / (pi * (power(sin(U(2, i)), 2) - power(sin(U(1, i)), 2))) ...
                        * obj.dx * obj.dphi;
                end
                %直接合并为一行,为每根接收光纤接收到的所有光通量
                fluxMatrix(H, :) = sum(tempMatrix, 1);
                %输出当前计算的进度
                sprintf(strcat("--> ", num2str(floor((count + H)* 1000 / points) / 10), "%% <--"))
            end
            %返回统计结果
            ic = obj.incluCount;
            nc = obj.ncaluCount;
            rc = obj.rcaluCount;
        end

    end
end