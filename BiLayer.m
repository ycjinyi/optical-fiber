classdef BiLayer < OptTool
%此类继承自OptCompute类,并进行双层介质情况下接收光通量的计算

    properties
        %下层介质的折射率实部和虚部
        NR1;
        NI1;
        %上层介质的折射率实部和虚部
        NR2;
        NI2;
        %外部介质的折射率实部
        NR3;

        %-----计算时需要使用到的参数-----

        %光源的波长(nm)
        lambda;
        %特定波长下的折射率数据
        nr1;
        ni1;
        nr2;
        ni2;
        nr3;
        %单个光源的光通量(lm)
        S;
        %最大出射角(rad)
        U;
        %接收光纤的半径(m)
        R;
        %dtheta(rad),dphi(rad)代表网格的大小
        dtheta;
        dphi;
        %加快计算的集合,rMap用于光通量求解
        rMap;
        %等效光源到光纤平面的距离d
        d;

        %-----计算时的统计数据-----
        %包含的点数
        incluCount;
        %需要计算的点数
        ncaluCount;
        %实际计算的点数
        rcaluCount;
    end

    methods
        function obj = BiLayer(NR1, NI1, NR2, NI2, NR3)
            obj = obj@OptTool();
            %在子类的构造函数中决定介质属性
            obj.NR1 = NR1;
            obj.NR2 = NR2;
            obj.NI1 = NI1;
            obj.NI2 = NI2;
            obj.NR3 = NR3;
        end

        %根据距离x、下层厚度H1、上层厚度H2
        %求上下层介质界面的入射角度和折射角度
        function [inTheta, outTheta] = thetaCompute(obj, x, H1, H2)
            %二分查找
            lower = 0;
            upper = pi / 2;
            %初始值
            inTheta = 0;
            outTheta = 0;
            while lower < upper
                inTheta = (lower + upper) / 2;
                outTheta = obj.snell(obj.nr1, obj.nr2, inTheta);
                nowDis = obj.disCompute(H1, H2, inTheta, outTheta);
                %误差限比网格精度小5倍
                errLim = obj.alphaCompute(H1, H2, inTheta, outTheta) * obj.dtheta / 5;
                if(nowDis - errLim > x)
                    upper = inTheta;
                elseif(nowDis + errLim < x)
                    lower = inTheta;
                else 
                    break;
                end
            end
        end

        %此函数根据H1,H2,inTheta,outTheta求距离x
        function x = disCompute(obj, H1, H2, inTheta, outTheta)
            x = (2 * H1 + obj.d) * tan(inTheta) + 2* H2 * tan(outTheta);
        end

        %此函数计算dx = alpha * dtheta 中的alpha
        function alpha = alphaCompute(obj, H1, H2, inTheta, outTheta)
            coff1 = (2 * H1 + obj.d) / power(cos(inTheta), 2);
            coff2 = 2 * H2 * obj.nr1 * cos(inTheta) / (obj.nr2 * power(cos(outTheta), 3));
            alpha = coff1 + coff2;
        end

        %该函数根据光源波长lambda(nm)得到介质的折射率实部和虚部
        function [nr1, ni1, nr2, ni2, nr3] = NCoffCompute(obj, lambda)
            nr1 = obj.findN(lambda, obj.NR1);
            ni1 = obj.findN(lambda, obj.NI1);
            nr2 = obj.findN(lambda, obj.NR2);
            ni2 = obj.findN(lambda, obj.NI2);
            nr3 = obj.findN(lambda, obj.NR3);
        end

        %该函数计算单点的光通量
        %该函数根据距离x(m), 网格大小dphi(rad),dtheta(rad), 波长lambda(nm)
        %上下层介质的折射率实部和虚部nr1, ni1, nr2, ni2, 外部介质折射率实部nr3
        %上层厚度H1(m), 下层厚度H2(m)
        %最大出射角度U(rad), 光源的光通量S(lm)
        %求得在dxdphi下的光通量dFlux(lm)
        function dFlux = dFluxCompute(obj, inTheta, outTheta, H1, H2)
            %计算光线进入外部介质中的折射角
            airTheta = obj.snell(obj.nr2, obj.nr3, outTheta);
            %根据上述角度计算3个界面的反射系数
            r1 = obj.ref(inTheta, outTheta);
            r2 = obj.ref(outTheta, airTheta);
            r3 = obj.ref(outTheta, inTheta);
            %计算下层和上层介质中的吸收系数
            k1 = obj.absCompute(obj.ni1, obj.lambda);
            k2 = obj.absCompute(obj.ni2, obj.lambda);
            %计算位置x处,dxdphi下的光通量
            coff1 = cos(inTheta) * sin(inTheta) * (1 - r1) * (1 - r3) * r2;
            coff2 = exp(-2 * k1 * H1 / cos(inTheta)) * ...
                exp(-2 * k2 * H2 / cos(outTheta));
            dFlux = coff1 * coff2;
        end

        %该函数计算单根光纤之间的光通量
        %返回以dphi(rad), dtheta(rad)为网格大小
        %在下层厚为H1(m), 上层厚为H2(m)
        %最大出射角度为U(rad), 光源的光通量为S(lm), 光源波长为lambda(nm)
        %接收光纤的中心位置距离光源x(m), 位于光源相角为phi(rad) 
        %接收半径为R(m)的条件下的接收光通量
        function flux = fluxCompute(obj, x, phi, minU, maxU, H1, H2)
            %计算可能包含的网格点范围
            %当前光纤相对于光源的半角
            theta = asin(obj.R / x);
            %方位角范围
            kphiLower = floor((phi - theta + 2 * pi) / obj.dphi);
            KphiUpper = ceil((phi + theta + 2 * pi) / obj.dphi);
            %天顶角范围
            [inTheta1, ~] = obj.thetaCompute(x - obj.R, H1, H2);
            [inTheta2, ~] = obj.thetaCompute(x + obj.R, H1, H2);
            lthetaLower = floor(inTheta1 / obj.dtheta);
            lthetaUpper = ceil(inTheta2 / obj.dtheta);
            %计算当前光纤的坐标位置
            x1 = x * cos(phi);
            y1 = x * sin(phi);
            %计算包含的网格点, 并将网格点中的dFlux累加
            flux = 0;
            %计算边界条件
            lmin = ceil(minU / obj.dtheta);
            lmax = floor(maxU / obj.dtheta);
            lthetaLower = max(lmin, lthetaLower);
            lthetaUpper = min(lmax, lthetaUpper);
            for i = lthetaLower: lthetaUpper
                %计算当前的入射角和折射角
                inTheta = i * obj.dtheta;
                outTheta = obj.snell(obj.nr1, obj.nr2, inTheta);
                %当前的距离
                nowX = obj.disCompute(H1, H2, inTheta, outTheta);
                %在x位置固定后,dphi的变化不会影响光通量的计算
                count = 0;
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
                    %只有该天顶角没有被计算过才需要计算
                    if ~isKey(obj.rMap, i)
                        obj.rMap(i) = obj.dFluxCompute(inTheta, outTheta, H1, H2);
                        obj.rcaluCount = obj.rcaluCount + 1;
                    end
                    dflux = obj.rMap(i);
                    flux = flux + count * dflux;
                end
            end
        end

        %该函数计算任意的发光光源和接收光纤的组合下,每个接收光纤最终接收到的光通量
        %posMatrix有3个维度[a, b, c]
        %a代表对应的发射光纤编号,b代表对应的接收光纤编号
        %c的取值是1和2,1代表距离之间的关系,2代表角度之间的关系,注意是弧度制
        %注意所有的光源都是统一参数的
        %lambda为光源的波长(nm)
        %S为1行n列,代表n个发射光纤的光通量(lm)
        %U代表n个光纤的最小入射角、最大入射角(rad),分别是第1、2行
        %n代表光源所在介质的折射率
        %ideal代表是否是理想情况，此时不考虑光源耦合，不需要转换角度
        %R代表接收光纤的半径(m)
        %dtheta(rad),dphi(rad)代表网格的大小
        %hPoints代表厚度点,,2行n列,第一行为介质1的厚度,第二行代表介质2的厚度
        %返回值fluxMatrix的列代表不同接收光纤的光通量,行代表不同介质1、2厚度组合下的情况
        %fluxMatrix的列数为接收光纤的个数,行数为不同介质厚度组合下的结果
        %idx表示当前计算的波段序号, points代表所有需要计算的波段数目*厚度数目,用于输出计算进度
        function [fluxMatrix, ic, nc, rc] = fluxMatrixCompute(obj, posMatrix, lambda,...
                hPoints, S, U, n, ideal, R, dtheta, dphi, idx, points)
            obj.incluCount = 0;
            obj.ncaluCount = 0;
            obj.rcaluCount = 0;
            %首先根据传入的参数设置计算时需要的参数
            obj.lambda = lambda;
            obj.R = R;
            obj.dtheta = dtheta;
            obj.dphi = dphi;
            %先根据lambda求出下层和上层介质的折射率实部和虚部
            [obj.nr1, obj.ni1, obj.nr2, obj.ni2, obj.nr3] = obj.NCoffCompute(obj.lambda);
            %计算光源参数
            sNumber = size(posMatrix, 1);
            %将入射角度数据转换为出射角度数据
            if ~ideal
                for i = 1: size(U, 1)
                    for j = 1: size(U, 2)
                        U(i, j) = obj.snell(n, obj.nr1, U(i, j));
                    end
                end
            end
            %分配返回值
            hNumber = size(hPoints, 2);
            rNumber = size(posMatrix, 2);
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
                    %每个发射光纤的等效光源位置需要更新
                    obj.d = obj.R / tan(U(2, i));
                    for j = 1: rNumber
                        tempMatrix(i, j) = obj.fluxCompute(posMatrix(i, j, 1), ...
                            posMatrix(i, j, 2), U(1, i), U(2, i), hPoints(1, H), hPoints(2, H));
                    end
                    tempMatrix(i, :) = tempMatrix(i, :)  * S(1, i) ...
                        / (pi * (power(sin(U(2, i)), 2) - power(sin(U(1, i)), 2))) ...
                        * obj.dtheta * obj.dphi;
                end
                %直接在按列方向合并为一行，为每根接收光纤接收到的所有光通量
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