classdef LightSource < OptTool
%此类根据光纤信号端面的排布信息计算每根光纤的光通量、最小入射角、最大入射角
%光源所在的介质是空气,为简化计算,认为其折射率实部不随波段变化

    properties
        %光源所在介质的折射率
        nr = 1;
        %光源的光通量(lm)
        S;
        %最大出射角(rad)
        U;
        %光纤的半径(m)
        R;
        %光纤纤芯的折射率
        NF;
        %光纤子午面内的数值孔径
        NA;
        %dx(m),dphi(rad)代表网格的大小
        dx;
        dphi;
        %加快计算的集合
        rMap;
    end

    methods
        function obj = LightSource()
            %构造父类即可
            obj = obj@OptTool();
        end

        %根据距离x、高度H, 求解入射角度
        function inTheta = thetaCompute(~, x, H)
            inTheta = atan(x / (2 * H));
        end

        %根据高度H, 出射角度theta, 计算距离x
        function x = disCompute(~, theta, H)
            x = 2 * H * tan(theta);
        end

        %此函数计算dx = alpha * dtheta 中的alpha
        function alpha = alphaCompute(~, H, inTheta)
            alpha = 2 * H / power(cos(inTheta), 2);
        end

        %该函数根据光源波长lambda(nm)得到光源所在介质的折射率, lambda作为保留参数
        function [nr] = NCoffCompute(obj, lambda)
            if lambda > 0
                nr = obj.nr;
            end
        end

        %该函数计算单点的光通量
        %该函数根据距离x(m),网格大小dphi(rad),dx(m)
        %介质的折射率实部nr
        %距离H(m),最大出射角度U(rad), 源的光通量S(lm)
        %求得在dxdphi下的光通量dFlux(lm)
        function dFlux = dFluxCompute(obj, inTheta)
            %计算在光纤中的出射角
            outTheta = obj.snell(obj.nr, obj.NF, inTheta);
            %根据上述角度计算界面的反射系数
            r = obj.ref(inTheta, outTheta);
            %计算位置x处,dxdphi下的光通量
            dFlux = cos(inTheta) * sin(inTheta) * (1 - r);
        end

        %该函数计算光源和单根光纤之间的光通量、最小入射角、最大入射角
        %返回以dphi(rad), dx(m)为网格大小, 在厚为H(m)
        %最大出射角度为U(rad), 光源的光通量为S(lm), 光源波长为lambda(nm)
        %信号光纤的中心位置距离光源x(m), 位于光源相角为phi(rad) 
        %半径为R(m)的条件下的接收光通量
        function [flux, minTheta, maxTheta] = fluxCompute(obj, x, phi, H)
            %计算可能包含的网格点范围
            %当前光纤相对于光源的半角
            theta = asin(obj.R / x);
            %角度范围
            kphiLower = floor((phi - theta + 2 * pi) / obj.dphi);
            KphiUpper = ceil((phi + theta + 2 * pi) / obj.dphi);
            %长度范围
            kxLower = floor((x - obj.R) / obj.dx);
            %计算当前光纤的坐标位置
            x1 = x * cos(phi);
            y1 = x * sin(phi);
            %计算包含的网格点, 并将网格点中的dFlux累加
            flux = 0;
            minTheta = pi / 2;
            maxTheta = 0;
            %计算边界条件
            maxD = obj.disCompute(obj.U, H);
            kxUpper = ceil(min(x + obj.R, maxD) / obj.dx);
            for i = kxLower: 1: kxUpper
                nowX = obj.dx * i;
                %计算当前距离下的入射角
                inTheta = obj.thetaCompute(nowX, H);
                %在x位置固定后,dphi的变化不会影响光通量的计算,但是可能会不在临界角范围内
                %因此根据入射角计算是否在临界角内,计算统计结果
                count = 0;
                for j = kphiLower: 1: KphiUpper
                    nowPhi = obj.dphi * j;
                    x2 = nowX * cos(nowPhi);
                    y2 = nowX * sin(nowPhi);
                    dis = power(x1 - x2, 2) + power(y1 - y2, 2);
                    %不在光纤接收范围内直接跳过
                    if dis > power(obj.R, 2)
                        continue;
                    end
                    if obj.CA == true
                        %计算临界角
                        theta = obj.criticalAngle(x, obj.R, nowPhi - phi, obj.nr, obj.NA);
                        %超过临界角直接跳过
                        if inTheta > theta
                            continue;
                        end
                    end
                    minTheta = min(minTheta, inTheta);
                    maxTheta = max(maxTheta, inTheta);   
                    count = count + 1;
                end
                if count > 0
                    %只有该距离没有被计算过才需要计算
                    if ~isKey(obj.rMap, nowX)
                        obj.rMap(nowX) = obj.dFluxCompute(inTheta);
                    end
                    dflux = obj.rMap(nowX);
                    flux = flux + count * dflux;
                end
            end
        end

        %该函数计算任意的发光光源和信号光纤的组合下
        %每个信号光纤最终接收到的光通量,最大入射角,最大出射角
        %posMatrix为光源和接收光纤中心的距离矩阵, [a, b, c] 3个维度
        %a为1代表距离关系,a为2代表相角关系
        %b代表对应的光源,c代表对应的信号光纤
        %注意所有的光源都是统一参数的
        %H代表光源到信号光纤平面的距离
        %S代表所有光源的光通量总和(lm), U代表最大出射角(rad)
        %R代表接收光纤的半径(m)
        %dx(m),dphi(rad)代表网格的大小
        %返回值: flux,返回光通量,1行n列,列号是不同的信号光纤(lm)
        %angle,2行n列,第1行代表最小入射角,第2行代表最大入射角(rad)
        function [fluxs, angle] = fluxMatrixCompute(obj, posMatrix, ...
                H, S, U, R, NA, NF, CA, dx, dphi)
            %首先根据传入的参数设置计算时需要的参数
            obj.U = U;
            obj.R = R;
            obj.NA = NA;
            obj.NF = NF;
            obj.CA = CA;
            obj.dx = dx;
            obj.dphi = dphi;
            %分配返回值
            rNumber = size(posMatrix, 3);
            fluxs = zeros(1, rNumber);
            angle = zeros(2, rNumber);
            %计算单个光源的光通量
            sNumber = size(posMatrix, 2);
            obj.S = S / sNumber;
            %初始化集合
            obj.rMap = containers.Map("KeyType", 'double', "ValueType", 'double');
            %临时保存结果的矩阵, 3个维度, 1代表光通量, 2代表最小入射角, 3代表最大入射角
            tempMatrix = zeros(sNumber, rNumber, 3);
            %计算发射光纤和接收光纤两两之间的响应
            for i = 1: sNumber
                for j = 1: rNumber
                     [flux, minTheta, maxTheta] = obj.fluxCompute(posMatrix(1, i, j), ...
                        posMatrix(2, i, j), H);
                     tempMatrix(i, j, :) = [flux, minTheta, maxTheta];
                end
            end
            %由于每个光源的参数都是一样的,因此将共同的项放到最后进行计算
            tempMatrix(:, :, 1) = tempMatrix(:, :, 1) * obj.S ...
                    / (pi * power(sin(obj.U), 2)) * obj.dx * obj.dphi;
            %光通量直接合并为一行,为每根信号光纤接收到的所有光通量
            fluxs(1, :) = sum(tempMatrix(:, :, 1), 1);
            %最小入射角和最大入射角直接按列取最大值和最小值
            angle(1, :) = min(tempMatrix(:, :, 2), [], 1);
            angle(2, :) = max(tempMatrix(:, :, 3), [], 1);
        end

    end
end