classdef OptCompute
%此类封装计算模型, 分别考虑耦合和不耦合光源的情况
    properties
    end

    methods
        function obj = OptCompute()
        end

        %此函数用于计算设定波段范围和厚度范围内的结果，是考虑了光源的模型
        %flux是接收光纤在各个波段的响应,hPoints是对应的介质厚度变化情况
        %sourcePosMatrix是光源和发射光纤之间的排布参数
        %posMatrix有3个维度[a, b, c]
        %a代表对应的发射光纤编号,b代表对应的接收光纤编号
        %c的取值是1和2,1代表距离之间的关系,2代表角度之间的关系,注意是弧度制
        %OC参数可以传入Bilayer和SingleLayer对象用于不同的计算
        %flux [a, b, c] 3维,a是波段,b是厚度点,c是接收光纤
        %hPoints是需要计算的厚度点数据
        %H是光源到光纤平面的高度(m)
        %V代表光源的视场角半角(rad)
        %S代表光源的光通量(lm)
        %R代表光纤半径(R)
        %NA代表光纤的数值孔径
        %NF代表光纤纤芯的折射率
        %NS代表光源所在介质的折射率, 主要用于临界角和反射比的计算
        %CA代表是否添加光线入射的临界角约束
        %dphi(rad), dtheta(rad)为网格大小
        function [fluxs, ic, nc, rc] = couplingCompute(~, OC, ...
                sourcePosMatrix, posMatrix, lambdas, hPoints, H, V, S,...
                R, NA, NF, NS, CA, dtheta, dphi)
            ic = 0; 
            nc = 0;
            rc = 0;
            %需要计算的总波段数目
            pNumber = size(lambdas, 2);
            %厚度数目
            hNumber = size(hPoints, 2);
            %数据集合
            rNumber = size(posMatrix, 2);
            fluxs = zeros(pNumber, hNumber, rNumber);
            %需要计算的所有波段数目*厚度数目,用于输出计算进度
            points = pNumber * hNumber;
            %计算光源参数,不需要区分波段
            LS = LightSource();
            [flux, angle] = LS.fluxMatrixCompute(sourcePosMatrix,...
                H, S, V, R, NA, NF, NS, CA, dtheta, dphi);
            %目前不考虑最小入射角的限制
            angle(1, :) = angle(1, :) * 0;
            %分波段进行计算
            for i = 1: pNumber
                %根据出射参数计算结果
                [fluxMatrix, ict, nct, rct] = ...
                    OC.fluxMatrixCompute(posMatrix, lambdas(1, i), ...
                    hPoints, flux, angle, NS, false, R,...
                    dtheta, dphi, i, points);
                %保存当前波段下计算结果
                fluxs(i, :, :) = fluxMatrix;
                %累加统计值
                ic = ic + ict;
                nc = nc + nct;
                rc = rc + rct;
            end
        end

        %此函数用于计算设定波段范围和厚度范围内的结果，不考虑光源耦合模型
        %flux是接收光纤在各个波段的响应,hPoints是对应的介质厚度变化情况
        %posMatrix有3个维度[a, b, c]
        %a代表对应的发射光纤编号,b代表对应的接收光纤编号
        %c的取值是1和2,1代表距离之间的关系,2代表角度之间的关系,注意是弧度制
        %OC参数可以传入Bilayer和SingleLayer对象用于不同的计算
        %flux [a, b, c] 3维,a是波段,b是厚度点,c是接收光纤
        %R代表接收光纤的半径(m)
        %dtheta(rad),dphi(rad)代表网格的大小
        function [fluxs, ic, nc, rc] = idealCompute(~, OC, posMatrix, ...
                lambdas, hPoints, U, S, R, dtheta, dphi)
            ic = 0; 
            nc = 0;
            rc = 0;
            %需要计算的总波段数目
            pNumber = size(lambdas, 2);
            %厚度数目
            hNumber = size(hPoints, 2);
            %数据集合
            rNumber = size(posMatrix, 2);
            fluxs = zeros(pNumber, hNumber, rNumber);
            %需要计算的所有波段数目*厚度数目,用于输出计算进度
            points = pNumber * hNumber;
            %理想光源, 所有发射光纤的最小出射角为0, 最大出射角统一
            sNumber = size(posMatrix, 1);
            angle = ones(2, sNumber) * U;
            angle(1, :) = angle(1, :) * 0;
            flux = ones(1, sNumber) * S / sNumber;
            %分波段进行计算
            for i = 1: pNumber
                %根据出射参数计算结果
                %由于不需要考虑光源耦合,因此不需要转换角度,因而光源所在介质的折射率设置为1即可
                [fluxMatrix, ict, nct, rct] = ...
                    OC.fluxMatrixCompute(posMatrix, lambdas(1, i), hPoints,...
                    flux, angle, 1, true, R,...
                    dtheta, dphi, i, points);
                %保存当前波段下计算结果
                fluxs(i, :, :) = fluxMatrix;
                %累加统计值
                ic = ic + ict;
                nc = nc + nct;
                rc = rc + rct;
            end
        end

    end
end