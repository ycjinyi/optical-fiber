classdef OptCompute
%此类封装光学计算参数，并进行计算
    properties
        %网格参数
        %dphi(rad), dtheta(rad)为网格大小
        dphi = pi / 2000; 
        dtheta = pi / 2000;
        % dphi = pi / 800; 
        % dtheta = pi / 800;
        %光源参数
        %最大出射平面孔径角为U(rad), 光源总的光通量为S(lm)
        U = 30 * pi / 180;
        S = 1e8;
        %光纤参数
        %接收光纤的接收半径为R(m)
        R = 1e-6 * 187.562 / 2;
        %光纤纤芯的折射率
        NF = 1.445;
        %光纤子午面内的数值孔径
        NA = 0.37;
        %光源到光纤平面的高度(m)
        H = 1.5e-3;
        %是否限制光源入射的临界角
        CAS = true;
    end

    methods
        function obj = OptCompute()
        end

        %此函数用于计算设定波段范围和厚度范围内的结果
        %flux是接收光纤在各个波段的响应,hPoints是对应的介质厚度变化情况
        %posMatrix是光纤的排布参数
        %OC参数可以传入Bilayer和SingleLayer对象用于不同的计算
        %flux [a, b, c] 3维,a是波段,b是厚度点,c是接收光纤
        function [fluxs, ic, nc, rc] = compute(obj, OC, SPM, posMatrix, lambdas, hPoints)
            tic;
            ic = 0; 
            nc = 0;
            rc = 0;
            %需要计算的总波段数目
            pNumber = size(lambdas, 2);
            %厚度数目
            hNumber = size(hPoints, 2);
            %数据集合
            rNumber = size(posMatrix, 3);
            fluxs = zeros(pNumber, hNumber, rNumber);
            %需要计算的所有波段数目*厚度数目,用于输出计算进度
            points = pNumber * hNumber;
            %不区分波段先计算光源参数
            LS = LightSource();
            [flux, angle] = LS.fluxMatrixCompute(SPM, obj.H, obj.S, obj.U, obj.R, ...
                obj.NA, obj.NF, obj.CAS, obj.dtheta, obj.dphi);
            % % angle(1, :) = angle(1, :) * 0;
            % % angle(2, :) = ones(1, size(angle, 2)) * 0.4887;
            % [~, idx] = sort(SPM(1, 1, :));
            % figure;
            % % plot(flux(idx));
            % plot(flux, '*');
            % grid on;
            % % xlabel("按距离排序后的光纤编号");
            % xlabel("光纤编号");
            % ylabel("光通量");
            % figure;
            % plot(angle(1, idx) * 180 / pi); hold on;
            % plot(angle(2, idx) * 180 / pi);
            % % plot(angle(1, :) * 180 / pi); hold on;
            % % plot(angle(2, :) * 180 / pi);
            % legend("最小入射角", "最大入射角");
            % xlabel("按距离排序后的光纤编号");
            % % xlabel("光纤编号");
            % ylabel("入射角");
            % grid on;
            %目前去不考虑最小入射角的限制
            angle(1, :) = angle(1, :) * 0;
            %分波段进行计算
            for i = 1: pNumber
                %根据出射参数计算结果
                [fluxMatrix, ict, nct, rct] = ...
                    OC.fluxMatrixCompute(posMatrix, lambdas(1, i), hPoints,...
                    flux, angle, LS.nr, obj.NF, obj.NA, obj.R,...
                    obj.dtheta, obj.dphi, i, points);
                %保存当前波段下计算结果
                fluxs(i, :, :) = fluxMatrix;
                %累加统计值
                ic = ic + ict;
                nc = nc + nct;
                rc = rc + rct;
            end
            toc;
        end
    end
end