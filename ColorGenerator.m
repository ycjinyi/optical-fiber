classdef ColorGenerator
%颜色生成器 

    properties
        %基准颜色
        %baseColor = [1 0 0; 0 1 0; 0 0 1; 0 1 1; 1 0 1; 1 1 0];
        baseColor = [0 0 1; 0 1 1; 1 0 1];
    end

    methods
        function obj = ColorGenerator()
            %构造函数什么也不需要做
        end

        %此函数根据目标数据个数产生对应的颜色表和标签
        function [colorTable, lambdaStr] = generate(obj, data)
            %颜色表
            colorNum = size(data, 2);
            baseNum = size(obj.baseColor, 1);
            colorStep = (baseNum - 1) / colorNum;
            colorTable = zeros(colorNum, 3);
            for i = 1: colorNum
                px = floor(i * colorStep) + 1;
                delta = i * colorStep + 1 - px;
                if(delta == 0)
                    colorTable(i, :) = obj.baseColor(px, :);
                    continue;
                end
                colorTable(i, :) = delta * (obj.baseColor(px + 1, :) - obj.baseColor(px, :))...
                    + obj.baseColor(px, :);
            end
            %标签
            lambdaStr = cell(1, colorNum);
            for i = 1: colorNum
                lambdaStr(1, i) = {mat2str(data(1, i))};
            end
        end
    end
end