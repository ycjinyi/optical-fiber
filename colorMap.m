function [colorIdx, colorTable] = colorMap(data, value)
    %数据的最大值
    dataMax = max(data);
    %数据的最小值
    dataMin = min(data);
    %数据范围
    range = dataMax - dataMin;
    %计算系数
    coff = 1 / range;
    %乘以系数
    data = (data - dataMin) * coff;
    %颜色映射
    colorIdx = round(data * (value - 1)) + 1;  
    colorIdx = min(colorIdx, 255);
    colorIdx = max(colorIdx, 1);
    colorIdx = 256 - colorIdx;
    %颜色表
    % CG = ColorGenerator([0 1 1; 1 0 1; 1 1 0]);
    CG = ColorGenerator([0 1 1; 1 0 1]);
    [color, ~] = CG.generate(1: 1: 255);
    colorTable = color(colorIdx, :);
    colorTable = colorTable * 255;
end
