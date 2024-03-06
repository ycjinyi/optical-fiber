function [fitresult, gof] = surfFit(X, Y, Z)

%创建一个拟合
[xData, yData, zData] = prepareSurfaceData(X, Y, Z);

% 设置 fittype 和选项。
ft = fittype( 'loess' );
opts = fitoptions( 'Method', 'LowessFit' );
opts.Normalize = 'on';
opts.Span = 0.35;

% 对数据进行模型拟合。
[fitresult, gof] = fit([xData, yData], zData, ft, opts);




