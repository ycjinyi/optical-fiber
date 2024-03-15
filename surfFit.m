function [fitresult, gof] = surfFit(X, Y, Z, span)

% 设置 fittype 和选项。
ft = fittype( 'loess' );
opts = fitoptions( 'Method', 'LowessFit' );
opts.Span = span;

% 对数据进行模型拟合。
[fitresult, gof] = fit([X, Y], Z, ft, opts);




