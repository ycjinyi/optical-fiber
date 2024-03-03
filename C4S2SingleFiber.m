clc;
clear;
close all;

%此脚本用于计算单光纤对情况下的多谱数据

%计算参数
%dphi(rad), dtheta(rad)为网格大小
dphi = pi / 1800; 
dtheta = pi / 1800;
%最大出射平面孔径角为U(rad), 光源总的光通量为S(lm)
U = 20 * pi / 180;
S = 1;
%接收光纤的接收半径为R(m)
R = 1e-6 * 90;
%波段数据(nm)
lambdas = 890;
%介质厚度数据(m)
H1 = (0.01: 0.005: 3) * 1e-3;

CG = ColorGenerator();
OC = OptCompute();
%设置介质层属性
OT = OptTool();
% BL = BiLayer(OT.INR, OT.INI, OT.WNR, OT.WNI, OT.ANR);
SL = SingleLayer(OT.INR, OT.INI, OT.ANR);

FG = FiberGenerator();

%---------------6、波段不变,<冰水比例><总厚度>变化---------------
lambda = [890, 1500];
coreDis = 1e-6 * 600;
[s, r1] = FG.singleFiberGen(coreDis);
[posMatrix] = FG.posConvert(s, r1);
H = (0.01: 0.01: 3) * 1e-3;
%厚度中水的占比
ratio = 0:0.1:1;
pNum = size(ratio, 2);
lNum = size(lambda, 2);
%创建保存结果的矩阵
%波段、厚度点、不同水占比
pic6Data = zeros(lNum, size(H, 2), pNum);
%考虑冰和水混合情况
BL1 = BiLayer(OT.INR, OT.INI, OT.WNR, OT.WNI, OT.ANR);
% SL1 = SingleLayer(OT.INR, OT.INI, OT.WNR);
tic;
for i = 1: pNum
    %水厚度
    Hw = H * ratio(1, i);
    %冰厚度
    Hi = H * (1 - ratio(1, i));
    %fluxs(波段,厚度,接收光纤)
    [fluxsB, ~] = OC.idealCompute(BL1, posMatrix, lambda, [Hi; Hw],...
                        U, S, R, dtheta, dphi);
    % [fluxsS, ~] = OC.idealCompute(SL1, posMatrix, lambda, Hi,...
    %                     U, S, R, dtheta, dphi);
    % pic6Data(:, :, i) = fluxsB(:, :, 1) + fluxsS(:, :, 1);
    pic6Data(:, :, i) = fluxsB(:, :, 1);
end
toc;
save pic6Data pic6Data;
load pic6Data.mat;
[colorTable, lambdaStr] = CG.generate(ratio);
%作图展示
figure(6);
for i = 1: pNum
    plot(H * 1e3, pic6Data(1, :, i), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=1); hold on;
end
grid on;
legend(lambdaStr);
xlabel("总厚度(mm)");
ylabel("光通量(lm)");
% set(gca, "YScale", "log");
title("890nm,冰水混合,不同水占比下的响应");

%作图展示
figure(7);
for i = 1: pNum
    plot(H * 1e3, pic6Data(2, :, i), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=1); hold on;
end
grid on;
legend(lambdaStr);
xlabel("总厚度(mm)");
ylabel("光通量(lm)");
% set(gca, "YScale", "log");
title("1500nm,冰水混合,不同水占比下的响应");
%---------------7、总厚度不变,<冰水比例><波段>变化---------------
lambdao = 1405:10:1565;
H = [0.8, 1.6] * 1e-3;
coreDis = 1e-6 * 600;
[s, r1] = FG.singleFiberGen(coreDis);
[posMatrix] = FG.posConvert(s, r1);
%厚度中水的占比
ratio = 0:0.05:1;
lNum = size(lambdao, 2);
hNum = size(H, 2);
rNum = size(ratio, 2);
%创建保存结果的矩阵
%总厚度、波段、不同水占比
pic7Data = zeros(hNum, lNum, rNum);
%考虑冰和水混合情况
BL1 = BiLayer(OT.INR, OT.INI, OT.WNR, OT.WNI, OT.ANR);
SL1 = SingleLayer(OT.INR, OT.INI, OT.WNR);
tic;
for i = 1: hNum
    %水厚度
    Hw = H(1, i) * ratio;
    %冰厚度
    Hi = H(1, i) * (1 - ratio);
    %fluxs(波段,厚度,接收光纤)
    [fluxsB, ~] = OC.idealCompute(BL1, posMatrix, lambdao, [Hi; Hw],...
                        U, S, R, dtheta, dphi);
    % [fluxsS, ~] = OC.idealCompute(SL1, posMatrix, lambdao, Hi,...
    %                     U, S, R, dtheta, dphi);
    % pic7Data(i, :, :) = fluxsB(:, :, 1) + fluxsS(:, :, 1);
    pic7Data(i, :, :) = fluxsB(:, :, 1);
end
toc;
save pic7Data pic7Data;
load pic7Data.mat;
[colorTable, lambdaStr] = CG.generate(lambdao);
%作图展示
figure(8);
for i = 1: lNum
    plot(ratio * 100, squeeze(pic7Data(1, i, :)), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=1); hold on;
end
grid on;
legend(lambdaStr);
xlabel("水占比(%)");
ylabel("光通量(lm)");
% set(gca, "YScale", "log");
title("总厚度0.8mm,冰水混合,不同波段下的响应");

%作图展示
figure(9);
for i = 1: lNum
    plot(ratio * 100, squeeze(pic7Data(2, i, :)), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=1); hold on;
end
grid on;
legend(lambdaStr);
xlabel("水占比(%)");
ylabel("光通量(lm)");
% set(gca, "YScale", "log");
title("总厚度1.6mm,冰水混合,不同波段下的响应");

