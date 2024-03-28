clc;
clear;
close all;

%此脚本用于计算单光纤对下的冰水混合以及波段数据展示

%<<<<<<<<<<<<<<<<<计算参数>>>>>>>>>>>>>>>>>>>>>
%天顶角和方位角的网格大小
dphi = pi / 1200; 
dtheta = pi / 1200;
%最大出射平面孔径角为U(rad), 光源总的光通量为S(lm)
U = 20 * pi / 180;
S = 1;
%接收光纤的接收半径为R(m)
R = 1e-6 * 90;

%需要使用的工具类
OT = OptTool();
OC = OptCompute();
FG = FiberGenerator();

%设置介质层属性
BL = BiLayer(OT.INR, OT.INI, OT.WNR, OT.WNI, OT.ANR);
SL = SingleLayer(OT.INR, OT.INI, OT.ANR);

%<<<<<<<<<<<<<<1、通过单层介质计算当前单光纤对的冰厚响应特性，并找到合适的厚度点>>>>>>>>
CG = ColorGenerator();
%介质厚度数据(m)
H1 = (0.01: 0.005: 3) * 1e-3;
%波段数据(nm)
lambdas = 890;
%光纤间距(um)
coreDisO = 300: 100: 1000;
pNum = size(coreDisO, 2);
%创建保存结果的矩阵
%行为不同的厚度点，列为不同距离下的结果
pic1Data = zeros(size(H1, 2), pNum);
% tic;
% for i = 1: pNum
%     [s, r1] = FG.singleFiberGen(coreDisO(1, i) * 1e-6);
%     %fluxs(波段,厚度,接收光纤)
%     [fluxs, ic, nc, rc] = OC.idealCompute(SL, FG.posConvert(s, r1), lambdas, H1,...
%                         U, S, R, dtheta, dphi);
%     pic1Data(:, i) = fluxs(1, :, :);
% end
% toc;
% save pic1Data pic1Data;
load pic1Data.mat;
[colorTable, lambdaStr] = CG.generate(coreDisO);
%作图展示
figure(1);
for i = 1: pNum
    plot(H1 * 1e3, pic1Data(:, i), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=1); hold on;
end
grid on;
legend(lambdaStr);
xlabel("冰厚度(mm)");
ylabel("光通量(lm)");
title("不同光纤距离下的响应");

%选择最终的纤芯间距(m)
finalCoreDis = 400 * 1e-6;
%峰值厚度点选择(m)
peakThick = 0.5 * 1e-3;

%<<<<<<<<<<<<<<<<<<2、总厚度不变,<冰水比例><波段>变化>>>>>>>>>>>>>>>>>>
CG = ColorGenerator([0 0 1; 0 1 1; 1 0 1; 1 1 0]);
lambdao = 780:10:1700;
H = peakThick;
coreDis = finalCoreDis;
[s, r1] = FG.singleFiberGen(coreDis);
[posMatrix] = FG.posConvert(s, r1);
%厚度中水的占比
ratio = 0.04:0.04:0.96;
lNum = size(lambdao, 2);
rNum = size(ratio, 2);
%创建保存结果的矩阵
%波段、不同水占比
pic2Data = zeros(lNum, rNum);
% %考虑冰和水混合情况
% tic;
% %水厚度
% Hw = H * ratio;
% %冰厚度
% Hi = H * (1 - ratio);
% %fluxs(波段,厚度,接收光纤)
% [fluxsB, ~] = OC.idealCompute(BL, posMatrix, lambdao, [Hi; Hw],...
%                     U, S, R, dtheta, dphi);
% pic2Data(:, :) = fluxsB(:, :, 1);
% toc;
% save pic2Data pic2Data;
load pic2Data.mat;
[colorTable, lambdaStr] = CG.generate(lambdao);
%作图展示
figure(2);
for i = 1: lNum
    plot(ratio * 100, pic2Data(i, :), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=1); hold on;
end
grid on;
legend(lambdaStr);
xlabel("水占比(%)");
ylabel("光通量(lm)");
xlim([0, 110]);
% set(gca, "YScale", "log");
title("冰水的总厚度为0.5mm");

ratio = round(ratio * 100);
a = pic2Data';

% <<<<<<<<<<<<<<<<<<<<<3、波段不变,<冰水比例><总厚度>变化>>>>>>>>>>>>>>>>>>>>
CG = ColorGenerator();
lambdao = [800, 1000, 1200, (1380: 5: 1450), (1460: 15: 1610) (1640: 30: 1700)];
%厚度中水的占比, 计算2个点
ratio = [0.1, 0.9];
pNum = size(ratio, 2);
lNum = size(lambdao, 2);
%厚度范围
H = (0.01: 0.01: 2) * 1e-3;
hNum = size(H, 2);
%创建保存结果的矩阵
%波段、厚度点、不同水占比
pic3Data = zeros(lNum, size(H, 2), pNum);
%单层模型，注意外层介质是水
SL = SingleLayer(OT.INR, OT.INI, OT.WNR);
%考虑冰和水混合情况
% tic;
% for i = 1: pNum
%     %水厚度
%     Hw = H * ratio(1, i);
%     %冰厚度
%     Hi = H * (1 - ratio(1, i));
%     %fluxs(波段,厚度,接收光纤)
%     [fluxsB, ~] = OC.idealCompute(BL, posMatrix, lambdao, [Hi; Hw],...
%                         U, S, R, dtheta, dphi);
%     % [fluxsS, ~] = OC.idealCompute(SL1, posMatrix, lambda, Hi,...
%     %                     U, S, R, dtheta, dphi);
%     % pic3Data(:, :, i) = fluxsB(:, :, 1) + fluxsS(:, :, 1);
%     pic3Data(:, :, i) = fluxsB(:, :, 1);
% end
% toc;
% save pic3Data pic3Data;
load pic3Data.mat;
diffData = pic3Data(:, :, 1) - pic3Data(:, :, 2);

[colorTable, lambdaStr] = CG.generate(lambdao);
%作图展示
%波段、厚度点、不同水占比
figure(3);
for i = 1: lNum
    plot(H * 1e3, diffData(i, :), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=1); hold on;
end
grid on;
legend(lambdaStr);
xlabel("总厚度(mm)");
xlim([0, 3]);
ylabel("光通量之差(lm)");
title("不同波段下的变化趋势");

H = H' * 1e3;
diffData = diffData';

% <<<<<<<<<<<<<<<<<<<<<4、波段不变,<冰水比例><总厚度>变化>>>>>>>>>>>>>>>>>>>>
lambdao = 800:1:1700;
%厚度中水的占比, 计算2个点
ratio = [0.1, 0.9];
pNum = size(ratio, 2);
lNum = size(lambdao, 2);
%厚度范围
H = [(0.4: 0.02: 0.58), (0.6: 0.1: 1.5)] * 1e-3;
hNum = size(H, 2);
%创建保存结果的矩阵
%波段、厚度点、不同水占比
pic4Data = zeros(lNum, size(H, 2), pNum);
%单层模型，注意外层介质是水
SL = SingleLayer(OT.INR, OT.INI, OT.WNR);
% %考虑冰和水混合情况
% tic;
% for i = 1: pNum
%     %水厚度
%     Hw = H * ratio(1, i);
%     %冰厚度
%     Hi = H * (1 - ratio(1, i));
%     %fluxs(波段,厚度,接收光纤)
%     [fluxsB, ~] = OC.idealCompute(BL, posMatrix, lambdao, [Hi; Hw],...
%                         U, S, R, dtheta, dphi);
%     % [fluxsS, ~] = OC.idealCompute(SL1, posMatrix, lambda, Hi,...
%     %                     U, S, R, dtheta, dphi);
%     % pic3Data(:, :, i) = fluxsB(:, :, 1) + fluxsS(:, :, 1);
%     pic4Data(:, :, i) = fluxsB(:, :, 1);
% end
% toc;
% save pic4Data pic4Data;
load pic4Data.mat;
diffData = pic4Data(:, :, 1) - pic4Data(:, :, 2);

[colorTable, lambdaStr] = CG.generate(H * 1e3);
figure(4);
%波段、厚度点、不同水占比
for i = 1: hNum
    plot(lambdao, diffData(:, i), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=1); hold on;
end
grid on;
legend(lambdaStr);
xlabel("波长(nm)");
ylabel("光通量之差(lm)");
title("不同厚度下的变化趋势");

lambdao = lambdao';
H = H* 1e3;