clc;
clear;
close all;

%此脚本只计算光纤探测模型, 不考虑光源耦合模型

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

%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<单层介质,单对光纤>>>>>>>>>>>>>>>>>>>>>>>>>>>
%---------------1、890nm/距离/出射角/介质不变, <半径>变化--------------
ro = 40: 10: 100;
r = ro * 1e-6;
coreDis = 1e-6 * 600;
[s, r1] = FG.singleFiberGen(coreDis);
[posMatrix] = FG.posConvert(s, r1);
pNum = size(ro, 2);
%创建保存结果的矩阵
%行为不同的厚度点，列为不同半径下的结果
pic1Data = zeros(size(H1, 2), pNum);
% tic;
% for i = 1: pNum
%     %fluxs(波段,厚度,接收光纤)
%     [fluxs, ic, nc, rc] = OC.idealCompute(SL, posMatrix, lambdas, H1,...
%                         U, S, r(1, i), dtheta, dphi);
%     pic1Data(:, i) = fluxs(1, :, 1);
% end
% toc;
% save pic1Data pic1Data;
load pic1Data.mat;
[colorTable, lambdaStr] = CG.generate(ro);
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
title("不同光纤半径下的响应");

%---------------2、890nm/半径/出射角/介质不变, <距离>变化--------------
coreDisO = 300: 100: 1000;
pNum = size(coreDisO, 2);
%创建保存结果的矩阵
%行为不同的厚度点，列为不同距离下的结果
pic2Data = zeros(size(H1, 2), pNum);
% tic;
% for i = 1: pNum
%     [s, r1] = FG.singleFiberGen(coreDisO(1, i) * 1e-6);
%     %fluxs(波段,厚度,接收光纤)
%     [fluxs, ic, nc, rc] = OC.idealCompute(SL, FG.posConvert(s, r1), lambdas, H1,...
%                         U, S, R, dtheta, dphi);
%     pic2Data(:, i) = fluxs(1, :, 1);
% end
% toc;
% save pic2Data pic2Data;
load pic2Data.mat;
[colorTable, lambdaStr] = CG.generate(coreDisO);
%作图展示
figure(2);
for i = 1: pNum
    plot(H1 * 1e3, pic2Data(:, i), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=1); hold on;
end
grid on;
legend(lambdaStr);
xlabel("冰厚度(mm)");
ylabel("光通量(lm)");
title("不同光纤距离下的响应");

%---------------3、890nm/距离/半径/介质不变, <出射角>变化--------------
uo = 12: 3: 42;
u = uo * pi / 180;
coreDis = 1e-6 * 600;
[s, r1] = FG.singleFiberGen(coreDis);
[posMatrix] = FG.posConvert(s, r1);
pNum = size(uo, 2);
%创建保存结果的矩阵
%行为不同的厚度点，列为不同出射角下的结果
pic3Data = zeros(size(H1, 2), pNum);
% tic;
% for i = 1: pNum
%     %fluxs(波段,厚度,接收光纤)
%     [fluxs, ic, nc, rc] = OC.idealCompute(SL, posMatrix, lambdas, H1,...
%                         u(1, i), S, R, dtheta, dphi);
%     pic3Data(:, i) = fluxs(1, :, 1);
% end
% toc;
% save pic3Data pic3Data;
load pic3Data.mat;
[colorTable, lambdaStr] = CG.generate(uo);
%作图展示
figure(3);
for i = 1: pNum
    plot(H1 * 1e3, pic3Data(:, i), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=1); hold on;
end
grid on;
legend(lambdaStr);
xlabel("冰厚度(mm)");
ylabel("光通量(lm)");
title("不同最大出射角下的响应");

%---------------4、距离/半径/出射角/介质不变, <波段>变化---------------
lambdao = 1465:10:1565;
coreDis = 1e-6 * 600;
[s, r1] = FG.singleFiberGen(coreDis);
[posMatrix] = FG.posConvert(s, r1);
pNum = size(lambdao, 2);
%创建保存结果的矩阵
%行为不同的厚度点，列为不同波段下的结果
pic4Data = zeros(size(H1, 2), pNum);
% tic;
% for i = 1: pNum
%     %fluxs(波段,厚度,接收光纤)
%     [fluxs, ic, nc, rc] = OC.idealCompute(SL, posMatrix, lambdao(1, i), H1,...
%                         U, S, R, dtheta, dphi);
%     pic4Data(:, i) = fluxs(1, :, 1);
% end
% toc;
% save pic4Data pic4Data;
load pic4Data.mat;
[colorTable, lambdaStr] = CG.generate(lambdao);
%作图展示
figure(4);
for i = 1: pNum
    plot(H1 * 1e3, pic4Data(:, i) + 1e-7, 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=1); hold on;
end
grid on;
legend(lambdaStr);
xlabel("冰厚度(mm)");
ylabel("光通量(lm)");
% set(gca, "YScale", "log");
title("不同波段下的响应");


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<单层介质,多对光纤>>>>>>>>>>>>>>>>>>>>>>>>>>>




