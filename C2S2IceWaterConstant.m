clc;
clear;
close all;

%生成论文中冰水的折射率实部、吸收系数以及反射比的关系

%加载数据
load BaseData.mat;

%波段范围，单位nm
b = 21;
iceNI = iceNI(b: end, :);
iceNR = iceNR(b: end, :);
waterNI = waterNI(b: end, :);
waterNR = waterNR(b: end, :);
airNR = airNR(b: end, :);

targetRange = waterNR(:, 1) * 1000;

%颜色表和标签
CG = ColorGenerator();
[colorTable, ~] = CG.generate(zeros(1, 17));

sidx = 3;
eidx = 16;

%----------------------------冰水的折射率实部-------------------------------
figure(1);
plot(targetRange, waterNR(:, 2), 'Color', ...
        [colorTable(sidx, :), 0.6], LineWidth=1); hold on;
plot(targetRange, iceNR(:, 2), 'Color', ...
        [colorTable(eidx, :), 0.6], LineWidth=1);
legend("水", "冰");
% set(gca, "YScale", "log");
xlabel("波长(um)");
ylabel("折射率实部");
title("冰水的折射率实部");
xlim([800, 1700]);
grid on;

%----------------------------冰水的折射率虚部-------------------------------
figure(2);
plot(targetRange, waterNI(:, 2), 'Color', ...
        [colorTable(sidx, :), 0.6], LineWidth=1); hold on;
plot(targetRange, iceNI(:, 2), 'Color', ...
        [colorTable(eidx, :), 0.6], LineWidth=1);
legend("水", "冰");
set(gca, "YScale", "log");
xlabel("波长(um)");
ylabel("折射率虚部");
title("冰水的折射率虚部");
xlim([800, 1700]);
grid on;

% %计算冰和水随波段变化的吸收距离和比值
% figure(3);
% Dwater = 1e-6 / (4 * pi) * waterNI(:, 1) ./ waterNI(:, 2);
% Dice = 1e-6 / (4 * pi) * iceNI(:, 1) ./ iceNI(:, 2);
% plot(targetRange, Dwater, 'Color', ...
%         [colorTable(sidx, :), 0.6], LineWidth=1); hold on;
% plot(targetRange, Dice, 'Color', ...
%         [colorTable(eidx, :), 0.6], LineWidth=1);
% legend("水", "冰");
% set(gca, "YScale", "log");
% xlabel("波长(um)");
% ylabel("吸收距离(m)");
% title("冰水的吸收距离");
% xlim([800, 1700]);
% grid on;

%首先需要确定波段
lambdas = 800: 100: 1700;
%直接获取对应的折射率实部
OT = OptTool();
lambdaNum = size(lambdas, 2);
iceNRList = zeros(1, lambdaNum);
waterNRList = zeros(1, lambdaNum);
airNRList = ones(1, lambdaNum);

for i = 1: lambdaNum
    iceNRList(1, i) = OT.findN(lambdas(1, i), OT.INR);
    waterNRList(1, i) = OT.findN(lambdas(1, i), OT.WNR);
end

%角度范围(rad)
inThetaO = 0:1:40;
inTheta = inThetaO * pi / 180;
thetaNum = size(inTheta, 2);
refI2W = zeros(lambdaNum, thetaNum);
refW2A = zeros(lambdaNum, thetaNum);

[colorTable, lambdaStr] = CG.generate(lambdas);

%-------------------不同波段从冰到水的反射比随入射角的变化--------------------
for j = 1: lambdaNum
    for i = 1: thetaNum
        refI2W(j, i) = OT.ref(inTheta(1, i), OT.snell(iceNRList(1, j), ...
            waterNRList(1, j), inTheta(1, i)));
    end
end

figure(3);
for i = 1: lambdaNum
    plot(inThetaO, refI2W(i, :), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=1); hold on;
end
legend(lambdaStr);
set(gca, "YScale", "log");
xlabel("入射角");
ylabel("反射比");
title("从冰进入水中");
grid on;

%-------------------不同波段从水到空气的反射比随入射角的变化-------------------
for j = 1: lambdaNum
    for i = 1: thetaNum
        refW2A(j, i) = OT.ref(inTheta(1, i), OT.snell(waterNRList(1, j), ...
            airNRList(1, j), inTheta(1, i)));
    end
end

figure(4);
for i = 1: lambdaNum
    plot(inThetaO, refW2A(i, :), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=1); hold on;
end
legend(lambdaStr);
set(gca, "YScale", "log");
xlabel("入射角");
ylabel("反射比");
title("从水进入空气中");
grid on;