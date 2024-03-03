clc;
clear;
close all;

%此脚本用于生成冰、水的折射率实部、虚部，并进行温度补偿

%找到目标文件路径
nowpath = pwd;
waterPath = strcat(nowpath, "\光学常数水\298.15k\");
icePath = strcat(nowpath, "\光学常数冰\266k\");
%读取对应的数据
waterNR = importdata(waterPath + "实部.txt");
waterNI = importdata(waterPath + "虚部.txt");
iceNR = importdata(icePath + "实部.txt");
iceNI = importdata(icePath + "虚部.txt");
%按照目标波段对数据进行插值统一,目标的波段范围,单位um
targetRange = 0.2: 0.01: 100;
targetRange = targetRange';
waterNR = linearInterpolation(waterNR, targetRange);
waterNI = linearInterpolation(waterNI, targetRange);
iceNR = linearInterpolation(iceNR, targetRange);
iceNI = linearInterpolation(iceNI, targetRange);

%颜色表和标签
CG = ColorGenerator();
%光学计算函数
OT = OptTool();

%<<<<<<<<<<<<<<<<<<<1、266k的冰和298.15k的水的复折射率<<<<<<<<<<<<<<<<
[colorTable, ~] = CG.generate(zeros(1, 17));
sidx = 7;
eidx = 16;
width = 1.2;
cof = 1;
%折射率实部
figure(1);
plot(targetRange, waterNR(:, 2), 'Color', ...
        [colorTable(sidx, :), cof], LineWidth=width); hold on;
plot(targetRange, iceNR(:, 2), 'Color', ...
        [colorTable(eidx, :), cof], LineWidth=width, LineStyle='-');
legend("water,298.15k", "ice,266k");
set(gca, "XScale", "log");
% set(gca, "YScale", "log");
xlabel("波长(um)");
ylabel("折射率实部");
xlim([0.2, 100]);
title("冰水的折射率实部");
grid on;

%折射率虚部
figure(2);
plot(targetRange, waterNI(:, 2), 'Color', ...
        [colorTable(sidx, :), cof], LineWidth=width); hold on;
plot(targetRange, iceNI(:, 2), 'Color', ...
        [colorTable(eidx, :), cof], LineWidth=width, LineStyle='-');
legend("water,298.15k", "ice,266k");
set(gca, "XScale", "log");
set(gca, "YScale", "log");
xlabel("波长(um)");
ylabel("折射率虚部");
xlim([0.2, 100]);
title("冰水的折射率虚部");
grid on;

%<<<<<<<<<<<<<<<<<<<2、266k的冰和298.15k的水的吸收系数和穿透深度<<<<<<<<<<<<<<
%将折射率虚部转换为吸收系数
waterAbs = 4 * pi * waterNI(:, 2) *  1e6 ./ waterNI(:, 1);
iceAbs = 4 * pi * iceNI(:, 2) *  1e6 ./ iceNI(:, 1);
waterPen = ones(size(waterAbs, 1), size(waterAbs, 2)) ./ waterAbs;
icePen = ones(size(iceAbs, 1), size(iceAbs, 2)) ./ iceAbs;
%吸收系数
figure(3);
plot(targetRange, waterAbs, 'Color', ...
        [colorTable(sidx, :), cof], LineWidth=width); hold on;
plot(targetRange, iceAbs, 'Color', ...
        [colorTable(eidx, :), cof], LineWidth=width, LineStyle='-');
legend("water,298.15k", "ice,266k");
set(gca, "XScale", "log");
set(gca, "YScale", "log");
xlabel("波长(um)");
ylabel("吸收系数");
xlim([0.2, 100]);
title("冰水的吸收系数");
grid on;

%穿透深度
figure(4);
plot(targetRange, waterPen, 'Color', ...
        [colorTable(sidx, :), cof], LineWidth=width); hold on;
plot(targetRange, icePen, 'Color', ...
        [colorTable(eidx, :), cof], LineWidth=width, LineStyle='-');
legend("water,298.15k", "ice,266k");
set(gca, "XScale", "log");
set(gca, "YScale", "log");
xlabel("波长(um)");
ylabel("穿透深度");
xlim([0.2, 100]);
title("冰水的穿透深度");
grid on;

%<<<<<<<<<<<<<<<<<<<3、266k的冰的温度依赖特性<<<<<<<<<<<<<<<<
width = 1;
%读取冰的温度依赖特性文件
iceTemp = importdata(strcat(nowpath, "\光学常数冰\冰的温度依赖特性.txt"));
tempRange = 160: 10: 270;
targetRange = (0.8: 0.005: 1.8)';
iceNITemp = [iceTemp(:, 1), iceTemp(:, 2 + size(tempRange, 2): end)];
[colorTable, legendStr] = CG.generate(tempRange);
%首先转换为目标范围
iceNITemp = linearInterpolation(iceNITemp, targetRange);
%将虚部转换为穿透深度
for i = 2: size(iceNITemp, 2)
    iceNITemp(:, i) = 4 * pi * iceNITemp(:, i) *  1e6 ./ iceNITemp(:, 1);
end
iceNITemp(:, 2: end) = ones(size(iceNITemp, 1), ...
    size(iceNITemp, 2) - 1) ./ iceNITemp(:, 2: end);
%作图
figure(5);
for i = 1: size(tempRange, 2)
    plot(targetRange, iceNITemp(:, i + 1), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=width); hold on;
end
legend(legendStr);
% set(gca, "XScale", "log");
set(gca, "YScale", "log");
xlabel("波长(um)");
ylabel("穿透深度");
% xlim([0.2, 100]);
title("冰穿透深度的温度依赖特性");
grid on;

%计算一下温度依赖特性
iceTempDep = zeros(size(iceNITemp, 1), size(iceNITemp, 2) - 2);
for i = 2: size(iceNITemp, 2)  - 1
    nowNI = iceNITemp(:, i);
    baseNI = iceNITemp(:, end);
    iceTempDep(:, i - 1) = (baseNI - nowNI) ./ baseNI /...
    (tempRange(1, end) - tempRange(1,i - 1));
end
%作图
figure(6);
for i = 1: size(tempRange, 2) - 1
    plot(targetRange, iceTempDep(:, i), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=width); hold on;
end
legend(legendStr(1, 1: end - 1));
% set(gca, "XScale", "log");
% set(gca, "YScale", "log");
xlabel("波长(um)");
ylabel("穿透深度");
% xlim([0.2, 100]);
title("冰穿透深度的温度依赖特性");
grid on;



%<<<<<<<<<<<<<<<<<<<4、水的温度补偿系数和由此计算的温度漂移特性<<<<<<<<<<<<<<<<
[colorTable, ~] = CG.generate(zeros(1, 17));
idx = 16;
%首先读取水的温度补偿系数,第一列为波长,第二列为波数,第三列为温度补偿系数,第四列为标准差
waterTempCoff = importdata(strcat(nowpath, "\光学常数水\水的温度补偿系数.txt"));
waterTempCoff = [waterTempCoff(:, 1) / 1e3, waterTempCoff(:, 3)];
targetRange = (0.8: 0.001: 1.8)';
%首先转换为目标范围
waterTempCoff = linearInterpolation(waterTempCoff, targetRange);
figure(7);
plot(waterTempCoff(:, 1), waterTempCoff(:, 2), 'Color', ...
        [colorTable(idx, :), 1], LineWidth=width);
legend("水");
xlabel("波长(um)");
ylabel("补偿系数(每米摄氏度)");
title("水的温度补偿系数");
grid on;
%读取水在293.15k的虚部数据
waterPath = strcat(nowpath, "\光学常数水\293.15k\");
waterNI = importdata(waterPath + "虚部.txt");
waterNI = linearInterpolation(waterNI, targetRange);
%以该温度为基准, 获得从274k - 294k温度变化下的数据
%首先转换为吸收系数
waterAbs1 = 4 * pi * waterNI(:, 2) *  1e6 ./ waterNI(:, 1);
tempRange1 = 274: 2: 294;
waterTemp = zeros(size(targetRange, 1), size(tempRange1, 2) + 1);
waterTemp(:, 1) = targetRange;
%然后计算不同温度下的数据
for i = 2: size(waterTemp, 2)
    waterTemp(:, i) = waterAbs1 + (tempRange1(1, i - 1) - 293.15) * waterTempCoff(:, 2);
end
waterTemp(:, 2: end) = ones(size(waterTemp, 1), ...
    size(waterTemp, 2) - 1) ./ waterTemp(:, 2: end);
[colorTable, legendStr] = CG.generate(tempRange1);
%作图展示
figure(8);
for i = 1: size(tempRange1, 2)
    plot(targetRange, waterTemp(:, i + 1), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=width); hold on;
end
legend(legendStr);
% set(gca, "XScale", "log");
set(gca, "YScale", "log");
xlabel("波长(um)");
ylabel("穿透深度");
% xlim([0.2, 100]);
title("水穿透深度的温度依赖特性");
grid on;

%<<<<<<<<<<<<<<<<<<<5、后续计算基于的水和冰的温度对应的光学常数<<<<<<<<<<<<<<<<


