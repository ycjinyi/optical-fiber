clc;
clear;

%获取光纤排布参数
FG = FiberGenerator();
[SPM, R1PM, R2PM, R3PM] = FG.realFiberGen();
% 不同的posMatrix之间通过第三个维度进行拼接
posMatrix = cat(3, R1PM, R2PM, R3PM);

%设置介质层属性
OT = OptTool();
% BL = BiLayer(OT.INR, OT.INI, OT.WNR, OT.WNI, OT.ANR);
% SL = SingleLayer(OT.INR, OT.INI, OT.ANR);
SL = SingleLayer(OT.WNR, OT.WNI, OT.ANR);

%设置波段数据
% lambdas = 890;
lambdas = [890, 1350, 1450, 1550];

%设置介质厚度数据
%接收光纤1的厚度点
H1 = [(0.05: 0.05: 1), (1.2: 0.1: 2.4), (2.6: 0.4: 8.4), (9: 0.5: 10)] * 1e-3;
%接收光纤2的厚度点
H2 = [(0.1: 0.3: 8.1), (8.5: 0.5: 10)] * 1e-3;
%接收光纤3的厚度点
H3 = (0.1: 0.5: 10.1) * 1e-3;

%传入不同模型的厚度数据和波段数据进行计算
OC = OptCompute();
%flux [a, b, c] 3维,a是波段,b是厚度点,c是接收光纤
[flux1, ~] = OC.compute(SL, SPM, R1PM, lambdas, H1, true);
[flux2, ~] = OC.compute(SL, SPM, R2PM, lambdas, H2, true);
[flux3, ~] = OC.compute(SL, SPM, R3PM, lambdas, H3, true);

%首先需要按照拼接时的关系对计算结果进行拆分,分别将对应接收光纤的光通量累计起来
flux1 = sum(flux1(:, :, :), 3);
flux2 = sum(flux2(:, :, :), 3);
flux3 = sum(flux3(:, :, :), 3);

save 2024031901.mat R1PM R2PM R3PM SPM SL lambdas H1 H2 H3 flux1 flux2 flux3;

%颜色表和标签
CG = ColorGenerator();
[~, lambdaStr] = CG.generate(lambdas);
[colorTable, ~] = CG.generate(ones(1, 38));
colorPen = 0.9;
coff = 9;

% load 2024030901.mat;
load 2024031901.mat;

figure(1);
for i = 1: size(lambdas, 2)
    plot(H1 * 1e3, flux1(i, :) + 1e-8, 'Color', ...
        [colorTable(i * coff, :), colorPen], LineWidth=1); hold on;
end
grid on;
legend(lambdaStr);
title("接收光纤1");
xlabel("水厚度(mm)");
ylabel("光通量(lm)");
xlim([0, 10]);
ylim([1e-8, 2* 1e-5]);
set(gca, "YScale", "log");

figure(2);
for i = 1: size(lambdas, 2)
    plot(H2 * 1e3, flux2(i, :)  + 1e-8, 'Color', ...
        [colorTable(i * coff, :), colorPen], LineWidth=1); hold on;
end
grid on;
legend(lambdaStr);
title("接收光纤2");
xlabel("水厚度(mm)");
ylabel("光通量(lm)");
xlim([0, 10]);
ylim([1e-8, 2* 1e-5]);
set(gca, "YScale", "log");

figure(3);
for i = 1: size(lambdas, 2)
    plot(H3 * 1e3, flux3(i, :)  + 1e-8, 'Color', ...
        [colorTable(i * coff, :), colorPen], LineWidth=1); hold on;
end
grid on;
legend(lambdaStr);
title("接收光纤3");
xlabel("水厚度(mm)");
ylabel("光通量(lm)");
xlim([0, 10]);
ylim([1e-8, 2* 1e-5]);
set(gca, "YScale", "log");