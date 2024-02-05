clc;
clear;
close all;

%测试并生成论文中光纤出射光通量的分布模型数据

%--------------光纤半径不变, 最大出射角不变, 变化轴向和径向距离--------------
%参数设置
%光纤半径(m)
R = 1e-6 * 187.562 / 2;
%光通量大小(lm)
S = 1;
%最大出射角(rad)
U = 30 * pi / 180;
%径向距离分辨率(m)
dx = 1e-6 * 1;
%径向距离点(m)
xPoints = 0: dx: 1e-6*200;
%轴向距离分辨率(m)
dh = 1e-6* 10;
%轴向距离点(m)
hPoints = 0: dh: 15*dh;

%计算求解
tic;
FD = FiberLmDistribution();
lmMatrix = FD.illuminanceCompute(R, U, S, xPoints, hPoints);
toc;

%颜色表和标签
CG = ColorGenerator();
[colorTable, lambdaStr] = CG.generate(hPoints * 1e6);

hNumber = size(hPoints, 2);
figure;
for i = 1: hNumber
    plot(xPoints * 1e6, lmMatrix(i, :), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=1); hold on;
end
legend(lambdaStr);
grid on;
% ylim([0, 7e7]);
% set(gca, "YScale", "log");
xlabel('径向距离(um)');
ylabel('光照度(lm/m^2)');

%--------------光纤半径不变, 轴向距离不变, 变化最大出射角--------------
%径向距离分辨率(m)
dx = 1e-6 * 0.1;
%径向距离点(m)
xPoints = 0: dx: 1e-6*160;
xNumber = size(xPoints, 2);
%轴向距离(m)
hPoints = 30 * 1e-6;
%最大出射角(rad)
du = 2;
u = 18: du: 50;
U = u * pi / 180;
uNumber = size(U, 2);
lmMatrix = zeros(uNumber, xNumber);
for i = 1: uNumber
    lmMatrix(i, :) = FD.illuminanceCompute(R, U(1, i), S, xPoints, hPoints);    
end

[colorTable, lambdaStr] = CG.generate(u);
figure;
for i = 1: uNumber
    plot(xPoints * 1e6, lmMatrix(i, :), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=1); hold on;
end
legend(lambdaStr);
grid on;
% set(gca, "YScale", "log");
% ylim([0, 7e7]);
xlabel('径向距离(um)');
ylabel('光照度(lm/m^2)');


%--------------最大出射角度不变, 轴向距离不变, 变化光纤半径--------------
U = 30 * pi / 180;
%径向距离分辨率(m)
dx = 1e-6 * 0.1;
%径向距离点(m)
xPoints = 0: dx: 1e-6*140;
xNumber = size(xPoints, 2);
%光纤半径(m)
dr = 5;
r = 40: dr: 100;
R = 1e-6 * r;
rNumber = size(R, 2);
lmMatrix = zeros(rNumber, xNumber);
for i = 1: rNumber
    lmMatrix(i, :) = FD.illuminanceCompute(R(1, i), U, S, xPoints, hPoints);    
end

[colorTable, lambdaStr] = CG.generate(r);
figure;
for i = 1: rNumber
    plot(xPoints * 1e6, lmMatrix(i, :), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=1); hold on;
end
legend(lambdaStr);
grid on;
% set(gca, "YScale", "log");
% ylim([0, 7e7]);
xlabel('径向距离(um)');
ylabel('光照度(lm/m^2)');









