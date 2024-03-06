clc;
clear;
close all;

%此脚本用于计算实际光纤排布下非耦合的响应

%天顶角和方位角的网格大小
dphi = pi / 900; 
dtheta = pi / 900;
%光源总的光通量为S(lm)
S = 1;
%接收光纤的接收半径为R(m)
R = 1e-6 * 187.562 / 2;
%最大出射平面孔径角(rad)
U = 10 * pi / 180;

%获取光纤排布参数
FG = FiberGenerator();
%posMatrix有3个维度[a, b, c]
%a代表对应的发射光纤编号,b代表对应的接收光纤编号
%c的取值是1和2,1代表距离之间的关系,2代表角度之间的关系,注意是弧度制
[SPM, R1PM, R2PM, R3PM] = FG.realFiberGen();
% posMatrix = [R1PM, R2PM, R3PM];
posMatrix = [R1PM, R2PM];

% check = R2PM;
% %检查矩阵是否计算正确
% figure;
% plot(sort(reshape(check(4, :, 2) * 180 / pi, [1, size(check, 2)])));
% figure;
% plot(sort(reshape(check(4, :, 1), [1, size(check, 2)])));

%设置介质层属性
OT = OptTool();
BL = BiLayer(OT.INR, OT.INI, OT.WNR, OT.WNI, OT.ANR);
SL = SingleLayer(OT.WNR, OT.WNI, OT.ANR);

%设置波段数据(nm)
lambdas = 890;
% lambdas = [800: 10: 980, 1250: 10: 1650];
% lambdas = [890, 1350, 1450, 1550];

%设置介质厚度数据(m)
H1 = (0.1: 0.1: 8) * 1e-3;
H2 = (0.02: 0.02: 8) * 1e-3;

%传入不同模型的厚度数据和波段数据进行计算
OC = OptCompute();

tic;
%flux[a, b, c] 3维,a是波段,b是厚度点,c是接收光纤
[flux, ~] = OC.idealCompute(SL, posMatrix, ...
               lambdas, H1, U, S, R, dtheta, dphi);
toc;

%按照拼接时的关系对计算结果进行拆分,将对应接收光纤的光通量累计起来
r1Num = size(R1PM, 2);
r2Num = size(R2PM, 2);
r3Num = size(R3PM, 2);
flux1 = sum(flux(:, :, 1: r1Num), 3);
flux2 = sum(flux(:, :, r1Num + 1: r1Num + r2Num), 3);
% flux3 = sum(flux(:, :, r1Num + r2Num + 1: r1Num + r2Num + r3Num), 3);
x = H1;

%颜色表和标签
CG = ColorGenerator();
[colorTable, lambdaStr] = CG.generate(lambdas);


% %归一化展示
% maxR = max(max(flux, [], 2));
% coff = maxR ./ max(flux, [], 2);
% flux = flux .* coff;

figure;
for i = 1: size(lambdas, 2)
    plot(x * 1e3, flux1(i, :), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=1); hold on;
end
grid on;
legend(lambdaStr);
xlabel("厚度");
ylabel("光通量");
title("接收光纤1");
%set(gca, "YScale", "log");
%xlim([0, 1e-3 * 1.2]);
%ylim([0, 3 * 1e-7]);

figure;
for i = 1: size(lambdas, 2)
    plot(x * 1e3, flux2(i, :), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=1); hold on;
end
grid on;
legend(lambdaStr);
xlabel("厚度");
ylabel("光通量");
title("接收光纤2");
%set(gca, "YScale", "log");
%xlim([0, 1e-3 * 1.2]);
%ylim([0, 3 * 1e-7]);

%save 2024010801.mat R1PM R2PM R3PM posMatrix SL lambdas H1 Sflux ic nc rc;