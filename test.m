clc;
clear;

%创建颜色
%颜色表和标签
CG = ColorGenerator();
[colorTable, lambdaStr] = CG.generate(zeros(1, 6));


%获取光纤排布参数
FG = FiberGenerator();
% [s, r1, r2] = FG.regFiberGen(1e-3 * 0.227, 6);
% posMatrix1 = FG.posConvert(s, r1);
% posMatrix2 = FG.posConvert(s, r2);
% 不同的posMatrix之间通过第三个维度进行拼接
% posMatrix = cat(3, posMatrix1, posMatrix2);
[SPM, R1PM, R2PM, R3PM] = FG.realFiberGen();
posMatrix = cat(3, R1PM, R2PM, R3PM);
%posMatrix = cat(3, R1PM, R2PM, R3PM);

% figure;
% plot(sort(reshape(SPM(2, 1, :) * 180 / pi, [1, size(SPM, 3)])));
% figure;
% plot(sort(reshape(SPM(1, 1, :), [1, size(SPM, 3)])));

%设置介质层属性
OT = OptTool();
% BL = BiLayer(OT.INR, OT.INI, OT.WNR, OT.WNI, OT.ANR);
SL = SingleLayer(OT.WNR, OT.WNI, OT.ANR);

%设置波段数据
%lambdas = 1200: 20: 1600;
lambdas = [890];
% a = 800: 10: 960;
% b = 1200: 10: 1700;
% lambdas = [a, b];

%设置介质厚度数据
H1 = (0.1: 0.1: 15) * 1e-3;
% H1 = (1: 1: 8) * 1e-3;
% H1 = ones(1, 20) * 1e-3 * 0.5;
% H2 = (0.1: 0.1: 2) * 1e-3 * 0.4;

%传入不同模型的厚度数据和波段数据进行计算
OC = OptCompute();
%[Bflux, ic, nc, rc] = OC.compute(BL, SPM, posMatrix, lambdas, [H1; H2]);
[Sflux, ic, nc, rc] = OC.compute(SL, SPM, posMatrix, lambdas, H1);
%flux [a, b, c] 3维,a是波段,b是厚度点,c是接收光纤

%首先需要按照拼接时的关系对计算结果进行拆分,分别将对应接收光纤的光通量累计起来
r1Num = size(R1PM, 3);
r2Num = size(R2PM, 3);
r3Num = size(R3PM, 3);
Sflux1 = sum(Sflux(:, :, 1: r1Num), 3);
Sflux2 = sum(Sflux(:, :, r1Num + 1: r1Num + r2Num), 3);
% Sflux3 = sum(Sflux(:, :, r1Num + r2Num + 1: r1Num + r2Num + r3Num), 3);
flux = Sflux2;
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
    plot(x, flux(i, :), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=1); hold on;
end
grid on;
legend(lambdaStr);
%legend(lambdaStr(tarIdx));
xlabel("水厚度");
ylabel("光通量");
%set(gca, "YScale", "log");
%xlim([0, 1e-3 * 1.2]);
%ylim([0, 3 * 1e-7]);

%save 2024010801.mat R1PM R2PM R3PM posMatrix SL lambdas H1 Sflux ic nc rc;