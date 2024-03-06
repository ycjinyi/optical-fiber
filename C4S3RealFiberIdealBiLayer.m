clc;
clear;

%获取光纤排布参数
FG = FiberGenerator();
[SPM, R1PM, R2PM, R3PM] = FG.realFiberGen();

%设置介质层属性
OT = OptTool();
BL = BiLayer(OT.INR, OT.INI, OT.WNR, OT.WNI, OT.ANR);

%设置波段数据
% lambdas = 890;
lambdas = [890, 1350, 1450, 1550];
pNum = size(lambdas, 2);

%设置介质厚度数据
%接收光纤1的厚度点
H1 = (0.05: 0.05: 0.95) * 1e-3;
H2 = 1e-3 - H1;
% H2 = (0.5: 0.5: 6) * 1e-3;

%转换为网格数据
%X的行和H1相同
%Y的列和H2相同
[X, Y] = meshgrid(H1, H2);
%H厚度点
H1Num = size(H1, 2);
H2Num = size(H2, 2);
H = zeros(2, H1Num * H2Num);
for i = 1: H1Num
    for j = 1: H2Num
        p = (i - 1) * H2Num + j;
        H(1, p) = H1(1, i);
        H(2, p) = H2(1, j);
    end
end

%传入不同模型的厚度数据和波段数据进行计算
OC = OptCompute();
%flux [a, b, c] 3维,a是波段,b是厚度点,c是接收光纤
[flux1, ~] = OC.compute(BL, SPM, R1PM, lambdas, [H1; H2], true);
flux1 = sum(flux1(:, :, :), 3);

CG = ColorGenerator();
[colorTable, lambdaStr] = CG.generate(lambdas);
colorPen = 0.85;
figure;
for i = 1: pNum
    plot(H1 * 1e3, flux1(i, :), 'Color', ...
        [colorTable(i, :), colorPen], LineWidth=1); hold on;
end
grid on;
legend(lambdaStr);


% save 2024030602.mat R1PM SPM BL lambdas H1 H2 H flux1;

%将数据再转换为网格形式
Z = zeros(pNum, size(X, 2), size(X, 1));
for i = 1: pNum
    for j = 1: H1Num
        for k = 1: H2Num
            Z(i, j, k) = flux1(i, (j - 1) * H2Num + k);
        end
    end
end

% load 2024030602.mat;

mesh(X, Y, squeeze(Z(1, :, :))); 
xlabel("冰厚");
ylabel("水厚");