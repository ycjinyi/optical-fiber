clc;
clear;
close all;

%此脚本计算不同光纤排布模式下的响应特点

%dphi(rad), dtheta(rad)为网格大小
dphi = pi / 720; 
dtheta = pi / 720;
%最大出射平面孔径角为U(rad), 光源总的光通量为S(lm)
U = 20 * pi / 180;
S = 1;
%接收光纤的接收半径为R(m)
R = 1e-6 * 40;
%光纤组中的光纤间距(m)
coreDis = 1e-6 * 90;
%波段数据(nm)
lambdas = 890;
%介质厚度数据(m)
H1 = (0: 0.01: 3) * 1e-3;
CG = ColorGenerator();
OC = OptCompute();
%设置介质层属性
OT = OptTool();
SL = SingleLayer(OT.INR, OT.INI, OT.ANR);
FG = FiberGenerator();
%发射光纤和接收光纤的纵向光纤个数
verticalNum = 6;
%光纤混合时光纤的根数
group = 1: 6;
%混合的模式数目
patternNum = size(group, 2);
%光纤混合的轮数
roundNum = 2;
%创建保存结果的矩阵
%不同的厚度点，不同路的接收光纤、不同光纤排布
% pic1Data = zeros(size(H1, 2), 2, patternNum);
% tic;
% for i = 1: patternNum
%     %首先根据不同模式产生发射光纤和接收光纤的横坐标
%     number = group(1, i);
%     allNum = roundNum * number;
%     Tpos = zeros(1, allNum);
%     R1pos = zeros(1, allNum);
%     R2pos = zeros(1, allNum);
%     idx = 0;
%     tidx = 1;
%     ridx = 1;
%     %发射光纤和接收光纤1
%     for j = 1: 2 *roundNum
%         for k = 1: number
%             if mod(j, 2) == 1
%                 Tpos(1, tidx) = idx * coreDis;
%                 tidx = tidx + 1;
%             else
%                 R1pos(1, ridx) = idx * coreDis;
%                 ridx = ridx + 1;
%             end
%             idx = idx + 1;
%         end
%     end
%     %接收光纤2
%     for j = 1: allNum
%         R2pos(1, j) = idx * coreDis;
%         idx = idx + 1;
%     end
%     %根据每组光纤的横坐标和纵向光纤数目进行拓展
%     TF = zeros(2, allNum * verticalNum);
%     R1F = zeros(2, allNum * verticalNum);
%     R2F = zeros(2, allNum * verticalNum);
%     for j = 1: allNum
%         for k = 1: verticalNum
%             idx = (j - 1) * verticalNum + k;
%             TF(1, idx) = Tpos(1, j);
%             TF(2, idx) = (k - 1) * coreDis;
%             R1F(1, idx) = R1pos(1, j);
%             R1F(2, idx) = (k - 1) * coreDis;
%             R2F(1, idx) = R2pos(1, j);
%             R2F(2, idx) = (k - 1) * coreDis;
%         end
%     end
%     %fluxs(波段,厚度,接收光纤)
%     [fluxs, ic, nc, rc] = OC.idealCompute(SL, ...
%                             FG.posConvert(TF, [R1F, R2F]),...
%                             lambdas, H1, U, S, R, dtheta, dphi);
%     %接收光纤
%     idx = size(R1F, 2);
%     fluxs1 = sum(fluxs(1, :, 1: idx), 3);
%     fluxs2 = sum(fluxs(1, :, idx + 1: end), 3);
%     pic1Data(:, 1, i) = fluxs1;
%     pic1Data(:, 2, i) = fluxs2;
% end
% toc;
% save pic1Data pic1Data;
load pic1Data.mat;
showI = [1, 2, 3, 4, 5, 6];
[colorTable, lambdaStr] = CG.generate(group(1, showI));
%作图展示
figure(1);
for i = 1: size(showI, 2)
    idx = showI(1, i);
    plot(H1 * 1e3, pic1Data(:, 1, idx), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=1); hold on;
end
grid on;
legend(lambdaStr);
xlabel("冰厚度(mm)");
ylabel("光通量(lm)");
ylim([0, 5e-3]);
title("不同光纤排布下接收光纤1的响应");
figure(2);
for i = 1: size(showI, 2)
    idx = showI(1, i);
    plot(H1 * 1e3, pic1Data(:, 2, idx), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=1); hold on;
end
grid on;
legend(lambdaStr);
xlabel("冰厚度(mm)");
ylabel("光通量(lm)");
ylim([0, 1.4e-3]);
title("不同光纤排布下接收光纤2的响应");