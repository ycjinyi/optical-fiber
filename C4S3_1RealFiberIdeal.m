clc;
clear;
close all;

%此脚本用于计算实际光纤排布下非耦合的响应, 计算0-1mm小范围

%天顶角和方位角的网格大小
dphi = pi / 1080; 
dtheta = pi / 1080;
%光源总的光通量为S(lm)
S = 1;
%接收光纤的接收半径为R(m)
R = 1e-6 * 187.562 / 2;
%最大出射平面孔径角(rad)
U = 12 * pi / 180;

%获取光纤排布参数
FG = FiberGenerator();
%posMatrix有3个维度[a, b, c]
%a代表对应的发射光纤编号,b代表对应的接收光纤编号
%c的取值是1和2,1代表距离之间的关系,2代表角度之间的关系,注意是弧度制
[SPM, R1PM, R2PM, R3PM] = FG.realFiberGen();
% posMatrix = [R1PM, R2PM, R3PM];
posMatrix = R1PM;

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
% lambdas = [890, 1400, 1650];
% lambdas = [800: 10: 980, 1250: 10: 1650];
lambdas = [890, 1350, 1450, 1550];
pNum = size(lambdas, 2);

%设置介质厚度数据(m)
H1 = (0.02: 0.02: 1) * 1e-3;
H2 = (0.02: 0.02: 1) * 1e-3;

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
% OC = OptCompute();
% 
% tic;
% %flux[a, b, c] 3维,a是波段,b是厚度点,c是接收光纤
% % [flux, ~] = OC.idealCompute(SL, posMatrix, ...
% %                lambdas, H1, U, S, R, dtheta, dphi);
% [flux, ~] = OC.idealCompute(BL, posMatrix, ...
%                lambdas, H, U, S, R, dtheta, dphi);
% toc;
% 
% 
% flux1 = sum(flux, 3);
% %将数据再转换为网格形式
% Z = zeros(pNum, size(X, 1), size(X, 2));
% for i = 1: pNum
%     for j = 1: H1Num
%         for k = 1: H2Num
%             Z(i, k, j) = flux1(i, (j - 1) * H2Num + k);
%         end
%     end
% end

% save 2024030601.mat R1PM BL lambdas H1 H2 H flux X Y Z;

load 2024030601.mat

% s = mesh(X, Y,  squeeze(Z(1, :, :))); 
% xlabel("冰厚");
% ylabel("水厚");
% s.FaceColor = 'flat';

z1 = squeeze(Z(1, :, :));
z2 = squeeze(Z(2, :, :));
z3 = squeeze(Z(3, :, :));
z4 = squeeze(Z(4, :, :));

level = 10;


%颜色表和标签
CG = ColorGenerator();
[colorTable, ~] = CG.generate(ones(1, 100));

X1 = X * 1e3;
Y1 = Y * 1e3;

figure;  
%对z1进行拟合
[fitresult, ~] = surfFit(X, Y, z1, 0.35);
[C,h] = contour(X1, Y1, 1e3 * fitresult(X, Y), level, 'LineWidth', 0.9, 'ShowText', 'on');
colormap(colorTable);
xlabel("冰厚(mm)");
ylabel("水厚(mm)");
title("890nm");
h.LevelList=round(h.LevelList,1);
clabel(C,h,'LabelSpacing',270);


figure;  
%对z2进行拟合
[fitresult, ~] = surfFit(X, Y, z2, 0.55);
[C,h] = contour(X1, Y1, 1e3 * fitresult(X, Y), level, 'LineWidth', 0.9, 'ShowText', 'on'); 
colormap(colorTable);
xlabel("冰厚(mm)");
ylabel("水厚(mm)");
title("1350nm");
h.LevelList=round(h.LevelList,1);
clabel(C,h,'LabelSpacing',270);
% colorbar;

figure;  
[C,h] = contour(X1, Y1, 1e3 * z3, level, 'LineWidth', 0.9, 'ShowText', 'on'); 
colormap(colorTable);
xlabel("冰厚(mm)");
ylabel("水厚(mm)");
title("1450nm");
h.LevelList=round(h.LevelList,1);
clabel(C,h,'LabelSpacing',270);
% colorbar;

figure;  
% 绘制等高线图
[C,h] = contour(X1, Y1, 1e3 * z4, level, 'LineWidth', 0.9, 'ShowText', 'on');
colormap(colorTable);
xlabel("冰厚(mm)");
ylabel("水厚(mm)");
title("1550nm");
h.LevelList=round(h.LevelList,1);
% clabel(C, h,'fontsize',9,'color','b','rotation',0)
clabel(C,h,'LabelSpacing',270);
% colorbar;
% grid on;


% set(gca, 'Box', 'off', ...                                        % 边框
%          'LineWidth', 1,...                                       % 线宽
%          'XGrid', 'off', 'YGrid', 'off', ...                      % 网格
%          'TickDir', 'out', 'TickLength', [.015 .015], ...         % 刻度
%          'XMinorTick', 'on', 'YMinorTick', 'on', ...              % 小刻度
%          'XColor', [.1 .1 .1],  'YColor', [.1 .1 .1])