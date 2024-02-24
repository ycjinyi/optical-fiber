clc;
clear;
close all;

%此脚本计算考虑光源耦合下的光纤探测模型

%计算参数
%dphi(rad), dtheta(rad)为网格大小
dphi = pi / 3000; 
dtheta = pi / 3000;
%光源视场角半角为V(rad), 光源总的光通量为S(lm)
V = 30 * pi / 180;
S = 1;
%光源和光纤端面的距离H(m)
H = 1e-3 * 2;
%接收光纤的接收半径为R(m)
R = 1e-6 * 90;
%波段数据(nm)
lambdas = 890;
%介质厚度数据(m)
hPoints = (0.01: 0.002: 4) * 1e-3;
%光纤纤芯的折射率
NF = 1.445;
%光纤子午面内的数值孔径
NA = 0.22;
%光源所在介质的折射率
NS = 1;

CG = ColorGenerator();
OC = OptCompute();
FG = FiberGenerator();

%设置介质层属性
OT = OptTool();
SL = SingleLayer(OT.INR, OT.INI, OT.ANR);

%-----------单光纤对,轴向距离、方位角不变,数值孔径、径向距离变化--------------
%生成光纤对的posMatrix
coreDis = 1e-6 * 600;
[s, r1] = FG.singleFiberGen(coreDis);
posMatrix = FG.posConvert(s, r1);
%不同的数值孔径
NAS = [0.2, 0.27];
%不同的径向距离
lightDisO = (0: 0.1: 1);
lightDis = 1e-3 * lightDisO;
disNum = size(lightDis, 2);
NANum = size(NAS, 2);
%数值孔径、表面介质厚度、光源距离光纤的径向距离
pic1Data = zeros(NANum, size(hPoints, 2), disNum);
tic;
for i = 1: NANum
    for j = 1: disNum
        [s1, r2] = FG.singleFiberGen(lightDis(1, j));
        [fluxs, ~] = OC.couplingCompute(SL, FG.posConvert(s1, r2),...
            posMatrix, lambdas, hPoints, H, V, S,...
                R, NAS(1, i), NF, NS, true, dtheta, dphi);
        pic1Data(i, :, j) = fluxs(1, :, 1);
    end
end
toc;
save pic1Data pic1Data;
load pic1Data.mat;
[colorTable, lambdaStr] = CG.generate(lightDisO);
%作图展示
figure(1);
for i = 1: disNum
    plot(hPoints * 1e3, pic1Data(1, :, i), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=1); hold on;
end
grid on;
legend(lambdaStr);
xlabel("冰厚度(mm)");
ylabel("光通量(lm)");
title("NA0.2 不同光源径向距离下的响应");
figure(2);
for i = 1: disNum
    plot(hPoints * 1e3, pic1Data(2, :, i), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=1); hold on;
end
grid on;
legend(lambdaStr);
xlabel("冰厚度(mm)");
ylabel("光通量(lm)");
title("NA0.27 不同光源径向距离下的响应");