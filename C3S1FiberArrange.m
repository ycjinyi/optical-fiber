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
H1 = (0: 0.1: 3) * 1e-3;
CG = ColorGenerator();
OC = OptCompute();
%设置介质层属性
OT = OptTool();
SL = SingleLayer(OT.INR, OT.INI, OT.ANR);
FG = FiberGenerator();
%发射光纤和接收光纤的横向光纤个数
fNum = 6;
%发射光纤和接收光纤的纵向光纤个数
verticalNum = 6;
%光纤组*1行x坐标 = 光纤数目,一发两收
group = [1, 2, 3, 4, ];
%排布模式的个数
pNum = size(group, 2);
Tpos = zeros(pNum, fNum);
Rpos = zeros(pNum, fNum);

%i遍历光纤组，代表多少个一组
for i = 1: pNum
    number = group(1, i);
    tidx = 1;
    ridx = 1;
    %直接遍历所有光纤
    for j = 0: 2 * fNum - 1
        %判断当前属于发射光纤还是接收光纤
        if mod(floor(j / number), 2) == 0 && tidx <= fNum
            Tpos(i, tidx) = j * coreDis;
            tidx = tidx + 1;
        else
            Rpos(i, ridx) = j * coreDis;
            ridx = ridx + 1;
        end
    end
end

%创建保存结果的矩阵
%行为不同的厚度点，列为不同光纤组下的结果
pic1Data = zeros(size(H1, 2), pNum);
tic;
for i = 1: pNum
    %首先根据每组光纤的横坐标和纵向光纤数目进行拓展
    TF = zeros(2, fNum * verticalNum);
    RF = zeros(2, fNum * verticalNum);
    for j = 1: fNum
        for k = 1: verticalNum
            idx = (j - 1) * verticalNum + k;
            TF(1, idx) = Tpos(i, j);
            TF(2, idx) = (k - 1) * coreDis;
            RF(1, idx) = Rpos(i, j);
            RF(2, idx) = (k - 1) * coreDis;
        end
    end
    %fluxs(波段,厚度,接收光纤)
    [fluxs, ic, nc, rc] = OC.idealCompute(SL, ...
                            FG.posConvert(TF, RF),...
                            lambdas, H1, U, S, R, dtheta, dphi);
    %划分为2路接收光纤
    idx1 = size(RF, 2);
    fluxs1 = sum(fluxs(1, :, 1: idx1), 3);
    pic1Data(:, i) = fluxs1;
end
toc;
save pic1Data pic1Data;
load pic1Data.mat;
showI = [1, 2, 3, 4, 5];
[colorTable, lambdaStr] = CG.generate(group(1, showI));
%作图展示
figure(1);
for i = 1: size(showI, 2)
    idx = showI(1, i);
    plot(H1 * 1e3, pic1Data(:, idx), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=1); hold on;
end
grid on;
legend(lambdaStr);
xlabel("冰厚度(mm)");
ylabel("光通量(lm)");
title("不同光纤排布下接收光纤的响应");