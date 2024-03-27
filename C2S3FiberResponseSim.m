clc;
clear;
close all;

%此脚本只计算光纤探测模型, 不考虑光源耦合模型

computeSwitch = false; 

%计算参数
%dphi(rad), dtheta(rad)为网格大小
dphi = pi / 1800; 
dtheta = pi / 1800;
%最大出射平面孔径角为U(rad), 光源总的光通量为S(lm)
U = 20 * pi / 180;
S = 1;
%接收光纤的接收半径为R(m)
R = 1e-6 * 90;
%波段数据(nm)
lambdas = 890;
%介质厚度数据(m)
H1 = (0.01: 0.005: 3) * 1e-3;

CG = ColorGenerator();
OC = OptCompute();
%设置介质层属性
OT = OptTool();
% BL = BiLayer(OT.INR, OT.INI, OT.WNR, OT.WNI, OT.ANR);
SL = SingleLayer(OT.INR, OT.INI, OT.ANR);

FG = FiberGenerator();

% %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<单层介质,单对光纤>>>>>>>>>>>>>>>>>>>>>>>>>>>
% %---------------1、890nm/距离/出射角/介质不变, <半径>变化--------------
% ro = 40: 10: 100;
% r = ro * 1e-6;
% coreDis = 1e-6 * 600;
% [s, r1] = FG.singleFiberGen(coreDis);
% [posMatrix] = FG.posConvert(s, r1);
% pNum = size(ro, 2);
% %创建保存结果的矩阵
% if computeSwitch
%     %行为不同的厚度点，列为不同半径下的结果
%     pic1Data = zeros(size(H1, 2), pNum);
%     tic;
%     for i = 1: pNum
%         %fluxs(波段,厚度,接收光纤)
%         [fluxs, ic, nc, rc] = OC.idealCompute(SL, posMatrix, lambdas, H1,...
%                             U, S, r(1, i), dtheta, dphi);
%         pic1Data(:, i) = fluxs(1, :, 1);
%     end
%     toc;
%     save pic1Data pic1Data;
% end
% load pic1Data.mat;
% [colorTable, lambdaStr] = CG.generate(ro);
% %作图展示
% figure(1);
% for i = 1: pNum
%     plot(H1 * 1e3, pic1Data(:, i), 'Color', ...
%         [colorTable(i, :), 0.6], LineWidth=1); hold on;
% end
% grid on;
% legend(lambdaStr);
% xlabel("冰厚度(mm)");
% ylabel("光通量(lm)");
% title("不同光纤半径下的响应");

% % H1 = H1' * 1e3;


% 
% %---------------2、890nm/半径/出射角/介质不变, <距离>变化--------------
% coreDisO = 300: 100: 1000;
% pNum = size(coreDisO, 2);
% %创建保存结果的矩阵
% %行为不同的厚度点，列为不同距离下的结果
% if computeSwitch
%     pic2Data = zeros(size(H1, 2), pNum);
%     tic;
%     for i = 1: pNum
%         [s, r1] = FG.singleFiberGen(coreDisO(1, i) * 1e-6);
%         %fluxs(波段,厚度,接收光纤)
%         [fluxs, ic, nc, rc] = OC.idealCompute(SL, FG.posConvert(s, r1), lambdas, H1,...
%                             U, S, R, dtheta, dphi);
%         pic2Data(:, i) = fluxs(1, :, 1);
%     end
%     toc;
%     save pic2Data pic2Data;
% end
% load pic2Data.mat;
% [colorTable, lambdaStr] = CG.generate(coreDisO);
% %作图展示
% figure(2);
% for i = 1: pNum
%     plot(H1 * 1e3, pic2Data(:, i), 'Color', ...
%         [colorTable(i, :), 0.6], LineWidth=1); hold on;
% end
% grid on;
% legend(lambdaStr);
% xlabel("冰厚度(mm)");
% ylabel("光通量(lm)");
% title("不同光纤距离下的响应");

% % H1 = H1' * 1e3;


% 
% %---------------3、890nm/距离/半径/介质不变, <出射角>变化--------------
% uo = 12: 3: 42;
% u = uo * pi / 180;
% coreDis = 1e-6 * 600;
% [s, r1] = FG.singleFiberGen(coreDis);
% [posMatrix] = FG.posConvert(s, r1);
% pNum = size(uo, 2);
% %创建保存结果的矩阵
% if computeSwitch
%     %行为不同的厚度点，列为不同出射角下的结果
%     pic3Data = zeros(size(H1, 2), pNum);
%     tic;
%     for i = 1: pNum
%         %fluxs(波段,厚度,接收光纤)
%         [fluxs, ic, nc, rc] = OC.idealCompute(SL, posMatrix, lambdas, H1,...
%                             u(1, i), S, R, dtheta, dphi);
%         pic3Data(:, i) = fluxs(1, :, 1);
%     end
%     toc;
%     save pic3Data pic3Data;
% end
% load pic3Data.mat;
% [colorTable, lambdaStr] = CG.generate(uo);
% %作图展示
% figure(3);
% for i = 1: pNum
%     plot(H1 * 1e3, pic3Data(:, i), 'Color', ...
%         [colorTable(i, :), 0.6], LineWidth=1); hold on;
% end
% grid on;
% legend(lambdaStr);
% xlabel("冰厚度(mm)");
% ylabel("光通量(lm)");
% title("不同最大出射角下的响应");
% 
% H1 = H1' * 1e3;

% 
% %---------------4、距离/半径/出射角/介质不变, <波段>变化---------------
% lambdao = 1465:10:1565;
% coreDis = 1e-6 * 600;
% [s, r1] = FG.singleFiberGen(coreDis);
% [posMatrix] = FG.posConvert(s, r1);
% pNum = size(lambdao, 2);
% %创建保存结果的矩阵
% if computeSwitch
%     %行为不同的厚度点，列为不同波段下的结果
%     pic4Data = zeros(size(H1, 2), pNum);
%     tic;
%     for i = 1: pNum
%         %fluxs(波段,厚度,接收光纤)
%         [fluxs, ic, nc, rc] = OC.idealCompute(SL, posMatrix, lambdao(1, i), H1,...
%                             U, S, R, dtheta, dphi);
%         pic4Data(:, i) = fluxs(1, :, 1);
%     end
%     toc;
%     save pic4Data pic4Data;
% end
% load pic4Data.mat;
% [colorTable, lambdaStr] = CG.generate(lambdao);
% %作图展示
% figure(4);
% for i = 1: pNum
%     plot(H1 * 1e3, pic4Data(:, i), 'Color', ...
%         [colorTable(i, :), 0.6], LineWidth=1); hold on;
% end
% grid on;
% legend(lambdaStr);
% xlabel("冰厚度(mm)");
% ylabel("光通量(lm)");
% % set(gca, "YScale", "log");
% title("不同波段下的响应");


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<5、单层介质,多对光纤>>>>>>>>>>>>>>>>>>>>>>>>>>
% %发射光纤和接收光纤的个数
% fNum = 8;
% %光纤组*2行坐标*光纤数目
% group = 1:1:fNum;
% Tpos = zeros(fNum, 2, fNum);
% Rpos = zeros(fNum, 2, fNum);
% %光纤组中的光纤间距(m)
% coreDis = 1e-6 * 200;
% %i遍历光纤组，代表多少个一组
% for i = 1: fNum
%     tidx = 1;
%     ridx = 1;
%     %直接遍历所有光纤
%     for j = 0: 2 * fNum - 1
%         %判断当前属于发射光纤还是接收光纤
%         if mod(floor(j / i), 2) == 0 && tidx <= fNum
%             Tpos(i, 1, tidx) = j * coreDis;
%             tidx = tidx + 1;
%         else
%             Rpos(i, 1, ridx) = j * coreDis;
%             ridx = ridx + 1;
%         end
%     end
% end
% %创建保存结果的矩阵
% if computeSwitch
%     %行为不同的厚度点，列为不同光纤组下的结果
%     pic5Data = zeros(size(H1, 2), fNum);
%     tic;
%     for i = 1: fNum
%         %fluxs(波段,厚度,接收光纤)
%         [fluxs, ic, nc, rc] = OC.idealCompute(SL, ...
%                                 FG.posConvert(squeeze(Tpos(i,:,:)), squeeze(Rpos(i,:,:))),...
%                                 lambdas, H1, U, S, R, dtheta, dphi);
%         pic5Data(:, i) = sum(fluxs, 3);
%     end
%     toc;
%     save pic5Data pic5Data;
% end
% load pic5Data.mat;
% [colorTable, lambdaStr] = CG.generate(group);
% %作图展示
% figure(5);
% for i = 1: fNum
%     plot(H1 * 1e3, pic5Data(:, i), 'Color', ...
%         [colorTable(i, :), 0.6], LineWidth=1); hold on;
% end
% grid on;
% legend(lambdaStr);
% xlabel("冰厚度(mm)");
% ylabel("光通量(lm)");
% title("多光纤对不同光纤组下的响应");

% %<<<<<<<<<<<<<<<<<<<<<<<<<<<<<双层介质,单对光纤>>>>>>>>>>>>>>>>>>>>>>>>>>>
% %---------------6、波段不变,<冰水比例><总厚度>变化---------------
% lambda = [890, 1500];
% coreDis = 1e-6 * 600;
% [s, r1] = FG.singleFiberGen(coreDis);
% [posMatrix] = FG.posConvert(s, r1);
% H = (0.01: 0.01: 3) * 1e-3;
% %厚度中水的占比
% ratio = 0:0.1:1;
% pNum = size(ratio, 2);
% lNum = size(lambda, 2);
% %创建保存结果的矩阵
% if computeSwitch
%     %波段、厚度点、不同水占比
%     pic6Data = zeros(lNum, size(H, 2), pNum);
%     %考虑冰和水混合情况
%     BL1 = BiLayer(OT.INR, OT.INI, OT.WNR, OT.WNI, OT.ANR);
%     % SL1 = SingleLayer(OT.INR, OT.INI, OT.WNR);
%     tic;
%     for i = 1: pNum
%         %水厚度
%         Hw = H * ratio(1, i);
%         %冰厚度
%         Hi = H * (1 - ratio(1, i));
%         %fluxs(波段,厚度,接收光纤)
%         [fluxsB, ~] = OC.idealCompute(BL1, posMatrix, lambda, [Hi; Hw],...
%                             U, S, R, dtheta, dphi);
%         % [fluxsS, ~] = OC.idealCompute(SL1, posMatrix, lambda, Hi,...
%         %                     U, S, R, dtheta, dphi);
%         % pic6Data(:, :, i) = fluxsB(:, :, 1) + fluxsS(:, :, 1);
%         pic6Data(:, :, i) = fluxsB(:, :, 1);
%     end
%     toc;
%     save pic6Data pic6Data;
% end
% load pic6Data.mat;
% [colorTable, lambdaStr] = CG.generate(ratio);
% %作图展示
% figure(6);
% for i = 1: pNum
%     plot(H * 1e3, pic6Data(1, :, i), 'Color', ...
%         [colorTable(i, :), 0.6], LineWidth=1); hold on;
% end
% grid on;
% legend(lambdaStr);
% xlabel("总厚度(mm)");
% ylabel("光通量(lm)");
% % set(gca, "YScale", "log");
% title("890nm,冰水混合,不同水占比下的响应");
% 
% %作图展示
% figure(7);
% for i = 1: pNum
%     plot(H * 1e3, pic6Data(2, :, i), 'Color', ...
%         [colorTable(i, :), 0.6], LineWidth=1); hold on;
% end
% grid on;
% legend(lambdaStr);
% xlabel("总厚度(mm)");
% ylabel("光通量(lm)");
% % set(gca, "YScale", "log");
% title("1500nm,冰水混合,不同水占比下的响应");

% % H = H' * 1e3;
% % ratio = round(ratio * 10) / 10;

% %---------------7、总厚度不变,<冰水比例><波段>变化---------------
lambdao = 1405:10:1565;
H = [0.8, 1.6] * 1e-3;
coreDis = 1e-6 * 600;
[s, r1] = FG.singleFiberGen(coreDis);
[posMatrix] = FG.posConvert(s, r1);
%厚度中水的占比
ratio = 0:0.05:1;
lNum = size(lambdao, 2);
hNum = size(H, 2);
rNum = size(ratio, 2);
%创建保存结果的矩阵
if computeSwitch
    %总厚度、波段、不同水占比
    pic7Data = zeros(hNum, lNum, rNum);
    %考虑冰和水混合情况
    BL1 = BiLayer(OT.INR, OT.INI, OT.WNR, OT.WNI, OT.ANR);
    SL1 = SingleLayer(OT.INR, OT.INI, OT.WNR);
    tic;
    for i = 1: hNum
        %水厚度
        Hw = H(1, i) * ratio;
        %冰厚度
        Hi = H(1, i) * (1 - ratio);
        %fluxs(波段,厚度,接收光纤)
        [fluxsB, ~] = OC.idealCompute(BL1, posMatrix, lambdao, [Hi; Hw],...
                            U, S, R, dtheta, dphi);
        % [fluxsS, ~] = OC.idealCompute(SL1, posMatrix, lambdao, Hi,...
        %                     U, S, R, dtheta, dphi);
        % pic7Data(i, :, :) = fluxsB(:, :, 1) + fluxsS(:, :, 1);
        pic7Data(i, :, :) = fluxsB(:, :, 1);
    end
    toc;
    save pic7Data pic7Data;
end
load pic7Data.mat;
[colorTable, lambdaStr] = CG.generate(lambdao);
%作图展示
figure(8);
for i = 1: lNum
    plot(ratio * 100, squeeze(pic7Data(1, i, :)), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=1); hold on;
end
grid on;
legend(lambdaStr);
xlabel("水占比(%)");
ylabel("光通量(lm)");
% set(gca, "YScale", "log");
title("总厚度0.8mm,冰水混合,不同波段下的响应");

%作图展示
figure(9);
for i = 1: lNum
    plot(ratio * 100, squeeze(pic7Data(2, i, :)), 'Color', ...
        [colorTable(i, :), 0.6], LineWidth=1); hold on;
end
grid on;
legend(lambdaStr);
xlabel("水占比(%)");
ylabel("光通量(lm)");
% set(gca, "YScale", "log");
title("总厚度1.6mm,冰水混合,不同波段下的响应");

ratio = round(ratio' * 100);
