clc;
clear;
close all;

%此脚本用于计算实际光纤排布下的响应

% % angle(1, :) = angle(1, :) * 0;
% % angle(2, :) = ones(1, size(angle, 2)) * 0.4887;
% [~, idx] = sort(SPM(1, :, 1));
% figure;
% % plot(flux(idx));
% plot(flux, '*');
% grid on;
% % xlabel("按距离排序后的光纤编号");
% xlabel("光纤编号");
% ylabel("光通量");
% figure;
% plot(angle(1, idx) * 180 / pi); hold on;
% plot(angle(2, idx) * 180 / pi);
% % plot(angle(1, :) * 180 / pi); hold on;
% % plot(angle(2, :) * 180 / pi);
% legend("最小入射角", "最大入射角");
% xlabel("按距离排序后的光纤编号");
% % xlabel("光纤编号");
% ylabel("入射角");
% grid on;