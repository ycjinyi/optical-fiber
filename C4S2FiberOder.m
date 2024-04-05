clc;
clear;
close all;

%颜色表和标签
CG = ColorGenerator([0 1 1; 1 0 1; 1 1 0]);
all = 256;
A = (1: 1: 161)';
% idx = floor(A * 255 / max(A));

[colorTable, lambdaStr] = CG.generate(A');

colorTable = floor(colorTable * 255);

writematrix(colorTable, "161根发射光纤排序颜色");