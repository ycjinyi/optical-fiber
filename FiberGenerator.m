classdef FiberGenerator
%此函数用于生成光纤相关的数据
    properties
       %实际光纤的半径数据,单位mm
       r = 1e-3 * 187.562 / 2;
       %光纤建模时像素和实际尺寸之间的关系,单位px/mm
       %光纤探测端面的系数
       coff1 = 108.64404;
       %光纤信号端面的系数
       coff2 = 141.7;
       %实际光纤数据的文件名
       ftName = '发射光纤探测端.txt';
       fr1Name = '接收光纤1.txt';
       fr2Name = '接收光纤2.txt';
       fr3Name = '接收光纤3.txt';
       fsName = '发射光纤信号端.txt'
    end

    methods

        function obj = FiberGenerator()
        end

        %------------------------产生光纤xy坐标的函数-----------------------
        
        %1束发射光纤2束接收光纤的排布
        %coreDis是纤芯之间的距离, number是横向光纤的个数(必须是偶数)
        %s是所有发射光纤的位置, r1是所有接收光纤1的位置, r2是所有接收光纤2的位置
        %一共有2行, 第一行是x坐标(横向), 第二行是y坐标(纵向)
        %发射光纤和接收光纤1混合排布,接收光纤2阶梯排布
        function [s, r1, r2] = regFiberGen(~, coreDis, number)
            s = zeros(2, number * number);
            r1 = zeros(2, number * number);
            r2 = zeros(2, number * number);
            %混合光纤分簇的个数
            num = number / 2;
            %发射光纤
            c = 1;
            for k = 1: 2
                d = (k - 1) * number;
                for i = num + 1 + d: num + num + d
                    for j = 1: number
                        s(1, c) = (j - 1) * coreDis;
                        s(2, c) = (i - 1) * coreDis;
                        c = c + 1;
                    end
                end
            end
            %接收光纤1
            c = 1;
            for k = 1: 2
                d = (k - 1) * number;
                for i = 1 + d: num + d
                    for j = 1: number
                        r1(1, c) = (j - 1) * coreDis;
                        r1(2, c) = (i - 1) * coreDis;
                        c = c + 1;
                    end
                end
            end
            %接收光纤2
            c = 1;
            for i = 2 * number + 1: 2 * number + number
                for j = 1: number
                    r2(1, c) = (j - 1) * coreDis;
                    r2(2, c) = (i - 1) * coreDis;
                    c = c + 1;
                end
            end
        end

        %单根发射光纤单根接收光纤的排布
        function [s, r] = singleFiberGen(~, coreDis)
            %发射光线作为原点
            s = [0; 0];
            %接收光线与发射光线y坐标相同
            r = [coreDis; 0];
        end
        
        %-----------------将xy坐标转换为实际的posMatrix的函数----------------
        
        %此函数将发射光纤和接收光纤的xy坐标转换为posMatrix
        %posMatrix有3个维度[a, b, c]
        %a的取值是1和2,1代表距离之间的关系,2代表角度之间的关系,注意是弧度制
        %b代表对应的发射光纤编号,c代表对应的接收光纤编号
        %s,r分别是光源和接收光纤的中心位置,2行,第1行是x坐标,第2行是y坐标
        function [posMatrix] = posConvert(obj, s, r)
            sNum = size(s, 2);
            rNum = size(r, 2);
            posMatrix = zeros(2, sNum, rNum);
            for i = 1: sNum
                x1 = s(1, i);
                y1 = s(2, i);
                for j = 1: rNum
                    x2 = r(1, j);
                    y2 = r(2, j);
                    dx = x2 - x1;
                    dy = y2 - y1;
                    dis = sqrt(power(dx, 2) + power(dy, 2));
                    posMatrix(1, i, j) = dis;
                    phi = obj.phiCompute(dx, dy);
                    posMatrix(2, i, j) = phi;
                end
            end
        end

        %根据当前点相对于目标原点的坐标dx, dy计算得到当前点的角度,注意是弧度
        %角度的范围是0-2*pi
        function [phi] = phiCompute(~, dx, dy)
            %考虑特殊情况
            if abs(dx) < 1e-10
                if dy >= 0
                    phi = 0.5 * pi;
                else
                    phi = 1.5 * pi;
                end
                return;
            end
            if abs(dy) < 1e-10
                if dx >= 0
                    phi = 0;
                else
                    phi = pi;
                end
                return;
            end
            %计算正常情况
            phi = atan(dy / dx);
            if dx > 0
                if dy > 0
                    return;
                end
                phi = phi + 2 * pi;
            else
                phi = phi + pi;
            end
        end

        %--------对实际的光纤端面建模,得到其所有接收光纤束的posMatrix----------
        
        function [SPM, R1PM, R2PM, R3PM] = realFiberGen(obj)
            %首先读取实际光纤的数据, 第一列是x坐标, 第二列是y坐标, 第三列是半径直接忽略
            ft = importdata(obj.ftName);
            fr1 = importdata(obj.fr1Name);
            fr2 = importdata(obj.fr2Name);
            fr3 = importdata(obj.fr3Name);
            fs = importdata(obj.fsName);
            %将数据转换为需要的格式(m)
            ft = 1e-3 * ft(:, 1: 2)' / obj.coff1;
            fr1 = 1e-3 * fr1(:, 1: 2)' / obj.coff1;
            fr2 = 1e-3 * fr2(:, 1: 2)' / obj.coff1;
            fr3 = 1e-3 * fr3(:, 1: 2)' / obj.coff1;
            fs = 1e-3 * fs(:, 1: 2)' / obj.coff2;
            %依次进行转换计算
            R1PM = obj.posConvert(ft, fr1);
            R2PM = obj.posConvert(ft, fr2);
            R3PM = obj.posConvert(ft, fr3);
            %对于信号光纤而言,将光源设置为所有光纤中心位置的均值
            %偏心光源
            %SPM = obj.posConvert(1e-3 * [506; 301] / obj.coff2, fs);
            SPM = obj.posConvert(mean(fs, 2), fs);
        end
       
    end
end