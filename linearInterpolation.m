function [data] = linearInterpolation(originData, dataPoints)
%originData为原始数据,至少2列,第一列为坐标,后续列为数据
%dataPoints为需要插值的数据点,1列,为坐标
%data为按照目标数据点插值得到的数据,至少2列,第一列和dataPoints一致,后续列为对应的数据
%这里面需要保证至少有2行数据，否则无法外推
    orgRow = size(originData, 1);
    dpRow = size(dataPoints, 1);
    data = zeros(dpRow, size(originData, 2));
    data(:, 1) = dataPoints;
    k = 1;
    for i = 1: dpRow
        %在originData中找到第一个大等于当前数的位置
        while(k <= orgRow && originData(k, 1) < data(i, 1))
            k = k + 1;
        end
        %需要外推的情况
        if(k == 1 || k > orgRow)
            p1 = 2;
            p2 = 1;
            if(k > orgRow)
                p1 = orgRow - 1;
                p2 = orgRow;
            end
            x1 = originData(p1, 1);
            y1 = originData(p1, 2: end);
            x2 = originData(p2, 1);
            y2 = originData(p2, 2: end);
            x3 = data(i, 1);
            data(i, 2: end) = y1 + (x3 - x1) * (y2 - y1) / (x2 - x1);
            continue;
        end
        % 相同就直接赋值了
        if (originData(k, 1) == data(i, 1))
            data(i, 2: end) = originData(k, 2: end);
            continue;
        end
        %内插
        x1 = originData(k - 1, 1);
        y1 = originData(k - 1, 2: end);
        x3 = originData(k, 1);
        y3 = originData(k, 2: end);
        x2 = data(i, 1);
        data(i, 2: end) = y1 + (x2 - x1) * (y3 - y1) / (x3 - x1);
    end
end