data = load('2dproj_raw.xvg');
d1 = data(1:10001, :);
d2 = data(10002:end, :);


xBininfo = 15;
yBininfo = 15;

d = d1;
Bins = getBin(d, xBininfo, yBininfo);
[xGridNum, yGridNum] = size(Bins);
x = d(:, 1);
y = d(:, 2);
X = linspace(min(x), max(x), xGridNum);
Y = linspace(min(y), max(y), yGridNum); 
contour(X, Y, Bins, 'Linecolor', 'r');
hold on;                        

d = d2;
Bins = getBin(d, xBininfo, yBininfo);
[xGridNum, yGridNum] = size(Bins);
x = d(:, 1);
y = d(:, 2);
X = linspace(min(x), max(x), xGridNum);
Y = linspace(min(y), max(y), yGridNum); 
contour(X, Y, Bins, 'Linecolor', 'b');
hold on;                                                                         
%{
d = d3;
Bins = getBin(d, xBininfo, yBininfo);
[xGridNum, yGridNum] = size(Bins);
x = d(:, 1);
y = d(:, 2);
X = linspace(min(x), max(x), xGridNum);
Y = linspace(min(y), max(y), yGridNum); 
contour(X, Y, Bins, 'Linecolor', 'g');
hold on;                                 

d = d4;
Bins = getBin(d, xBininfo, yBininfo);
[xGridNum, yGridNum] = size(Bins);
x = d(:, 1);
y = d(:, 2);
X = linspace(min(x), max(x), xGridNum);
Y = linspace(min(y), max(y), yGridNum); 
contour(X, Y, Bins, 'Linecolor', 'k');
hold on; 
%}
