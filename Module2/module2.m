[x, y] = meshgrid(0:0.05:5);
[x1, y1] = meshgrid(0:0.05:5);
[x2, y2] = meshgrid(0:0.05:5);
[x3, y3] = meshgrid(0:0.05:5);
[x4, y4] = meshgrid(0:0.05:5);
txs.x1 = [0.5, 2.5, 4.5];
txs.y1 = [1.5, 3.5, 1.5];
txs.x2 = [0.5, 2.5, 4.5];
txs.y2 = [3.5, 1.5, 3.5];
txs.x3 = [1.5, 3.5];
txs.y3 = [0.5, 4.5];
txs.x4 = [1.5, 3.5];
txs.y4 = [4.5, 0.5];

avg_txs = zeros(10, 10000);
cnt_avg = zeros(10, 1);
max_tx = zeros(101);

tx_h = 50;
rx_h = 1;
fq = 28000;
N = 0.001;
bw = 10e6;
rx_power = 43.01;

%%PathLoss
d1=zeros(101);
d2=zeros(101);
d3=zeros(101);
d4=zeros(101);
for i = 1:3
    d1(:, :, i)= sqrt((x1-txs.x1(i)).^2 + (y1-txs.y1(i)).^2);
    d2(:, :, i)= sqrt((x2-txs.x2(i)).^2 + (y2-txs.y2(i)).^2);
end
for i = 1:2
    d3(:, :, i)= sqrt((x3-txs.x3(i)).^2 + (y3-txs.y3(i)).^2);
    d4(:, :, i)= sqrt((x4-txs.x4(i)).^2 + (y4-txs.y4(i)).^2);
end
% PL_sui1 = 46.3 + 33.9 * log10(fq) - 13.82*log10(tx_h) - 3.2 * (log10(11.75 * rx_h))^2 - 4.97 + (44.9 - 6.55 * log10(tx_h)) * log10(d1) + 3;
PL_sui1 = 20*log10(4*pi*100*fq)-147.55+10*(3.6-0.005*tx_h+20/tx_h)*log10(d1/100)+6*log10(fq/2000)-20*log10(rx_h/2);
% PL_sui2 = 46.3 + 33.9 * log10(fq) - 13.82*log10(tx_h) - 3.2 * (log10(11.75 * rx_h))^2 - 4.97 + (44.9 - 6.55 * log10(tx_h)) * log10(d2) + 3;
PL_sui2 = 20*log10(4*pi*100*fq)-147.55+10*(3.6-0.005*tx_h+20/tx_h)*log10(d2/100)+6*log10(fq/2000)-20*log10(rx_h/2);
% PL_sui3 = 46.3 + 33.9 * log10(fq) - 13.82*log10(tx_h) - 3.2 * (log10(11.75 * rx_h))^2 - 4.97 + (44.9 - 6.55 * log10(tx_h)) * log10(d3) + 3;
PL_sui3 = 20*log10(4*pi*100*fq)-147.55+10*(3.6-0.005*tx_h+20/tx_h)*log10(d3/100)+6*log10(fq/2000)-20*log10(rx_h/2);
% PL_sui4 = 46.3 + 33.9 * log10(fq) - 13.82*log10(tx_h) - 3.2 * (log10(11.75 * rx_h))^2 - 4.97 + (44.9 - 6.55 * log10(tx_h)) * log10(d4) + 3;
PL_sui4 = 20*log10(4*pi*100*fq)-147.55+10*(3.6-0.005*tx_h+20/tx_h)*log10(d4/100)+6*log10(fq/2000)-20*log10(rx_h/2);

PL_max1=zeros(101);
PL_max2=zeros(101);
for c = 1:1:101
    for r = 1:1:101
        PL_max1(c, r) = 10000;
        PL_max2(c, r) = 10000;
        for i = 1:3
            if PL_max1(c, r) > PL_sui1(c, r, i)
                PL_max1(c, r) = PL_sui1(c, r, i);
                max_tx(c, r, 1) = i;
            end
            if PL_max2(c, r) > PL_sui2(c, r, i)
                PL_max2(c, r) = PL_sui2(c, r, i);
                max_tx(c, r, 2) = i;
            end
        end
    end
end
PL_max3=zeros(101);
PL_max4=zeros(101);
for c = 1:1:101
    for r = 1:1:101
        PL_max3(c, r) = 10000;
        PL_max4(c, r) = 10000;
        for i = 1:2
            if PL_max3(c, r) > PL_sui3(c, r, i)
                PL_max3(c, r) = PL_sui3(c, r, i);
                max_tx(c, r, 3) = i;
            end
            if PL_max4(c, r) > PL_sui4(c, r, i)
                PL_max4(c, r) = PL_sui4(c, r, i);
                max_tx(c, r, 4) = i;
            end
        end
    end
end

PL_else1=zeros(101, 101, 2);
PL_else2=zeros(101, 101, 2);
for c = 1:1:101
    for r = 1:1:101
        j=1;
        k=1;
        for i = 1:3
            if PL_sui1(c, r, i) ~= PL_max1(c, r)
                PL_else1(c, r, j) = PL_sui1(c, r, i);
                j=j+1;
            end
            if PL_sui2(c, r, i) ~= PL_max2(c, r)
                PL_else2(c, r, k) = PL_sui2(c, r, i);
                k=k+1;
            end
        end
        if j==1
            PL_else1(c, r, 1) = PL_max1(c, r);
            PL_else1(c, r, 2) = PL_max1(c, r);
        end
        if k==1
            PL_else2(c, r, 1) = PL_max2(c, r);
            PL_else2(c, r, 2) = PL_max2(c, r);
        end
        if j==2
            PL_else1(c, r, 2) = PL_max1(c, r);
        end
        if k==2
            PL_else2(c, r, 2) = PL_max2(c, r);
        end
    end
end
PL_else3=zeros(101, 101, 1);
PL_else4=zeros(101, 101, 1);
for c = 1:1:101
    for r = 1:1:101
        j=1;
        k=1;
        for i = 1:2
            if PL_sui3(c, r, i) ~= PL_max3(c, r)
                PL_else3(c, r, j) = PL_sui3(c, r, i);
                j=j+1;
            end
            if PL_sui4(c, r, i) ~= PL_max4(c, r)
                PL_else4(c, r, k) = PL_sui4(c, r, i);
                k=k+1;
            end
        end
        if j==1
            PL_else3(c, r, 1) = PL_max3(c, r);
        end
        if k==1
            PL_else4(c, r, 1) = PL_max4(c, r);
        end
    end
end

signal_max1 = rx_power - PL_max1;
signal_else1 = rx_power - PL_else1;
signal_max2 = rx_power - PL_max2;
signal_else2 = rx_power - PL_else2;
signal_max3 = rx_power - PL_max3;
signal_else3 = rx_power - PL_else3;
signal_max4 = rx_power - PL_max4;
signal_else4 = rx_power - PL_else4;

%%SINR
SINR_w1 = zeros(101);
SINR_w2 = zeros(101);
SINR_w3 = zeros(101);
SINR_w4 = zeros(101);
for c = 1:1:101
    for r = 1:1:101
        itf1 = 0;
        itf2 = 0;
        for i=1:2
            itf1 = itf1 + 10^(signal_else1(c, r, i)/10);
            itf2 = itf2 + 10^(signal_else2(c, r, i)/10);
        end
        SINR_w1(c, r) = 10^(signal_max1(c, r)/10)/itf1 + N;
        SINR_w2(c, r) = 10^(signal_max2(c, r)/10)/itf2 + N;
        SINR_w3(c, r) = 10^(signal_max3(c, r)/10)/10^(signal_else3(c, r, 1)/10) + N;
        SINR_w4(c, r) = 10^(signal_max4(c, r)/10)/10^(signal_else4(c, r, 1)/10) + N;
    end
end

SINR_w = zeros(101);
for c = 1:1:101
    for r = 1:1:101
        if SINR_w1(c, r)>=SINR_w2(c, r) && SINR_w1(c, r)>=SINR_w3(c, r) && SINR_w1(c, r)>=SINR_w4(c, r)
            if isinf(SINR_w1(c, r))
                SINR_w1(c, r) = SINR_w1(c-1, r-1);
            end
            SINR_w(c, r) = SINR_w1(c, r);

            if max_tx(c, r, 1) == 1
                cnt_avg(3) = cnt_avg(3) + 1;
                avg_txs(3, cnt_avg(3)) = bw*log2(1+SINR_w(c, r));
            elseif max_tx(c, r, 1) == 2
                cnt_avg(7) = cnt_avg(7) + 1;
                avg_txs(7, cnt_avg(7)) = bw*log2(1+SINR_w(c, r));
            else
                cnt_avg(5) = cnt_avg(5) + 1;
                avg_txs(5, cnt_avg(5)) = bw*log2(1+SINR_w(c, r));
            end

        elseif SINR_w2(c, r)>=SINR_w1(c, r) && SINR_w2(c, r)>=SINR_w3(c, r) && SINR_w2(c, r)>=SINR_w4(c, r)
            if isinf(SINR_w2(c, r))
                SINR_w2(c, r) = SINR_w2(c-1, r-1);
            end
            SINR_w(c, r) = SINR_w2(c, r);

            if max_tx(c, r, 2) == 1
                cnt_avg(6) = cnt_avg(6) + 1;
                avg_txs(6, cnt_avg(6)) = bw*log2(1+SINR_w(c, r));
            elseif max_tx(c, r, 2) == 2
                cnt_avg(4) = cnt_avg(4) + 1;
                avg_txs(4, cnt_avg(4)) = bw*log2(1+SINR_w(c, r));
            else
                cnt_avg(8) = cnt_avg(8) + 1;
                avg_txs(8, cnt_avg(8)) = bw*log2(1+SINR_w(c, r));
            end

        elseif SINR_w3(c, r)>=SINR_w1(c, r) && SINR_w3(c, r)>=SINR_w2(c, r) && SINR_w3(c, r)>=SINR_w4(c, r)
            if isinf(SINR_w3(c, r))
                SINR_w3(c, r) = SINR_w3(c-1, r-1);
            end
            SINR_w(c, r) = SINR_w3(c, r);

            if max_tx(c, r, 3) == 1
                cnt_avg(1) = cnt_avg(1) + 1;
                avg_txs(1, cnt_avg(1)) = bw*log2(1+SINR_w(c, r));
            else
                cnt_avg(10) = cnt_avg(10) + 1;
                avg_txs(10, cnt_avg(10)) = bw*log2(1+SINR_w(c, r));
            end

        else
            if isinf(SINR_w4(c, r))
                SINR_w4(c, r) = SINR_w4(c-1, r-1);
            end
            SINR_w(c, r) = SINR_w4(c, r);

            if max_tx(c, r, 4) == 1
                cnt_avg(2) = cnt_avg(2) + 1;
                avg_txs(2, cnt_avg(2)) = bw*log2(1+SINR_w(c, r));
            else
                cnt_avg(9) = cnt_avg(9) + 1;
                avg_txs(9, cnt_avg(9)) = bw*log2(1+SINR_w(c, r));
            end

        end
    end
end
SINR = 10*log10(SINR_w4);

%%DataRate
DataRate = bw*log2(1+SINR_w);

for i = 1:1:10
    minDR_tx = mink(nonzeros(avg_txs(i, :)), round(cnt_avg(i)/10));
    maxDR_tx = maxk(nonzeros(avg_txs(i, :)), round(cnt_avg(i)/10));
    totalDR_tx_min = sum(minDR_tx);
    totalDR_tx_max = sum(maxDR_tx);
    totalDR_tx_average = sum(avg_txs(i, :));
    fprintf("%d_min : %f\n", i, totalDR_tx_min/round(cnt_avg(i)/10));
    fprintf("%d_max : %f\n", i, totalDR_tx_max/round(cnt_avg(i)/10));
    fprintf("%d_average : %f\n", i, totalDR_tx_average/cnt_avg(i));
end

dum = reshape(DataRate, [1, 10201]);
minDR = mink(dum, 1020);
maxDR = maxk(dum, 1020);
totalDR_min = sum(minDR);
totalDR_max = sum(maxDR);
totalDR = sum(dum);
fprintf("total_min : %f\n", totalDR_min/1020);
fprintf("total_max : %f\n", totalDR_max/1020);
fprintf("total_average : %f\n", totalDR/10201);

%%3D
% surf(x, y, DataRate)
% xlabel('X')
% ylabel('Y')
% zlabel('DataRate')

%%등고선
% contour(x, y, DataRate, 50)

%%등고선 특정 위치
% contour(x, y, DataRate, [0 3 50], "ShowText", true, "LabelFormat", "%d dB")