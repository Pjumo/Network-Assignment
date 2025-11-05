[x, y] = meshgrid(0:0.05:5);
txs.x = rand(1, 10)*5;
txs.y = rand(1, 10)*5;
tx_h = 50;
rx_h = 1;
fq = 28000;
N = 0.001;
bw = 40e6;

%%PathLoss
d=zeros(101);
for i = 1:10
    d(:, :, i)= sqrt((x-txs.x(i)).^2 + (y-txs.y(i)).^2);
end
%%PL_cost = 46.3 + 33.9 * log10(fq) - 13.82*log10(tx_h) - 3.2 * (log10(11.75 * rx_h))^2 - 4.97 + (44.9 - 6.55 * log10(tx_h)) * log10(d) + 3;
PL_sui = 20*log10(4*pi*100*fq)-147.55+10*(3.6-0.005*tx_h+20/tx_h)*log10(d/100)+6*log10(fq/2000)-20*log10(rx_h/2);

PL_max=zeros(101);
for c = 1:1:101
    for r = 1:1:101
        PL_max(c, r) = 10000;
        for i = 1:10
            if PL_max(c, r) > PL_sui(c, r, i)
                PL_max(c, r) = PL_sui(c, r, i);
            end
        end
    end
end

PL_else=zeros(101, 101, 9);
for c = 1:1:101
    for r = 1:1:101
        j=1;
        for i = 1:10
            if PL_sui(c, r, i) ~= PL_max(c, r)
                PL_else(c, r, j) = PL_sui(c, r, i);
                j=j+1;
            end
        end
    end
end

signal_max = 46.021 - PL_max;
signal_else = 46.021 - PL_else;

%%SINR
SINR_w = zeros(101);
for c = 1:1:101
    for r = 1:1:101
        itf = 0;
        for i=1:9
            itf = itf + 10^(signal_else(c, r, i)/10);
        end
        SINR_w(c, r) = 10^(signal_max(c, r)/10)/itf + N;
    end
end
SINR = 10*log10(SINR_w);

%%DataRate
DataRate = bw*log10(1+SINR_w);

%%3D
surf(x, y, DataRate)
xlabel('X')
ylabel('Y')
zlabel('DataRate')

%%등고선
% contour(x, y, DataRate, 10)

%%등고선 특정 위치
% contour(x, y, DataRate, [0 3 50], "ShowText", true, "LabelFormat", "%d dB")