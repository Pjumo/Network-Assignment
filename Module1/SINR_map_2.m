tx1.x = 0;
tx1.y = 0;
tx2.x = 2.7*sqrt(2);
tx2.y = 2.7*sqrt(2);
tx_h = 50;
rx_h = 1;
fq=800;
[x, y] = meshgrid(-5:0.1:10);
d_1 = sqrt((x-tx1.x).^2 + (y-tx1.y).^2);
d_2 = sqrt((x-tx2.x).^2 + (y-tx2.y).^2);
d_close = d_1;
d_far = d_2;
for c = 1:1:151
    for r = 1:1:151
        if d_1(c, r)>d_2(c, r)
            d_close(c, r) = d_2(c, r);
            d_far(c, r) = d_1(c, r);
        end
    end
end
PL_hata_1 = 69.55 + 26.16*log10(fq) - 13.82*log10(tx_h) - 3.2 * (log10(11.75 * rx_h))^2 - 4.97 + (44.9 - 6.55 * log10(tx_h)) * log10(d_close);
PL_hata_2 = 69.55 + 26.16*log10(fq) - 13.82*log10(tx_h) - 3.2 * (log10(11.75 * rx_h))^2 - 4.97 + (44.9 - 6.55 * log10(tx_h)) * log10(d_far);
PL_cost_1 = 46.3 + 33.9 * log10(fq) - 13.82*log10(tx_h) - 3.2 * (log10(11.75 * rx_h))^2 - 4.97 + (44.9 - 6.55 * log10(tx_h)) * log10(d_close) + 3;
PL_cost_2 = 46.3 + 33.9 * log10(fq) - 13.82*log10(tx_h) - 3.2 * (log10(11.75 * rx_h))^2 - 4.97 + (44.9 - 6.55 * log10(tx_h)) * log10(d_far) + 3;
PL_sui_1 = 20*log10(4*pi*100*fq)-147.55+10*(3.6-0.005*tx_h+20/tx_h)*log10(d_close/100)+6*log10(fq/2000)-20*log10(rx_h/2);
PL_sui_2 = 20*log10(4*pi*100*fq)-147.55+10*(3.6-0.005*tx_h+20/tx_h)*log10(d_far/100)+6*log10(fq/2000)-20*log10(rx_h/2);
rx_tx1 = 23.01 - PL_hata_1;
rx_tx2 = 23.01 - PL_hata_2;
rx_tx1(isinf(rx_tx1)) = 125;
rx_tx2(isinf(rx_tx2)) = 125;
z = rx_tx1 - rx_tx2;
contour(x, y, rx_tx1, [0 -105], "ShowText", true, "LabelFormat", "%d dB")
hold on
contour(x, y, rx_tx2, [0 -105], "ShowText", true, "LabelFormat", "%d dB")
hold on
contour(x, y, z, [0 3], "ShowText", true, "LabelFormat", "%d dB")
hold off