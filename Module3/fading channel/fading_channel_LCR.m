c0 = 108e7;
f0 = 2000e6;
v = 300;
N = [5, 10, 15, 20];
f_max = f0*v/c0;
t = 0:0.01666:100;

m1 = zeros(size(t));
m2 = zeros(size(t));
for i=1:4
    for n=1:N(i)-1
        c1_in = (2/sqrt(N(i)-1/2))*sin(pi*n/N(i)-1);
        c2_in = (2/sqrt(N(i)-1/2))*cos(pi*n/N(i)-1);
        f_in = f_max*cos(pi*n/(2*(N(i)-1/2)));
        r =(c1_in*cos(2*pi*f_in*(t/1000)));
        m1 = m1 + r;
        m2 = m2 + (c2_in*cos(2*pi*f_in*(t/1000)));
    end
    m1 = m1 + (2/sqrt(N(i)-1/2))*cos(2*pi*f_max*(t/1000));
    m2 = m2 + (2/sqrt(N(i)-1/2))*cos(2*pi*f_max*(t/1000));
    res = sqrt(m1.^2 + m2.^2);
    res = 20*log10(res);
    figure(1)
    plot(t, res)
    hold on
end
xlabel('t(ms)')
ylabel('received signal(dB)')
ylim([-40 30])
legend('N=5', 'N=10', 'N=15', 'N=20')
hold off