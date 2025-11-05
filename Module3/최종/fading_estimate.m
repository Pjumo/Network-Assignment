c0 = 108e7;
f0 = 2400e6;
v = 100;
N = 7;
f_max = f0*v/c0;
t = 0:0.06:100;

m1 = zeros(size(t));
m2 = zeros(size(t));
for n=1:N-1
    c1_in = (2/sqrt(N-1/2))*sin(pi*n/N-1);
    c2_in = (2/sqrt(N-1/2))*cos(pi*n/N-1);
    f_in = f_max*cos(pi*n/(2*(N-1/2)));
    m1 = m1 + (c1_in*cos(2*pi*f_in*(t/1000)));
    m2 = m2 + (c2_in*cos(2*pi*f_in*(t/1000)));
end
m1 = m1 + (2/sqrt(N-1/2))*cos(2*pi*f_max*(t/1000));
m2 = m2 + (2/sqrt(N-1/2))*cos(2*pi*f_max*(t/1000));

m = sqrt(m1.^2 + m2.^2);

es_t = zeros(1, fix(length(t)/8)+1);
for i=1:length(t)
    if mod(i, 8)==1
        es_t(fix((i-1)/8)+1) = m(i);
    end
end

es_m = zeros(1, length(t));
for i=1:length(es_t)-1
    pos = (i-1)*8+1;
    rat = (es_t(i+1) - es_t(i))/8;
    es_m(pos) = es_t(i);
    for j=1:7
        es_m(pos+j) = es_t(i) + rat*j;
    end
end

es_m = 20*log10(es_m);
m = 20*log10(m);


figure(1)
plot(t, m)
xlabel('t(ms)')
ylabel('received signal(dB)')
ylim([-40 30])
% hold on
% plot(t, es_m)
% xlabel('t(ms)')
% ylabel('received signal(dB)')
% ylim([-40 30])
% legend('simulated', 'estimated')
% hold off
