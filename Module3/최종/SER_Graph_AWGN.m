clearvars; clc;
%% mapping BPSK using fadding channel
SNR_db = 0:0.5:20;
SER_bpsk = zeros(size(SNR_db));
SER_qpsk = zeros(size(SNR_db));

t = 0:0.06:6000;

nSym = 40;

PSK_bpsk = zeros(1, nSym);
PSK1_qpsk = zeros(1, nSym);
PSK2_qpsk = zeros(1, nSym);
for i=1:length(SNR_db)
    cnt_error_bpsk = 0;
    cnt_error_qpsk = 0;

    for j=1:40:length(t)-1
        code = randomCodeGenerator();

        %% BPSK
        scale_PSK = zeros(1, nSym);
        for k=1:nSym
            if code(k)==0
                scale_PSK(k) = -1;
            else
                scale_PSK(k) = 1;
            end
        end
    
        PSK_a = scale_PSK.*sqrt(10^(SNR_db(i)/10));
        PSK_bpsk = noise_generator(PSK_a);

        recieved_data = zeros(1, nSym);
        for k=1:nSym
            if PSK_bpsk(k) < 0
                recieved_data(k) = 0;
            else
                recieved_data(k) = 1;
            end
        end
        
        for k=1:nSym
            if recieved_data(k) ~= code(k)
                cnt_error_bpsk = cnt_error_bpsk + 1;
            end
        end

        %% QPSK
        scale_PSK1 = zeros(1, nSym/2);
        scale_PSK2 = zeros(1, nSym/2);
        for k=1:nSym
            if mod(k, 2) == 1
                if code(k)==0
                    scale_PSK1((k+1)/2) = -1;
                else
                    scale_PSK1((k+1)/2) = 1;
                end
            else
                if code(k)==0
                    scale_PSK2(k/2) = -1;
                else
                    scale_PSK2(k/2) = 1;
                end
            end
        end
    
        PSK1_a = scale_PSK1.*(sqrt(10^(SNR_db(i)/10))*sqrt(1/2));
        PSK2_a = scale_PSK2.*(sqrt(10^(SNR_db(i)/10))*sqrt(1/2));
        PSK1_qpsk = noise_generator(PSK1_a);
        PSK2_qpsk = noise_generator(PSK2_a);

        recieved_data = zeros(1, nSym);
        for k=1:nSym/2
            if PSK1_qpsk(k) < 0
                recieved_data(k*2-1) = 0;
            else
                recieved_data(k*2-1) = 1;
            end
            
            if PSK2_qpsk(k) < 0
                recieved_data(k*2) = 0;
            else
                recieved_data(k*2) = 1;
            end
        end
        
        for k=1:nSym
            if mod(k, 2) == 1
                eSym = 0;
            end
            if recieved_data(k) ~= code(k)
                if eSym == 0
                    cnt_error_qpsk = cnt_error_qpsk + 1;
                    eSym = 1;
                end
            end
        end

    end
    fprintf('%d    cnt_error : %d\n', i, cnt_error_bpsk);
    fprintf('          full : %d\n', length(t)-1);
    SER_bpsk(i) = cnt_error_bpsk/(length(t));
    SER_bpsk(i) = 10*log10(SER_bpsk(i));
    SER_qpsk(i) = cnt_error_qpsk/(length(t)/2);
    SER_qpsk(i) = 10*log10(SER_qpsk(i));
end

%% plot SER graph
figure(1)
plot(SNR_db, SER_bpsk)
hold on
plot(SNR_db, SER_qpsk)
xlabel('SNR(dB)')
ylabel('SER')
xlim([0 12])
ylim([-40 0])
legend('BPSK', 'QPSK')
hold off




%% add noise
function n = noise_generator(s)
    noise = randn(1, length(s))*(1/sqrt(2));
    n = s + noise;
end


function hamming_code = randomCodeGenerator()
    %% random data 1 Frame
    random_data = randi([0 1], 1, 40);
    hamming_code = random_data;
end