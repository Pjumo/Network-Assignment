clearvars; clc;
%% mapping BPSK using fadding channel
SNR_db = 0:0.5:20;
FER = zeros(size(SNR_db));

[m1, m2] = fadingChannel();

nSym = 98;
nPilot = 7;
nSum = nSym+nSym/nPilot+1;

for i=1:length(SNR_db)
    cnt_error = 0;
    frame_count = 0;

    %% estimated channel
    PSK = zeros(1, nSum);
    for j=1:nSum:length(m1)
        if length(m1)-j<nSum
            break;
        end
        frame_count=frame_count+1;
        code = randomCodeGenerator();
        % code = reshape(code, 7, 14).';
        % code = reshape(code, 1, 98);
        scale_PSK = zeros(1, nSum);
        scale_PSK(1) = 1;
        pilot_count = 1;
        for k=1:nSym
            if code(k)==0
                scale_PSK(k+pilot_count) = -1;
            else
                scale_PSK(k+pilot_count) = 1;
            end

            if mod(k, nPilot) == 0
                pilot_count = pilot_count + 1;
                scale_PSK(k+pilot_count) = 1;
            end
        end

        PSK1_a = scale_PSK.*sqrt(10^(SNR_db(i)/10));
        m = complex(m1, m2);

        %% Generate noise
        for k=1:nSum
            PSK(k) = PSK1_a(k) * m(j+k-1) + (1/sqrt(2))*randn;
        end
        
        %% channel estimate
        estimate_pilot = zeros(1, nSym/nPilot+1);
        uncorrect_data = zeros(1, nSym);
        data_count = 1;
        for k=1:nSum
            if mod(k, nPilot+1) == 1
                estimate_pilot(fix(k/(nPilot+1))+1) = PSK(k) / sqrt(10^(SNR_db(i)/10));
            else
                uncorrect_data(data_count) = PSK(k);
                data_count = data_count + 1;
            end
        end

        inter_es = zeros(1, nSym/nPilot);
        for k=1:nSym/nPilot
            inter_es(k) = (estimate_pilot(k+1)-estimate_pilot(k))/8;
        end
        % fprintf('inter_es %s\n', mat2str(inter_es));

        recieved_data = zeros(1, nSym);
        for k=1:nSym 
            m_pos = fix((k-1)/7)+1;
            % fprintf('m_pos %d\n', m_pos);
            m_es = estimate_pilot(m_pos) + (mod(k-1, 7)+1) * inter_es(m_pos);
            % fprintf('m %f\n', m);
            if uncorrect_data(k)/m_es < 0
                recieved_data(k) = 0;
            else
                recieved_data(k) = 1;
            end
        end

        % recieved_data = reshape(recieved_data, 14, 7).';
        % recieved_data = reshape(recieved_data, 1, 98);

        cnt_error=cnt_error+check_frame_error(recieved_data);
    end
    fprintf('%d    cnt_error : %d\n', i, cnt_error);
    fprintf('          full : %d\n', frame_count);
    FER(i) = cnt_error/frame_count;
    FER(i) = log10(FER(i));
end

%% plot SER graph
figure(1)
plot(SNR_db, FER)
yticks([-3 -2 -1 0])
ytickformat('10e%d');
xlabel('SNR(dB)')
ylabel('SER(log10)')
ylim([-3 0])
xlim([0 20])
hold on





%% fadding channel
function [m1, m2] = fadingChannel()
    c0 = 108e7;
    f0 = 2400e6;
    v = 60;
    N = 13;
    f_max = f0*v/c0;
    t = 600:0.06:6600;
    
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
end


%% check error hamming and crc
function error = check_frame_error(check_code)
    error = 0;
    code = reshape(check_code, 7, 14);
    code_str = num2str(code);
    str_c = '';
    for i=1:40
        for j=1:7
            if code_str(j, i) ~= ' '
                str_c = strcat(str_c, code_str(j, i));
            end
        end
    end

    code = reshape(str_c, 7, 14);

    crc_code = zeros(4, 14);
    for i=1:14
        hd = code(:, i).';
        e4 = xor(xor(str2double(hd(4)), str2double(hd(5))), xor(str2double(hd(6)), str2double(hd(7))));
        e2 = xor(xor(str2double(hd(2)), str2double(hd(3))), xor(str2double(hd(6)), str2double(hd(7))));
        e1 = xor(xor(str2double(hd(1)), str2double(hd(3))), xor(str2double(hd(5)), str2double(hd(7))));
        if e4==0&&e2==0&&e1==0
            crc_code(:, i) = [str2double(hd(3)), str2double(hd(5)), str2double(hd(6)), str2double(hd(7))];
        elseif e4==0&&e2==0&&e1==1
            crc_code(:, i) = [str2double(hd(3)), str2double(hd(5)), str2double(hd(6)), str2double(hd(7))];
        elseif e4==0&&e2==1&&e1==0
            crc_code(:, i) = [str2double(hd(3)), str2double(hd(5)), str2double(hd(6)), str2double(hd(7))];
        elseif e4==0&&e2==1&&e1==1
            crc_code(:, i) = [not(str2double(hd(3))), str2double(hd(5)), str2double(hd(6)), str2double(hd(7))];
        elseif e4==1&&e2==0&&e1==0
            crc_code(:, i) = [str2double(hd(3)), str2double(hd(5)), str2double(hd(6)), str2double(hd(7))];
        elseif e4==1&&e2==0&&e1==1
            crc_code(:, i) = [str2double(hd(3)), not(str2double(hd(5))), str2double(hd(6)), str2double(hd(7))];
        elseif e4==1&&e2==1&&e1==0
            crc_code(:, i) = [str2double(hd(3)), str2double(hd(5)), not(str2double(hd(6))), str2double(hd(7))];
        elseif e4==1&&e2==1&&e1==1
            crc_code(:, i) = [str2double(hd(3)), str2double(hd(5)), str2double(hd(6)), not(str2double(hd(7)))];
        end
    end


    crc_code = reshape(crc_code, 28, 2);
    remainder = bin2dec(num2str(crc_code(:, 1).'));
    data_length = 28;

    divisor_list = [1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
    divisor = bin2dec(num2str(divisor_list));
    divisor_degree = 16;

    divisor = bitshift(divisor, data_length-divisor_degree-1);

    for i=1:data_length-divisor_degree
        if bitget(remainder, data_length)
            remainder = bitxor(remainder, divisor);
        end
        remainder = bitshift(remainder, 1);
        % fprintf('%2dth Remainder : ', i);
        % 
        % for j = 1:(data_length - i - ...
        %     strlength(dec2bin(bitshift(remainder, -i))))
        %     fprintf('0');
        % end
        % fprintf('%s\n', dec2bin(bitshift(remainder, -i)));
    end

    remainder = bitshift(remainder, 16) + bin2dec(num2str(crc_code(:, 2).'));
    data_length = 44;

    divisor_list = [1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
    divisor = bin2dec(num2str(divisor_list));

    divisor = bitshift(divisor, data_length-divisor_degree-1);

    for i=1:data_length-divisor_degree
        if bitget(remainder, data_length)
            remainder = bitxor(remainder, divisor);
        end
        remainder = bitshift(remainder, 1);
        % fprintf('%2dth Remainder : ', i);
        % 
        % for j = 1:(data_length - i - ...
        %     strlength(dec2bin(bitshift(remainder, -i))))
        %     fprintf('0');
        % end
        % fprintf('%s\n', dec2bin(bitshift(remainder, -i)));
    end

    if remainder ~= 0
        error = 1;
    end
end


function hamming_code = randomCodeGenerator()
    %% random data 1 Frame
    random_data = randi([0 1], 1, 40);
    data = bin2dec(num2str(random_data));
    data_length = 40;
    
    %% CRC for 1Frame data
    divisor_list = [1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
    divisor = bin2dec(num2str(divisor_list));
    divisor_degree = 16;
    
    % fprintf('     Real Data : %s\n', dec2bin(data));
    % fprintf('       Divisor : %s\n\n', dec2bin(divisor));
    
    remainder = bitshift(data, divisor_degree);
    divisor = bitshift(divisor, data_length-1);
    % fprintf(' 0th Remainder : %s\n', dec2bin(remainder));
    
    for i = 1:data_length
        if bitget(remainder, data_length+divisor_degree)
            remainder = bitxor(remainder, divisor);
        end
        remainder = bitshift(remainder, 1);
    
        % fprintf('%2dth Remainder : ', i);
        for j = 1:(data_length + divisor_degree - i - ...
            strlength(dec2bin(bitshift(remainder, -i))))
            % fprintf('0');
        end
        % fprintf('%s\n', dec2bin(bitshift(remainder, -i)));
    end
    
    CRC_value = bitshift(remainder,-data_length);
    CRC_string = num2str(dec2bin(CRC_value));
    if strlength(dec2bin(CRC_value))<divisor_degree
        dummy = '';
        for j = 1:(divisor_degree - strlength(dec2bin(CRC_value)))
            dummy = strcat(dummy, '0');
        end
        CRC_string = strcat(dummy, CRC_string);
    end
    % fprintf('\n      CRC code : %s\n', CRC_string);
    
    send_data = strcat(num2str(dec2bin(data)), CRC_string);
    if strlength(send_data) < 56
        dummy = '';
        for i = 1:56-strlength(send_data)
            dummy = strcat(dummy, '0');
        end
        send_data = strcat(dummy, send_data);
    end
    send_data = splitGraphemes(send_data);
    send_data = reshape(send_data, 4, 14);
    
    code = cell2mat(send_data);
    
    %% hamming code for CRC composed data
    hamming_code = zeros(7, 14);
    for i = 1:14
        switch code(:, i).'
            case ['0', '0', '0', '0']
                hamming_code(:, i) = [0, 0, 0, 0, 0, 0, 0];
            case ['0', '0', '0', '1']
                hamming_code(:, i) = [1, 1, 0, 1, 0, 0, 1];
            case ['0', '0', '1', '0']
                hamming_code(:, i) = [0, 1, 0, 1, 0, 1, 0];
            case ['0', '0', '1', '1']
                hamming_code(:, i) = [1, 0, 0, 0, 0, 1, 1];
            case ['0', '1', '0', '0']
                hamming_code(:, i) = [1, 0, 0, 1, 1, 0, 0];
            case ['0', '1', '0', '1']
                hamming_code(:, i) = [0, 1, 0, 0, 1, 0, 1];
            case ['0', '1', '1', '0']
                hamming_code(:, i) = [1, 1, 0, 0, 1, 1, 0];
            case ['0', '1', '1', '1']
                hamming_code(:, i) = [0, 0, 0, 1, 1, 1, 1];
            case ['1', '0', '0', '0']
                hamming_code(:, i) = [1, 1, 1, 0, 0, 0, 0];
            case ['1', '0', '0', '1']
                hamming_code(:, i) = [0, 0, 1, 1, 0, 0, 1];
            case ['1', '0', '1', '0']
                hamming_code(:, i) = [1, 0, 1, 1, 0, 1, 0];
            case ['1', '0', '1', '1']
                hamming_code(:, i) = [0, 1, 1, 0, 0, 1, 1];
            case ['1', '1', '0', '0']
                hamming_code(:, i) = [0, 1, 1, 1, 1, 0, 0];
            case ['1', '1', '0', '1']
                hamming_code(:, i) = [1, 0, 1, 0, 1, 0, 1];
            case ['1', '1', '1', '0']
                hamming_code(:, i) = [0, 0, 1, 0, 1, 1, 0];
            case ['1', '1', '1', '1']
                hamming_code(:, i) = [1, 1, 1, 1, 1, 1, 1];
        end
    end
    
    hamming_code = reshape(hamming_code, 1, 98);
    
    % fprintf('  Hamming Code : %s\n', mat2str(hamming_code));
end