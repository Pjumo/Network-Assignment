data = 0b11010011101100u32;
data_length = 14;
divisor = 0b1011u32;
divisor_degree = 3;

fprintf('     Real Data : %s\n', dec2bin(data));
fprintf('       Divisor : %s\n\n', dec2bin(divisor));

remainder = bitshift(data, divisor_degree);
divisor = bitshift(divisor, data_length-1);
fprintf(' 0th Remainder : %s\n', dec2bin(remainder));

for i = 1:data_length
    if bitget(remainder, data_length+divisor_degree)
        remainder = bitxor(remainder, divisor);
    end
    remainder = bitshift(remainder, 1);

    fprintf('%2dth Remainder : ', i);
    for j = 1:(data_length + divisor_degree - i - ...
            strlength(dec2bin(bitshift(remainder, -i))))
        fprintf('0');
    end
    fprintf('%s\n', dec2bin(bitshift(remainder, -i)));
end

CRC_value = bitshift(remainder,-data_length);
fprintf('\n      CRC code : %s\n', dec2bin(CRC_value));

send_data = bitshift(data, divisor_degree) + CRC_value;
fprintf('     Send Data : %s\n\n\n', dec2bin(send_data));


recieved_data = 0b10110000101100000u32;
fprintf(' Recieved Data : %s\n\n', dec2bin(recieved_data));

remainder = recieved_data;

for i = 1:data_length
    if bitget(remainder,data_length+divisor_degree)
        remainder = bitxor(remainder,divisor);
    end
    remainder = bitshift(remainder,1);

    fprintf('%2dth Remainder : ', i);
    for j = 1:(data_length + divisor_degree - i - ...
            strlength(dec2bin(bitshift(remainder, -i))))
        fprintf('0');
    end
    fprintf('%s\n', dec2bin(bitshift(remainder, -i)));
end
if remainder == 0
    fprintf('\n---Error Free---\n');
else
    fprintf('\n---Error Detected---\n');
end




