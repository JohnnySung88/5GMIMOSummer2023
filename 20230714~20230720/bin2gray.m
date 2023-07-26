function gray_dec = bin2gray(data_dec, method, QAM)

if strcmp(method, 'qam') && QAM == 4
    gray_dec = zeros(size(data_dec));
    for idx = 1:numel(data_dec)
        switch data_dec(idx)
            case 0
                gray_dec(idx) = bin2dec('00');
            case 1
                gray_dec(idx) = bin2dec('01');
            case 2
                gray_dec(idx) = bin2dec('11');
            case 3
                gray_dec(idx) = bin2dec('10');
            otherwise
                error('Input data_dec is out of range 0-3');
        end
    end
else
    error('This function currently only supports QAM4');
end

end
