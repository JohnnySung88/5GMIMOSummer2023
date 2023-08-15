function X_Interpolation = InterpolationDetector(Y, H_frame, Tx, W)
    [N, S, ~, ~] = size(H_frame); % N - subcarriers, S - symbols
    X_Interpolation = zeros(N, S, Tx);

    for tx = 1:Tx
        for s = 1:S
            known_channels = H_frame(2:2:end, 3:14:s, :, tx); % 使用LS估計的通道作為已知點
            unknown_channels = H_frame(1:2:end, 3:14:s, :, tx);
            
            % 進行線性插值
            interpolated_channels = interp1(2:2:N, known_channels, 1:N, 'linear', 'extrap');
            
            % 使用插值的通道進行LMMSE估計
            H_interpolated = reshape(interpolated_channels, [N, 1, 1, 1]);
            H_interpolated = repmat(H_interpolated, [1, 1, size(H_frame, 3), 1]);
            X_LMMSE = LMMSE(Y, H_interpolated, tx);
            
            % 將結果存儲到輸出變量中
            X_Interpolation(:, s, tx) = X_LMMSE(:, s, tx) * W(s);
        end
    end
end