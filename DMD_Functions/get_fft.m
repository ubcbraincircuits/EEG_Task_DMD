function [S, f] = get_fft(X, fs, N)
    S = abs(fft(X, N));
    S = S(1:N/2,:);
    if rem(N,2) == 0
        f = (0:N/2-1)*fs/(N);
    else
        f = (0:(N-1)/2)*fs/(N);
    end
end
