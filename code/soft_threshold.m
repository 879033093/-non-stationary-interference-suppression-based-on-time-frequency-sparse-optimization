function x_soft = soft_threshold(x, threshold)
    x_soft = sign(x) .* max(abs(x) - threshold, 0);
end
