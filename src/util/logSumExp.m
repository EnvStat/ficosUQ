function [y] = logSumExp(x)
    %LOGSUMEXP Numerically stable way to calculate y = log(sum(exp(x)));
    %   y = logSumExp(x)
    c = max(x);
    y = c + log(sum(exp(x-c)));
end

