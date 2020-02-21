function [estim_params] = fit_sim(mt, dist, init_params)
fcn = @(params)loss_fcn(mt, run_sim(dist, params));
estim_params = fminsearch(fcn, init_params);
end

function [mt_hat] = run_sim(dist, params)
a=params(1); b=params(2);
mt_hat = a + b * (log(2 * (dist / 0.2))); % fitts' law 
end

function [loss] = loss_fcn(y, y_hat)
n = numel(y); loss = zeros(n, 1);
for i = 1:n; loss(i, 1) = (y(i) * log(y_hat(i)) + (1-y(i)) * log(1-y_hat(i))); end
loss = -(nansum(loss));
end