function [param, param_mat, lambda_vec] = est_mag_bias_scale_and_rotation_RLS_only_mag(mag, lambda_min, p0)

N = size(mag, 1);
param_mat = zeros(N, 10);
lambda_vec = zeros(N, 1);

theta = zeros(9, 1);
U = eye(9);
D = p0 * eye(9);
lambda = lambda_min;

alpha = zeros(9, 1);
k0 = zeros(9, 1);

for i = 1:N

    x = mag(i,1);
    y = mag(i,2);
    z = mag(i,3);

    phi = [x.*x, y.*y, z.*z, x.*y, x.*z, y.*z, x, y, z];

    y_hat = phi * theta;
    y     = 1;
    e = y - y_hat;

    % U D U^T
    f = U.' * phi.';
    v = D * f;
    
    % j = 1
    alpha(1) = lambda + v(1) * f(1);
    D(1,1) = D(1,1) / alpha(1);
    k = k0;
    k(1) = v(1);
    % j = 2:9
    Upast = U; % to calculate k
    for j = 2:9
        alpha(j) = alpha(j-1) + v(j) * f(j);
        D(j,j) = D(j,j) * alpha(j-1) / ( alpha(j) * lambda );
        U(:,j) = U(:,j) - ( f(j) / alpha(j-1) ) * k;
        k = k + v(j) * Upast(:,j);
    end
        
    % parameter-update
    gamma = k ./ alpha(9);
    theta = theta + gamma * e;

    % variable lambda
    zn = lambda / ( lambda + phi * U * v );
    lambda = lambda_min + ( 1 - lambda_min ) * zn^2;
    
    param_mat(i,:) = [theta.', -1];
    lambda_vec(i,:) = lambda;

end

param = [theta.', -1];

end

