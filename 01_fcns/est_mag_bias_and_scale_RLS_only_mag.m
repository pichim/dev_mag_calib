function [theta, theta_mat, lambda_vec] = est_mag_bias_and_scale_RLS_only_mag(mag, lambda_min, p0)

N = size(mag, 1);
theta_mat = zeros(N, 6);
lambda_vec = zeros(N, 1);

theta = zeros(6, 1);
U = eye(6);
D = p0 * eye(6);
lambda = lambda_min;

alpha = zeros(6, 1);
k0 = zeros(6, 1);

for i = 1:N

    x = mag(i,1);
    y = mag(i,2);
    z = mag(i,3);

    phi = [x.*x, y.*y, z.*z, x, y, z];

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
    % j = 2:6
    Upast = U; % to calculate k
    for j = 2:6
        alpha(j) = alpha(j-1) + v(j) * f(j);
        D(j,j) = D(j,j) * alpha(j-1) / ( alpha(j) * lambda );
        U(:,j) = U(:,j) - ( f(j) / alpha(j-1) ) * k;
        k = k + v(j) * Upast(:,j);
    end
        
    % parameter-update
    gamma = k ./ alpha(6);
    theta = theta + gamma * e;

    % variable lambda
    zn = lambda / ( lambda + phi * U * v );
    lambda = lambda_min + ( 1 - lambda_min ) * zn^2;
    
    theta_mat(i,:) = theta.';
    lambda_vec(i,:) = lambda;

end

end

