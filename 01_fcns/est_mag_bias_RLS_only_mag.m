function [b, b_mat, lambda_vec] = est_mag_bias_RLS_only_mag(mag, lambda_min, p0, scale_mag)

mag = mag / scale_mag;

N = size(mag, 1);
b_mat = zeros(N, 3);
lambda_vec = zeros(N, 1);

theta = zeros(4, 1);
U = eye(4);
D = p0 * ones(4,1);
lambda = lambda_min;

for i = 1:N
    phi = [sum(mag(i,:).^2, 2), mag(i,:)].';

    y_hat = phi.' * theta;
    y     = 1;
    e = y - y_hat;

    % U D U^T
    f = U.' * phi;
    v = D .* f;
    
    % j = 1
    alpha = zeros(4, 1);
    alpha(1) = lambda + v(1) * f(1);
    D(1) = D(1) / alpha(1);
    k = zeros(4, 1);
    k(1) = v(1);
    % j = 2:4
    for j = 2:4
        alpha(j) = alpha(j-1) + v(j) * f(j);
        D(j) = D(j) * alpha(j-1) / ( alpha(j) * lambda );
        for l = 1:j-1
            dU = - (f(j) / alpha(j-1) ) * k(l);
            k(l) = k(l) +  v(j) * U(l,j);
            U(l,j) = dU + U(l,j);
        end
        k(j) = k(j) + v(j);
    end
        
    % parameter-update
    gamma = k ./ alpha(4);
    theta = theta + gamma * e;

    % variable lambda
    zn = lambda / ( lambda + phi.' * U * v );
    lambda = lambda_min + ( 1 - lambda_min ) * zn^2;
    
    b_mat(i,:) = -0.5 * theta(2:4).' ./ theta(1);
    lambda_vec(i,:) = lambda;

end

b = -0.5 * theta(2:4) ./ theta(1);

b_mat = b_mat * scale_mag;
b     = b     * scale_mag;

end