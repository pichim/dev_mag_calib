function [b, b_mat, lambda_vec] = est_mag_bias_RLS(mag, dmag, gyro, lambda_min, p0)

N = size(mag, 1);
b_mat = zeros(N, 3);
lambda_vec = zeros(N, 1);

b = zeros(3, 1);
P = p0 * eye(3);
lambda = lambda_min;

zn    = zeros(3, 1);
Gamma = zeros(3, 3);
for i = 1:N
    Sw = getSkew( gyro(i,:) );

    y_hat = Sw * b;
    y     = dmag(i,:).' + Sw * mag(i,:).';
    e = y - y_hat;
    
    for j = 1:3
        zn(j) = 1 / ( Sw(j,:) * P * Sw(j,:).' + lambda );
        Gamma(:,j) = P * Sw(j,:).' * zn(j);
        b = b + Gamma(:,j) * e(j);
        P = ( eye(3) - Gamma(:,j) * Sw(j,:) ) * P;
    end

    zn(j) = zn(j) * lambda;
    P = P / lambda;

    lambda = lambda_min + (1 - lambda_min) * ( zn.' * zn ) / 3;

%     zn = svd( ( lambda * eye(3) ) / ( Sw * P * Sw.' + lambda * eye(3) ) ).';
%     % zn = sort( eig( ( lambda * eye(3) ) / ( Sw * P * Sw.' + lambda * eye(3) ) ), 'descend' );
%     % - 3 parameters -> 2 singularvalues not equal to 1 and <= 1
%     % - since the matrix is symmetrix svd = sort(eig, 'descend')
%     lambda = lambda_min + (1 - lambda_min) * ( zn(2)^2 + zn(3)^2 ) / 2;
    
    b_mat(i,:) = b.';
    lambda_vec(i,:) = lambda;
end

end

