function [b, scale, b_mat, scale_mat, lambda_vec] = est_mag_bias_and_scale_RLS_only_mag_c_implementation(mag, lambda_min, p0, scale_mag)

N = size(mag, 1);
b_mat = zeros(N, 3);
scale_mat = zeros(N, 3);
lambda_vec = zeros(N, 1);

cBE = [];
cBE = magBiasEstimatorInit(cBE, lambda_min, p0);

sq = @(x) x*x;

for x = 1:N

    mag_scaled = mag(x,:) / scale_mag;

    % phi = [x.*x, y.*y, z.*z, x, y, z];
    phi = [sq(mag_scaled(1)), sq(mag_scaled(2)), sq(mag_scaled(3)), ...
           mag_scaled(1), mag_scaled(2), mag_scaled(3)].';

    % e = 1.0 - phi.' * cBE.theta;
    e = 1.0;
    for i = 1:6
        e =  e - phi(i) * cBE.theta(i);
    end

    % U D U^T
    % f = cBE.U.' * phi;
    % v = cBE.D .* f;
    f = rand(6,1);
    v = rand(6,1);
    for i = 1:6
        f(i) = 0.0;
        for j = 1:i
            f(i) = f(i) + cBE.U(j,i) * phi(j);
        end
        v(i) = cBE.D(i) * f(i);
    end
    
    % first iteration
    % alpha = zeros(6, 1);
    % k = zeros(6, 1);
    % alpha(1) = cBE.lambda + v(1) * f(1);
    % cBE.D(1) = cBE.D(1) / alpha(1);
    % k(1) = v(1);
    alpha = rand(6,1);
    k = zeros(6,1);
    alpha(1) = cBE.lambda + v(1) * f(1);
    cBE.D(1) = cBE.D(1) / alpha(1);
    k(1) = v(1);


    % rest of the iterations
    % for j = 2:6
    %     alpha(j) = alpha(j-1) + v(j) * f(j);
    %     cBE.D(j) = cBE.D(j) * alpha(j-1) / ( alpha(j) * cBE.lambda );
    %     for l = 1:j-1
    %         dU = - (f(j) / alpha(j-1) ) * k(l);
    %         k(l) = k(l) +  v(j) * cBE.U(l,j);
    %         cBE.U(l,j) = dU + cBE.U(l,j);
    %     end
    %     k(j) = k(j) + v(j);
    % end
    for i = 2:6
        alpha(i) = alpha(i - 1) + v(i) * f(i);
        cBE.D(i) = cBE.D(i) * alpha(i - 1) / (alpha(i) * cBE.lambda);
        for j = 1:i-1
            dU = -(f(i) / alpha(i - 1)) * k(j);
            k(j) = k(j) + v(i) * cBE.U(j,i);
            cBE.U(j,i) = cBE.U(j,i) + dU;
        end
        k(i) = k(i) + v(i);
    end
        
    % parameter-update
    % gamma = k ./ alpha(4);
    % cBE.theta = cBE.theta + gamma * e;
    for i = 1:6
        cBE.theta(i) = cBE.theta(i) + (k(i) / alpha(6)) * e;
    end

    % bias update
    % cBE.b = (-0.5 * cBE.theta(4:6) ./ cBE.theta(1:3)) * scale_mag;
    % cBE.scale = sqrt(cBE.theta(1:3));
    % cBE.scale = cBE.scale ./ mean(cBE.scale);
    scale_sum = 0.0;
    for i = 1:3
        cBE.b(i) = (-0.5 * cBE.theta(i+3) / cBE.theta(i)) * scale_mag;
        cBE.scale(i) = sqrt(cBE.theta(i));
        scale_sum = scale_sum + cBE.scale(i);
    end
    for i = 1:3
        cBE.scale(i) = cBE.scale(i) * 3.0 / scale_sum;
    end

    % compute zn
    % zn = cBE.lambda / ( cBE.lambda + phi.' * (cBE.U * v) );
    U_v = rand(1);
    phiTrans_U_v = 0.0;
    for i = 1:6
        U_v = 0.0;
        for j = i:6
            U_v = U_v + cBE.U(i,j) * v(j);
        end
        phiTrans_U_v = phiTrans_U_v + phi(i) * U_v;
    end
    zn = cBE.lambda / (cBE.lambda + phiTrans_U_v);

    % update lambda
    % cBE.lambda = cBE.lambda_min + ( 1 - cBE.lambda_min ) * zn^2;
    cBE.lambda = cBE.lambda_min + (1.0 - cBE.lambda_min) * sq(zn);
    
    b_mat(x,:) = cBE.b.';
    scale_mat(x,:) = cBE.scale.';
    lambda_vec(x,:) = cBE.lambda;

end

b = cBE.b;
scale = cBE.scale;

end


function cBE = magBiasEstimatorInit(cBE, lambda_min, p0)

    % typedef struct compassBiasEstimator_s {
    %     float lambda_min, lambda;
    %     float b[3];
    %     float scale[3];
    %     float theta[6];
    %     float U[6][6];
    %     float D[6];
    % } compassBiasEstimator_t;

    % memset(cBE, 0, sizeof(*cBE)); // zero contained IEEE754 floats
    cBE.lambda_min = 0.0;
    cBE.lambda = 0.0;
    cBE.b = zeros(3,1);
    cBE.scale = zeros(3,1);
    cBE.theta = zeros(6,1); %cBE.theta(1:3) = 1/9;
    cBE.U = zeros(6,6);
    cBE.D = zeros(6,1);
    % create identity matrix
    for i = 1:6
        cBE.U(i,i) = 1.0;
    end

    % compassBiasEstimatorUpdate(cBE, lambda_min, p0); 
    cBE.lambda_min = lambda_min;
    % update diagonal entries for faster convergence
    for i = 1:3
        cBE.D(i) = p0;
    end
    for i = 4:6
        cBE.D(i) = p0;
    end

    cBE.lambda = lambda_min;
end
