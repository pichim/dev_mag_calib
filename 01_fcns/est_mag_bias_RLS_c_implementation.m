function [b, b_mat, lambda_vec] = est_mag_bias_RLS_c_implementation(mag, dmag, gyro, lambda_min, p0)

N = size(mag, 1);
b_mat = zeros(N, 3);
lambda_vec = zeros(N, 1);

mBE = [];
mBE = magBiasEstimatorInit(mBE, lambda_min, p0);

for i = 1:N

    e = [dmag(i,1) + gyro(i,3) * ( mBE.b(2) - mag(i,2) ) - gyro(i,2) * ( mBE.b(3) - mag(i,3) ); ...
         dmag(i,2) - gyro(i,3) * ( mBE.b(1) - mag(i,1) ) + gyro(i,1) * ( mBE.b(3) - mag(i,3) ); ...
         dmag(i,3) + gyro(i,2) * ( mBE.b(1) - mag(i,1) ) - gyro(i,1) * ( mBE.b(2) - mag(i,2) )];

    zn = zeros(3,1);

    % iteration 1 : k = 1; i = 2; j = 3; sign =  1.0;
    [mBE, zn] = magBiasEstimatorSolveRecursively(mBE, zn, e, gyro(i,:), 1, 2, 3,  1.0);
    % iteration 2 : k = 2; i = 1; j = 3; sign = -1.0;
    [mBE, zn] = magBiasEstimatorSolveRecursively(mBE, zn, e, gyro(i,:), 2, 1, 3, -1.0);
    % iteration 3 : k = 3; i = 1; j = 2; sign =  1.0;
    [mBE, zn] = magBiasEstimatorSolveRecursively(mBE, zn, e, gyro(i,:), 3, 1, 2,  1.0);

    for l = 1:3
        zn(l) = zn(l) * mBE.lambda;
        for m = 1:3
            mBE.P(l,m) = mBE.P(l,m) / mBE.lambda;
        end
    end
    
    mBE.lambda = mBE.lambda_min - ( mBE.lambda_min - 1.0 ) * ( zn(1) * zn(1) + zn(2) * zn(2) + zn(3) * zn(3) ) / 3.0;

    b_mat(i,:) = mBE.b.';
    lambda_vec(i,:) = mBE.lambda;

end

    b = mBE.b;

end

function mBE = magBiasEstimatorInit(mBE, lambda_min, p0)

    mBE.lambda_min = lambda_min;
    mBE.p0 = p0;

    mBE.lambda     = lambda_min;
    for i = 1:3
        mBE.b(i) = 0.0;
        for j = 1:3
            if (i == j)
                mBE.P(i, j) = p0;
            else
                mBE.P(i, j) = 0.0;
            end
        end
    end
end

function [mBE, zn] =  magBiasEstimatorSolveRecursively(mBE, zn, e, gyro, k, i, j, sign)

    dP = sign * [mBE.P(1,i) * gyro(j) - mBE.P(1,j) * gyro(i); ...
                 mBE.P(2,i) * gyro(j) - mBE.P(2,j) * gyro(i); ...
                 mBE.P(3,i) * gyro(j) - mBE.P(3,j) * gyro(i)];

    zn(k) = 1.0 / ( mBE.lambda + sign * ( dP(i) * gyro(j) - dP(j) * gyro(i) ) );

    g(1) = zn(k) * dP(1);
    g(2) = zn(k) * dP(2);
    g(3) = zn(k) * dP(3);

    mBE.b(1) = mBE.b(1) - e(k) * g(1);
    mBE.b(2) = mBE.b(2) - e(k) * g(2);
    mBE.b(3) = mBE.b(3) - e(k) * g(3);
    
    mBE.P(1,1) = mBE.P(1,1) - g(1) * dP(1);
    mBE.P(2,1) = mBE.P(2,1) - g(2) * dP(1);
    mBE.P(2,2) = mBE.P(2,2) - g(2) * dP(2);
    mBE.P(3,1) = mBE.P(3,1) - g(3) * dP(1);
    mBE.P(3,2) = mBE.P(3,2) - g(3) * dP(2);
    mBE.P(3,3) = mBE.P(3,3) - g(3) * dP(3);
    mBE.P(1,2) = mBE.P(2,1);
    mBE.P(1,3) = mBE.P(3,1);
    mBE.P(2,3) = mBE.P(3,2);
    
end
