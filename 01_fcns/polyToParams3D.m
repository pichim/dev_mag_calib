function [center, axes, rotationMatrix] = polyToParams3D(vec, printMe)

    % Convert the polynomial form of the 3D-ellipsoid to parameters
    % Center, axes, and transformation matrix
    % Vec is the vector whose elements are the polynomial coefficients A..J
    % Returns (center, axes, rotation matrix)

    % Algebraic form: X.T * Amat * X --> polynomial form

    if nargin == 1
        printMe = false;
    end

    if printMe
        fprintf('\nPolynomial\n');
        disp(vec);
    end

    Amat = [
        vec(1),     vec(4)/2.0, vec(5)/2.0, vec(7)/2.0;
        vec(4)/2.0, vec(2),     vec(6)/2.0, vec(8)/2.0;
        vec(5)/2.0, vec(6)/2.0, vec(3),     vec(9)/2.0;
        vec(7)/2.0, vec(8)/2.0, vec(9)/2.0, vec(10) 
    ];

    if printMe
        fprintf('\nAlgebraic form of polynomial\n');
        disp(Amat);
    end

    % See B.Bartoni, Preprint SMU-HEP-10-14 Multi-dimensional Ellipsoidal Fitting
    % Equation 20 for the following method for finding the center
    A3 = Amat(1:3, 1:3);
    A3inv = inv(A3);
    ofs = vec(7:9)/2.0;
    center = -A3inv * ofs;
    
    if printMe
        fprintf('\nCenter at:\n');
        disp(center);
    end

    % Center the ellipsoid at the origin
    Tofs = eye(4);
    Tofs(4, 1:3) = center';
    R = Tofs * Amat * Tofs';
    
    if printMe
        fprintf('\nAlgebraic form translated to center\n');
        disp(R);
        fprintf('\n');
    end

    R3 = R(1:3, 1:3);
    R3test = R3 / R3(1, 1);
    if printMe
        fprintf('Normed \n');
        disp(R3test);
        fprintf('\n');
    end
    s1 = -R(4, 4);
    R3S = R3 / s1;
    [ec, el] = eig(R3S);

    recip = 1.0 ./ abs(diag(el));
    axes = sqrt(recip);
    
    if printMe
        fprintf('\nAxes are\n');
        disp(axes);
        fprintf('\n');
    end

    rotationMatrix = inv(ec)'; % Inverse is actually the transpose here
    
    if printMe
        fprintf('\nRotation matrix\n');
        disp(rotationMatrix);
    end
end