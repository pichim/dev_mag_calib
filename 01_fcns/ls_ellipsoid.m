function eansa = ls_ellipsoid(mag)

    % Change xx from vector of length N to Nx1 matrix so we can use hstack
    x = mag(:,1);
    y = mag(:,2);
    z = mag(:,3);

    % Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz = 1
    J = [x.*x, y.*y, z.*z, x.*y, x.*z, y.*z, x, y, z];
    K = ones(size(x)); % Column of ones

    % Transpose J
    JT = J.';
    % Calculate JTJ
    JTJ = JT * J;
    % Calculate inverse of JTJ
    InvJTJ = inv(JTJ);
    % Calculate ABC
    ABC = InvJTJ * (JT * K);

    % Rearrange, move the 1 to the other side
    % Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz - 1 = 0
    % or
    % Ax^2 + By^2 + Cz^2 +  Dxy +  Exz +  Fyz +  Gx +  Hy +  Iz + J = 0
    % where J = -1
    eansa = [ABC; -1];
end