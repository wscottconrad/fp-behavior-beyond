function basis = createRaisedCosineBasis(n_lags, n_basis)
% creates raised cosine filters for  glm
% n_lags (number of frames dependent signal lags from ind signal
% n_basis is number of filters/basis?
% scott conrad 2025

% Linearly spaced centers

    c = linspace(0, n_lags - 1, n_basis);
    w = (c(2) - c(1));  % Width of each cosine

    basis = zeros(n_lags, n_basis);
    for i = 1:n_basis
        for j = 1:n_lags
            z = (j - 1 - c(i)) * pi / w;
            if abs(z) <= pi
                basis(j, i) = 0.5 * (cos(z) + 1);
            else
                basis(j, i) = 0;
            end
        end
    end
end