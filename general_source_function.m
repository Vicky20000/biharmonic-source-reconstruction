function S = general_source_function(y)
    % Evaluate the four overlapping polynomial bumps at N-by-2 points.

    x = y(:,1);
    z = y(:,2);

    S = zeros(size(y, 1), 1);

    % Each row gives the center of one bump.
    centers = [
        -0.65,  0.45;
         0.35,  0.45;
        -0.35, -0.45;
         0.65, -0.40
    ];

    % Support radii.
    radii = [1.05, 0.98, 0.95, 0.90];

    % Amplitudes.
    amps = [1.00, 0.80, 0.65, 0.50];

    % The fourth power gives a C^3 compactly supported bump.
    m = 4;

    % Add the contributions where the supports overlap.
    for j = 1:length(amps)
        x0 = centers(j,1);
        z0 = centers(j,2);
        R  = radii(j);

        r2 = ((x - x0).^2 + (z - z0).^2) / R^2;

        in_region = r2 < 1;

        S(in_region) = S(in_region) + amps(j) .* (1 - r2(in_region)).^m;
    end
end
