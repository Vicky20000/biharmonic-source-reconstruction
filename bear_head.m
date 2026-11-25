function z = bear_head(points)
% Evaluate the piecewise-constant source at N-by-2 Cartesian points.
% Features assigned later in the loop overwrite the head value.

    n = size(points, 1);
    z = zeros(n, 1);

    for i = 1:n
        x = points(i, 1);
        y = points(i, 2);

        % Head.
        head_radius = 1.2;
        head_center = [0, 0];
        head_dist = sqrt((x - head_center(1)).^2 + (y - head_center(2)).^2);
        head_mask = head_dist <= head_radius;
        if head_mask
            z(i) = 0.3;
        end

        % Outer ears.
        ear_radius = 0.5;
        left_ear_center = [-0.8, 0.8];
        right_ear_center = [0.8, 0.8];

        left_ear_dist = sqrt((x - left_ear_center(1)).^2 + (y - left_ear_center(2)).^2);
        right_ear_dist = sqrt((x - right_ear_center(1)).^2 + (y - right_ear_center(2)).^2);

        left_ear_mask = left_ear_dist <= ear_radius;
        right_ear_mask = right_ear_dist <= ear_radius;

        if left_ear_mask || right_ear_mask
            z(i) = 0.4;
        end

        % Inner ears.
        inner_ear_radius = 0.3;
        left_inner_ear_mask = left_ear_dist <= inner_ear_radius;
        right_inner_ear_mask = right_ear_dist <= inner_ear_radius;

        if left_inner_ear_mask || right_inner_ear_mask
            z(i) = 0.2;
        end

        % Eyes.
        eye_x_radius = 0.15;
        eye_y_radius = 0.2;

        left_eye_mask = ((x + 0.4).^2/eye_x_radius^2 + (y - 0.2).^2/eye_y_radius^2) <= 1;
        right_eye_mask = ((x - 0.4).^2/eye_x_radius^2 + (y - 0.2).^2/eye_y_radius^2) <= 1;

        if left_eye_mask || right_eye_mask
            z(i) = 0.1;
        end

        % Pupils.
        pupil_radius = 0.05;
        left_pupil_mask = sqrt((x + 0.4).^2 + (y - 0.2).^2) <= pupil_radius;
        right_pupil_mask = sqrt((x - 0.4).^2 + (y - 0.2).^2) <= pupil_radius;

        if left_pupil_mask || right_pupil_mask
            z(i) = 0.9;
        end

        % Mouth interior (trapezoid).
        nose_mask = (y <= -0.3) & (y >= -0.6) & (abs(x) <= 0.3 - 0.5*(y + 0.3));
        if nose_mask
            z(i) = 0.7;
        end

        % Lips (annular sector).
        [theta, rho] = cart2pol(x, y+0.15);
        mouth_mask = (rho <= 0.4) & (rho >= 0.3) & (theta >= -3*pi/4) & (theta <= -1*pi/4);
        if mouth_mask
            z(i) = 0.1;
        end

        % Cheeks.
        cheek_radius = 0.25;
        left_cheek_mask = sqrt((x + 0.6).^2 + (y + 0.2).^2) <= cheek_radius;
        right_cheek_mask = sqrt((x - 0.6).^2 + (y + 0.2).^2) <= cheek_radius;

        if left_cheek_mask || right_cheek_mask
            z(i) = 0.6;
        end
    end
end
