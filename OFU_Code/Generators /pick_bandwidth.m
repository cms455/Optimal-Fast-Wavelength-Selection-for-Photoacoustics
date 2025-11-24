%% ========================================================================
% Helper Function: pick_bandwidth
% Crops a spectral curve to a wavelength window and interpolates smoothly.
% ========================================================================
function selected_curve_interp = pick_bandwidth(curve, min_w, max_w, num_points)
% curve is 2×N: row 1 = wavelength, row 2 = absorption

    % Find exact indices for min/max wavelengths (assumes data sampled exactly)
    start_idx = find(curve(:,1) == min_w);
    end_idx   = find(curve(:,1) == max_w);

    % Original wavelength indices in that range
    wavelengths = start_idx:end_idx;

    % Dense wavelength sampling for interpolation
    dense_wavelengths = linspace(start_idx, end_idx, num_points);

    % Extract absorption values in the chosen band
    selected_curve = curve(start_idx:end_idx, 2);

    % Interpolate absorption onto dense grid (smooth spline interpolation)
    selected_curve_interp = interp1(wavelengths, ...
                                    selected_curve, ...
                                    dense_wavelengths, ...
                                    'Spline');

    % Plotting code (disabled)
    %{
    figure; 
    hold on;
    plot(wavelengths, selected_curve, 'ro');          % original samples
    plot(dense_wavelengths, selected_curve_interp,'b'); % interpolated
    %}
end
