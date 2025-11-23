function A = build_absorption_matrix_v2(min_w, max_w, species_bool, num_points)
% BUILD_ABSORPTION_MATRIX_V2
%   Constructs an absorption matrix A by selecting a subset of species
%   (Hb, HbO, ICG, Lipid, Methylene Blue), cropping each spectrum to a
%   wavelength window [min_w, max_w], and interpolating to a fixed number
%   of points.
%
% INPUTS:
%   min_w, max_w  : wavelength bounds (must match entries in input spectra)
%   species_bool  : logical vector selecting which species to include
%                   Format: [Hb, HbO, Lipid, ICG, MethyleneBlue]
%   num_points     : number of interpolated wavelength samples
%
% OUTPUT:
%   A : num_species × num_points absorption matrix
%       each row = selected + interpolated absorption spectrum

%% -----------------------------------------------------------
% Load combined absorption spectra (HbO, Hb, ICG, Lipid)
% ------------------------------------------------------------
load_A = load('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/combined_spectra_v2.mat');
load_A = load_A.combined_spectra;

% First column is wavelength; remaining columns are absorption spectra
wavelengths = load_A(:,1);

% Transpose so each species becomes a separate row:
% Row 1 = wavelengths, rows 2–5 = species spectra
full_A = load_A(:,:)';   

% Extract individual species spectra.
% Each is 2×N: row 1 = wavelengths, row 2 = absorption values.
HbO_data   = full_A([1,2], :);   % HbO spectrum
Hb_data    = full_A([1,3], :);   % Hb spectrum
ICG_data   = full_A([1,4], :);   % ICG spectrum
lipid_data = full_A([1,5], :);   % Lipid spectrum

%% -----------------------------------------------------------
% Load Methylene Blue separately and rescale it
% ------------------------------------------------------------
load_methyl_blue = readmatrix('/Users/calvinsmith/Bouma_lab/Analytical_Spectral_Unmixing/Spectrum_Data/MethyleneBlue_spectrum.csv');

% Skip header rows
methyl_blue = load_methyl_blue(3:end,:);

% Scale amplitude to match magnitude of other species (empirical step)
methyl_blue(:,2) = 0.1 * methyl_blue(:,2);

% Bundle all species together in standard order
species_data = {Hb_data, HbO_data, lipid_data, ICG_data, methyl_blue};

% Count how many species are selected
num_species = sum(species_bool);

%% -----------------------------------------------------------
% Initialize output absorption matrix
% ------------------------------------------------------------
A = zeros(num_species, num_points);

% species_bool selects which species to include:
% Order is: [Hb, HbO, Lipid, ICG, MethyleneBlue]

count = 1;
for n = 1:length(species_bool)

    % Only include species where species_bool(n) == 1
    if species_bool(n)

        % Extract raw (wavelength, absorption) curve
        curve = species_data{n};

        % Crop to [min_w, max_w] and interpolate to num_points uniformly
        selected_curve_interp = pick_bandwidth(curve, min_w, max_w, num_points);

        % Store into the absorption matrix
        A(count, :) = selected_curve_interp;

        count = count + 1;
    end
end

end  % end main function


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
