function A = generate_spectrum_curve(num_wavelengths, num_species, min_wavelength, max_wavelength, num_peaks)
% GENERATE_SPECTRUM_CURVE
%   Generates a synthetic absorption matrix where each row is a simulated
%   absorption spectrum for a different species. Each spectrum is built
%   using the helper function BUILD_CURVE, which creates smooth curves
%   with a specified number of peaks.
%
% INPUTS:
%   num_wavelengths : number of wavelength samples (columns)
%   num_species     : number of simulated species (rows)
%   min_wavelength  : lower bound of wavelength range
%   max_wavelength  : upper bound of wavelength range
%   num_peaks       : number of Gaussian-like peaks per species
%
% OUTPUT:
%   A : num_species × num_wavelengths matrix
%       each row is a synthetic absorption spectrum

    % If true, build_curve will plot each synthetic spectrum.
    % Disabled by default for speed and cleanliness.
    plot_flag = false;

    % Preallocate the absorption matrix:
    % Each row i will hold the ith species' synthetic spectrum.
    A = zeros(num_species, num_wavelengths);

    % Generate one absorption curve per species
    for i = 1:num_species

        % Build a smooth, random multi-peak absorption profile
        absorption_curve = build_curve( ...
            num_wavelengths, ...
            min_wavelength, ...
            max_wavelength, ...
            num_peaks, ...
            plot_flag);

        % Store the generated spectrum in row i
        A(i, :) = absorption_curve;
    end
end
