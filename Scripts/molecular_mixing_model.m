function [frac, CNHO, rmse, r_squared, new_mat] = molecular_mixing_model(NMR_data)
% MOLECULAR_MIXING_MODEL Calculates fractional compositions of different five
% molecular classes using the molecular mixing model.
% 
% [frac, CNHO, rmse, r_squared, new_mat] = molecular_mixing_model(NMR_data)
% calculates the fractional compositions of different components based on 
% NMR data using the molecular mixing model.
%
% Input Arguments:
% - NMR_data: A structure containing NMR data.
%
% Output Arguments:
% - frac: Fractional compositions of different components.
% - CNHO: Table containing C, N, H, O compositions and molecular formulas.
% - rmse: Root Mean Square Error of the simulted NMR and obs NMR data
% - r_squared: R-squared value of the simulted NMR and obs NMR data
% - new_mat: Updated mixing matrix.
%
% Input Data Structure (NMR_data):
% NMR_data should be a structure with the following fields:
% - source: Source of the data.
% - C_N: Carbon to Nitrogen ratio.
% - C_conc: Carbon concentration.
% - N_conc: Nitrogen concentration.
% - AALKYL0_45Ppm to GCARBOX165_210Ppm: NMR data for different chemical shifts.
% - a to f: Availability spectra.
% - system: Type of system (terrestrial or aquatic).
% - Csource: Carbon source (e.g., plantC, MicrobialC, pure, soil, etc.).
%
% Example:
% [frac, CNHO, rmse, r_squared, new_mat] = molecular_mixing_model(NMR_data);
% see NMR Data Analysis Script for example
% Author: Arjun Chakrawal
% Affiliation: Stockhom Univeristy
% Contact: arjun.chakrawal@natgeo.su.se
%

SRC = NMR_data.source{1};
% Calculating NCobs based on C_N ratio or C_conc and N_conc
if ~isnan(NMR_data.C_N)
    NCobs = (12 / 14) / NMR_data.C_N; %Molar
else
    NCobs = (12 / 14) / (NMR_data.C_conc/NMR_data.N_conc); %Molar
end

% Extracting NMR data and availability spectra
NMR = [NMR_data.AALKYL0_45Ppm, NMR_data.BMETHOX45_60Ppm, NMR_data.CO_ALKYL60_95Ppm, ...
    NMR_data.DDI_O_ALK95_110ppm, NMR_data.EAROM110_145Ppm, NMR_data.FPHEN145_165Ppm, ...
    NMR_data.GCARBOX165_210Ppm]';
NMR = NMR(~isnan(NMR));
avail_spectra = [NMR_data.a, NMR_data.b, NMR_data.c, NMR_data.d, ...
    NMR_data.e, NMR_data.f, NMR_data.g];
avail_spectra = avail_spectra(~cellfun(@isempty, avail_spectra));

% Mixing matrices and constraints based on system and carbon source
% Assignment below is from Table 1 of Nelson and Baldock
% Carbohydrate	Protein	Lignin	Lipid	Carbonyl	Char
Mixing_terrestrial = [; ...
    0, 39.6, 10.5, 75.6, 0, 0; ... %0-45
    4.3, 21.9, 13.8, 4.5, 0, 0; ... % 45-60
    79.0, 2.1, 12.5, 9, 0, 0; ... %60-95
    15.7, 0, 8.6, 0, 0, 4.3; ... %95-110
    1, 7.5, 30.6, 3.6, 0, 73.9; ... %110-145
    0, 2.5, 19.5, 0.7, 0, 16.1; ... %145-165
    0, 26.4, 4.6, 6.6, 100, 5.6, ... % 165-215
    ] ./ 100;

Mixing_aquatic = [; ...
    0, 36.6, 10.5, 94.4, 0, 1.3; ... %0-45
    0, 24.7, 13.8, 0, 0, 1.2; ... % 45-60
    83.3, 2.9, 12.5, 0, 0, 1.3; ... %60-95
    16.7, 0, 8.6, 0, 0, 6.3; ... %95-110
    0, 4.5, 30.6, 0, 0, 64.9; ... %110-145
    0, 1.0, 19.5, 0, 0, 17.5; ... %145-165
    0, 30.4, 4.6, 5.6, 100, 7.7, ... % 165-215
    ] ./ 100;

CNHO_terrestrial = [1, 1, 1, 1, 1, 1; ...
    0, 0.27, 0, 0, 0, 0; ...
    1.67, 1.10, 1.24, 1.94, 1, 0.45; ...
    0.83, 0.16, 0.43, 0.24, 2, 0.41];

CNHO_aquatic = [1, 1, 1, 1, 1, 1; ...
    0, 0.27, 0, 0, 0, 0; ...
    1.67, 1.17, 1.24, 1.89, 1, 0.45; ...
    0.83, 0.12, 0.43, 0.11, 2, 0.41];


% Define Mixing matrices
if (strcmp(NMR_data.system, 'terrestrial'))
    MixingMatrix = Mixing_terrestrial; % or Mixing_aquatic
    switch SRC
        case {'plantC', 'MicrobialC', 'pure'} % no char
            Assignment = MixingMatrix(:, 1:5);
            CNHO_Assignment = CNHO_terrestrial(:, 1:5);
        case {'soil', 'compostANDmanure', 'peat', 'wastewater'}
            Assignment = MixingMatrix;
            CNHO_Assignment = CNHO_terrestrial;
    end
else
    MixingMatrix = Mixing_aquatic; % or Mixing_aquatic
    switch SRC
        case 'plantC' % no char
            Assignment = MixingMatrix(:, 1:5);
            CNHO_Assignment = CNHO_aquatic(:, 1:5);
        case 'MicrobialC'
            Assignment = MixingMatrix(:, [1, 2, 4, 5]);
            CNHO_Assignment = CNHO_aquatic(:, [1, 2, 4, 5]);
        case 'sea_floor'
            Assignment = MixingMatrix;
            CNHO_Assignment = CNHO_aquatic;
    end
end

if (strcmp(NMR_data.Csource, 'biochar'))
    Assignment = MixingMatrix;
    CNHO_Assignment = CNHO_terrestrial;
end


% Calculate new_mat
str = ["a", "b", "c", "d", "e", "f", "g"];
new_mat = zeros(length(avail_spectra), size(Assignment, 2));
idx = 0;
for i = 1:length(avail_spectra)
    temp = avail_spectra{i};
    idx = idx + numel(temp);
    if numel(temp) == 1
        new_mat(i, :) = Assignment(idx, :);
    else
        id = zeros(numel(temp), 1);
        for j = 1:numel(temp)
            id(j) = find(strcmp(temp(j), str));
        end
        new_mat(i, :) = sum(Assignment(id, :));
    end
end


% Define constraints Aeq and beq
prot_fraction=nan;
if isnan(NCobs)
    Aeq = ones(1, size(Assignment, 2)); % sum of fractions = 1
    beq = 1;
else
    prot_fraction = NCobs/0.27;
    new_mat(:,2) = [];
    Aeq = [ones(1, size(new_mat, 2))]; % sum of fractions = 1, N:C ratio
    beq = 1-prot_fraction;
end

% Define optimization options
options = optimoptions('fmincon', 'Display', 'none', 'Algorithm', 'interior-point', ...
    'EnableFeasibilityMode', true, 'SubproblemAlgorithm', 'cg');

% Define the objective function
fun = @(x) sqrt(mean((NMR - new_mat * x').^2));

% Define lower and upper bounds for x
lb = zeros(1, size(new_mat, 2));
ub = ones(1, size(new_mat, 2));

% Initial guess for x
guess_xopt = ones(1, size(new_mat, 2)) * 0.01;

% Create the optimization problem
problem = createOptimProblem('fmincon', 'objective', fun, 'x0', guess_xopt, ...
    'lb', lb, 'ub', ub, 'Aeq', Aeq, 'beq', beq, 'options', options);

% Solve the optimization problem
x = fmincon(problem);
if ~isnan(prot_fraction)
    frac = [x(1), prot_fraction,x(2:end)];
else
    frac =x;
end

% Calculate CNHO compositions and molecular formulas
CNHO = array2table((CNHO_Assignment * frac')', VariableNames = {'C', 'N', 'H', 'O'});
CNHO.molecularFormula = strcat('C', num2str(CNHO.C, 2), 'H', num2str(CNHO.H, 2), ...
    'N', num2str(CNHO.N, 2), 'O', num2str(CNHO.O, 2));
CNHO.Cox = 4 - (4 * CNHO.C + CNHO.H - 2 * CNHO.O - 3 * CNHO.N);

% Calculate RMSE and R-squared values
est_nmr =  new_mat * x';
rmse = sqrt(mean((NMR - est_nmr).^2));
tss = sum((NMR - mean(NMR)).^2);
rss = sum((NMR - est_nmr).^2);
r_squared = 1 - (rss / tss);
end
