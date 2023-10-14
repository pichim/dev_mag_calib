clc, clear all % variables

addpath 01_fcns\
addpath ..\bf_function_libary\
%%

% parameters
do_i_use_matlab        = false;
do_show_MgnCalibration = false;
do_show_RLS_results    = false;
do_compare_RLS_results = false;
do_use_filtered_data   = false;
f_cut = 20; filter_type = 'pt3';


linewidth = 1.2;
set(0, 'defaultAxesColorOrder', get_my_colors);


if ~do_i_use_matlab
    try
        pkg load control
        pkg load signal
    catch exception
        % nothing
    end
end

fs_mag = 200;        % assumed sampling frequency of mag unit
lambda_min = 0.95;   % minimal adaptive forgetting factor, range: [0.90, 0.99]
p0 = 1e2;            % value to initialize P(0) = diag([P0, P0, P0]), typically in range: (1, 1000)

% these parameters do not matter in current c implementation
lambda = 1.00;       % adaptive forgetting factor, range: [0.9, 1.0]
scale_mag = 1.0e0;   % unnescessary
                     % in the range (1000, 2000)
lambda_min_bias_and_scale = 0.97;
lambda_min_bias_scale_and_rotation = 0.98;

% measurements

% quad armed, probs off
% - set blackbox_mode = NORMAL
% file_name = '20231008_apex5_mag_on_tpu_00.bbl.csv';
% T_eval = [4.0, 136.0];

% online calibration using stick commands, blackboxmode
% - set blackbox_mode = ALWAYS
% - result from fc: mag_calibration = 1010,505,549
% file_name = '20231009_apex5_mag_on_tpu_00.bbl.csv';
% T_eval = [25.7226, inf];
% - result from fc: mag_calibration = 1011,503,559
% file_name = '20231009_apex5_mag_on_tpu_01.bbl.csv';
% T_eval = [20.7976, inf];

% online calibration using stick commands
% - set blackbox_mode = ALWAYS
% - result from fc: mag_calibration = 892,398,466
% file_name = '20231014_apex5_mag_on_tpu_00.bbl.csv';
% T_eval = [16.0919 + 0.0*(+ 0.2 - 1.0), inf];

% - result from fc: mag_calibration = 945,417,503
% file_name = '20231014_apex5_mag_on_tpu_01.bbl.csv';
% T_eval = [26.4394, inf];
% - result from fc: mag_calibration = 958,420,506
% file_name = '20231014_apex5_mag_on_tpu_02.bbl.csv';
% T_eval = [11.8071, inf];
% - result from fc: mag_calibration = 957,417,502
file_name = '20231014_apex5_mag_on_tpu_03.bbl.csv';
T_eval = [15.6461, inf];


% extract header information
file_path = ['00_data/', file_name];
[para, Nheader, ind] = extract_header_information(file_path);


% read the data
tic
try
   load([file_path(1:end-8), '.mat'])
catch exception
   % data = readmatrix(file_path, 'NumHeaderLines', Nheader);
   import_data = importdata(file_path, ',', Nheader);
   data = import_data.data;
   save([file_path(1:end-8), '.mat'], "data");
end
[Ndata, Nsig] = size(data) %#ok
toc


% convert and evaluate time
time = (data(:,ind.time) - data(1,ind.time)) * 1.0e-6;
dtime_meas_mus = diff(time) * 1.0e6;


% unscale highResolutionGain
if para.blackbox_high_resolution
    blackbox_high_resolution_scale = 10.0;
    ind_bb_high_res = [ind.gyroADC, ind.gyroUnfilt, ind.rcCommand, ind.setpoint(1:3)];
    data(:, ind_bb_high_res) = 1.0 / blackbox_high_resolution_scale * data(:, ind_bb_high_res);
end


% create different sampling times
Ts      = para.looptime * 1.0e-6;             % gyro
Ts_cntr = para.pid_process_denom * Ts;        % cntrl
Ts_log  = para.frameIntervalPDenom * Ts_cntr; % logging


%% helper functions

draw_line = @() fprintf(' ------------------------------------\n');
draw_matrix = @(M) fprintf('%10.4f, %10.4f, %10.4f\n', M.'); % have to transpose that is is shown corret


%% signal processing

% filter data
if do_use_filtered_data
    [~, Bf, Af] = get_filter(filter_type, f_cut, Ts_log); %#ok
    ind_f = [ind.magADC, ind.gyroADC, ind.accSmooth];
    % data(:,ind_f) = filter(Bf, Af, data(:,ind_f));
    data(:,ind_f) = filtfilt(Bf, Af, data(:,ind_f));
end


% downasmple data
n_ds = (1/Ts_log) / fs_mag; % sample from (1/Ts_log) Hz to fs_mag Hz
data = data(1:n_ds:end,:);
time = time(1:n_ds:end);
Ts = Ts_log * n_ds;


% extract and scale relevant data
mag  = data(:,ind.magADC );
gyro = data(:,ind.gyroADC) * pi / 180;      % rad/sec
acc  = data(:,ind.accSmooth) / 2000 * 9.81; % m/s^2

figure(1)
ax(1) = subplot(221);
plot(time, mag), grid on, xlim([0 time(end)]), ylabel('magADC')
title('mag')
ax(2) = subplot(223);
plot(time, sqrt(sum(mag.^2, 2))), grid on, xlim([0 time(end)]), ylabel('|magADC|'), xlabel('Time (sec)')
% plot(time, data(:, ind.debug(:,8))), grid on, ylabel('lambda'), xlabel('Time (sec)')
ax(3) = subplot(222);
plot(time, gyro), grid on, xlim([0 time(end)]), ylabel('gyroADC')
title('gyro')
ax(4) = subplot(224);
plot(time, sqrt(sum(gyro.^2, 2))), grid on, xlim([0 time(end)]), ylabel('gyroADC'), xlabel('Time (sec)')
linkaxes(ax, 'x'), clear ax


% use only the part of the measurement where copter was lifted from ground
ind_eval = time >= T_eval(1) & time < T_eval(2);
mag  = mag(ind_eval,:);
gyro = gyro(ind_eval,:);
time = time(ind_eval); time = time - time(1);
N = size(gyro, 1);


% calculate mag derivative
dmag = [zeros(1,3); diff(mag) / Ts];


% write text file that can be run in mag_calib dev project: https://github.com/pichim/mag_calib
fileID = fopen(['03_input_for_c_project/', file_name(1:end-8), '.txt'], 'w');
for i = 1:N
    fprintf(fileID, '%0.12f, %0.12f, %0.12f, %0.12f, %0.12f, %0.12f\n', ...
            mag(i,1), mag(i,2), mag(i,3), gyro(i,1), gyro(i,2), gyro(i,3));
end
fclose(fileID);


%% Matlab File-Exchange function MgnCalibration

% carefull, return values are correspond to model A * (mag - b)

if (do_show_MgnCalibration)
    % MgnCalibration full calibration
    C = [[0 0 1]; ...
         [0 1 0]; ...
         [1 0 0]]; %#ok
    [A, b] = MgnCalibration( mag * C.' );
    A = C.' * A * C;
    expgfs0 = 1 / mean( eig( A ) );
    A = A * expgfs0;
    b = C.' * b;
    draw_line()
    fprintf(" MgnCalibration\n")
    draw_matrix(A)
    draw_matrix(b)
end


%% Matlab functions magcal

if (do_i_use_matlab)
    % Algorithm 1: magcal only bias
    [~, b] = magcal(mag, 'eye'); %#ok
    b = b.';
    draw_line()
    fprintf(" magcal only bias\n")
    draw_matrix(b)


    % Algorithm 2: magcal bias and scaling
    [A, b] = magcal(mag, 'diag');
    b = b.';
    draw_line()
    fprintf(" magcal bias and scaling\n")
    draw_matrix(A)
    draw_matrix(b)


    % Algorithm 3: magcal full calibration
    [A, b] = magcal(mag, 'sym');
    b = b.';
    draw_line()
    fprintf(" magcal full calibration\n")
    draw_matrix(A)
    draw_matrix(b)
end


%% Least-Squares Solutions

% https://github.com/pronenewbits/Arduino_AHRS_System
% - almost identical to magcal with option 'eye' (<1% difference)
% mag_scaled = mag ./ 2000; cond([mag_scaled, ones(N, 1)])
b = [eye(3), zeros(3,1)] * ([mag, ones(N, 1)] \ (0.5 * sum(mag.^2, 2)));
draw_line()
fprintf(" alternative LS solution only bias\n")
draw_matrix(b)


% http://www.juddzone.com/ALGORITHMS/least_squares_3D_ellipsoid.html
% mag_scaled = mag ./ 2000; cond([sum(mag_scaled.^2, 2), mag_scaled])
theta = [sum(mag.^2, 2), mag] \ ones(N,1);
b = -0.5 * theta(2:4) ./ theta(1);
draw_line()
fprintf(" Algorithm 1: LS solution only bias\n")
draw_matrix(b)

calib.eye.A = eye(3);
calib.eye.b = b;


% http://www.juddzone.com/ALGORITHMS/least_squares_3D_ellipsoid.html
% - almost identical to magcal with option 'diag' (<1% difference)
% scale_mag_test = 1e2; cond([(mag/scale_mag_test).^2, mag/scale_mag_test])
% mag_scaled = mag ./ 2000; cond([mag_scaled.^2, mag_scaled])
theta = [mag.^2, mag] \ ones(N,1);
b = -0.5 * theta(4:6) ./ theta(1:3);
a = sqrt(theta(1:3));
A = diag(a);
A = A ./ mean( diag(A) );
draw_line()
fprintf(" Algorithm 2: LS solution bias and scaling\n")
draw_matrix(A)
draw_matrix(b)

calib.diag.A = A;
calib.diag.b = b;


% http://www.juddzone.com/ALGORITHMS/least_squares_3D_ellipsoid.html
% - almost identical to magcal with option 'sym' (<1% difference)
[b, axes, R] = polyToParams3D( ls_ellipsoid(mag) );
A = R * diag(1./axes) * R.';
A = A ./ mean( eig(A) );
draw_line()
fprintf(" Algorithm 3: LS solution full calibration\n")
draw_matrix(A)
draw_matrix(b)

calib.sym.A = A;
calib.sym.b = b;


% https://www.roboticsproceedings.org/rss09/p50.pdf
M = [];
y = [];
for i = 1:N
    Sw = getSkew( gyro(i,:) );
    M(end+1:end+3,:) = Sw;
    y(end+1:end+3,:) = dmag(i, :).' + Sw * mag(i, :).';
end
b = M \ y;
draw_line()
fprintf(" LS solution only bias using gyro and mag data\n")
draw_matrix(b)


%% Recursive-Least-Squares Solutions

if (do_show_RLS_results)
    [b, b_mat, lambda_vec] = ...
        est_mag_bias_RLS_only_mag(mag, lambda_min, p0, scale_mag); %#ok
    draw_line()
    fprintf(" Algorithm 1: adaptive RLS solution only bias\n")
    draw_matrix(b)

    est_rls.alg1_0.b_mat = b_mat;
    est_rls.alg1_0.lambda_vec = lambda_vec;


    % same as above
    [b, b_mat, lambda_vec] = ...
        est_mag_bias_RLS_only_mag_c_implementation(mag, lambda_min, p0, scale_mag);
    draw_line()
    fprintf(" Algorithm 1: adaptive RLS solution only bias c implementation\n")
    draw_matrix(b)

    est_rls.alg1_1.b_mat = b_mat;
    est_rls.alg1_1.lambda_vec = lambda_vec;


    [theta, theta_mat, lambda_vec] = ...
        est_mag_bias_and_scale_RLS_only_mag(mag, lambda_min_bias_and_scale, p0, scale_mag); %#ok
    b = (-0.5 * theta(4:6) ./ theta(1:3)) * scale_mag;
    a = sqrt(theta(1:3));
    A = diag(a);
    A = A ./ mean( diag(A) );
    draw_line()
    fprintf(" Algorithm 2: adaptive RLS solution bias and scaling\n")
    draw_matrix(A)
    draw_matrix(b)

    b_mat = zeros(size(b_mat));
    scale_mat = zeros(size(b_mat));
    for i = 1:size(b_mat, 1)
        b_mat(i,:) = (-0.5 * theta_mat(i,4:6) ./ theta_mat(i,1:3)) * scale_mag;
        scale_mat(i,:) = sqrt(theta_mat(i,1:3));
        scale_mat(i,:) = scale_mat(i,:) ./ mean(scale_mat(i,:));
    end
    est_rls.alg2_0.b_mat = b_mat;
    est_rls.alg2_0.scale_mat = scale_mat;
    est_rls.alg2_0.lambda_vec = lambda_vec;


    [b, scale, b_mat, scale_mat, lambda_vec] = ...
        est_mag_bias_and_scale_RLS_only_mag_c_implementation(mag, lambda_min_bias_and_scale, p0, scale_mag); %#ok
    draw_line()
    fprintf(" Algorithm 2: adaptive RLS solution bias and scaling c implementation\n")
    draw_matrix(diag(scale))
    draw_matrix(b)

    est_rls.alg2_1.b_mat = b_mat;
    est_rls.alg2_1.scale_mat = scale_mat;
    est_rls.alg2_1.lambda_vec = lambda_vec;
    

    [param, param_mat, lambda_vec] = ...
        est_mag_bias_scale_and_rotation_RLS_only_mag(mag, lambda_min_bias_scale_and_rotation, p0); %#ok
    [b, axes, R] = polyToParams3D( param.' );
    A = R * diag(1./axes) * R.';
    A = A ./ mean( eig(A) );
    draw_line()
    fprintf(" Algorithm 3: adaptive RLS solution full calibration\n")
    draw_matrix(A)
    draw_matrix(b)


    [b, b_mat, lambda_vec] = ...
        est_mag_bias_RLS(mag, dmag, gyro, lambda_min, p0);
    draw_line()
    fprintf(" adaptive RLS solution only bias using gyro and mag data\n")
    draw_matrix(b)


    % same as above
    [b, b_mat, lambda_vec] = ...
        est_mag_bias_RLS_c_implementation(mag, dmag, gyro, lambda_min, p0);
    draw_line()
    fprintf(" adaptive RLS solution only bias using gyro and mag data c implementation\n")
    draw_matrix(b)
end


%% visualization

mag_calib = (mag - calib.eye.b.') * calib.eye.A.';
mag_calib_norm(:,1) = sqrt( sum( mag_calib.^2  , 2) );
mag_calib = (mag - calib.diag.b.') * calib.diag.A.';
mag_calib_norm(:,2) = sqrt( sum( mag_calib.^2  , 2) );
mag_calib = (mag - calib.sym.b.') * calib.sym.A.';
mag_calib_norm(:,3) = sqrt( sum( mag_calib.^2  , 2) );

figure(2), clf
ax(1) = subplot(221);
plot(time, mag), grid on
title('uncalibrated')
ax(2) = subplot(222);
plot(time, mag_calib), grid on
title('calibrated')
linkaxes(ax, 'xy'), clear ax
xlim([0 time(end)])
ax(1) = subplot(223);
plot(time, sqrt( sum( mag.^2, 2) )), grid on
ax(2) = subplot(224);
plot(time, mag_calib_norm), grid on
legend('bias', 'bias and scaling', 'bias, scaling and rotation', 'location', 'northeast')
linkaxes(ax, 'xy'), clear ax
xlim([0 time(end)])


% only show datapoints every Ts_show sec
Ts_show = 0.05;
ind_show = (1:ceil(Ts_show / Ts):length(time)).';

figure(3)
plot3(mag(ind_show,1), mag(ind_show,2), mag(ind_show,3), 'b.', 'MarkerSize', 10), grid on, hold on
plot3(mag_calib(ind_show,1), mag_calib(ind_show,2), mag_calib(ind_show,3), '.', 'color', [0 0.5 0], 'MarkerSize', 10), hold off
axis equal
legend('uncalibrated', 'calibrated', 'location', 'northeast')

%%

dx_bin = 5;
x_bins = floor(min(min(mag_calib_norm))/dx_bin)*dx_bin:dx_bin:ceil(max(max(mag_calib_norm))/dx_bin)*dx_bin;

% figure(6), clf
% ax(1) = subplot(311);
% bin_count = histc(mag_calib_norm(:,1), x_bins); %#ok
% bar(x_bins, bin_count, 'hist'), grid on
% set(findobj(gca, 'Type', 'patch'), 'FaceColor', [0 0 1])
% title('bias')
% ax(2) = subplot(312);
% bin_count = histc(mag_calib_norm(:,2), x_bins); %#ok
% bar(x_bins, bin_count, 'hist'), grid on
% set(findobj(gca, 'Type', 'patch'), 'FaceColor', [0 0.5 0])
% title('bias and scaling')
% ax(3) = subplot(313);
% bin_count = histc(mag_calib_norm(:,3), x_bins); %#ok
% bar(x_bins, bin_count, 'hist'), grid on
% set(findobj(gca, 'Type', 'patch'), 'FaceColor', [1 0 0])
% title('bias, scaling and rotation')
% linkaxes(ax, 'xy'), clear ax

figure(6), clf
bin_count = histc(mag_calib_norm, x_bins); %#ok
h = bar(x_bins, bin_count(:,1), 'hist'); grid on, hold on
set (h, 'FaceColor', [0 0 1], "FaceAlpha", 0.5);
h = bar(x_bins, bin_count(:,2), 'hist');
set (h, 'FaceColor', [0 0.5 0], "FaceAlpha", 0.5);
h = bar(x_bins, bin_count(:,3), 'hist'); hold off
set (h, 'FaceColor', [1 0 0], "FaceAlpha", 0.5);
legend('bias', 'bias and scaling', 'bias, scaling and rotation', 'location', 'northeast')


%% not important for other users

if (do_compare_RLS_results)

    y_lim_upper = ceil (max(calib.eye.b) / 100) * 100 + 100;
    y_lim_lower = floor(min(calib.eye.b) / 100) * 100 - 100;

    % y_lim_upper = 1200;
    % y_lim_lower = 0;

    figure(7)
    subplot(321)
    plot(time, [est_rls.alg1_0.b_mat, est_rls.alg1_1.b_mat]), grid on, hold on, title('Alg.1: bias')
    plot([time(1); time(end)], [calib.eye.b.'; calib.eye.b.']), hold off
    xlim([0 time(end)]), ylim([y_lim_lower, y_lim_upper])
    subplot(322)
    plot(time, [est_rls.alg1_0.lambda_vec, est_rls.alg1_1.lambda_vec]), grid on, title('Alg.1: lamda')
    xlim([0 time(end)]), ylim([lambda_min, 1.0])
    subplot(323)
    plot(time, [est_rls.alg2_0.b_mat, est_rls.alg2_1.b_mat]), grid on, hold on, title('Alg.2: bias')
    plot([time(1); time(end)], [calib.diag.b.'; calib.diag.b.']), hold off
    xlim([0 time(end)]), ylim([y_lim_lower, y_lim_upper])
    subplot(324)
    plot(time, [est_rls.alg2_0.lambda_vec, est_rls.alg2_1.lambda_vec]), grid on, title('Alg.2: lamda')
    xlim([0 time(end)]), ylim([lambda_min, 1.0])
    subplot(325)
    est_rls.alg2_0.scale_mat(imag(est_rls.alg2_0.scale_mat) ~= 0) = 0;
    est_rls.alg2_1.scale_mat(imag(est_rls.alg2_1.scale_mat) ~= 0) = 0;
    plot(time, [est_rls.alg2_0.scale_mat, est_rls.alg2_1.scale_mat]), grid on, hold on, title('Alg.2: scale')
    plot([time(1); time(end)], [diag(calib.diag.A).'; diag(calib.diag.A).']), hold off
    xlim([0 time(end)]), ylim([0.95, 1.05])
    subplot(326)
    plot(time, [est_rls.alg2_0.lambda_vec, est_rls.alg2_1.lambda_vec]), grid on%, title('Alg.2: lamda')
    xlim([0 time(end)]), ylim([lambda_min, 1.0])

    if (false)
        % compare it to c project
        data_from_c_project = readmatrix(['05_output_from_c_project/', file_name(1:end-8), '.txt']); %#ok
        figure(8)
        subplot(121)
        plot(time, [est_rls.alg1_1.b_mat, data_from_c_project(:,1:3)]), grid on, title('bias')
        xlim([0 time(end)]), ylim([y_lim_lower, y_lim_upper])
        subplot(122)
        plot(time, [est_rls.alg1_1.lambda_vec, data_from_c_project(:,4)]), grid on, title('lamda')
        xlim([0 time(end)]), ylim([lambda_min, 1.0])
    end
end

