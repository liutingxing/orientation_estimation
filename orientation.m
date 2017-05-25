clear all;
clc;
addpath('.\functions\basic');
addpath('.\functions\mag_calibration');

%% read data from file

ahrs_data = load('.\data\rotate.txt');

Time = ahrs_data(:, 1);                      % ( ms )
Roll = ahrs_data(:, 3);                      % ( degree )
Pitch = ahrs_data(:, 4);                     % ( degree )
Yaw = ahrs_data(:, 5);                       % ( degree )
Gyro = ahrs_data(:, 6:8) * pi / 180;         % ( rad/s )
Acc = ahrs_data(:, 9:11) * 9.8;              % ( m/s2 )
Mag = ahrs_data(:, 12:14);                   % ( count )
M = length(Time);                            % the number of measurement

%% configure

SampleRate = 100;   % Hz
AlignmentTime = 2;  % second

%% mag calibration 

% (mag_calibrated = inv(W)*(mag_raw - V)
[Cali_B, Cali_V, Cali_W_inv, Cali_Error] = cali7eig(Mag);
% calibrated mag data
Mag_c = zeros(M, 3);
for i = 1:M
    Bp = Mag(i, :);
    Bc = Cali_W_inv*(Bp - Cali_V)';
    Mag_c(i, :) = Bc;
end

%% initial alignment

N = SampleRate*AlignmentTime;
initial_alignment_flag = 0;
acc_alignment = Acc(1:N, :);
mag_alignment = Mag(1:N, :);

% static check
acc_x_var = var(acc_alignment(1:N, 1));
acc_y_var = var(acc_alignment(1:N, 2));
acc_z_var = var(acc_alignment(1:N, 3));
mag_x_var = var(mag_alignment(1:N, 1));
mag_y_var = var(mag_alignment(1:N, 2));
mag_z_var = var(mag_alignment(1:N, 3));
acc_var_threshold = 0.01;
mag_var_threshold = 1;
if acc_x_var < acc_var_threshold && acc_y_var < acc_var_threshold && acc_z_var < acc_var_threshold && ...
   mag_x_var < mag_var_threshold && mag_y_var < mag_var_threshold && mag_z_var < mag_var_threshold
   initial_alignment_flag = 1;
end

if initial_alignment_flag == 1
    g_x_mean = -mean(acc_alignment(1:N, 1));
    g_y_mean = -mean(acc_alignment(1:N, 2));
    g_z_mean = -mean(acc_alignment(1:N, 3));
    mag_x_mean = mean(mag_alignment(1:N, 1));
    mag_y_mean = mean(mag_alignment(1:N, 2));
    mag_z_mean = mean(mag_alignment(1:N, 3));
    Cnb_initial = ecompass_ned([g_x_mean, g_y_mean, g_z_mean], [mag_x_mean, mag_y_mean, mag_z_mean]);
    Cbn_initial = Cnb_initial';
    [yaw_initial, pitch_initial, roll_initial] = dcm2euler(Cbn_initial);
    % compute the geomagnetic inclination angle
    geoB = norm([mag_x_mean, mag_y_mean, mag_z_mean]);
    gmod = norm([g_x_mean, g_y_mean, g_z_mean]);
    geo_inclination = asin(dot([mag_x_mean, mag_y_mean, mag_z_mean], [g_x_mean, g_y_mean, g_z_mean])/geoB/gmod); % rad
%     Mag_vector = [geoB*cos(geo_inclination), 0, geoB*sin(geo_inclination)]';
    Mag_vector = [Cali_B*cos(geo_inclination), 0, Cali_B*sin(geo_inclination)]';
end

%% sensor fusion variable

Ge = 9.80665;
G_vector = [0, 0, Ge]';
yaw = zeros(N, 1);
pitch = zeros(N, 1);
roll = zeros(N, 1);
gyro_bias = zeros(3, 1);
acc_bias = zeros(3, 1);
mag_jamming = zeros(3, 1);
Cnb = eye(3, 3);
Cbn = eye(3, 3);

StateNum = 12;
MeasNumG = 3;
MeasNumM = 3;
x = zeros(StateNum, 1); % roll, pitch, yaw, gyro_bias_x, gyro_bias_y, gyro_bias_z, acc_bias_x, acc_bias_y, acc_bias_z, mag_jamming_x, mag_jamming_y, mag_jamming_z
Corr_time_gyro = 0.01;
Corr_time_acc = 0.01;
Corr_time_mag = 0.1;
sigma_Win = 1.0e-6;
sigma_gyro_1 = (2*pi/180/3600)^2;   % Markov process
sigma_gyro_2 = (1*pi/180/3600)^2;   % Random walk
sigma_acc = ((5.0e-2)*Ge)^2;  % Markov process
sigma_mag = 0.1*0.1;  % Markov process

error_phim_e = 1.0*pi/180;
error_phim_n = 1.0*pi/180;
error_phim_u = 1.0*pi/180;
error_gyro = 1000*pi/180/3600;
error_acc = 0.3;
error_mag_jamming = 5;
P = zeros(StateNum, StateNum);
P(1, 1) = error_phim_e^2;
P(2, 2) = error_phim_n^2;
P(3, 3) = error_phim_u^2;
P(4, 4) = error_gyro^2;
P(5, 5) = error_gyro^2;
P(6, 6) = error_gyro^2;
P(7, 7) = error_acc^2;
P(8, 8) = error_acc^2;
P(9, 9) = error_acc^2;
P(10, 10) = error_mag_jamming^2;
P(11, 11) = error_mag_jamming^2;
P(12, 12) = error_mag_jamming^2;

%% AHRS estimation

q = euler2q(yaw_initial, pitch_initial, roll_initial);
for i = 1 : N
    yaw(i) = yaw_initial;
    pitch(i) = pitch_initial;
    roll(i) = roll_initial;
end

for i = N+1 : M
    dt = 1/SampleRate;
    
    % strapdown mechanization
    Cbn = q2dcm(q);
    Cnb = Cbn';
    
    Wepp = zeros(3, 1); % no latitude information in computing latitude and longitude rate
    Wiep = zeros(3, 1); % no latitude information in computing earth rate in the navigation frame
    Wipp = Wiep + Wepp;
    Wipb = Cnb * Wipp;
    Wpbb = Gyro(i, :)' - gyro_bias - Wipb;
    
    dq = zeros(4, 1);
    dq(1) = -(Wpbb(1)*q(2) + Wpbb(2)*q(3) + Wpbb(3)*q(4))/2;
    dq(2) = (Wpbb(1)*q(1) + Wpbb(3)*q(3) - Wpbb(2)*q(4))/2;
    dq(3) = (Wpbb(2)*q(1) - Wpbb(3)*q(2) + Wpbb(1)*q(4))/2;
    dq(4) = (Wpbb(3)*q(1) + Wpbb(2)*q(2) - Wpbb(1)*q(3))/2;
    
    q = q + dq*dt;
    q = q_norm(q);
    
    [yaw(i), pitch(i), roll(i)] = dcm2euler(Cbn);   % pure gyro estimation

    % kalman sensor fusion
    Cbn = q2dcm(q);
    
	Fg = diag([-1/Corr_time_gyro, -1/Corr_time_gyro, -1/Corr_time_gyro]);
    Fa = diag([-1/Corr_time_acc, -1/Corr_time_acc, -1/Corr_time_acc]);
    Fm = diag([-1/Corr_time_mag, -1/Corr_time_mag, -1/Corr_time_mag]);
    Fnull = zeros(3, 3);
    
    F = [Fnull, -Cbn, Fnull, Fnull;
         Fnull, Fg,   Fnull, Fnull;
         Fnull, Fnull,Fa,    Fnull;
         Fnull, Fnull,Fnull, Fm  ];
     
    qdt = diag([sigma_Win, sigma_Win, sigma_Win, sigma_gyro_1, sigma_gyro_1, sigma_gyro_1, sigma_gyro_2,...
                sigma_gyro_2, sigma_gyro_2, sigma_acc, sigma_acc, sigma_acc, sigma_mag, sigma_mag, sigma_mag]);
    I = eye(3, 3);
    Gnull = zeros(3, 3);
    G = [-Cbn,    Gnull,  -Cbn,  Gnull,  Gnull;
          Gnull,      I,  Gnull, Gnull,  Gnull;
          Gnull,  Gnull,  Gnull,     I,  Gnull;
          Gnull,  Gnull,  Gnull, Gnull,      I];
	% Q matrix discretization-2 order
    Q_basic = G*qdt*G';
    M1 = Q_basic;
    M2 = Q_basic*F'+F*Q_basic;
    Q = dt*M1 + 1/2*dt*dt*M2;
    
    % PHIM matrix discretization-2 order
    I = eye(StateNum, StateNum);
    PHIM = I + dt*F + 1/2*dt*dt*F*F;
    
    % predict
    x = PHIM*x;
    P = PHIM*P*PHIM' + Q;
    
    % update from acc
    H = zeros(3, StateNum);
    H(1, 2) = G_vector(3);
    H(2, 1) = -G_vector(3);
    H(1, 7) = Cbn(1, 1);
    H(1, 8) = Cbn(1, 2);
    H(1, 9) = Cbn(1, 3);
    H(2, 7) = Cbn(2, 1);
    H(2, 8) = Cbn(2, 2);
    H(2, 9) = Cbn(2, 3);
    H(3, 7) = Cbn(3, 1);
    H(3, 8) = Cbn(3, 2);
    H(3, 9) = Cbn(3, 3);
    
    R = eye(MeasNumG, MeasNumG);
    R(1, 1) = 0.5^2;
    R(2, 2) = 0.5^2;
    R(3, 3) = 0.5^2;
    
    acc_liner_b = zeros(3, 1);
% only for acc bias estimate test
%     Acc(i, 1) = Acc(i, 1) + 1;
%     Acc(i, 2) = Acc(i, 2) + 1;
%     Acc(i, 3) = Acc(i, 3) + 1;
% only for acc bias estimate test
    g_estimate = Cbn*(acc_bias - Acc(i, :)' + acc_liner_b);
    Z = G_vector - g_estimate;
    K = P*H'*((H*P*H'+R)^-1);
    x = x + K*(Z - H*x);
    P = (I - K*H)*P;
    
    [deltCbn] = euler2dcm (x(3), x(2), x(1)); % (I+P)Cbn
    Cbn = deltCbn*Cbn;
    [yaw(i), pitch(i), roll(i)] = dcm2euler(Cbn);   % gyro + acc estimation
    q = euler2q(yaw(i), pitch(i), roll(i));
    q = q_norm(q);
    
    % update from mag
    H = zeros(3, StateNum);
    H(1, 2) = Mag_vector(3);
    H(2, 1) = -Mag_vector(3);
    H(2, 3) = Mag_vector(1);
    H(3, 2) = -Mag_vector(1);
    H(1, 10) = -Cbn(1, 1);
    H(1, 11) = -Cbn(1, 2);
    H(1, 12) = -Cbn(1, 3);
    H(2, 10) = -Cbn(2, 1);
    H(2, 11) = -Cbn(2, 2);
    H(2, 12) = -Cbn(2, 3);
    H(3, 10) = -Cbn(3, 1);
    H(3, 11) = -Cbn(3, 2);
    H(3, 12) = -Cbn(3, 3);
    
    R = eye(3, 3);
    R(1, 1) = 5^2;
    R(2, 2) = 5^2;
    R(3, 3) = 5^2;
    
    mag_estimate = Cbn*Mag(i, :)';
    Z = Mag_vector - mag_estimate;
    
    K = P*H'*((H*P*H'+R)^-1);
    x = x + K*(Z - H*x);
    
    % mag jamming check
    mag_jamming_b = x(10:12);
    if norm(mag_jamming_b) > 2*Cali_B^2
        disp('mag jamming exist at %d', Time(i));
    else
        P = (I - K*H)*P;
        [deltCbn] = euler2dcm (x(3), x(2), x(1)); % (I+P)Cbn
        Cbn = deltCbn*Cbn;
        [yaw(i), pitch(i), roll(i)] = dcm2euler(Cbn);   % gyro + acc + mag estimation
        q = euler2q(yaw(i), pitch(i), roll(i));
        q = q_norm(q);
        
        % mag vector correction
        mag_jamming_n = Cbn*mag_jamming_b;
        mag_vector_temp = Mag_vector - mag_jamming_n;
        geoB_temp = sqrt(mag_vector_temp(1)^2 + mag_vector_temp(3)^2);
        sindelta = mag_vector_temp(3)/geoB_temp;
        cosdelta = mag_vector_temp(1)/geoB_temp;
        SINDELTAMAX = 0.9063078; % sin of max +ve geomagnetic inclination angle: here 65.0 deg
        COSDELTAMAX = 0.4226183; % cos of max +ve geomagnetic inclination angle: here 65.0 deg
        if sindelta > SINDELTAMAX || sindelta < -SINDELTAMAX
            sindelta = SINDELTAMAX;
            cosdelta = COSDELTAMAX;
        end
        geo_inclination = asin(sindelta);
        Mag_vector = [Cali_B*cos(geo_inclination), 0, Cali_B*sin(geo_inclination)]';
    end
    
    % feedback
    acc_bias = acc_bias + x(7:9);
    gyro_bias = gyro_bias + x(4:6);
    x(1:StateNum) = 0;
end

%% display result

% yaw
figure;
plot(yaw*180/pi, 'r');
hold on;
plot(Yaw, 'b');
legend('unicore', 'daoyuan');
title('yaw comparison');
xlabel('sample point');
ylabel('yaw (degree)');

% pitch
figure;
plot(pitch*180/pi, 'r');
hold on;
plot(Pitch, 'b');
legend('unicore', 'daoyuan');
title('pitch comparison');
xlabel('sample point');
ylabel('pitch (degree)');

% roll
figure;
plot(roll*180/pi, 'r');
hold on;
plot(Roll, 'b');
legend('unicore', 'daoyuan');
title('roll comparison');
xlabel('sample point');
ylabel('roll (degree)');









































