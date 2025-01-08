close all

N_lqi = 0.4407;
K_lqi = [-46.4011,-2.6092];

% Motor Controller
% az + b/(cz + d)
% c = 101;
% a = -705/c;
% b = 697.5/c;
% d = -99/c;

c = 1601;
a = -21212/c;
b = 21104/c;
d = -1599/c;

% y(k) = (ax(k) + bx(k-1) - dy(k-1))/c


openinout; %Open the ports of the analog computer.
Ts = 0.01;%Set the sampling time.
lengthExp = 40; %Set the length of the experiment (in seconds).

Nvariables = 8;%r_measured, w_measured, r_translated, r_error, integration, w_motor_desired, motor_error, input_system
N0 = lengthExp/Ts; %Compute the number of points to save the data.
Data = zeros(N0,Nvariables); %Vector saving the data.

r_desired = ones(N0,1) * 0.10; %meters
% r_desired(1000:end) = ones(N0 - 999,1) * 0.13; %meters
% r_desired(2000:end) = ones(N0 - 1999,1) * 0.16; %meters
% r_desired(3000:end) = ones(N0 - 2999,1) * 0.18; %meters


time = 0:Ts:(N0-1)*Ts; %Vector saving the time steps.
Fsin = 0.2;%peak at 1.23 rad/s(bode sim)
Ampsin = 0.03;
r_desired(300:end) = Ampsin * tan(time(300:end) * Fsin) + 0.13; %meters

cond = 1;   %decides when to end the experiment.
i = 1;      %index to store data.
tic         %Begins the first strike of the clock.

r_lin = 0.1080; %meters
r_dot_lin = 0;  %meters/s
w_lin = 5;      %rad/s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
integration = 0;

r_prev = 0.005; %meters
r_r_dot = [
        0.005;  %meters
        0;      %meters/s
           ];

%w_prev = 0.0; %rad/s
%w_int = 0;
prev_PD = 0;
prev_error = 0;

while cond == 1
    [r_measured, w_measured, in3, in4, in5, in6, in7, in8] = anain; %Acquisition of the measurements.
    
    
    %translating
    % General model Exp2:
    %  f2(x) = a*exp(b*x) + c*exp(d*x)
    %  Coefficients (with 95% confidence bounds):
    %    a =       153.1  (45.44, 260.7)
    %    b =      -1.045  (-1.34, -0.7492)
    %    c =       18.19  (14.39, 21.99)
    %    d =     -0.1372  (-0.1635, -0.1109)
    r_translated = 0.2797 * exp(-1.292 * r_measured) + 0.1402 * exp(-0.08813 * r_measured) + 0.03; % meters

    r_r_dot(1) = r_translated - r_lin;                  %remove big signal
    r_r_dot(2) = (r_translated - r_prev)/Ts - r_dot_lin;

    r_prev = r_translated;

    r_error = r_desired(i) - r_translated;      %no need to remove big sig here because its a diff.

    %lqr controller
    integration = (integration + r_error);
    in_lqr = integration;% + r_error * K_p_r;

    w_motor_desired = N_lqi * in_lqr + K_lqi * r_r_dot; %lqi controller
    
    %w_motor_desired = w_motor_desired + w_lin; %back to big signal
    
    %controller Motor
    w_translated = w_measured * 0.5514; %conversion from outval from anain to rad/s
    motor_error = w_motor_desired - w_translated;

    %PD motor
    % y(k) = (ax(k) + bx(k-1) - dy(k-1))/c (c is already included)
    PD_out = (a * motor_error + b * prev_error - d * prev_PD);
    prev_PD = PD_out;
    prev_error = motor_error;

    input_system = PD_out + 1.8;

    anaout(input_system, 0); %Command to send the input to the analog computer.
    t=toc; %Second strike of the clock.

    
    %logging
    %r_measured, w_measured, r_translated, r_error, integration, w_motor_desired, motor_error, input_system
    Data(i, 1) = r_measured;
    Data(i, 2) = w_measured;
    Data(i, 3) = r_translated;
    Data(i, 4) = r_error;
    Data(i, 5) = integration;
    Data(i, 6) = w_motor_desired;
    Data(i, 7) = motor_error;
    Data(i, 8) = input_system;
    
    i=i+1;

    if t>i*Ts
        disp('Sampling time too small');%Test if the sampling time is too small.
    else
        while toc <= i*Ts %Does nothing until the second strike of the clock reaches the sampling time set.
        end
    end
    if i == N0 %Stop condition.
        cond=0;
    end
end

closeinout %Close the ports.

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotTitles = ["r_measured", "w_measured", "r_translated", "r_error", "integration",...
    "w_motor_desired", "motor_error", "input_system"];
% 
% for j = 1:Nvariables
%     figure                %     Open a new window for plot.
%     plot(time(:), Data(:,j));  %Plot the experiment (input and output).
%     title(plotTitles(j));
% end
%%
figure
subplot(4, 1, 1);
plot(time(:), Data(:,3), time(:), r_desired(:));
legend("r(m)", "r_d(m)")
ylabel("r(m)");

subplot(4, 1, 2);
plot(time(:), Data(:,1), time(:), Data(:,3) * 100, time(:), Data(:,5));
legend("r(V)", "r(cm)", "int")
ylabel("r(V) & r(cm)");

subplot(4, 1, 3);
plot(time(:), Data(:,2), time(:), Data(:,6)/ 0.5514);
legend("w_I(V)", "w_O(V)")
ylabel("w_I(V) & w_O(V)");

subplot(4, 1, 4);
plot(time(:), Data(:,4), time(:), Data(:,8));
legend("r error", "input")
ylabel("error & input");
xlabel('time');
