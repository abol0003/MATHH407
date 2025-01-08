close all

K_i_w = 0.401;      %0
K_d_w = 0.00221;    %0.11
K_p_w = 0.444;      %22


%az + b/(cz + d)
c = 401;
a = -6208/c;
b = 6184/c;
d = 399/c;
%y(k) = (ax(k) + bx(k-1) - dy(k-1))/c


openinout; %Open the ports of the analog computer.
Ts = 0.01;%Set the sampling time.
lengthExp = 20; %Set the length of the experiment (in seconds).

Nvariables = 5;%r_measured, w_measured, w_motor_desired, motor_error, input_system
N0 = lengthExp/Ts; %Compute the number of points to save the data.
Data = zeros(N0,Nvariables); %Vector saving the data.

w_motor_desired = ones(N0,1) * 4; %should become rad/s but at the moment just stabalize for output

cond = 1;   %decides when to end the experiment.
i = 1;      %index to store data.
tic         %Begins the first strike of the clock.
time = 0:Ts:(N0-1)*Ts; %Vector saving the time steps.

r_lin = 0.1080; %meters
r_dot_lin = 0;  %meters/s
w_lin = 5;      %rad/s

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%w_prev = 0.0; %rad/s
%w_int = 0;
prev_PID = 0;

while cond == 1
    [r_measured, w_measured, in3, in4, in5, in6, in7, in8] = anain; %Acquisition of the measurements.

    %controller Motor
    w_translated = w_measured;
    motor_error = w_motor_desired(i) - w_translated;

    %PID motor
    %w_int = w_int + motor_error * K_p_w * Ts;
    %w_dot = (w_translated - w_prev) * K_d_w/Ts;
    %PID_out = K_p_w * motor_error + w_int + w_dot;
    
    % y(k) = (ax(k) + bx(k-1) - dy(k-1))/c (c is already included)
    PID_out = (a * motor_error + b * Data(i - 1, 4) - d * prev_PID);
    prev_PID = PID_out;

    %translate back to input (if necessary)
    input_system = PID_out;

    anaout(input_system, 0); %Command to send the input to the analog computer.
    t=toc; %Second strike of the clock.

    
    %logging
    %r_measured, w_measured, r_translated, r_error, integration, w_motor_desired, motor_error, input_system
    Data(i, 1) = r_measured;
    Data(i, 2) = w_measured;
    Data(i, 3) = w_motor_desired;
    Data(i, 4) = motor_error;
    Data(i, 5) = input_system;
    
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
plotTitles = ["r_measured", "w_measured", "w_motor_desired", "motor_error", "input_system"];

for j = 1:Nvariables
    figure                %     Open a new window for plot.
    plot(time(:), Data(:,j));  %Plot the experiment (input and output).
    title(plotTitles(j));
end