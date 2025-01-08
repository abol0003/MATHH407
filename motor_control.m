close all

% lower: -1.75 1.7
% sat: -8.42 8.30
% no integrator because lqi

Ts = 0.01;
DataCommand = load("motorIdent\DataCommands.mat").DataCommands;
MotorSpeed = load("motorIdent\Data.mat").Data(:, 2);

DataCommand_steady = DataCommand(1001:length(DataCommand)) - DataCommand(1001);
MotorSpeed_steady = MotorSpeed(1001:length(MotorSpeed)) - MotorSpeed(1001);

scale = DataCommand_steady(2);

time = 0:Ts:(length(MotorSpeed) - 1)*Ts;

time_steady = 0:Ts:(length(MotorSpeed) - 1001)*Ts;
%%
figure
plot(time, DataCommand);

figure
plot(time, MotorSpeed);

figure
plot(time_steady, DataCommand_steady);

figure
plot(time_steady, MotorSpeed_steady);

H = tf([3.3286], [2 1]);

Hd = c2d(H, Ts, "tustin");

%C = pidtune(H, "PD")

figure
stepplot(Hd*scale);

hold on
plot(time_steady, MotorSpeed_steady);

figure
lsim(Hd, DataCommand_steady, time_steady)

hold on
plot(time_steady, MotorSpeed_steady);

%%
close all
Ts = 0.01;

Md = -tf([0,0,0.000138421208189068], [0.0172840168098617, -0.0244481169001018, 0.00720478979897442], Ts)

Mc = d2c(Md, 'tustin');
%%
figure
lsim(Md , DataCommand_steady, time_steady)
hold on
plot(time_steady, -MotorSpeed_steady);
%%
close all
%use lead compensator to remove effect of noise on the D action

Ts = 0.01;
Td = 8.5;
Tdl = 8;
Kp = 12;

%Cc = tf([Kp], [1]) + tf([Td 0], [Tdl 1]);

Cd = -tf([(Td + Kp*Tdl + Kp*Ts/2)*2/Ts Kp*Tdl - (Td + Tdl*Kp)*2/Ts], [(1 + 2*Tdl/Ts) (1 - 2*Tdl/Ts)], Ts);

%az + b/(cz + d)
%y(k) = (ax(k) + bx(k-1) - dy(k-1))/c
c = 1601;
a = -21212/c;
b = 21104/c;
d = -1599/c;

Cd = tf([-21212 21104], [1601 -1599], Ts)

ZOH = tf([1 1], [2 0], Ts);

% figure
% pzmap(Cd)
% 
% figure
% bode(Cd)
% 
% figure
% pzmap(Md)
% 
% figure
% rlocus(Cd * Md * ZOH);

figure
pzmap(Md * Cd * ZOH / (1 + (Md * Cd * ZOH)))

figure
bode(Md * Cd * ZOH / (1 + Md * Cd * ZOH))

figure
step(Md * Cd * ZOH / (1 + Md * Cd * ZOH))
ylabel("Rotational Speed (/)")

time = 0:Ts:100*Ts;
figure
lsim(Md * Cd * ZOH / (1 + Md * Cd * ZOH), time, time)%talud
title("Talud Response");
ylabel("Rotational speed (/)")

%%
close all

outW_to_rotPer10s = [
    2 5.5;
    2.5 7;
    2.5 7.2;
    3 8.8;
    3 8.5;
    3 8.5;
    3 8.2;
    3.3 9.2;
    3.3 9.5;
    3.5 10;
    3.5 10.5;
    3.5 10.8;
    3.5 11;
    3.5 10.5;
    3.5 10;
    4 11;
    4 12;
    4 11;
    4.5 13;
    4.5 13;
    4.5 13;
    5 14;
];

outW_to_radPersec = outW_to_rotPer10s;
outW_to_radPersec(:, 2) = outW_to_radPersec(:, 2)* 2 * pi/10; %2pi radians per rev

flin = fit(outW_to_radPersec(:, 1), outW_to_radPersec(:, 2), "a*x + b" )
fzero = fit(outW_to_radPersec(:, 1), outW_to_radPersec(:, 2), "a*x" )

figure
plot(flin, outW_to_radPersec(:, 1), outW_to_radPersec(:, 2));
xlim([0 6])
ylim([0 13])

figure
plot(fzero, outW_to_radPersec(:, 1), outW_to_radPersec(:, 2));
xlim([0 6])
ylim([0 13])

flin_inv = fit(outW_to_radPersec(:, 2), outW_to_radPersec(:, 1), "a*x + b" )
fzero_inv = fit(outW_to_radPersec(:, 2), outW_to_radPersec(:, 1), "a*x" )

figure
plot(flin_inv, outW_to_radPersec(:, 2), outW_to_radPersec(:, 1));
ylim([0 6])
xlim([0 13])

figure
plot(fzero_inv, outW_to_radPersec(:, 2), outW_to_radPersec(:, 1));
ylim([0 6])
xlim([0 13])

%-> a = 1.811 (rad/s)/outW
%=> a = 0.5514 outW/(rad/s)

%desired out, real out, rot/10s
%  2                2     5.5
%  2.5            2.5       7
%  2.5            2.5       7.2

% 3               3.028     8.8
% 3               3.028     8.5
% 3               3.028     8.5

% 3               3.028     8.2
%3.3                        9.2
%3.3                        9.5

% 3.5             3.495     10
% 3.5             3.495     10.5
% 3.5             3.495     10.8

% 3.5             3.495     11
% 3.5             3.495     10.5
% 3.5             3.499     10

% 4               3.996     11
% 4               3.996     12
%4                          11

%4.5              4.529     13
%4.5              4.529     13
%4.5              4.5       13
%5                5.064     14