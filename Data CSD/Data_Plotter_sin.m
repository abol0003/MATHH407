
Data1 = load("Sine Wave Data\Fsin1_Data.mat").Data(:,3);
r_desired1 = load("Sine Wave Data\Fsin1_Desired.mat").r_desired;

Ts = 0.01;
N0 = length(r_desired1);
time1 = 0:Ts:(N0-1)*Ts;

figure
plot(time1(1:4000), Data1(1:4000), time1(1:4000), r_desired1(1:4000));
legend("r(m)", "r_d(m)")
ylabel("r(m)");
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data02 = load("Sine Wave Data\Fsin02_Data.mat").Data(:,3);
r_desired02 = load("Sine Wave Data\Fsin02_Desired.mat").r_desired;

Ts = 0.01;
N0 = length(r_desired02);
time02 = 0:Ts:(N0-1)*Ts;

figure
plot(time02(:), Data02(:), time02(:), r_desired02(:));
legend("r(m)", "r_d(m)")
ylabel("r(m)");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data2 = load("Sine Wave Data\Fsin2_Data.mat").Data(:,3);
r_desired2 = load("Sine Wave Data\Fsin2_Desired.mat").r_desired;

Ts = 0.01;
N0 = length(r_desired2);
time2 = 0:Ts:(N0-1)*Ts;

figure
plot(time2(:), Data2(:), time2(:), r_desired2(:));
legend("r(m)", "r_d(m)")
ylabel("r(m)");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data10 = load("Sine Wave Data\Fsin10_Data.mat").Data(:,3);
r_desired10 = load("Sine Wave Data\Fsin10_Desired.mat").r_desired;

Ts = 0.01;
N0 = length(r_desired10);
time10 = 0:Ts:(N0-1)*Ts;

figure
plot(time10(:), Data10(:), time10(:), r_desired10(:));
legend("r(m)", "r_d(m)")
ylabel("r(m)");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data15 = load("Sine Wave Data\Fsin15_Data.mat").Data(:,3);
r_desired15 = load("Sine Wave Data\Fsin15_Desired.mat").r_desired;

Ts = 0.01;
N0 = length(r_desired15);
time15 = 0:Ts:(N0-1)*Ts;

figure
plot(time15(:), Data15(:), time15(:), r_desired15(:));
legend("r(m)", "r_d(m)")
ylabel("r(m)");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data123 = load("Sine Wave Data\Fsin123_Data.mat").Data(:,3);
r_desired123 = load("Sine Wave Data\Fsin123_Desired.mat").r_desired;

Ts = 0.01;
N0 = length(r_desired123);
time123 = 0:Ts:(N0-1)*Ts;

figure
plot(time123(:), Data123(:), time123(:), r_desired123(:));
legend("r(m)", "r_d(m)")
ylabel("r(m)");
%%
figure
subplot(3, 1, 1);
plot(time02(:), Data02(:), time02(:), r_desired02(:));
legend("r(m)", "r_d(m)", "FontSize", 13, "Location", "southeast")
ylabel("r (m)", "FontSize", 13);
xlabel("Time (s)", "FontSize", 13)
title("Sine Wave (F = 0.2 rad/s)", "FontSize", 14)
ylim([0 0.2]);

subplot(3, 1, 2);
plot(time1(1:4000), Data1(1:4000), time1(1:4000), r_desired1(1:4000));
legend("r(m)", "r_d(m)", "FontSize", 13, "Location", "southeast")
ylabel("r(m)", "FontSize", 13);
xlabel("Time (s)", "FontSize", 13)
title("Sine Wave (F = 1 rad/s)", "FontSize", 14)
ylim([0 0.2]);

subplot(3, 1, 3);
plot(time15(1:4000), Data15(1:4000), time15(1:4000), r_desired15(1:4000));
legend("r(m)", "r_d(m)", "FontSize", 13, "Location", "southeast")
ylabel("r(m)", "FontSize", 13);
xlabel("Time (s)", "FontSize", 13)
title("Sine Wave (F = 1.5 rad/s)", "FontSize", 14)
ylim([0 0.2]);

%%
freq_resp = [
    1.23 0.03 0.023
    1.5 0.03 0.021
    10 0.03 0.00006
    2 0.03 0.019
    0.2 0.03 0.03
    1 0.03 0.023
];

freq_resp = sortrows(freq_resp);

figure
plot(freq_resp(:,1), freq_resp(:,2), freq_resp(:,1), freq_resp(:,3), "Marker","o");
ylim([0 0.05])
%xlim(0.01, 100);

%xticks([0.08 0.2 1 2 1.23 1.5 10 12])
%xticks([0.01 0.1 1 10 100])
%xlim(0.01, 100);
ylim([0 0.05])
xscale("log");
yscale("log");
title("Frequency Response");
ylabel("Amplitude (m)")
xlabel("Frequency (rad/s)")
legend("r_d", "r")


%%
Data = zeros(2000, 5);
Data(:, 1) = Data1(2001:4000);
Data(:, 2) = Data10(2001:4000);
Data(:, 3) = Data123(2001:4000);
Data(:, 4) = Data15(2001:4000);
Data(:, 5) = Data2(2001:4000);

figure                %     Open a new window for plot.
t = time123(2001:4000);
for j = 1:5
    plot(t(:), Data(:,j));  %Plot the experiment (input and output).
    hold on
end

% L = 2000;
% Fs = 1/Ts;
% t = (0:L-1)*Ts;        % Time vector
% 
% 
% plot(Fs/L*(-L/2:L/2-1), abs(fftshift(Data(:, 1))), "LineWidth", 3)
% title("fft Spectrum in the Positive and Negative Frequencies")
% xlabel("f (Hz)")
% ylabel("|fft(X)|")