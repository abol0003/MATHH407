
distTranslate2 = [
                 9.103 6;
                 8.85 6.5;
                 8.21 7;
                 6.5 8;
                 4.381 9.5;
                 3.24 11;
                 1.264 18;
                 ];


f5 = fit(distTranslate2(:, 1), distTranslate2(:, 2), 'exp2')


plot(f5,distTranslate2(:, 1), distTranslate2(:, 2))

%%
%old dist translate

distTranslate = [9.35 5;
                 8.2 6;
                 7 7;
                 6.25 8;
                 5.6 9;
                 4.9 10;
                 4.6 11;
                 4.2 12;
                 4 13;
                 3.7 14;
                 3.55 15;
                 3.4 16;
                 3.2 17;
                 3.1 17.9];

global f1
global f2
f2 = fit(distTranslate(:, 1), distTranslate(:, 2), 'exp2')

f1 = fit(distTranslate(:, 1), distTranslate(:, 2), 'exp1')

fp = fit(distTranslate(:, 1), distTranslate(:, 2),'poly2')

finv = fit(distTranslate(:, 1), distTranslate(:, 2).^-1,'poly2')


%General model Exp2:
%f2(x) = a*exp(b*x) + c*exp(d*x)
%Coefficients (with 95% confidence bounds):
%   a =       153.1  (45.44, 260.7)
%   b =      -1.045  (-1.34, -0.7492)
%   c =       18.19  (14.39, 21.99)
%   d =     -0.1372  (-0.1635, -0.1109)

function y = f3(x)
    global f1;
    global f2;
    if x > 3.56
        y = f2(x);
    else
        y = f1(x);
    end
end

figure
plot(f2, distTranslate(:, 1), distTranslate(:, 2))
figure
plot(f1, distTranslate(:, 1), distTranslate(:, 2))
figure
plot(fp, distTranslate(:, 1), distTranslate(:, 2))

correctedDist = zeros(14, 1);
for i = 1:numel(distTranslate(:, 1))
  correctedDist(i) = 1/finv(distTranslate(i));
end

figure
plot(distTranslate(:, 1), correctedDist(:, 1), "blue", distTranslate(:, 1), distTranslate(:, 2), "red")

%%

Ts = 0.01;
DataCommand = load("measurement\dataComA.mat").DataCommands;
Dist = load("measurement\dataInA.mat").Data(:, 1);
MotorSpeed = load("measurement\dataInA.mat").Data(:, 2);

time = 0:Ts:(length(Dist) - 1)*Ts;

figure
plot(time, DataCommand(:, 1), time, Dist(:, 1), time, MotorSpeed(:, 1));

figure
plot(time, DataCommand(:, 1), time, f2(Dist(:, 1)), time, MotorSpeed(:, 1));

figure
plot(time, DataCommand(:, 1), time, f1(Dist(:, 1)), time, MotorSpeed(:, 1));

correctedDist = zeros(length(Dist), 1);
for i = 1:numel(Dist)
  correctedDist(i) = f3(Dist(i));
end

figure
plot(time, DataCommand(:, 1), time, correctedDist(:, 1), time, MotorSpeed(:, 1));


figure
plot(time, DataCommand(:, 1), time, fp(Dist(:, 1)), time, MotorSpeed(:, 1));

figure
plot(time, DataCommand(:, 1), time, finv(Dist(:, 1)).^-1, time, MotorSpeed(:, 1));