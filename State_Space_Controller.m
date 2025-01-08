alfa = 20 * pi/180;    %rad
m = 0.0165;         %kg
n = 0.3;            %kg/s
r_lin = 0.108;      %m
w_lin = 5;          %rad/s

A = [0 1;
     w_lin^2*cos(alfa) -n/m];

B = [0
     2*w_lin*r_lin*cos(alfa)];

C = [1 0];

D = [0];

sys = ss(A, B, C, D);

dsys = c2d(sys, 0.01, "zoh");


step(dsys);

%%

H = tf(sys);

R = [B A*B];
O = [C ; C*A];
rank(O);%=2-> observable


%%
%optimal controller using p86 controller
A_ext = [
        A(1, 1) A(1, 2) 0;
        A(2, 1) A(2, 2) 0;
        -C(1, 1) -C(1, 2) 1
        ];

B_ext = [
         B(1, 1)
         B(2, 1)
         0];

Q = [10 0 0; %r
     0 1 0; %dr/dt
     0 0 2; % N*integrator(e(t))
     ];
R = [1];

[F S P] = lqr(A_ext, B_ext, Q, R);

P
F

K = F(1:2) %rad/s.m rad/m
N = F(3) %rad/s²


%% same but using matlab func

Q = [10 0 0; %r
     0 100 0; %dr/dt
     0 0 2; % integrator(e(t))
     ];
R = [10];

[F,S,e] = lqi(dsys, Q, R);

K = -F(1:2)
N = -F(3)

% closed loop with LQI


% energy of error
% Q = [100 0 0; %r
%      0 20 0; %dr/dt
%      0 0 10; % N*integrator(e(t))
%      ];
% 
% % energy of input
% R = [2];
% 
% % cross-energy between input and error
% T = [1; 10; 10];
% 
% [F,S,e] = lqi(sys, Q, R, T);
% 
% K = -F(1:2)
% N = -F(3)

close all

Ts = 0.01;

dA = dsys.A;
dB = dsys.B;
dC = dsys.C;
dD = dsys.D;

dA_ext = [
        dA(1, 1) dA(1, 2) 0;
        dA(2, 1) dA(2, 2) 0;
        -C(1, 1) -C(1, 2) 1
        ];

dB_ext = [
         dB(1, 1)
         dB(2, 1)
         0];

dC_ext = [dC 0];

integrator = tf([N], [1 -1], Ts);

closed_loop_sys = ss(dA_ext - dB_ext*F, dB_ext, dC_ext, dD, Ts);

figure
margin(integrator*closed_loop_sys);

[Gm, Pm, Wcg, Wcp] = margin(integrator*closed_loop_sys);

figure
pzmap(integrator*closed_loop_sys)
% 
% figure
% bode(integrator*closed_loop_sys);