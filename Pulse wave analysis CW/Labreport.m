% Toby Thomas, 2018
clear;close all;
%% Load data
foot = load('foot_to_foot.mat');
single = load('single_point.mat');

figure();
subplot(2,2,1); plot(single.p1);title('Single Pressure 1');xlabel('Index');ylabel('Pressure (mm/hg)');
subplot(2,2,2); plot(single.Q);title('Single Velocity');xlabel('Index');ylabel('Velocity (m/min)');
subplot(2,2,3); plot(foot.p1);title('Foot-Foot Pressure 1');xlabel('Index');ylabel('Pressure (mm/hg)');
subplot(2,2,4); plot(foot.p2);title('Foot-Foot Pressure 2');xlabel('Index');ylabel('Pressure (mm/hg)');

% Systole (4th):
%   t = 2899:3191 --> 8.957:9.2490 seconds
%   t = 2897:3191 --> 8.9550:9.249 seconds

% Define window of 4th pulse
win_start = 2899;
win_stop = 3899;
s_pulse_p_4 = single.p1(win_start:win_stop);
s_pulse_q_4 = single.Q(win_start:win_stop);
f_pulse_p1_4 = foot.p1(win_start:win_stop);
f_pulse_p2_4 = foot.p2(win_start:win_stop);
s_pulse_time_4 = single.t(win_start:win_stop);
f_pulse_time_4 = foot.t(win_start:win_stop);

figure(); 
subplot(2,2,1);plot(s_pulse_time_4,single.p1(win_start : win_stop));title('Single pressure (4th)');xlabel('Time (Seconds)');ylabel('Pressure (mm/hg)');
subplot(2,2,2);plot(s_pulse_time_4,single.Q(win_start : win_stop));title('Single velocity (4th)');xlabel('Time (Seconds)');ylabel('Velocity (m/min)');
subplot(2,2,3);plot(f_pulse_time_4,foot.p1(win_start : win_stop));title('Foot pressure 1(4th)');xlabel('Time (Seconds)');ylabel('Pressure (mm/hg)');
subplot(2,2,4);plot(f_pulse_time_4,foot.p2(win_start : win_stop));title('Foot pressure 2(4th)');xlabel('Time (Seconds)');ylabel('Pressure (mm/hg)');

%% (i) Local pulse wave velocity
% c = 1/r * dP/dU
r = 1000;
% find peak pressure
peak = 3191;
trough = 3146;

% find dP and dQ
dP = single.p1(peak) - single.p1(peak-1);
dQ = single.Q(peak) - single.Q(peak-1);
dP = dP*133.322; % Pa
dQ = dQ/60000; 

% Find area of Lumen
dia = 0.02;
th = 0.001;
Lumen_rad = (dia-2*th)/2;
Lumen_area = pi*Lumen_rad^2;

% find dU
dU = dQ./Lumen_area;
c = (1/r)*dP/dU;
fprintf('c at pressure peak =           %f m/s\n', c);

for i = 1:(peak-trough)
    dP2(i) = single.p1(i+1) - single.p1(i);
    dQ2(i) = single.Q(i+1) - single.Q(i);
end

dP2 = dP2 .*133.322; % Pa
dQ2 = dQ2 ./60000; 
dU2 = dQ2 ./Lumen_area;

c2 = ((1/r).*dP2./dU2);
c2_avg = mean(c2);
fprintf('Average systolic c =           %f m/s \n', c2_avg);


%% ii) Foot to Foot pulse velocity c_ff

% eq 11: c_ff = (Dx) / (Dt)
Dx = 0.25; %m

% find fifference in peak times
m1 = max(foot.p1(win_start : win_stop));
m2 = max(foot.p2(win_start : win_stop));

peak_1 = find(foot.p1(win_start:win_stop) == m1);
peak_2 = find(foot.p2(win_start:win_stop) == m2);

Dt = abs(foot.t(peak_2)-foot.t(peak_1));

cff = Dx / Dt;

fprintf('C_ff =                         %f \n', cff);


%% iii) Seporate P(t) and U(t) into P_f(t), P_b(t), U_f(t) and U_b(t)
% use eq 6 to 9

P = single.p1(win_start:win_start+1000) .*133.322;
U = single.Q(win_start: win_start+1000) ./60 ./(Lumen_area*r);

for i  = 1:length(P)-1
    dP(i) = P(i+1) - P(i);
    dU(i) = U(i+1) - U(i);
end

for t = 1:length(dP)
    dP_f(t) = 0.5 .* ( dP(t) + r.*c.*dU(t) ); 
    dP_b(t) = 0.5 .* ( dP(t) - r.*c.*dU(t) );
    dU_f(t) = 0.5 * ( dU(t) + ( dP(t)./(r.*c) ) );
    dU_b(t) = 0.5 * ( dU(t) - ( dP(t)./(r.*c) ) );
end

for t = 1:length(dP)
    P_f(t) = sum(dP_f) + P(t)/2 ;
    P_b(t) = sum(dP_b) + P(t)/2 ;
    U_f(t) = sum(dU_f) + U(t)/2 ;
    U_b(t) = sum(dU_b) + U(t)/2 ;
end

time = single.t(win_start:win_stop);

figure(); 
plot(time,P); hold on
plot(time(1:end-1),P_f);
plot(time(1:end-1),P_b);
legend('Relative pressure','Forward component', 'Backward component');
xlabel('time (Seconds)'); ylabel('pressure (Pa)'); title('Pressure flow');
hold off

figure(); 
plot(time,U); hold on
plot(time(1:end-1),U_f);
plot(time(1:end-1),U_b);
legend('velocity','Forward component', 'Backward component');
xlabel('time (Seconds)'); ylabel('velocity (m/s)'); title('Velocity flow');
hold off


%% Describe results

% i) The local pulse velocity varys through out systole due the the noise
% in the data provded. This is shown in the figure below showing Q vs time
% and c vs time

figure(); 
subplot(1,2,1);plot(single.t(trough:peak),single.Q(trough:peak));
title('Flow rate velocity - systole');
xlabel('Time (seconds)'); ylabel('m^3 /min');

subplot(1,2,2); plot(single.t(trough:peak-1), c2);
title('local pulse velocity (c)- systole')
xlabel('Time (seconds)'); ylabel('Local pulse velocity m/s')

% iii)  The backward component of pressure was consistantly very slightly higher than
% the forward comonent during systole. 


