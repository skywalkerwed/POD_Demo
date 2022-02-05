clear;clc;close all

R=5e4;
C=1e-5;
dt= 1e-5;

t=0:dt:10;
u=exp(-0.1*t);
u1 = [0];
u2 = [1/2];
u3 = [1];
%numeric
for i = 2:length(t)
    u3(i) = u(i);
    u2(i) = 1/2*(u3(i)+u1(i-1));
    u1(i) = dt/R/C*(u2(i-1)-u1(i-1))+u1(i-1);
end
%excat
u1t = 10/9*exp(-0.1*t)-10/9*exp(-t);
u2t = u1t+5e-1*(-1/9*exp(-0.1*t)+10/9*exp(-t));
u3t = u;

A = [C 0 0; 0 0 0; 0 0 0];
B = [1/R -1/R 0; -1/R 2/R -1/R; 0 0 1];
f = [0 0 1]';
%% case1, take snapshots at T=[1 1/5*lenth(t) 4/5*length(t)]
T = [1 ceil(1/5*length(t)) ceil(4/5*length(t))];
u1T = u1(T);
u2T = u2(T);
u3T = u3(T);

W=[u1T;u2T;u3T];
[U, S, V]=svd(W);
Eng = diag(S).^2/norm(diag(S),2)^2;
N = find(Eng < 1e-4);N=N(1)-1;%cut-of engery lower than 0.01%
phi = U(:, 1:N);
A_r = phi'*A*phi;
B_r = phi'*B*phi;
f_r = phi'*f;

%In the following, we exactly know there are two modes.
%If not, the diff. eq.s should be solved in other ways.
u_r0 = pinv(phi)*[0 1/2 1]';%ini. cond.
u_r1 = u_r0(1);
u_r2 = u_r0(2);
A_r_inv = inv(A_r+B_r*dt);
for i = 2:length(t)
    u_rsol = A_r_inv*[dt*f_r(1)*u(i)+A_r(1,1)*u_r1(i-1)+A_r(1,2)*u_r2(i-1);
        dt*f_r(2)*u(i)+A_r(2,1)*u_r1(i-1)+A_r(2,2)*u_r2(i-1)];
    u_r1(i) = u_rsol(1);
    u_r2(i) = u_rsol(2);
end
u_pod = phi*[u_r1;u_r2];
u_pod1 = u_pod(1, :);
u_pod2 = u_pod(2, :);
u_pod3 = u_pod(3, :);
figure(1)
subplot(2,2,1)
plot(t, u1t, 'linewidth', [2]);
hold on
plot(t, u2t, 'linewidth', [2]);
hold on
plot(t, u3t, 'linewidth', [2]);
hold on
plot(t(T), u1t(T), 'ro', 'linewidth', [1]);
hold on
plot(t(T), u2t(T), 'ro', 'linewidth', [1]);
hold on
plot(t(T), u3t(T), 'ro', 'linewidth', [1]);
title('exact')

subplot(2,2,2)
plot(t, u_pod1, 'linewidth', [2]);
hold on
plot(t, u_pod2, 'linewidth', [2]);
hold on
plot(t, u_pod3, 'linewidth', [2])
title('POD')

e1 = abs(u_pod1-u1t);
e2 = abs(u_pod2-u2t);
e3 = abs(u_pod3-u3t);

subplot(2,2,[3 4])
plot(t, e1, 'linewidth', [2]);
hold on
plot(t, e2, 'linewidth', [2]);
hold on
plot(t, e3, 'linewidth', [2]);
title('error')
suptitle('snapshots at T=[0 1/5 4/5], u=e\^-t')

%for case3
phi2 = phi;
%% case2, take u = cos(t), and snapshots at T=[1 1/5*length(t) 4/5*length(t)]
u=cos(t);

u1 = [0];
u2 = [1/2];
u3 = [1];

for i = 2:length(t)
    u3(i) = u(i);
    u2(i) = 1/2*(u3(i)+u1(i-1));
    u1(i) = dt/R/C*(u2(i-1)-u1(i-1))+u1(i-1);
end

u1t = 0.5*(sin(t)+cos(t)-exp(-t));
u2t = u1t+0.25*(cos(t)-sin(t)+exp(-t));
u3t = u;

T = [1 ceil(1/5*length(t)) ceil(4/5*length(t))];
u1T = u1(T);
u2T = u2(T);
u3T = u3(T);

W=[u1T;u2T;u3T];
[U, S, V]=svd(W);
Eng = diag(S).^2/norm(diag(S),2)^2;
N = find(Eng < 1e-4);N=N(1)-1;%cut-of engery lower than 0.01%
phi = U(:, 1:N);
A_r = phi'*A*phi;
B_r = phi'*B*phi;
f_r = phi'*f;

u_r0 = pinv(phi)*[0 1/2 1]';%ini. cond.
u_r1 = u_r0(1);
u_r2 = u_r0(2);
A_r_inv = inv(A_r+B_r*dt);
for i = 2:length(t)
    u_rsol = A_r_inv*[dt*f_r(1)*u(i)+A_r(1,1)*u_r1(i-1)+A_r(1,2)*u_r2(i-1);
        dt*f_r(2)*u(i)+A_r(2,1)*u_r1(i-1)+A_r(2,2)*u_r2(i-1)];
    u_r1(i) = u_rsol(1);
    u_r2(i) = u_rsol(2);
end
u_pod = phi*[u_r1;u_r2];
u_pod1 = u_pod(1, :);
u_pod2 = u_pod(2, :);
u_pod3 = u_pod(3, :);

figure(2)
subplot(2,2,1)
plot(t, u1t, 'linewidth', [2]);
hold on
plot(t, u2t, 'linewidth', [2]);
hold on
plot(t, u3t, 'linewidth', [2]);
hold on
plot(t(T), u1t(T), 'ro', 'linewidth', [1]);
hold on
plot(t(T), u2t(T), 'ro', 'linewidth', [1]);
hold on
plot(t(T), u3t(T), 'ro', 'linewidth', [1]);
title('exact')

subplot(2,2,2)
plot(t, u_pod1, 'linewidth', [2]);
hold on
plot(t, u_pod2, 'linewidth', [2]);
hold on
plot(t, u_pod3, 'linewidth', [2])
title('POD')

e1 = abs(u_pod1-u1t);
e2 = abs(u_pod2-u2t);
e3 = abs(u_pod3-u3t);

subplot(2,2,[3 4])
plot(t, e1, 'linewidth', [2]);
hold on
plot(t, e2, 'linewidth', [2]);
hold on
plot(t, e3, 'linewidth', [2]);
title('error')
suptitle('snapshots at T=[0 1/5 4/5], u=cos(t)')
%% case3, take u = cos(t), use phi obtained by former simulation
A_r = phi2'*A*phi2;
B_r = phi2'*B*phi2;
f_r = phi2'*f;

u_r0 = pinv(phi2)*[0 1/2 1]';%ini. cond.
u_r1 = u_r0(1);
u_r2 = u_r0(2);
A_r_inv = inv(A_r+B_r*dt);
for i = 2:length(t)
    u_rsol = A_r_inv*[dt*f_r(1)*u(i)+A_r(1,1)*u_r1(i-1)+A_r(1,2)*u_r2(i-1);
        dt*f_r(2)*u(i)+A_r(2,1)*u_r1(i-1)+A_r(2,2)*u_r2(i-1)];
    u_r1(i) = u_rsol(1);
    u_r2(i) = u_rsol(2);
end
u_pod = phi2*[u_r1;u_r2];
u_pod1 = u_pod(1, :);
u_pod2 = u_pod(2, :);
u_pod3 = u_pod(3, :);

figure(3)
subplot(2,2,1)
plot(t, u1t, 'linewidth', [2]);
hold on
plot(t, u2t, 'linewidth', [2]);
hold on
plot(t, u3t, 'linewidth', [2]);
title('exact')

subplot(2,2,2)
plot(t, u_pod1, 'linewidth', [2]);
hold on
plot(t, u_pod2, 'linewidth', [2]);
hold on
plot(t, u_pod3, 'linewidth', [2])
title('POD')

e1 = abs(u_pod1-u1t);
e2 = abs(u_pod2-u2t);
e3 = abs(u_pod3-u3t);

subplot(2,2,[3 4])
plot(t, e1, 'linewidth', [2]);
hold on
plot(t, e2, 'linewidth', [2]);
hold on
plot(t, e3, 'linewidth', [2]);
title('error')
suptitle('u=cos(t), POD obtained by u=e\^-t')
%% compare two phi via gauss elimination
phi = phi/phi(1, 1);
phi(:, 2) = phi(:, 2)-phi(:, 1)*phi(1,2);
phi(:, 2) = phi(:, 2)/phi(2, 2);
phi(:, 1) = phi(:, 1)-phi(:, 2)*phi(2, 1)

phi2 = phi2/phi2(1, 1);
phi2(:, 2) = phi2(:, 2)-phi2(:, 1)*phi2(1,2);
phi2(:, 2) = phi2(:, 2)/phi2(2, 2);
phi2(:, 1) = phi2(:, 1)-phi2(:, 2)*phi2(2, 1)