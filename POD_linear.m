clear;clc;close all
A = [1 0 0;-1 5/2 -1; -1 -1 3];
t = 0:0.02:5;
u = exp(-t);
usol = inv(A)*[u; zeros(1, length(t)); zeros(1, length(t))];
u1 = usol(1, :);
u2 = usol(2, :);
u3 = usol(3, :);

%% case1, take snapshots at T=[1 20]
T = [1 20];
u1T = u1(T);
u2T = u2(T);
u3T = u3(T);

W = [u1T;u2T;u3T];
[U, S, V]=svd(W); 
Eng = diag(S).^2/norm(diag(S),2)^2;
N = find(Eng < 1e-4);N=N(1)-1;%cut-of engery lower than 0.01%
phi = U(:, 1:N)
asol = inv(phi.'*A*phi)*(phi.'*[u; zeros(1, length(t)); zeros(1, length(t))]);
u_pod = phi*asol;
u_pod1 = u_pod(1, :);
u_pod2 = u_pod(2, :);
u_pod3 = u_pod(3, :);

figure(1)

subplot(2,2,1)
plot(t, u1, 'linewidth', [2]);
hold on
plot(t, u2, 'linewidth', [2]);
hold on
plot(t, u3, 'linewidth', [2]);
hold on
plot(t(T), u1(T), 'ro', 'linewidth', [1])
hold on
plot(t(T), u2(T), 'ro', 'linewidth', [1])
hold on
plot(t(T), u3(T), 'ro', 'linewidth', [1])
title("full order")
subplot(2,2,2)
plot(t, u_pod1, 'linewidth', [2]);
hold on
plot(t, u_pod2, 'linewidth', [2]);
hold on
plot(t, u_pod3, 'linewidth', [2]);
title("reduced order")
subplot(2, 2, [3 4])
e1 = abs(u1-u_pod1);
e2 = abs(u2-u_pod2);
e3 = abs(u3-u_pod3);
plot(t,e1, 'linewidth', [2]);
hold on
plot(t, e2, 'linewidth', [2]);
hold on
plot(t, e3, 'linewidth', [2]);
title("error")
suptitle("snapshots at T=[1 20], u=e\^-t")
%% case2, take snapshots at T=[20 200];
T = [20 200];
u1T = u1(T);
u2T = u2(T);
u3T = u3(T);

W = [u1T;u2T;u3T];
[U, S, V]=svd(W);
Eng = diag(S).^2/norm(diag(S),2)^2;
N = find(Eng < 1e-4);N=N(1)-1;%cut-of engery lower than 0.01%
phi = U(:, 1:N)
asol = inv(phi.'*A*phi)*(phi.'*[u; zeros(1, length(t)); zeros(1, length(t))]);
u_pod = phi*asol;
u_pod1 = u_pod(1, :);
u_pod2 = u_pod(2, :);
u_pod3 = u_pod(3, :);

figure(2)

subplot(2,2,1)
plot(t, u1, 'linewidth', [2]);
hold on
plot(t, u2, 'linewidth', [2]);
hold on
plot(t, u3, 'linewidth', [2]);
hold on
plot(t(T), u1(T), 'ro', 'linewidth', [1])
hold on
plot(t(T), u2(T), 'ro', 'linewidth', [1])
hold on
plot(t(T), u3(T), 'ro', 'linewidth', [1])
title("full order")
subplot(2,2,2)
plot(t, u_pod1, 'linewidth', [2]);
hold on
plot(t, u_pod2, 'linewidth', [2]);
hold on
plot(t, u_pod3, 'linewidth', [2]);
title("reduced order")
subplot(2, 2, [3 4])
e1 = abs(u1-u_pod1);
e2 = abs(u2-u_pod2);
e3 = abs(u3-u_pod3);
plot(t,e1, 'linewidth', [2]);
hold on
plot(t, e2, 'linewidth', [2]);
hold on
plot(t, e3, 'linewidth', [2]);
title("error")
suptitle("snapshots at T=[20 200], u=e\^-t")
%% case3, take u = sin(t), and snapshots at T=[1 20]
u = sin(t);

usol = inv(A)*[u; zeros(1, length(t)); zeros(1, length(t))];
u1 = usol(1, :);
u2 = usol(2, :);
u3 = usol(3, :);

T = [1 20];
u1T = u1(T);
u2T = u2(T);
u3T = u3(T);

W = [u1T;u2T;u3T];
[U, S, V]=svd(W);
Eng = diag(S).^2/norm(diag(S),2)^2;
N = find(Eng < 1e-4);N=N(1)-1;%cut-of engery lower than 0.01%
phi = U(:, 1:N)
asol = inv(phi.'*A*phi)*(phi.'*[u; zeros(1, length(t)); zeros(1, length(t))]);
u_pod = phi*asol;
u_pod1 = u_pod(1, :);
u_pod2 = u_pod(2, :);
u_pod3 = u_pod(3, :);

figure(3)

subplot(2,2,1)
plot(t, u1, 'linewidth', [2]);
hold on
plot(t, u2, 'linewidth', [2]);
hold on
plot(t, u3, 'linewidth', [2]);
hold on
plot(t(T), u1(T), 'ro', 'linewidth', [1])
hold on
plot(t(T), u2(T), 'ro', 'linewidth', [1])
hold on
plot(t(T), u3(T), 'ro', 'linewidth', [1])
title("full order")
subplot(2,2,2)
plot(t, u_pod1, 'linewidth', [2]);
hold on
plot(t, u_pod2, 'linewidth', [2]);
hold on
plot(t, u_pod3, 'linewidth', [2]);
title("reduced order")
subplot(2, 2, [3 4])
e1 = abs(u1-u_pod1);
e2 = abs(u2-u_pod2);
e3 = abs(u3-u_pod3);
plot(t,e1, 'linewidth', [2]);
hold on
plot(t, e2, 'linewidth', [2]);
hold on
plot(t, e3, 'linewidth', [2]);
title("error")
suptitle("snapshots at T=[1 20], u=sin(t)")