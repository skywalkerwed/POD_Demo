function u_t = myTranseq(t, u, dummy, a, dx, A)
u_t = a^2/dx^2*A*u;
end