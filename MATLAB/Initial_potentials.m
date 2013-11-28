% Initial potentials
V = zero; 
dV = zero;
for i=1:3
    A{i} = zero;
    dA{i} = zero;
end

TMAX = 1;
clear M; 
n = 0;
t = 0.0;

t_hist = []; 
e_err_hist = []; 
p_err_hist = []; 
mde_hist = {}; 
mdp_hist = {}; 
spe_hist = {}; 
spp_hist = {};

while t < TMAX
    n = n + 1;
    t = t + dt;
    Vex = 5 * sin(pi*x).^4 .* sin(pi*y).^4;
    Aex = {0, 0, 0};
    [psi, V, dV, A, dA] = md_step(psi, V, dV, A, dA, space, fspace, dt, epsilon, delta, Vex, Aex);
    [phi_e, phi_p] = sp_step(phi_e, phi_p, Vex, dt, space, fspace);
    abspsi = md_abs2(psi);
    psinorm = sum(abspsi(:)) * volume / N^3;
    clf;
    hold on;
    wvplot(x, y, N, {psi{1},psi{2}}, psinorm, 2, 2, 1, 'MD, electronic (psi_{1,2})');
    wvplot(x, y, N, {psi{3},psi{4}}, psinorm, 2, 2, 2, 'MD, positronic (psi_{3,4})');
    wvplot(x, y, N, {phi_e}, psinorm, 2, 2, 3, 'SP, electronic (phi_e)');
    wvplot(x, y, N, {phi_p}, psinorm, 2, 2, 4, 'SP, positronic (phi_p)');
    drawnow;
end