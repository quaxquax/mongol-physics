epsilon = 1.0;
delta = 0.01;
dt = 1./128;
N = 32;
[space, fspace] = md_domain([-0.5,0.5],[-0.5,0.5],[-0.5,0.5],N)
x = space{1};
y = space{2};
z = space{3};
zero = x .*0;
xyz2 = x .^2 + y .^2 + z .^2; 
volume = 1;


% Initial wavefunction
wv_e = cos(pi*(x+0.1)).^4 .* cos(pi*(y+0.1)) .^ 4; wv_p = cos(pi*(x-0.1)).^4 .* cos(pi*(y-0.1)) .^ 4;
psi{1} = wv_e; psi{2} = zero; psi{3} = wv_p; psi{4} = zero;
phi_e = wv_e; phi_p = wv_p;