clc
clearvars
close all
format long
global rbf_type Npoly do_scaling
global rbfw_type scaling_size
global X Xc h_y rcov
global h v rbfscale
a = 0.25; b = 0.5; c = 0.25; d = 0.5;
v = 1/250;
h = 0.03/2;
t_F = 1;
rbf_type = 'MQ';
rbfscale = 2;
Npoly = 1;
do_scaling = 2;
rbfw_type = 'Wendland32';
h_y = 5*h;
rcov = sqrt(2)*h_y;
scaling_size = rcov;
dist_type = 'Halton'; 
domain_type = 1;
dt_values = [0.001/2,0.001/4,0.001/8];
num_dt = length(dt_values);
eerp_values = zeros(num_dt, 1);
eeru_values = zeros(num_dt, 1);
eerv_values = zeros(num_dt, 1);
[Xc, Xci, Xcb] = points_sq_unity(a, b, c, d, h_y, 'R');
Points = Omega(a, b, c, d, h, dist_type, domain_type);
Xi = Points{1, 1};
Xb = Points{2, 1};
X = [Xi; Xb];
n = size(X, 1);
for k = 1:num_dt
    tic
    dt = dt_values(k);
    niT = ceil(t_F / dt);
    AG1 = PUmat('x', X, Xi);
    AG2 = PUmat('y', X, Xi);
    delAi = PUmat('L', X, Xi);
    Ai = PUmat('0', X, Xi);
    Ab = PUmat('0', X, Xb);
    Uf{1} = f(X, 0, 'u1');
    Vf{1} = f(X, 0, 'u2');
    for i = 1:niT
        MassB1m = Ai - 0.5 * dt * v * delAi + dt * psi1(Ai * Uf{i}, AG1) + dt * psi1(Ai * Vf{i}, AG2);
        MassB2m = Ab;
        Bm = [MassB1m; MassB2m];
        [Lm, Um, Pm] = lu(Bm);
        Inv_Massm = Um \ (Lm \ Pm) * eye(n);
        RhsI1m = Ai * Uf{i};
        RhsI2m = Ai * Vf{i};
        RhsB1m = f(Xb, i * dt, 'u1');
        RhsB2m = f(Xb, i * dt, 'u2');
        Rhsm1 = [RhsI1m; RhsB1m];
        Rhsm2 = [RhsI2m; RhsB2m];
        Ufm(:, i) = Inv_Massm * Rhsm1;
        Vfm(:, i) = Inv_Massm * Rhsm2;
        MassB1p = -dt * delAi;
        MassB2p = Ab;
        Bp = [MassB1p; MassB2p];
        [Lp, Up, Pp] = lu(Bp);
        Inv_Massp = Up \ (Lp \ Pp) * eye(n);
        RhsIp = -AG1 * Ufm(:, i) - AG2 * Vfm(:, i) - dt * AG1 * g(X, i * dt, 'u1') - dt * AG2 * g(X, i * dt, 'u2');
        RhsBp = p(Xb, i * dt);
        Rhsp = [RhsIp; RhsBp];
        P(:, i) = Inv_Massp * Rhsp;
        MassB1 = Ai - 0.5 * dt * v * delAi;
        MassB2 = Ab;
        B1 = [MassB1; MassB2];
        [L1, U1, P1] = lu(B1);
        Inv_Mass = U1 \ (L1 \ P1) * eye(n);
        Rhs1I = Ai * Ufm(:, i) - dt * AG1 * P(:, i) + dt * g(Xi, i * dt, 'u1');
        Rhs2I = Ai * Vfm(:, i) - dt * AG2 * P(:, i) + dt * g(Xi, i * dt, 'u2');
        Rhs1B = f(Xb, i * dt, 'u1');
        Rhs2B = f(Xb, i * dt, 'u2');
        Rhs1 = [Rhs1I; Rhs1B];
        Rhs2 = [Rhs2I; Rhs2B];
        
        Uf{i+1} = Inv_Mass * Rhs1;
        Vf{i+1} = Inv_Mass * Rhs2;
        Uf{i} = Uf{i+1};
        Vf{i} = Vf{i+1};
    end
    eerp_values(k) = norm(P(:, end) - p(X, t_F), 2) / sqrt(n)
    eeru_values(k) = norm(Uf{end} - f(X, t_F, 'u1'), 2) / sqrt(n)
    eerv_values(k) = norm(Vf{end} - f(X, t_F, 'u2'), 2) / sqrt(n)
    toc
end
