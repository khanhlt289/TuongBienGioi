function[r]=mtctoadoxop(beta, alpha, Cx, Cy, Cz)

cos_beta = cos(beta);
sin_beta = sin(beta);

cos_gamma = (Cx^2 + Cz^2);
sin_gamma = Cy;

% alpha = 0;

T_alpha = [1,0,0;...
    0,cos(alpha),sin(alpha);...
    0,-sin(alpha),cos(alpha)];

T_gamma = [cos_gamma, sin_gamma, 0;...
    -sin_gamma, cos_gamma, 0;...
    0,0,1];

T_beta = [1, 0, 0; ...
    0,1,0;...
    0, 0, 1];

zeros_T = zeros(3,3);

T = T_alpha * T_gamma * T_beta;

r = [T,   zeros_T,    zeros_T,    zeros_T; ...
    zeros_T,    T,  zeros_T,    zeros_T; ...
    zeros_T,    zeros_T,    T,  zeros_T; ...
    zeros_T,    zeros_T,    zeros_T,    T];