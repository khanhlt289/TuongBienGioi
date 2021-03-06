% Chuong trinh PCG Diagonal co ma tran tien dieu kien la duong cheo
% Chuong trinh PCG tien dieu kien tach A = CC', hoac A = C^(1/2)*C^(1/2)'
% C la phan tich Cholesky cua ma tran A
% C'*A*C*C^(-1)*x = C'*b;
% C^(1/2)'*A*C^(1/2)*C^(-1/2)*x = C^(1/2)'*b;

function [x, converged, iter_cnt, res_norm] = PCG_Diagonal(a_val, a_row_ptr, a_col_idx, b, res_tol, max_iter, m_diagonal_val)
% (Left) Preconditioned Conjugate Gradient method
	n = size(a_row_ptr, 1) - 1;
	
	if (nargin < 3)	res_tol  = 1e-9; end
	if (nargin < 4)	max_iter = 1000; end
	if (nargin < 5) M = eye(n);      end
	
	x = zeros(n, 1);
	r = b - dotproduct(a_val, a_row_ptr, a_col_idx, x);
	z = r ./m_diagonal_val;        % Left preconditioning; MATLAB does not suggest saving inv(M)
	p = z;
	rho = r' * z;
	rn_stop = norm(r, 2) * res_tol;
	
	iter_cnt = 1;
	res_norm(iter_cnt) = norm(r);
	
	converged = 0;
	while ((iter_cnt < max_iter) && (res_norm(iter_cnt) > rn_stop))
		s     = dotproduct(a_val, a_row_ptr, a_col_idx, p);
		alpha = rho / (p' * s);

		x = x + alpha * p;
		r = r - alpha * s;
		
		rho_0 = rho;
		z     = r ./ m_diagonal_val;
		rho   = r' * z;
		
		beta = rho / rho_0;
		p    = z + beta * p;
		
		iter_cnt = iter_cnt + 1;
		res_norm(iter_cnt) = norm(r, 2);
	end
	if (res_norm(iter_cnt) <= rn_stop) converged = 1; end
end