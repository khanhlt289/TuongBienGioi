function y_val = dotproduct(a_val, a_row_ptr, a_col_idx, x_val)
n = size(a_row_ptr,1) - 1; iter_cnt = 0;
y_val = zeros(n,1);
    for i = 1:n
        j1 = a_row_ptr(i);
        j2 = a_row_ptr(i+1) - 1;
        for j = j1:j2
            k = a_col_idx(j);
            y_val(i) = y_val(i) + a_val(j) * x_val(k);
            iter_cnt = iter_cnt + 1;
        end
    end
% end