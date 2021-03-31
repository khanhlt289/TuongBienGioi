clear all; clc;
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Muc dich: Bao cao Seminar Khoa hoc tai truong Dai hoc KHTN - DHQGHN

% De tai: Tinh toan noi luc ket cau tuong rong ruot theo phuong phap PTHH
% voi tai trong gio tinh, giai he phuong trinh dai so tuyen tinh A * x = b
% bang xu ly tien dieu kien PCG theo phan tich ma tran duong cheo kk

% Nguoi thuc hien: Luu Truong Khanh
% Huong dan khoa hoc phan Dai so tuyen tinh: TS. Nguyen Trung Hieu

format long % Dinh dang so
nel=84;   % So luong phan tu cua he
nnel=2;     % So nut cua 1 phan tu
ndof=6;     % So chuyen vi nut dau phan tu
nnode=49;   % So nut cua he
sdof=nnode*ndof;    % So chuyen vi nut cua he
lx = 1.0; % 
lz = 1.0; % 
L = 1.0; % 
% Kich thuoc dam 
bb = 0.2; % Chieu rong dam la 30 cm.
hh = 0.6; % Chieu cao dam la 60 cm.
E = 2.65*10^6;   % Mo dun dan hoi vat lieu T/m2
Iz =  bb*hh^3/12;
Iy = hh*bb^3/12;
Ix = Iz*0.75;
G = E/2; % Mo dun dan hoi truot vat lieu T/m2       
Area = hh*bb;    % Dien tich tiet dien dam
n = 36; % So luong o san la 36 o (dang don gian)

for i = 1:7
    x(i) = (i - 1) * 1.0;
    z(i) = 0;
    y(i) = 0;
end

for i = 8:14
    x(i) = (i - 1) * 1.0;
    z(i) = 1.0;
    y(i) = 0;
end

for i = 15:21
    x(i) = (i - 1) * 1.0;
    z(i) = 2.0;
    y(i) = 0;
end

for i = 22:28
    x(i) = (i - 1) * 1.0;
    z(i) = 3.0;
    y(i) = 0;
end

for i = 29:35
    x(i) = (i - 1) * 1.0;
    z(i) = 4.0;
    y(i) = 0;
end

for i = 36:42
    x(i) = (i - 1) * 1.0;
    z(i) = 5.0;
    y(i) = 0;
end

for i = 43:49
    x(i) = (i - 1) * 1.0;
    z(i) = 6.0;
    y(i) = 0;
end

for i = 1:6
    node(i,1) = i;
    node(i,2) = i + 1;
end
for i = 7:13
    node(i,1) = i - 6;
    node(i,2) = i + 1;
end
for i = 14:19
    node(i,1) = i - 6;
    node(i,2) = i - 5;
end
for i = 20:26
    node(i,1) = i - 12;
    node(i,2) = i - 5;
end
for i = 27:32
    node(i,1) = i - 12;
    node(i,2) = i - 11;
end
for i = 33:39
    node(i,1) = i - 18;
    node(i,2) = i - 11;
end
for i = 40:45
    node(i,1) = i - 18;
    node(i,2) = i - 17;
end
for i = 46:52
    node(i,1) = i - 24;
    node(i,2) = i - 17;
end
for i = 53:58
    node(i,1) = i - 24;
    node(i,2) = i - 23;
end
for i = 59:65
    node(i,1) = i - 30;
    node(i,2) = i - 23;
end
for i = 66:71
    node(i,1) = i - 30;
    node(i,2) = i - 29;
end
for i = 72:78
    node(i,1) = i - 36;
    node(i,2) = i - 29;
end
for i = 79:84
    node(i,1) = i - 36;
    node(i,2) = i - 35;
end

nbc=12;  % So luong chuyen vi nut bi rang buoc
% Nut 1. Chuyen vi nut bi rang buoc la q1, q2, q3, q4, q5, q6
nq1 = 1;
nq2 = 2;
nq3 = 3;
nq4 = 4;
nq5 = 5;
nq6 = 6;
i = 1;
bcdof(1)=(i-1)*6 + nq1; bcval(1)=0;
bcdof(2)=(i-1)*6 + nq2; bcval(2)=0;
bcdof(3)=(i-1)*6 + nq3; bcval(3)=0;
bcdof(4)=(i-1)*6 + nq4; bcval(4)=0;
bcdof(5)=(i-1)*6 + nq5; bcval(5)=0;
bcdof(6)=(i-1)*6 + nq6; bcval(6)=0;
% Nut 8
i = 8;
bcdof(7)=(i-1)*6 + nq1; bcval(7)=0;
bcdof(8)=(i-1)*6 + nq2; bcval(8)=0;
bcdof(9)=(i-1)*6 + nq3; bcval(9)=0;
bcdof(10)=(i-1)*6 + nq4; bcval(10)=0;
bcdof(11)=(i-1)*6 + nq5; bcval(11)=0;
bcdof(12)=(i-1)*6 + nq6; bcval(12)=0;
% Nut 15
i = 15;
bcdof(13)=(i-1)*6 + nq1; bcval(13)=0;
bcdof(14)=(i-1)*6 + nq2; bcval(14)=0;
bcdof(15)=(i-1)*6 + nq3; bcval(15)=0;
bcdof(16)=(i-1)*6 + nq4; bcval(16)=0;
bcdof(17)=(i-1)*6 + nq5; bcval(17)=0;
bcdof(18)=(i-1)*6 + nq6; bcval(18)=0;
% Nut 22
i = 22;
bcdof(19)=(i-1)*6 + nq1; bcval(19)=0;
bcdof(20)=(i-1)*6 + nq2; bcval(20)=0;
bcdof(21)=(i-1)*6 + nq3; bcval(21)=0;
bcdof(22)=(i-1)*6 + nq4; bcval(22)=0;
bcdof(23)=(i-1)*6 + nq5; bcval(23)=0;
bcdof(24)=(i-1)*6 + nq6; bcval(24)=0;
% Nut 29
i = 29;
bcdof(25)=(i-1)*6 + nq1; bcval(25)=0;
bcdof(26)=(i-1)*6 + nq2; bcval(26)=0;
bcdof(27)=(i-1)*6 + nq3; bcval(27)=0;
bcdof(28)=(i-1)*6 + nq4; bcval(28)=0;
bcdof(29)=(i-1)*6 + nq5; bcval(29)=0;
bcdof(30)=(i-1)*6 + nq6; bcval(30)=0;
% Nut 36
i = 36;
bcdof(31)=(i-1)*6 + nq1; bcval(31)=0;
bcdof(32)=(i-1)*6 + nq2; bcval(32)=0;
bcdof(33)=(i-1)*6 + nq3; bcval(33)=0;
bcdof(34)=(i-1)*6 + nq4; bcval(34)=0;
bcdof(35)=(i-1)*6 + nq5; bcval(35)=0;
bcdof(36)=(i-1)*6 + nq6; bcval(36)=0;
% Nut 43
i = 43;
bcdof(37)=(i-1)*6 + nq1; bcval(37)=0;
bcdof(38)=(i-1)*6 + nq2; bcval(38)=0;
bcdof(39)=(i-1)*6 + nq3; bcval(39)=0;
bcdof(40)=(i-1)*6 + nq4; bcval(40)=0;
bcdof(41)=(i-1)*6 + nq5; bcval(41)=0;
bcdof(42)=(i-1)*6 + nq6; bcval(42)=0;

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Bac tu do, tong so nut, tong so bac tu do, tong so doan thanh (phan tu)
% q1 - Chuyen vi doc truc tai nut 1 (u_1). *
% q2 - Chuyen vi doc truc tai nut 2 (u_2).
% q3 - Chuyen vi ngang theo phuong y tai nut 1 (v_1). *
% q4 - Chuyen vi ngang theo phuong z tai nut 2 (v_2).
% q5 - Chuyen vi ngang theo phuong z tai nut 1 (w_1).
% q6 - Chuyen vi ngang theo phuong y tai nut 2 (w_2). *
% q7 - Goc xoan tai nut 1 (teta_1x).
% q8 - Goc xoan tai nut 2 (teta_2x).
% q9 - Goc xoay quan truc y tai dau 1 (teta_1y).
% q10 - Goc xoay quan truc y tai dau 2 (teta_2y).
% q11 - Goc xoay quan truc z tai dau 1 (teta_1z).
% q12 - Goc xoay quan truc z tai dau 2 (teta_2z).

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Thanh phan noi luc tai 2 nut cua (phan tu) doan thanh.
% N1 - luc doc tai nut 1 (N_1).
% N2 - luc doc tai nut 2 (N_2).
% Q1y - luc cat tren phuong y tai nut 1.
% Q2y - luc cat tren phuong y tai nut 2.
% Q1z - luc cat tren phuong z tai nut 1.
% Q2z - luc cat tren phuong z tai nut 2.
% M1x - mo men xoan tai nut 1.
% M2x - mo men xoan tai nut 2.
% M1y - mo men quanh truc y tai nut 1.
% M2y - mo men quanh truc y tai nut 2.
% M1z - mo men quanh truc z tai nut 1.
% M2z - mo men quanh truc z tai nut 2.

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Xay dung ma tran do cung cua thanh hai dau han trong khong gian 3 chieu
% Ma tran do cung dia phuong (cuc bo) co ma tran 12 x 12
S1 = E*Area/L;
S2 = 12*E*Iz/L^3;
S3 = 6*E*Iz/L^2;
S4 = 12*E*Iy/L^3;
S5 = 6*E*Iy/L^2;
S6 = G*Ix/L;
S7 = 4*E*Iy/L;
S8 = 2*E*Iy/L;
S9 = 4*E*Iz/L;
S10 = 2*E*Iz/L;

k = [S1, 0, 0, 0, 0, 0, -S1, 0, 0, 0, 0, 0; ... 
           0, S2, 0, 0, 0, S3, 0, -S2, 0, 0, 0, S3; ... 
           0, 0, S4, 0, -S5, 0, 0, 0, -S4, 0, -S5, 0; ... 
           0, 0, 0, S6, 0, 0, 0, 0, 0, -S6, 0, 0; ... 
           0, 0, -S5, 0, S7, 0, 0, 0, S5, 0, S8, 0; ... 
           0, S3, 0, 0, 0, S9, 0, -S3, 0, 0, 0, S10; ... 
           -S1, 0, 0, 0, 0, 0, S1, 0, 0, 0, 0, 0; ... 
           0, -S2, 0, 0, 0, -S3, 0, S2, 0, 0, 0, -S3; ... 
           0, 0, -S4, 0, S5, 0, 0, 0, S4, 0, S5, 0; ... 
           0, 0, 0, -S6, 0, 0, 0, 0, 0, S6, 0, 0; ... 
           0, 0, -S5, 0, S8, 0, 0, 0, S5, 0, S7, 0; ... 
           0, S3, 0, 0, 0, S10, 0, -S3, 0, 0, 0, S9]; 
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Ghep noi ma tran do cung dia phuong thanh ma tran do cung tong the

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Gan gia tri 0
ff=zeros(sdof,1);
kk=zeros(sdof,sdof);
rf=zeros(sdof,1);
index=zeros(nnel*ndof,1);
noiluc=zeros(12,nel);
cvnutpt=zeros(12,1);
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Tai trong quy nut dam
% Gia su tai trong dinh tuong co wmax = 100 kg/m2
wmax = 100 * 1e-3;

w1 = wmax;
W1 = wmax * L * L/2;

w2 = wmax - wmax / 6;
W2 = w2 * L * L;

w3 = wmax - wmax / 6 * 2;
W3 = w3 * L * L;

w4 = wmax - wmax / 6 * 3;
W4 = w4 * L * L;

w5 = wmax - wmax / 6 * 4;
W5 = w5 * L * L;

w6 = wmax - wmax / 6 * 5;
W6 = w6 * L * L * 3/2;

w7 = 0;
W7 = 0;

% Gan tai trong nut cho cac nut cung do cao
for i = 7:7:49
    ff((i-1)*ndof+nq2) = W1;
end
for i = 6:7:48
    ff((i-1)*ndof+nq2) = W2;
end
for i = 5:7:47
    ff((i-1)*ndof+nq2) = W3;
end
for i = 4:7:46
    ff((i-1)*ndof+nq2) = W4;
end
for i = 3:7:45
    ff((i-1)*ndof+nq2) = W5;
end
for i = 2:7:44
    ff((i-1)*ndof+nq2) = W6;
end
for i = 1:7:43
    ff((i-1)*ndof+nq2) = W7;
end

% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
%Xac dinh ma tran do cung cua he
for iel=1:nel
    index=mahoacvinut(iel,node,nnel,ndof);
    
    nd(1)=node(iel,1);
    nd(2)=node(iel,2);
    
    x1=x(nd(1));z1=z(nd(1));y1=y(nd(1));
    x2=x(nd(2));z2=z(nd(2));y2=y(nd(2));
    
    leng=sqrt((x2-x1)^2+(y2-y1)^2+(z2-z1)^2);
    if (x2-x1)==0
        beta=pi/2;
    else
        beta=atan((y2-y1)/(x2-x1));
    end

%     L = ((x2-x1)^2 + (y2-y1)^2 + (z2-z1)^2)^0.5

    % Chuyen ma tran tu he toa do dia phong ve he toa do tong the

    Cx = (x2 - x1)/L;
    Cy = (y2 - y1)/L;
    Cz = (z2 - z1)/L;
    alpha = 0; beta = 0;
    r=mtctoadoxop(beta, alpha, Cx, Cy, Cz);
    [k] = mtdcdam(E,Area,Iz,Iy,Ix,L);
%     [k]=mtdcdamUBoot(el,area,xi,b,ko,leng);
    k=r'*k*r;
    [kk]=mtdche(kk,k,index);
%     [kk]=mtdche(kk,k,index);
end
% Khu dieu kien bien
for i = 1:sdof
    ff(i)=ff(i)+rf(i);
end
[kk,ff]=khudkbien(kk,ff,bcdof,bcval);
% [kk,ff]=khudkbien(kk,ff,bcdof,bcval)

%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Giai he phuong trinh dai so tuyen tinh
% kk * xop = ff
% kk la ma tran do cung co so dieu kien 3.0e+09.
% ff la vec to tai quy nut cua gio tinh co wmax = 100 kg/m2.
% xop la vec to chuyen vi tai nut - don vi m.
% bang phuong phap tien xu ly dieu kien PCG duong cheo cua kk
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

format long;
res_tol = 1e-10;
max_iter = 3000;

A = kk;
b = ff;
n = size(A,1); xop = zeros(n,1);

A = sparse(A);

[a_val, a_row_ptr, a_col_idx] = sparse2csr(kk,1);
m_diagonal_val = diag(kk);
b = ff;
[xop, converged, iter_cnt, res_norm] = ...
    PCG_Diagonal(a_val, a_row_ptr, a_col_idx, b, res_tol, max_iter, m_diagonal_val);

disp = xop;
for iel=1:nel
    nd(1)=node(iel,1);
    nd(2)=node(iel,2);
    
    x1=x(nd(1));z1=z(nd(1));y1=y(nd(1));
    x2=x(nd(2));z2=z(nd(2));y2=y(nd(2));
    k = mtdcdam(E,Area,Iz,Iy,Ix,L);

    cvinutpt(1,1)=disp(6*nd(1)-5,1);
    cvinutpt(2,1)=disp(6*nd(1)-4,1);
    cvinutpt(3,1)=disp(6*nd(1)-3,1);
    cvinutpt(4,1)=disp(6*nd(1)-2,1);
    cvinutpt(5,1)=disp(6*nd(1)-1,1);
    cvinutpt(6,1)=disp(6*nd(1),1);
    cvinutpt(7,1)=disp(6*nd(2)-5,1);
    cvinutpt(8,1)=disp(6*nd(2)-4,1);
    cvinutpt(9,1)=disp(6*nd(2)-3,1);
    cvinutpt(10,1)=disp(6*nd(2)-2,1);
    cvinutpt(11,1)=disp(6*nd(2)-1,1);
    cvinutpt(12,1)=disp(6*nd(2),1);
    cvinutpt=r*cvinutpt; % chuyen ve Hetoado diaphuong
    noiluc(:,iel)=k*cvinutpt;
end
