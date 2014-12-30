c_ = sparse(23,1);
c_(1) = 1;
b_ = sparse(17,1);
b_(3:3) = a;
b_(6:6) = b;
b_(12:12) = 1.0;
b_(13:13) = -0.5;
b_(14:14) = -0.5;
A_ = sparse(17, 23);
A_(1:1, 3:3) = 1.0; A_(1:1, 2:2) = 1.0; A_(1:1, 14:14) = -1.0;
A_(2:2, 15:15) = -1.0; A_(2:2, 3:3) = 1.0; A_(2:2, 4:4) = -1.0;
A_(3:3, 14:14) = 1.0;
A_(4:4, 6:6) = 1.0; A_(4:4, 5:5) = 1.0; A_(4:4, 16:16) = -1.0;
A_(5:5, 17:17) = 1.0; A_(5:5, 6:6) = 1.0; A_(5:5, 7:7) = -1.0;
A_(6:6, 16:16) = 1.0;
A_(7:7, 18:18) = 1.0; A_(7:7, 19:19) = 1.0; A_(7:7, 1:1) = -1.0;
A_(8:8, 20:20) = 1.0; A_(8:8, 17:17) = 1.0; A_(8:8, 18:18) = -1.0;
A_(9:9, 15:15) = -1.0; A_(9:9, 20:20) = -1.0;
A_(10:10, 11:11) = 0.5; A_(10:10, 19:19) = 0.5; A_(10:10, 8:8) = -1.0;
A_(11:11, 11:11) = 0.5; A_(11:11, 19:19) = -0.5; A_(11:11, 9:9) = -1.0;
A_(12:12, 10:10) = 1.0;
A_(13:13, 21:21) = 0.5; A_(13:13, 12:12) = -1.0;
A_(14:14, 21:21) = -0.5; A_(14:14, 13:13) = -1.0;
A_(15:15, 22:22) = 1.0; A_(15:15, 23:23) = -1.0; A_(15:15, 21:21) = -1.0;
A_(16:16, 15:15) = a; A_(16:16, 22:22) = -1.0;
A_(17:17, 17:17) = b; A_(17:17, 23:23) = -1.0;
G_ = sparse(13, 23);
G_(1:1:1, 2:2) = -speye(1, 1);
G_(2:1:2, 3:3) = -speye(1, 1);
G_(3:1:3, 4:4) = -speye(1, 1);
G_(4:1:4, 5:5) = -speye(1, 1);
G_(5:1:5, 6:6) = -speye(1, 1);
G_(6:1:6, 7:7) = -speye(1, 1);
G_(7:1:7, 11:11) = -speye(1, 1);
G_(8:1:8, 8:8) = -speye(1, 1);
G_(9:1:9, 9:9) = -speye(1, 1);
G_(10:1:10, 10:10) = -speye(1, 1);
G_(11:1:11, 12:12) = -speye(1, 1);
G_(12:1:12, 13:13) = -speye(1, 1);
G_(13:1:13, 11:11) = -speye(1, 1);
h_ = zeros(13, 1);
dims.q = [3,3];
dims.l = 7;
[x_codegen, y_, info_] = conelp(full(c_), G_, h_, dims, A_, full(b_));
t2 = x_codegen(1:1);
t9 = x_codegen(2:2);
t8 = x_codegen(3:3);
t8z0 = x_codegen(4:4);
t11 = x_codegen(5:5);
t10 = x_codegen(6:6);
t10z0 = x_codegen(7:7);
t7z0 = x_codegen(8:8);
t7z1 = x_codegen(9:9);
t7z2 = x_codegen(10:10);
t6 = x_codegen(11:11);
t6z0 = x_codegen(12:12);
t6z1 = x_codegen(13:13);
pa = x_codegen(14:14);
x1 = x_codegen(15:15);
pb = x_codegen(16:16);
x2 = x_codegen(17:17);
t1 = x_codegen(18:18);
t7 = x_codegen(19:19);
t0 = x_codegen(20:20);
t4 = x_codegen(21:21);
t3 = x_codegen(22:22);
t5 = x_codegen(23:23);
ecos_optval = 1*info_.pcost;
