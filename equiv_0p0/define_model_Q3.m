function [call_model,func_struct] = define_model_Q3()

% General functions
L = @(widthL,centerL,x) (widthL.^2)./((x-centerL).^2 + widthL.^2);
G = @(widthG,centerG,x) exp(-((x - centerG)./widthG).^2);
S = @(heightS,slopeS,centerS,x) heightS./( 1 + exp( -2.*slopeS.*(x-centerS) ) );

% First exciton
first_ex_L = @(a,w1,E1,x) a.*L(w1,E1,x);
first_ex_GL = @(a,wL,wG,Ec,x) a.*piecewiseGL(wL,wG,Ec,-1,x);
undoped_1s = @(a,w1,E1,wL,wG,Ec,x) first_ex_L(a,w1,E1,x) + first_ex_GL(a,wL,wG,Ec,x);

% Band structure
BAND_step = @(h,k,Eg,x) S(h,k,Eg,x);
BAND_fine1 = @(c1,bw,B,x) c1.*G(bw,B,x);
BAND_fine2 = @(c2,bw,B,beep,x) c2.*G(bw,B + beep,x);
BAND_fine3 = @(c3,bw,B,beep,x) c3.*G(bw,B + 2*beep,x);
BAND = @(Eg,h,k,c1,c2,c3,B,bw,beep,x) BAND_step(h,k,Eg,x) + BAND_fine1(c1,bw,B,x) + BAND_fine2(c2,bw,B,beep,x) + BAND_fine3(c3,bw,B,beep,x);

% Higher order excitons
Hx2 = @(u,a,w2,E2,x) u.*(a./(1.5^2)).*G(w2,E2,x);
Hx3 = @(u,a,w3,E3,x) u.*(a./(2.5^2)).*G(w3,E3,x);
Hx4 = @(u,a,w4,E4,x) u.*(a./(3.5^2)).*G(w4,E4,x);
Hx = @(u,a,w2,E2,w3,E3,w4,E4,x) Hx2(u,a,w2,E2,x) + Hx3(u,a,w3,E3,x) + Hx4(u,a,w4,E4,x);

% Combined function
TotalAbs = @(Eg,a,w1,E1,wL,wG,Ec,u,w2,E2,w3,E3,w4,E4,h,k,c1,c2,c3,B,bw,beep,x) ...
    undoped_1s(a,w1,E1,wL,wG,Ec,x) + ... 
    Hx(u,a,w2,E2,w3,E3,w4,E4,x) + ... 
    BAND(Eg,h,k,c1,c2,c3,B,bw,beep,x) ;

% Call functions again to only require vector of parameters and data x
call_model = @(W,x) TotalAbs(W(1),W(2),W(3),W(4),W(5),W(6),W(7),W(8),W(9),W(10),W(11),W(12),W(13),W(14),W(15),W(16),W(17),W(18),W(19),W(20),W(21),W(22),x);
call_first_ex_L = @(W,x) first_ex_L(W(2),W(3),W(4),x);
call_first_ex_GL = @(W,x) first_ex_GL(W(2),W(5),W(6),W(7),x);
call_undoped_1s = @(W,x) undoped_1s(W(2),W(3),W(4),W(5),W(6),W(7),x);

call_BAND_step = @(W,x) BAND_step(W(15),W(16),W(1),x);
call_BAND_fine1 = @(W,x) BAND_fine1(W(17),W(21),W(20),x);
call_BAND_fine2 = @(W,x) BAND_fine2(W(18),W(21),W(20),W(22),x);
call_BAND_fine3 = @(W,x) BAND_fine3(W(19),W(21),W(20),W(22),x);
call_BAND = @(W,x) BAND(W(1),W(15),W(16),W(17),W(18),W(19),W(20),W(21),W(22),x);

call_Hx2 = @(W,x) Hx2(W(8),W(2),W(9),W(10),x);
call_Hx3 = @(W,x) Hx3(W(8),W(2),W(11),W(12),x);
call_Hx4 = @(W,x) Hx4(W(8),W(2),W(13),W(14),x);
call_Hx = @(W,x) Hx(W(8),W(2),W(9),W(10),W(11),W(12),W(13),W(14),x);

% Compact into a struct
func_struct = struct();
func_struct.model = call_model;
func_struct.first_ex_L = call_first_ex_L;
func_struct.first_ex_GL = call_first_ex_GL;
func_struct.undoped_1s = call_undoped_1s;
func_struct.BAND_step = call_BAND_step;
func_struct.BAND_fine1 = call_BAND_fine1;
func_struct.BAND_fine2 = call_BAND_fine2;
func_struct.BAND_fine3 = call_BAND_fine3;
func_struct.BAND = call_BAND;
func_struct.Hx2 = call_Hx2;
func_struct.Hx3 = call_Hx3;
func_struct.Hx4 = call_Hx4;
func_struct.Hx = call_Hx;

end

