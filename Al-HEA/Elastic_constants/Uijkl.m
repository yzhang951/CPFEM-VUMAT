function [a b c] = Uijkl(C11, C12, C44, K, G)
% calculate the diffraction elastic constants

L1 = [C11+2*C12   C11-C12   2*C44];
L0 = [3*K   2*G   2*G];
S  = [3*K/(3*K+4*G) 2*(3*K+6*G)/(15*K+20*G) 2*(3*K+6*G)/(15*K+20*G)];
I  = [1 1 1];
T  = (I+S.*(L0.^(-1)).*(L1-L0)).^(-1);
U  = T.*(L0.^(-1));

a  = U(1)/3;
b  = U(2)/2;
c  = U(3)/2;
