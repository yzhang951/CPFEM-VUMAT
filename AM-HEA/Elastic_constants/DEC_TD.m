function E_hkl = DEC_TD(C11, C12, C44, K, G, hkl)
% calculate the diffraction elastic constants
[a b c] = Uijkl(C11,C12,C44,K,G);

h=hkl(1);k=hkl(2);l=hkl(3);
GAMMA = (h*h*k*k+l*l*k*k+h*h*l*l)/(h*h+k*k+l*l)^2;

S_hkl = (3*a-2*b)/3 + 2*(b-c)*GAMMA;
E_hkl = 1/S_hkl;
