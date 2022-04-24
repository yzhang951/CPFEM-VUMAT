clearvars;
step = 1e-3;
% Find the elastic constants using gradient descent method
DEC_target = [	136.53;	 % FCC 200 LD
				208.64;  % FCC 220 LD
				238.72;  % FCC 111 LD
				175.14;  % FCC 311 LD
				222.78;  % FCC 331 LD
				-399.34; % FCC 200 TD
				-579.90; % FCC 220 TD
				-784.18; % FCC 111 TD
				-525.59; % FCC 311 TD
				-646.36; % FCC 331 TD
				107.75;  % BCC 200 LD
				197.31;  % BCC 110 LD
				192.44;  % BCC 211 LD
				192.13;  % BCC 321 LD
				-274.29; % BCC 200 TD
				-542.89; % BCC 110 TD
				-511.56; % BCC 211 TD
				-491.94; % BCC 321 TD
				 0;		 % bulk modulus self consistancy
				 0;		 % shear modulus self consistancy
				 171];	 % sample Young's modulus
DEC_target = DEC_target';
Weight = eye(21);
Weight(1,1) = 8;
Weight(1:5,1:5) = 4*Weight(1:5,1:5);
Weight(11:14,11:14) = 4*Weight(11:14,11:14);
Weight(19,19) = 3; Weight(20,20) = 5; Weight(21,21) = 3;
%	  FCC_C11 FCC_C12 FCC_C44 BCC_C11 BCC_C12 BCC_44 K       G	
%x0 = [246.5   147.3   124.7   239.8   58.6    99.0   180.37  66.90];
%x0 = [262.3207  175.9534  125.5797  286.0391  257.4653  129.9958  225.2889   62.2489];
x0 = [262.3207  175.9534  125.5797  486.0391  57.4653  129.9958  225.2889   62.2489];
x1 = x0;
res0  = f(x0)-DEC_target;
res1  = res0;
F0    = res0*Weight*res0';
F1    = F0;
delta_F = 100;
t  = 0;
while (delta_F > 1e-6)
	x0 = x1;
	F0 = F1;
	res0 = res1;
	gradF = x0;
	for i=1:length(x0)
		x = x0;
		x(i) = x(i) + 1e-4;
		res = f(x)-DEC_target;
		F   = res*Weight*res';
		gradF(i) = (F-F0)/1e-4;
	end
	x1 = x0 - step*gradF;
	res1  = f(x1)-DEC_target;
	F1    = res1*Weight*res1';
	delta_F = F0-F1;
	delta_F;
	t = t+1;
end

x1
K = x1(7);G = x1(8);
E = 9*K*G/(3*K+G)
K_FCC = (x1(1)+2*x1(2))/3
K_BCC = (x1(4)+2*x1(5))/3
K
0.67*K_FCC + 0.33*K_BCC
