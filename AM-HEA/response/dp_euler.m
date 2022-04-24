clear all;

load 'FCC.inp';
load 'BCC.inp';

noel = 1551;
BCC_euler = [0,0,0];

% FCC (1 -1  1) // BCC (1  1  0)
% FCC [1  1  0] // BCC [1 -1  1]
% FCC [1 -1 -2] // BCC [1 -1 -2]

e1 = [-1 1 0]/sqrt(2);
e2 = [1 1 1]/sqrt(3);
e3 = [1 1 -2]/sqrt(6);
Q_FCC = [e1;e2;e3]
Q_BCC = [e2;-e1;e3]

R = Q_FCC

th = mod(-acosd(R(3,3)),180);
ph = mod(atan(-R(3,1)/R(3,2))*180/pi,360);
om = mod(atan(R(1,3)/R(2,3))*180/pi,360);

Q = [	cosd(ph)*cosd(om)-sind(ph)*sind(om)*cosd(th)	sind(ph)*cosd(om)+cosd(ph)*sind(om)*cosd(th)	sind(om)*sind(th);
		-cosd(ph)*sind(om)-sind(ph)*cosd(om)*cosd(th)	-sind(ph)*sind(om)+cosd(ph)*cosd(om)*cosd(th)	cosd(om)*sind(th);
		sind(ph)*sind(th)								-cosd(ph)*sind(th)								cosd(th);];

th_FCC = th; ph_FCC = ph; om_FCC = om;


R = Q_BCC

th = mod(-acosd(R(3,3)),180);
ph = mod(atan(-R(3,1)/R(3,2))*180/pi,360);
om = mod(atan(R(1,3)/R(2,3))*180/pi,360);

Q = [	cosd(ph)*cosd(om)-sind(ph)*sind(om)*cosd(th)	sind(ph)*cosd(om)+cosd(ph)*sind(om)*cosd(th)	sind(om)*sind(th);
		-cosd(ph)*sind(om)-sind(ph)*cosd(om)*cosd(th)	-sind(ph)*sind(om)+cosd(ph)*cosd(om)*cosd(th)	cosd(om)*sind(th);
		sind(ph)*sind(th)								-cosd(ph)*sind(th)								cosd(th);];

th_BCC = th; ph_BCC = ph; om_BCC = om;


fp = fopen('aeuler','w');
fprintf(fp,'\n\n2\n\n\n\n\n');

for i=1:noel
	if ismember(i, FCC) == 1
		fprintf(fp,'%f %f %f 0 0 0 1\n',ph_FCC,th_FCC,om_FCC);
%		fprintf(fp,'%f %f %f 0 0 0 2\n',ph_BCC,th_BCC,om_BCC);
	elseif ismember(i, BCC) == 1
%		fprintf(fp,'%f %f %f 0 0 0 1\n',ph_FCC,th_FCC,om_FCC);
		fprintf(fp,'%f %f %f 0 0 0 2\n',ph_BCC,th_BCC,om_BCC);
	end
end

length(FCC)
length(BCC)
length(FCC)+length(BCC)
