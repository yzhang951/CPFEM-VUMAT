clear all;
load 'lam.rpt';
load 'FCC.rpt';
load 'BCC.rpt';

%plot(lam(:,1),lam(:,3),'o',FCC(:,1),FCC(:,2), lam(:,1),lam(:,2), 'o', BCC(:,1),BCC(:,2));


fcc = 644; bcc = 346;
c = fcc/(fcc+bcc);
c = 0.65;

lam_ave = lam(:,3)*c + lam(:,2)*(1-c);
ave_lam = FCC(:,2)*c + BCC(:,2)*(1-c);

plot(lam(:,1),lam_ave, 'r', lam(:,1),ave_lam, 'b');

output = [lam(:,1) lam_ave ave_lam];
save stress.dat output -ascii

