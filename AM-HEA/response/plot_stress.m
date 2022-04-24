clear all;
clf;
load 'fcc.rpt';
load 'bcc.rpt';
load 'lam.rpt';
load 'lam_ph.rpt';


c = 0.65;
mix_stress = fcc(:,2)*c + bcc(:,2)*(1-c);
plot(bcc(:,1),bcc(:,2), fcc(:,1),fcc(:,2), lam(:,1),lam_ph(:,2),lam(:,1),lam_ph(:,3));
legend('bcc_{mix}','fcc_{mix}', 'bcc_{lam}', 'fcc_{lam}');

plot(lam(:,1),mix_stress, lam(:,1),lam(:,2));
legend('mix','lam');
