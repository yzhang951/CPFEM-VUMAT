clearvars;

noel = 8000;
volume_fraction = 0.33;

theta_threshold = 5;

hkl_fcc = [ 2 0 0;
			2 2 0;
			1 1 1;
			3 1 1;
			3 3 1;];

for i=1:length(hkl_fcc(:,1))
	FCC_LD(i).num = 0;
	FCC_LD(i).ID = 0;
	FCC_TD(i).num = 0;
	FCC_TD(i).ID = 0;
end


hkl_bcc = [ 2 0 0;
			1 1 0;
			2 1 1;
			3 2 1;];


for i=1:length(hkl_bcc(:,1))
	BCC_LD(i).num = 0;
	BCC_LD(i).ID = 0;
	BCC_TD(i).num = 0;
	BCC_TD(i).ID = 0;
end


fp = fopen('aeuler','w');
fprintf(fp,'\n\n2\n\n\n\n\n');
for i=1:noel
	a = rand*360;
	b = rand*360;
	c = rand*360;
	phase = rand;
	if(phase>volume_fraction)
		p = 1;
		for j=1:length(hkl_fcc(:,1))
			theta = check_angle(a,b,c,hkl_fcc(j,1),hkl_fcc(j,2),hkl_fcc(j,3),0,1,0);
			if(theta < theta_threshold)
				FCC_LD(j).num = FCC_LD(j).num + 1;
				FCC_LD(j).ID(FCC_LD(j).num) = i;
			end

			theta = check_angle(a,b,c,hkl_fcc(j,1),hkl_fcc(j,2),hkl_fcc(j,3),1,0,0);
			if(theta < theta_threshold)
				FCC_TD(j).num = FCC_TD(j).num + 1;
				FCC_TD(j).ID(FCC_TD(j).num) = i;
			end

		end
	else
		p = 2;
		for j=1:length(hkl_bcc(:,1))
			theta = check_angle(a,b,c,hkl_bcc(j,1),hkl_bcc(j,2),hkl_bcc(j,3),0,1,0);
			if(theta < theta_threshold)
				BCC_LD(j).num = BCC_LD(j).num + 1;
				BCC_LD(j).ID(BCC_LD(j).num) = i;
			end

			theta = check_angle(a,b,c,hkl_bcc(j,1),hkl_bcc(j,2),hkl_bcc(j,3),1,0,0);
			if(theta < theta_threshold)
				BCC_TD(j).num = BCC_TD(j).num + 1;
				BCC_TD(j).ID(BCC_TD(j).num) = i;
			end

		end
	end
	fprintf(fp,'%f %f %f 0 0 0 %d\n',a,b,c,p);

end
fclose(fp);

fp = fopen('grain_family.inp','w');

for i=1:length(hkl_fcc(:,1))
	fprintf(fp,'*Elset, elset=FCC-LD%d%d%d, instance=Part-1-1\n',hkl_fcc(i,1),hkl_fcc(i,2),hkl_fcc(i,3));
	for j=1:FCC_LD(i).num
		fprintf(fp,'%7d,',FCC_LD(i).ID(j));
		if(mod(j,4)==0)
			fprintf(fp,'\n');
		end
	end
	fprintf(fp,'\n');

	fprintf(fp,'*Elset, elset=FCC-TD%d%d%d, instance=Part-1-1\n',hkl_fcc(i,1),hkl_fcc(i,2),hkl_fcc(i,3));
	for j=1:FCC_TD(i).num
		fprintf(fp,'%7d,',FCC_TD(i).ID(j));
		if(mod(j,4)==0)
			fprintf(fp,'\n');
		end
	end
	fprintf(fp,'\n');

end

for i=1:length(hkl_bcc(:,1))
	fprintf(fp,'*Elset, elset=BCC-LD%d%d%d, instance=Part-1-1\n',hkl_bcc(i,1),hkl_bcc(i,2),hkl_bcc(i,3));	
	for j=1:BCC_LD(i).num
		fprintf(fp,'%7d,',BCC_LD(i).ID(j));
		if(mod(j,4)==0)
			fprintf(fp,'\n');
		end
	end
	fprintf(fp,'\n');

	fprintf(fp,'*Elset, elset=BCC-TD%d%d%d, instance=Part-1-1\n',hkl_bcc(i,1),hkl_bcc(i,2),hkl_bcc(i,3));	
	for j=1:BCC_TD(i).num
		fprintf(fp,'%7d,',BCC_TD(i).ID(j));
		if(mod(j,4)==0)
			fprintf(fp,'\n');
		end
	end
	fprintf(fp,'\n');
end

fclose(fp);
