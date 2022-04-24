function theta = check_angle(a, b, c, h, k, l, m1, m2, m3) 
        phi = deg2rad(a);
        theta = deg2rad(b);
        omega = deg2rad(c);
        Q = zeros(3,3);
        Q(1,1) = cos(phi)*cos(omega) - sin(phi)*sin(omega)*cos(theta);
        Q(1,2) = sin(phi)*cos(omega) + cos(phi)*sin(omega)*cos(theta);
        Q(1,3) = sin(omega)*sin(theta);
        Q(2,1) = -cos(phi)*sin(omega)- sin(phi)*cos(omega)*cos(theta);
        Q(2,2) = -sin(phi)*sin(omega)+ cos(phi)*cos(omega)*cos(theta);
        Q(2,3) = cos(omega)*sin(theta);
        Q(3,1) = sin(phi)*sin(theta);
        Q(3,2) = -cos(phi)*sin(theta);
        Q(3,3) = cos(theta);                       % Rotation Tensor
		m = [m1 m2 m3];
		m = m/norm(m);
		p=zeros(24,3);
		p(1:6,:) = perms([h k l]);
		p(7:12,:) = perms([-h k l]);
		p(13:18,:) = perms([-h -k l]);
		p(19:24,:) = perms([h -k -l]);

		theta = 90;

		for i=1:24
			x=p(i,:)*Q';
			x=x/norm(x);
			alpha = acosd(abs(x*m'));
			if(alpha<theta)
				theta = alpha;
			end
		end
end
