function [ dm e_m e_n single_phase] = calc_mass_transfer(k_mt,dmdt_max,m01,m02,P,T,Pc,Tc,w,tk,MW,dt,ASindex)

N = numel(MW);
invMW = zeros(N,1);
for i=1:N
    invMW(i,1) = 1/MW(i,1);
end

ntot01 = m01'*invMW;
n01 = m01.*invMW;
x01 = n01/ntot01; 

ntot02 = m02'*invMW;
n02 = m02.*invMW;
x02 = n02/ntot02; 

ntot = ntot01 + ntot02;
c0 = ntot01/(ntot);

xi1 = x01;
xi2 = x02;

single_phase = 0;
failed_LLE = 0;

[x1 x2 c] = Matlab_mLLE_new(c0,P,T,Pc,Tc,w,tk,x01,x02);

if (c==NaN | c==Inf)
    failed_LLE = 1;
end
for i=1:N
    if (x1(i)==NaN | x1(i)==Inf)
        failed_LLE = 1;
    end
    if (x2(i)==NaN | x2(i)==Inf)
        failed_LLE = 1;
    end
end

if (failed_LLE)
    n1 = n01;
    n2 = n02;
else
    if (c==0 | (abs(x1(ASindex,1)-x2(ASindex,1))<1e-6))
        n2 = n01 + n02;
        n1 = n01 - n01;

        x2 = n2/ntot;
        x1 = x2;

        single_phase = 1;
    else
        ntot1 = c*ntot;
        ntot2 = ntot - ntot1;

        n1 = x1*ntot1;
        n2 = x2*ntot2;
    end
end

m1 = n1.*MW;
m2 = n2.*MW;

%conservation errors for debug
e_m = abs(m1+m2-m01-m02);
e_n = abs(n1+n2-n01-n02);

if (single_phase == 1)
    dm = (m2 - m02);
else
    dm = k_mt*(m2 - m02);
    for i=1:numel(dm)
        dm(i,1) = sign(dm(i,1))*min(abs(dm(i,1)),dmdt_max(i,1)*dt);
    end
end