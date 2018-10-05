function [dm1 dm2 rr1 rr2 err_R1 err_R2 nChemSteps] = calc_chem_source(m1,m2,V1,V2,k1,k2,dt,single_phase,chem_tol)

N=numel(m1);
dm1 = zeros(N,1);
dm2 = zeros(N,1);
rr1 = zeros(N+1,1);
rr2 = zeros(N+1,1);
err_R1 = zeros(N,1);
err_R2 = zeros(N,1);
nChemSteps = 0;
dtTotal = 0;
dt1 = dt;

while (dtTotal < dt)
    if (single_phase==1)
        [dm2] = calc_chem_RR_Wiehe4(V2,m2,k2,dt1);
        dm1 = zeros(N,1);
    else              
        [dm1] = calc_chem_RR_Wiehe4(V1,m1,k1,dt1);
        [dm2] = calc_chem_RR_Wiehe4(V2,m2,k2,dt1);
    end
    dtFactorMin = 1;
    for i=1:numel(m2)
        if ((m2(i,1) + dm2(i,1)) < -chem_tol)
            dtFactor = m2(i,1)/(abs(dm2(i,1)) + 1e-12);
            dtFactorMin = min(dtFactorMin,dtFactor);
        end
    end
    for i=1:numel(m1)
        if ((m1(i,1) + dm1(i,1)) < -chem_tol)
            dtFactor = m1(i,1)/(abs(dm1(i,1)) + 1e-12);
            dtFactorMin = min(dtFactorMin,dtFactor);
        end
    end
    dt1 = dt1*dtFactorMin;
    dm1 = dm1*dtFactorMin;
    dm2 = dm2*dtFactorMin;
    m1 = m1 + dm1;
    m2 = m2 + dm2;          
    rr1(1:N) = rr1(1:N) + dm1;
    rr2(1:N) = rr2(1:N) + dm2;

    tmp_m01 = m1;
    tmp_m02 = m2;
    for i=1:N
        if (m1(i,1) < 0)
            m1(i,1) = 0;
        end
        if (m2(i,1) < 0)
            m2(i,1) = 0;
        end
    end
    tmp_m1 = m1;
    tmp_m2 = m2;
    err_R1 = err_R1 + abs(tmp_m1 - tmp_m01);
    err_R2 = err_R2 + abs(tmp_m2 - tmp_m02);

    dtTotal = dtTotal + dt1;
    dt1 = dt - dtTotal;
    nChemSteps = nChemSteps + 1;
end
