# FUNCTION TO CALCULATE ADJOINT OPTIMAL POWER FLOW OPTION - RECTANGULAR
function adj_rect(crit_node,Vm0,θ0,QgY0)

    # Initialize inputs
    ygen1= ygen[2:end,:]; ngy1 = ngy_inv;
    ngy_nph1 = Int(sum(bus_count[Int.(ygen1[:,1])]))
    max_iter = 20; global err_vec = zeros(50)
    Cbg1 = zeros(nb_nph,ngy_nph1); count_c = 1
    for i=1:ngy1
        count_r = 1;
        for j=1:nb
            if bus[j,1] == ygen1[i,1]
                if bus_count[Int(bus[j,1])] == 3
                    Cbg1[count_r:count_r+2,count_c:count_c+2] = Matrix{Float64}(I, 3, 3)
                    count_r = count_r+3; count_c = count_c+3;
                elseif bus_count[Int(bus[j,1])] == 2
                    Cbg1[count_r:count_r+1,count_c:count_c+1] = Matrix{Float64}(I, 2, 2)
                    count_r = count_r+2; count_c = count_c+2
                elseif bus_count[Int(bus[j,1])] == 1
                    Cbg1[count_r,count_c] = 1
                    count_r = count_r+1; count_c = count_c+1
                end
            else
                count_r = count_r+Int(bus_count[Int(bus[j,1])])
            end
        end
    end
    global Cbg_mod = Cbg1[idx_noref,:]

    # Provide a start point and initialize model
    (Vm0,θ0,PgY0,QgY1,iter) = FP_pf(Vm0,θ0,QgY0); QgY0 = QgY1[4:end];
    # start_point_pol(4); global Vm0 = s0[1:nb_nph]; global θ0 = s0[nb_nph+1:2*nb_nph]; global QgY0 = zeros(ngy_nph-3);
    V0 = Vm0[idx_ref].*exp.(im*θ0[idx_ref]); global w = -YLL_inv*YL0*V0;
    Vd0 = Vm0.*cos.(θ0); Vq0 = Vm0.*sin.(θ0);
    err = 1e4; c_iter = 0; stat = "LOCALLY_SOLVED";
    global PgY = zeros(ngy_nph1);
    global DloadP = zeros(nld_nph); global DloadQ = zeros(nld_nph);

    # Start iterations with specified maximum iterations
    while err>1e-4 && c_iter<max_iter && (stat == "LOCALLY_SOLVED" || stat == "ALMOST_LOCALLY_SOLVED")
        c_iter = c_iter + 1
        global m = Model(with_optimizer(Ipopt.Optimizer,warm_start_init_point="yes",print_level=0, tol=1e-6));
        JuMP.register(m, :FP_pf_adj, 2*nb_nph+ngy_inv, FP_pf_adj, autodiff=true)

        # Define variables
        global Vd = @variable(m, [i=1:3] )
        global Vq = @variable(m,  [i=1:3] )
        global QgY = @variable(m, [i=1:ngy_nph1], start = QgY0[i] )

        # Define bounds on active and reactive power
        Pg_idx = [2 6 10];  Qg_idx = [4 8 12];
        count = 1;
        for i=1:ngy1
            j = Int(ygen1[i,1]); count1 = 1
            for k = 3:2:7
                if bus[j,k] != 0
                    PgY[count] = ygen1[i,Pg_idx[count1]]/baseMVA
                    @constraint(m, ygen1[i,Qg_idx[count1]]/baseMVA <= QgY[count] <= ygen1[i,Qg_idx[count1]+1]/baseMVA )
                    count = count + 1;
                end
                count1 = count1 + 1;
            end
        end
        Cbg_mod = Cbg1[idx_noref,:];
        Ygen_d = Cbg_mod*PgY; Ygen_q = Cbg_mod*QgY;

        # Calculate overall wye-connected loads in system
        global Sy=zeros(nly_nph,7); count = 1; Py_idx = [2 3 4];
        for i=1:nly
            j = Int(yload[i,1]); count1 = 1;  count2 = 1
            for k = 3:2:7
                if bus[j,k] != 0
                    Sy[count,1:6] = yload[i,Py_idx[count1]:3:Py_idx[count1]+15]/baseMVA
                    Sy[count,7] = sum(bus_count[1:j-1]) + count2
                    count = count + 1; count2 = count2 + 1;
                end
            count1 = count1 + 1;
            end
        end

        YloadP = zeros(nly_nph); YloadQ = zeros(nly_nph);
        for i=1:nly_nph
            YloadP[i] = Sy[i,1]+Sy[i,3]*Vm0[Int(Sy[i,7])]+Sy[i,5]*Vm0[Int(Sy[i,7])]^2
            YloadQ[i] = Sy[i,2]+Sy[i,4]*Vm0[Int(Sy[i,7])]+Sy[i,6]*Vm0[Int(Sy[i,7])]^2
        end
        Yload = CblY_mod*(-YloadP - im*YloadQ);

        # Calculate overall delta-connected loads in system
        global Sd = zeros(nld_nph,7); count = 1; Pd_idx = [2 3 4]; Pd_idx1 = 2:3:17;
        for i=1:nld
            j = Int(dload[i,1]); count1 = 1;  count2 = 1
            for k = 1:nld_count[j]
                Sd[count,1:6] = dload[i,Pd_idx[count1]:3:Pd_idx[count1]+15]/baseMVA
                if nld_count[j] == 1
                    for kk=1:6
                        if !isempty(findall(idx->idx!=0,dload[i,Pd_idx1[kk]:Pd_idx1[kk]+2]))
                            dum1 = findall(idx->idx!=0,dload[i,Pd_idx1[kk]:Pd_idx1[kk]+2])[1]
                            Sd[count,kk] = dload[i,dum1+Pd_idx1[kk]-1]/baseMVA
                        end
                    end
                end
                Sd[count,7] = sum(bus_count[1:j-1]) + count2
                count = count + 1; count2 = count2 + 1;
                count1 = count1 + 1;
            end
        end
        for i=1:nld_nph
            DloadP[i] = Sd[i,1]+Sd[i,3]*Vm0[Int(Sd[i,7])]+Sd[i,5]*Vm0[Int(Sd[i,7])]^2
            DloadQ[i] = Sd[i,2]+Sd[i,4]*Vm0[Int(Sd[i,7])]+Sd[i,6]*Vm0[Int(Sd[i,7])]^2
        end
        Dload = CblD_mod*(-DloadP - im*DloadQ)

        # Power balance constraints
        Vm_old =  Vm0[idx_noref]; θ_old = θ0[idx_noref]
        V = Vm_old.*exp.(im*θ_old);
        V_cdY = diagm(0=>conj(V)); V_cdD = diagm(0=>H*conj(V));
        dum1 = YLL_inv*pinv(V_cdY); dum1_d = real(dum1)*Cbg_mod; dum1_q = imag(dum1)*Cbg_mod;
        V_new1 = w + dum1*conj(Yload) + YLL_inv*H'*pinv(V_cdD)*conj(Dload)
        Vd_dep =  real(V_new1) + dum1_d*PgY + dum1_q*QgY;
        Vq_dep = imag(V_new1) + dum1_q*PgY - dum1_d*QgY;
        Vd_all = [Vd_dep[1:length(idx_noref1)]; Vd; Vd_dep[length(idx_noref1)+1:end] ];
        Vq_all = [Vq_dep[1:length(idx_noref1)]; Vq; Vq_dep[length(idx_noref1)+1:end] ];

        # Define bounds on voltage magnitude and set slack bus angle to 0
        @constraint(m, Vd .== [cos(0); cos(-2*pi/3); cos(2*pi/3)] )
        @constraint(m, Vq .== [sin(0); sin(-2*pi/3); sin(2*pi/3)] )
        global v_soft_min =  @variable(m,[i=1:nb_nph1], start = 0 )
        global v_soft_max =  @variable(m,[i=1:nb_nph1],start = 0 )
        @constraint(m,  v_soft_min .>= 0 ); @constraint(m,  v_soft_max .>= 0 )
        global soft_Vlim = @expression(m, sum(v_soft_min[i] for i=1:nb_nph1) + sum(v_soft_max[i] for i=1:nb_nph1) )
        # Vlim_lin = Vd0.^2+Vq0.^2+2*Vd0.*(Vd-Vd0)+2*Vq0.*(Vq-Vq0);
        count=1;
        for i=1:nb
            for j = 3:2:7
                if bus[i,j] != 0  && bus[i,2] != 3
                    # @constraint(m, bus[i,j]^2 <= Vd_dep[count]^2+Vq_dep[count]^2 <= bus[i,j+1]^2 )
                    # @NLconstraint(m, bus[i,j]^2 <= FP_pf_adj(Vm0,θ0,QgY)[count] <= bus[i,j+1]^2 )
                    @constraint(m, v_soft_min[count] >= bus[i,j]^2- Vd_dep[count]^2-Vq_dep[count]^2 )
                    @constraint(m, v_soft_max[count] >= Vd_dep[count]^2+Vq_dep[count]^2 - bus[i,j+1]^2 )
                    # @constraint(m, bus[i,j]^2 <= Vlim_lin[count] <= bus[i,j+1]^2 )
                    count = count + 1
                end
            end
        end

        # Cost function to minimize losses
        if obj_op == 1 || obj_op == 6
            @objective(m, Min, 0 )
        end

        if obj_op == 2
            if crit_node == "All" || crit_node == "all" || crit_node == "ALL"
                crit_node = bus_ID[idx_bus_3ph]
            end
            n_VU = length(crit_node); M = 1e3;
            global bus_target = zeros(n_VU)
            for i=1:length(crit_node)
                dum1 = findfirst(idx -> idx == crit_node[i], bus_ID)
                bus_target[i] = dum1
            end
            global idx_bus_3ph2 = zeros(n_VU*3)
            for i=1:n_VU
                j = Int(bus_target[i])
                idx_bus_3ph2[i*3-2:i*3] = sum(bus_count[1:j-1])+1:sum(bus_count[1:j])
            end
            ad = cos(2*pi/3); aq = sin(2*pi/3)
            Vn = @variable(m, [i=1:n_VU], start = 0 );
            Vp = @variable(m, [i=1:n_VU], start = 1 );
            count = 1;
            for i=1:n_VU
                j=Int.(idx_bus_3ph2[i*3-2:i*3]);
                @constraint(m, Vn[count] == (Vd_all[j[1]]+ad*(Vd_all[j[2]]+Vd_all[j[3]])+aq*(Vq_all[j[2]]-Vq_all[j[3]]))^2 + (Vq_all[j[1]]+ad*(Vq_all[j[2]]+Vq_all[j[3]])-aq*(Vd_all[j[2]]-Vd_all[j[3]]))^2 )
                @constraint(m, Vp[count] == (Vd_all[j[1]]+ad*(Vd_all[j[2]]+Vd_all[j[3]])-aq*(Vq_all[j[2]]-Vq_all[j[3]]))^2 + (Vq_all[j[1]]+ad*(Vq_all[j[2]]+Vq_all[j[3]])+aq*(Vd_all[j[2]]-Vd_all[j[3]]))^2 )
                count = count + 1;
            end
            # global cost_VUF = @NLexpression(m, sum(Vn[i]/Vp[i] for i=1:n_VU) )
            global cost_VUF = @NLexpression(m, sum(Vn[i] for i=1:n_VU) )
            @NLobjective(m, Min, M*cost_VUF)# + 1e1*t )
            # @NLobjective(m, Min, M*cost_VUF + 1e-1*V_pen )
            # @NLobjective(m, Min, M*cost_VUF + 1e-5*sum(Qg[i]^2 for i=4:ngy_nph))
        end
        if obj_op == 3
            @objective(m, Min, 1e3*sum((QgY[i]-QgY0[i])^2 for i=1:ngy_nph1) + 1e2*soft_Vlim )
        end

        # Solve and update solution
        @time status = optimize!(m); stat = string(termination_status(m))
        QgY0 = JuMP.value.(QgY);
        Vd0_dep=  real(V_new1) + dum1_d*PgY + dum1_q*QgY0;
        Vd1 = [Vd0_dep[1:length(idx_noref1)]; JuMP.value.(Vd); Vd0_dep[length(idx_noref1)+1:end]];
        Vq0_dep = imag(V_new1) + dum1_q*PgY - dum1_d*QgY0;
        Vq1 = [Vq0_dep[1:length(idx_noref1)]; JuMP.value.(Vq); Vq0_dep[length(idx_noref1)+1:end]];
        V0 = Vm0.*exp.(im*θ0); V1 = Vd1+im*Vq1;
        diff = abs.(V1-V0); err = sum(diff);
        Vm0 = abs.(V1); θ0 = angle.(V1);
        print("OPF- (",c_iter,") Error = ",err,", ",stat,"\n"); err_vec[c_iter] = err;
    end

    # Calculate missing solutions
    Pg_SS = zeros(3); Qg_SS = zeros(3); count = 1;
    for i=1:nb_nph
        if count<4 && i == idx_ref[count]
            Pg_SS[count] = sum(Vm0[i]*Vm0[j]*(G[i,j]*cos(θ0[i]-θ0[j])+B[i,j]*sin(θ0[i]-θ0[j])) for j=1:nb_nph)
            Qg_SS[count] = sum(Vm0[i]*Vm0[j]*(G[i,j]*sin(θ0[i]-θ0[j])-B[i,j]*cos(θ0[i]-θ0[j])) for j=1:nb_nph)
            count = count + 1
        end
    end
    PgY0 = [Pg_SS; PgY];
    Ill1 = zeros(nld_nph)*(1+im*0); count = 1
    for i=1:nld
        j = Int(dload[i,1])
        for k=Int.(sum(bus_count[1:j-1])+1:sum(bus_count[1:j]))
            kk = Int(Tl[k])
            if k != kk
                Vll1 = Vm0[k]*(cos(θ0[k])+im*sin(θ0[k]))-Vm0[kk]*(cos(θ0[kk])+im*sin(θ0[kk]))
                Ill1[count] =  (DloadP[count]-im*DloadQ[count])/conj(Vll1)
                count = count + 1
            end
        end
    end
    Ip1 = Tlp*Ill1; count = 1; P_d2y = zeros(nldy_nph)
    for i=1:nld
        j = Int(dload[i,1])
        for k=Int.(sum(bus_count[1:j-1])+1:sum(bus_count[1:j]))
            P_d2y[count] =  Vm0[k]*cos(θ0[k])*real(Ip1[count]) + Vm0[k]*sin(θ0[k])*imag(Ip1[count])
            count = count + 1
        end
    end
    Vll0 = zeros(nb_nph)
    for i=1:nb_nph
        k = Int(Tl[i])
        Vll0[i] =  sqrt(Vm0[i]^2 + Vm0[k]^2 - 2*Vm0[i]*Vm0[k]*cos(θ0[i]-θ0[k]) )
    end
    print_vuf_pol(PgY0,P_d2y,Vm0,θ0,Vll0);
    return Vm0,θ0,QgY0,stat
end

#----------------------------------------------------------------------------------------#
# FUNCTION TO CALCULATE POWER FLOW FOR ADJOINT METHOD
function FP_pf_adj(Vm0,θ0,QgY0)

    # Other inputs
    CblY_mod = CblY[idx_noref,:];

    # Power balance constraints
    V_new = zeros(nb_nph1); iter = 0;
    Vm_old =  Vm0[idx_noref]; θ_old = θ0[idx_noref]
    YloadP = zeros(nly_nph); YloadQ = zeros(nly_nph);
    DloadP = zeros(nld_nph); DloadQ = zeros(nld_nph);
    while iter<50
        V = Vm_old.*exp.(im*θ_old);
        V_cdY = diagm(0=>conj(V)); V_cdD = diagm(0=>H*conj(V));
        for i=1:nly_nph
            YloadP[i] = Sy[i,1]+Sy[i,3]*Vm0[Int(Sy[i,7])]+Sy[i,5]*Vm0[Int(Sy[i,7])]^2
            YloadQ[i] = Sy[i,2]+Sy[i,4]*Vm0[Int(Sy[i,7])]+Sy[i,6]*Vm0[Int(Sy[i,7])]^2
        end
        Yload = CblY_mod*(-YloadP - im*YloadQ);
        for i=1:nld_nph
            DloadP[i] = Sd[i,1]+Sd[i,3]*Vm0[Int(Sd[i,7])]+Sd[i,5]*Vm0[Int(Sd[i,7])]^2
            DloadQ[i] = Sd[i,2]+Sd[i,4]*Vm0[Int(Sd[i,7])]+Sd[i,6]*Vm0[Int(Sd[i,7])]^2
        end
        Dload = CblD_mod*(-DloadP - im*DloadQ)
        dum1 = YLL_inv*pinv(V_cdY); dum1_d = real(dum1)*Cbg_mod; dum1_q = imag(dum1)*Cbg_mod;
        V_new1 = w + dum1*conj(Yload) + YLL_inv*H'*pinv(V_cdD)*conj(Dload)
        Vd_dep =  real(V_new1) + dum1_d*PgY + dum1_q*QgY0;
        Vq_dep = imag(V_new1) + dum1_q*PgY - dum1_d*QgY0;
        Vm_old = @NLexpression(m, [i=1:ngy_nph-3], sqrt(Vd_dep[i]^2+Vq_dep[i]^2));
        θ_old = atan.(Vq_dep./Vd_dep);
        Vm0 = [Vm_old[1:length(idx_noref1)]; Vm0[idx_ref]; Vm_old[length(idx_noref1)+1:end] ]
        θ0 = [θ_old[1:length(idx_noref1)]; θ0[idx_ref]; θ_old[length(idx_noref1)+1:end] ]
        iter = iter + 1;
    end
    return Vm_old.^2
end

#----------------------------------------------------------------------------------------#
