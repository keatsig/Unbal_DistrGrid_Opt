### FUNCTION TO CALCULATE OPTIMAL POWER FLOW OPTION-2 (FP_OPF2) - POLAR
function FP_opf_pol(crit_node,Vm0,θ0,QgY0,dnc)

    ## Initialize inputs
    # (Vm0,θ0,QgY0,stat) = TP_OPF(crit_node,Vm0,θ0,QgY0)
    # (Vm0,θ0,PgY0,QgY1,iter) = FP_pf(Vm0,θ0,QgY0,20); QgY0 = QgY1[4:end];
    # start_point(4); Vm0 = s0[1:nb_nph]; θ0 = s0[nb_nph+1:2*nb_nph];

    ## Define active power of inverters
    Pg_idx = [2 6 10];  Qg_idx = [4 8 12]; 
    PgY = spzeros(ngy_inv); num_count= 1;
    if PV_en == 1
        for i=1:ngy_inv
            j = Int(ygen1[i,1]); count1 = 1
            for k = 3:2:7
                if bus[j,k] != 0
                    PgY[num_count]= ygen1[i,Pg_idx[count1]]/baseMVA
                    num_count= num_count+ 1;
                end
                count1 = count1 + 1;
            end
        end
    end

    ## Generate load matrix
    ZIP_load()

    ## Start OPF iterations
    max_iter = 30; global obj_add = 1;
    opf_err = 1; opf_iter = 0; stat = "LOCALLY_SOLVED";
    global err_vec = zeros(max_iter);
    DloadP = zeros(nld_nph); DloadQ = zeros(nld_nph); Vll0 = zeros(nb_nph)
    YloadP = zeros(nly_nph); YloadQ = zeros(nly_nph);
    while opf_err>1e-4 && opf_iter<max_iter && (stat == "LOCALLY_SOLVED" || stat == "ALMOST_LOCALLY_SOLVED")

        ## Setup model and check for slow convergence
        opf_iter = opf_iter + 1;
        if (opf_iter == 11 && pf_iter >= 4)
            opf_iter = max_iter;
        elseif opf_iter>4 && ((err_vec[opf_iter-2:opf_iter-1] == err_vec[opf_iter-4:opf_iter-3]) || (err_vec[opf_iter-1] == err_vec[opf_iter-2] && err_vec[opf_iter-1] == err_vec[opf_iter-3]) || (err_vec[opf_iter-1] >= err_vec[opf_iter-2] >= err_vec[opf_iter-3]))
            opf_iter = max_iter;
        else
            global m = Model(optimizer_with_attributes(Ipopt.Optimizer,"warm_start_init_point"=>"yes","print_level"=>0, "tol"=>1e-7));
            (Vm0,θ0,QgY0,pf_iter) = FP_pf(Vm0,θ0,QgY0,30); print("No of PF iterations: ", pf_iter,"\n")
            for i=1:nb_nph
                k = Int(Tl[i])
                Vll0[i] =  sqrt(Vm0[i]^2 + Vm0[k]^2 - 2*Vm0[i]*Vm0[k]*cos(θ0[i]-θ0[k]) )
            end

            ## Define variables
            global Vm = @variable(m, [i=1:nb_nph], start = Vm0[i] )
            global θ = @variable(m,  [i=1:nb_nph], start = θ0[i] )
            global Vll = @variable(m, [i=1:nb_nph], start = Vll0[i] )
            global QgY = @variable(m, [i=1:ngy_inv], start = 0 )

            ## Define bounds on voltage magnitude and set slack bus angle to 0
            @constraint(m, Vm[idx_ref] .== 1 )
            @constraint(m, θ[idx_ref] .== [0;-2*pi/3;2*pi/3] )
            # global v_soft_min =  @variable(m,[i=1:nb_nph1], start = 0 )
            # global v_soft_max =  @variable(m,[i=1:nb_nph1],start = 0 )
            # @constraint(m,  v_soft_min .>= 0 ); @constraint(m,  v_soft_max .>= 0 )
            # global soft_Vlim = @expression(m, sum(v_soft_min[i] for i=1:nb_nph1) + sum(v_soft_max[i] for i=1:nb_nph1) )
            num_count = 1
            for i=1:nb
                if  bus[i,2] == 3
                    num_count = num_count + 3
                end
                for j = 3:2:7
                    if bus[i,j] != 0  && bus[i,2] != 3
                        @constraint(m, bus[i,j] <= Vm[num_count] <= bus[i,j+1] )
                        # @constraint(m, v_soft_min[num_count-3] >= Vm[num_count] - bus[i,j] )
                        # @constraint(m, v_soft_max[num_count-3] >= bus[i,j+1] - Vm[num_count] )
                        num_count = num_count + 1
                    end
                end
            end
            for i=1:nb_nph
                k = Int(Tl[i])
                @NLconstraint(m, Vm[i]^2 + Vm[k]^2 - 2*Vm[i]*Vm[k]*cos(θ[i]-θ[k]) == Vll[i]^2 )
            end

            ## Define bounds on reactive power
            if PV_en == 1
                num_count= 1;
                for i=1:ngy_inv
                    j = Int(ygen1[i,1]); count1 = 1
                    for k = 3:2:7
                        if bus[j,k] != 0
                            @constraint(m, ygen1[i,Qg_idx[count1]]/baseMVA <= QgY[num_count]<= ygen1[i,Qg_idx[count1]+1]/baseMVA )
                            num_count= num_count+ 1;
                        end
                        count1 = count1 + 1;
                    end
                end
            end

            ## Calculate overall star-connected loads in system
            for i=1:nl_y
                dum1 = (Sy[i,3] + im*Sy[i,4])/(1*exp(im*Sy[i,8]))*(Vm0[Int(Sy[i,7])]*exp(im*θ0[Int(Sy[i,7])]))
                YloadP[i] = Sy[i,1]+real(dum1)+Sy[i,5]*Vm0[Int(Sy[i,7])]^2
                YloadQ[i] = Sy[i,2]+imag(dum1)+Sy[i,6]*Vm0[Int(Sy[i,7])]^2
            end
            for i=nl_y+1:nly_nph
                YloadP[i] = Sy[i,1]+Sy[i,3]*Vm0[Int(Sy[i,7])]+Sy[i,5]*Vm0[Int(Sy[i,7])]^2
                YloadQ[i] = Sy[i,2]+Sy[i,4]*Vm0[Int(Sy[i,7])]+Sy[i,6]*Vm0[Int(Sy[i,7])]^2
            end
            Yload = CblY_mod*(-YloadP - im*YloadQ);

            ## Calculate overall delta-connected loads in system
            for i=1:nld_nph
                j = Int(Tl[Int(Sd[i,7])]);
                dum1 = Vm0[Int(Sd[i,7])]*exp(im*θ0[Int(Sd[i,7])])-Vm0[j]*exp(im*θ0[j])
                dum2 = (Sd[i,3] + im*Sd[i,4])/(sqrt(3)*exp(im*Sd[i,8]))*dum1
                DloadP[i] = Sd[i,1]+real(dum2)+Sd[i,5]*Vll0[Int(Sd[i,7])]^2/3
                DloadQ[i] = Sd[i,2]+imag(dum2)+Sd[i,6]*Vll0[Int(Sd[i,7])]^2/3
            end
            Dload = CblD_FP*(-DloadP - im*DloadQ)

            ## Power balance constraints
            V0 = Vm0[idx_ref].*exp.(im*θ0[idx_ref]);
            V = Vm0[idx_noref].*exp.(im*θ0[idx_noref]);
            V_cdY = diagm(0=>conj(V)); V_cdD = diagm(0=>H*conj(V));
            dum1 = YLL_inv*pinv(V_cdY); dum1_d = real(dum1)*Cbg_mod; dum1_q = imag(dum1)*Cbg_mod;
            V_new1 = w + dum1*conj(Yload) + YLL_inv*H'*pinv(V_cdD)*conj(Dload)
            V_new_d = @NLexpression(m, [i=1:nb_nph1], real(V_new1)[i] + sum(dum1_d[i,k]*PgY[k] for k=1:ngy_inv) + sum(dum1_q[i,k]*QgY[k] for k=1:ngy_inv) )
            V_new_q = @NLexpression(m, [i=1:nb_nph1], imag(V_new1)[i] + sum(dum1_q[i,k]*PgY[k] for k=1:ngy_inv) - sum(dum1_d[i,k]*QgY[k] for k=1:ngy_inv) )
            count=1
            for i in idx_noref
                @NLconstraint(m, V_new_d[count] == Vm[i]*cos(θ[i]) )
                @NLconstraint(m, V_new_q[count] == Vm[i]*sin(θ[i]) )
                count = count+1
            end

            ## Define cost
            cost_polar(Vll);

            ## Solve and update solution
            @time status = optimize!(m); stat = string(termination_status(m))
            V0 = Vm0.*exp.(im*θ0); Vm0 = JuMP.value.(Vm); θ0 = JuMP.value.(θ); QgY0 = JuMP.value.(QgY);
            V1 = Vm0.*exp.(im*θ0); opf_err = sum(abs.(V1-V0));
            err_vec[opf_iter] = round(opf_err, RoundNearest, digits=4);
            print("FP_OPF_pol: (",opf_iter,") Error = ",opf_err,", ",stat,"\n");
            # (Vm0,θ0,QgY0,pf_iter) = FP_pf(Vm0,θ0,QgY0,50); print("No of PF iterations: ", pf_iter,"\n");
        end
    end

    ## Run power flow to convergence if OPF did not converge
    if opf_iter == max_iter
        println("\n****************OPF DID NOT CONVERGE!!*******************\n")
        (Vm0,θ0,QgY0,pf_iter) = FP_pf(Vm0,θ0,QgY0,50); print("No of PF iterations: ", pf_iter,"\n");
        if dnc == 0
            stat = "DID NOT CONVERGE";
        end
    end

    ## Calculate missing solutions
    Vph = Vm0.*exp.(im*θ0);
    Iph = (G+im*B)*Vph; S_inj = Vph.*conj(Iph); P_inj = real(S_inj);
    for i=1:nb_nph
        k = Int(Tl[i])
        Vll0[i] =  sqrt(Vm0[i]^2 + Vm0[k]^2 - 2*Vm0[i]*Vm0[k]*cos(θ0[i]-θ0[k]) )
    end
    print_VU(Vm0,θ0,Vll0);
    global P_loss = sum(P_inj)*baseMVA*1e3
    println("TOTAL TIME:")
    return Vm0,θ0,QgY0,stat
end

#----------------------------------------------------------------------------------------#

# FUNCTION TO CALCULATE OPTIMAL POWER FLOW (FP_OPF) - RECTANGULAR
function FP_opf_rect(crit_node,Vm0,θ0,QgY0,dnc)

    ## Initialize inputs
    # (Vm0,θ0,QgY0,stat) = TP_OPF(crit_node,Vm0,θ0,QgY0)
    # start_point(4); Vm0 = s0[1:nb_nph]; θ0 = s0[nb_nph+1:2*nb_nph];

    ## Define active power of inverters
    Pg_idx = [2 6 10];  Qg_idx = [4 8 12];
    PgY = spzeros(ngy_inv); num_count= 1;
    if PV_en == 1
        for i=1:ngy_inv
            j = Int(ygen1[i,1]); count1 = 1
            for k = 3:2:7
                if bus[j,k] != 0
                    PgY[num_count]= ygen1[i,Pg_idx[count1]]/baseMVA
                    num_count= num_count+ 1;
                end
                count1 = count1 + 1;
            end
        end
    end

    ## Generate load matrix
    ZIP_load()

    ## Start OPF iterations
    max_iter = 25; global obj_add = 1;
    opf_err = 1e4; opf_iter = 0; stat = "LOCALLY_SOLVED";
    global err_vec = zeros(max_iter); Vll0 = zeros(nb_nph)
    YloadP = zeros(nly_nph); YloadQ = zeros(nly_nph);
    DloadP = zeros(nld_nph); DloadQ = zeros(nld_nph);
    while opf_err>1e-4 && opf_iter<max_iter && (stat == "LOCALLY_SOLVED" || stat == "ALMOST_LOCALLY_SOLVED")

        ## Setup model and check for slow convergence
        opf_iter = opf_iter + 1;
        if (opf_iter == 11 && pf_iter >= 4)
            opf_iter = max_iter;
        elseif opf_iter>4 && ((err_vec[opf_iter-2:opf_iter-1] == err_vec[opf_iter-4:opf_iter-3]) || (err_vec[opf_iter-1] == err_vec[opf_iter-2] && err_vec[opf_iter-1] == err_vec[opf_iter-3]) || (err_vec[opf_iter-1] >= err_vec[opf_iter-2] >= err_vec[opf_iter-3]))
            opf_iter = max_iter;
        else
            global m = Model(optimizer_with_attributes(Ipopt.Optimizer,"warm_start_init_point"=>"yes","print_level"=>0, "tol"=>1e-7));
            # (Vm0,θ0,PgY0,QgY1,iter) = FP_pf(Vm0,θ0,QgY0,50); print("No of PF iterations: ", iter,"\n")
            Vd0 = Vm0.*cos.(θ0); Vq0 = Vm0.*sin.(θ0);
            for i=1:nb_nph
                k = Int(Tl[i])
                Vll0[i] =  sqrt(Vm0[i]^2 + Vm0[k]^2 - 2*Vm0[i]*Vm0[k]*cos(θ0[i]-θ0[k]) )
            end

            ## Define variables
            global Vd = @variable(m, [i=1:nb_nph], start = Vd0[i] )
            global Vq = @variable(m, [i=1:nb_nph], start = Vq0[i] )
            global Vll = @variable(m, [i=1:nb_nph], start = Vll0[i] )
            global QgY = @variable(m, [i=1:ngy_inv], start = 0 )

            ## Define bounds on voltage magnitude and set slack bus angle to 0
            @constraint(m, Vd[idx_ref] .== [cos(0); cos(-2*pi/3); cos(2*pi/3)] )
            @constraint(m, Vq[idx_ref] .== [sin(0); sin(-2*pi/3); sin(2*pi/3)] )
            # global v_soft_min =  @variable(m,[i=1:nb_nph], start = 0 )
            # global v_soft_max =  @variable(m,[i=1:nb_nph],start = 0 )
            # @constraint(m,  v_soft_min .>= 0 ); @constraint(m,  v_soft_max .>= 0 )
            # global soft_Vlim = @expression(m, sum(v_soft_min[i] for i=1:nb_nph) + sum(v_soft_max[i] for i=1:nb_nph) )
            # Vlim_lin = -Vd0.^2-Vq0.^2+2*Vd0.*Vd+2*Vq0.*Vq;
            num_count = 1;
            for i=1:nb
                if  bus[i,2] == 3
                    num_count = num_count + 3
                end
                for j = 3:2:7
                    if bus[i,j] != 0  && bus[i,2] != 3
                        @constraint(m, bus[i,j]^2 <= Vd[num_count]^2+Vq[num_count]^2 <= bus[i,j+1]^2 )
                        # @constraint(m, v_soft_min[num_count]>= bus[i,j]^2- Vd[num_count]^2-Vq[num_count]^2 )
                        # @constraint(m, v_soft_max[num_count]>= Vd[num_count]^2+Vq[num_count]^2 - bus[i,j+1]^2 )
                        # @constraint(m, bus[i,j]^2 <= Vlim_lin[num_count]<= bus[i,j+1]^2 )
                        num_count = num_count + 1
                    end
                end
            end
            for i=1:nb_nph
                k = Int(Tl[i])
                @constraint(m, (Vd[i]-Vd[k])^2 + (Vq[i]-Vq[k])^2 == Vll[i]^2 )
            end

            ## Define bounds on reactive power
            if PV_en == 1
                num_count= 1;
                for i=1:ngy_inv
                    j = Int(ygen1[i,1]); count1 = 1
                    for k = 3:2:7
                        if bus[j,k] != 0
                            @constraint(m, ygen1[i,Qg_idx[count1]]/baseMVA <= QgY[num_count]<= ygen1[i,Qg_idx[count1]+1]/baseMVA )
                            num_count= num_count+ 1;
                        end
                        count1 = count1 + 1;
                    end
                end
            end

            ## Calculate overall star-connected loads in system
            for i=1:nl_y
                dum1 = (Sy[i,3] + im*Sy[i,4])/(1*exp(im*Sy[i,8]))*(Vm0[Int(Sy[i,7])]*exp(im*θ0[Int(Sy[i,7])]))
                YloadP[i] = Sy[i,1]+real(dum1)+Sy[i,5]*Vm0[Int(Sy[i,7])]^2
                YloadQ[i] = Sy[i,2]+imag(dum1)+Sy[i,6]*Vm0[Int(Sy[i,7])]^2
            end
            for i=nl_y+1:nly_nph
                YloadP[i] = Sy[i,1]+Sy[i,3]*Vm0[Int(Sy[i,7])]+Sy[i,5]*Vm0[Int(Sy[i,7])]^2
                YloadQ[i] = Sy[i,2]+Sy[i,4]*Vm0[Int(Sy[i,7])]+Sy[i,6]*Vm0[Int(Sy[i,7])]^2
            end
            Yload = CblY_mod*(-YloadP - im*YloadQ);

            ## Calculate overall delta-connected loads in system
            for i=1:nld_nph
                j = Int(Tl[Int(Sd[i,7])]);
                dum1 = Vm0[Int(Sd[i,7])]*exp(im*θ0[Int(Sd[i,7])])-Vm0[j]*exp(im*θ0[j])
                dum2 = (Sd[i,3] + im*Sd[i,4])/(sqrt(3)*exp(im*Sd[i,8]))*dum1
                DloadP[i] = Sd[i,1]+real(dum2)+Sd[i,5]*Vll0[Int(Sd[i,7])]^2/3
                DloadQ[i] = Sd[i,2]+imag(dum2)+Sd[i,6]*Vll0[Int(Sd[i,7])]^2/3
            end
            Dload = CblD_FP*(-DloadP - im*DloadQ)

            ## Power balance constraints
            V = Vm0[idx_noref].*exp.(im*θ0[idx_noref]);
            V_cdY = diagm(0=>conj(V)); V_cdD = diagm(0=>H*conj(V));
            dum1 = YLL_inv*pinv(V_cdY); dum1_d = real(dum1)*Cbg_mod; dum1_q = imag(dum1)*Cbg_mod;
            V_new1 = w + dum1*conj(Yload) + YLL_inv*H'*pinv(V_cdD)*conj(Dload)
            @constraint(m, real(V_new1) + dum1_d*PgY + dum1_q*QgY .== Vd[idx_noref] )
            @constraint(m, imag(V_new1) + dum1_q*PgY - dum1_d*QgY .== Vq[idx_noref] )

            ## Define cost
            cost_rect(Vll);

            ## Solve and update solution
            @time status = optimize!(m); stat = string(termination_status(m))
            V0 = Vm0.*exp.(im*θ0); Vd0 = JuMP.value.(Vd); Vq0 = JuMP.value.(Vq);
            V1 = Vd0+im*Vq0; opf_err = sum(abs.(V1-V0));
            err_vec[opf_iter] = round(opf_err, RoundNearest, digits=4);
            Vm0 = abs.(V1); θ0 = angle.(V1); QgY0 = JuMP.value.(QgY);
            print("FP_OPF_rect: (",opf_iter,") Error = ",opf_err,", ",stat,"\n");

            ## Check if regulator taps need to change
            # j = Int(tf_branch[vr_idx,1]); jj = Int(sum(bus_count[1:j-1])+1):Int(sum(bus_count[1:j]))
            # k = Int(tf_branch[vr_idx,2]); kk = Int(sum(bus_count[1:k-1])+1):Int(sum(bus_count[1:k]))
            # tap_new = [tf_branch[vr_idx,6], tf_branch[vr_idx,7], tf_branch[vr_idx,8]];
            # Vm_act = Vm0[kk]*bus_Vnom[j];
            # if VR_config[1,2] != 0
            #     Vm_act = Vm0[kk]*bus_Vnom[k]/VR_config[1,2];
            # end
            # if any(x->x>VR_config[1,4] || x<VR_config[1,3], Vm_act)
            #     for kkk = 1:length(kk)
            #         if Vm_act[kkk] < VR_config[1,3]
            #             tap_new[kkk] = tap_new[kkk]+1;
            #         elseif Vm_act[kkk] > VR_config[1,4]
            #             tap_new[kkk] = tap_new[kkk]-1;
            #         end
            #     end
            #     tap_old = tf_branch[vr_idx,6:8]*0.00625 .+ 1;
            #     yt = 1e5*(1+0*im)/(findmax(tap_old)[1])^2
            #     Y1_old = yt*Matrix{Float64}(I, 3, 3); tap_old = Diagonal(tap_old);
            #     tf_branch[vr_idx,6:8] = tap_new
            #     tap = tap_new*0.00625 .+ 1; yt = 1e5*(1+0*im)/(findmax(tap)[1])^2
            #     Y1_new = yt*Matrix{Float64}(I, 3, 3); global tap_new = Diagonal(tap);
            #     G[jj,jj] = G[jj,jj]+real(Y1_new).*tap_new^2-real(Y1_old).*tap_old^2;
            #     B[jj,jj] = B[jj,jj]+imag(Y1_new).*tap_new^2-imag(Y1_old).*tap_old^2;
            #     G[kk,kk] = G[kk,kk]+real(Y1_new)-real(Y1_old);
            #     B[kk,kk] = B[kk,kk]+imag(Y1_new)-imag(Y1_old);
            #     G[jj,kk] = real(-Y1_new).*tap_new;
            #     B[jj,kk] = imag(-Y1_new).*tap_new;
            #     G[kk,jj] = real(-Y1_new).*tap_new;
            #     B[kk,jj] = imag(-Y1_new).*tap_new;
            #     Y = G+im*B; YLL = [Y[idx_noref1,idx_noref1] Y[idx_noref1,idx_noref2]; Y[idx_noref2,idx_noref1] Y[idx_noref2,idx_noref2] ];
            #     global YLL_inv = pinv(YLL); global YL0 = [Y[idx_ref,idx_noref1] Y[idx_ref,idx_noref2]]';
            #     global w = -YLL_inv*YL0*(Vm0[idx_ref].*exp.(im*θ0[idx_ref]));
            # end
        end
    end

    ## Run power flow to convergence if OPF did not converge
    if opf_iter == max_iter
        println("\n****************OPF DID NOT CONVERGE!!*******************\n")
        (Vm0,θ0,QgY0,pf_iter) = FP_pf(Vm0,θ0,QgY0,50); print("No of PF iterations: ", pf_iter,"\n");
        if dnc == 0
            # stat = "DID NOT CONVERGE";
        end
    end

    ## Calculate missing solutions
    Vph = Vm0.*exp.(im*θ0);
    Iph = (G+im*B)*Vph; S_inj = Vph.*conj(Iph); P_inj = real(S_inj);
    for i=1:nb_nph
        k = Int(Tl[i])
        Vll0[i] =  sqrt(Vm0[i]^2 + Vm0[k]^2 - 2*Vm0[i]*Vm0[k]*cos(θ0[i]-θ0[k]) )
    end
    print_VU(Vm0,θ0,Vll0);
    global P_loss = sum(P_inj)*baseMVA*1e3
    println("TOTAL TIME:")
    return Vm0,θ0,QgY0,stat
end
