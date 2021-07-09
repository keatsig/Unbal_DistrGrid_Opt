### FUNCTION TO CALCULATE OPTIMAL POWER FLOW (POLAR-CONVENTIONAL)
function TPOPF_act_constr(crit_node,Vm0,θ0,QgY0)

    ## Initialize model
    start_point(4); Vm0 = s0[1:nb_nph]; θ0 = s0[1+nb_nph:2*nb_nph];
    (Vm0,θ0,QgY0,pf_iter) = FP_pf(Vm0,θ0,QgY0,30);
    Vll0 = zeros(nb_nph)
    for i=1:nb_nph
        k = Int(Tl[i])
        Vll0[i] =  sqrt(Vm0[i]^2 + Vm0[k]^2 - 2*Vm0[i]*Vm0[k]*cos(θ0[i]-θ0[k]) )
    end

    ## Define active power of inverters
    Pg_idx = [2 6 10];  Qg_idx = [4 8 12]; SgY = yinv_ratedS/baseMVA;
    PgY = spzeros(ngy_inv); Qmax = spzeros(ngy_inv);  num_count= 1;
    if PV_en == 1
        for i=1:ngy_inv
            j = Int(ygen1[i,1]); count1 = 1
            for k = 3:2:7
                if bus[j,k] != 0
                    PgY[num_count]= ygen1[i,Pg_idx[count1]]/baseMVA
                    Qmax[num_count]= ygen1[i,Qg_idx[count1]+1]/baseMVA
                    num_count= num_count+ 1;
                end
                count1 = count1 + 1;
            end
        end
    end

    ## Generate load matrix
    ZIP_load()

    ## Start iterative linearized OPF
    max_iter = 20; opf_iter = 0; opf_err = 1; stat = "LOCALLY_SOLVED";
    while opf_err != 0 && opf_iter<max_iter && (stat == "LOCALLY_SOLVED" || stat == "ALMOST_LOCALLY_SOLVED")

        opf_iter = opf_iter + 1;

        ## Check for constraint violations
        global Vm_idx = findall(idx -> idx < 0.9 || idx > 1.1, Vm0)
        QgY_idx1 = zeros(ngy_inv);
        for i=1:ngy_inv
            if PgY[i]^2 + QgY0[i]^2 > SgY[i]^2
                QgY_idx1[i] = i;
            end
        end
        global QgY_idx = Int.(QgY_idx1[QgY_idx1 .!= 0])

        ## Define variables
        global m = Model(optimizer_with_attributes(Ipopt.Optimizer,"warm_start_init_point"=>"yes","print_level"=>0, "tol"=>1e-7));
        global Vm = @variable(m, [i=1:nb_nph], start = Vm0[i] )
        global θ = @variable(m,  [i=1:nb_nph], start = θ0[i] )
        global Vll = @variable(m, [i=1:nb_nph], start = Vll0[i] )
        global QgY = @variable(m, [i=1:ngy_inv], start = QgY0[i] )

        ## Define bounds on voltage magnitude and set slack bus angle to 0
        @constraint(m, Vm[idx_ref] .== 1 )
        @constraint(m, θ[idx_ref] .== [0;-2*pi/3;2*pi/3] )
        for i=1:nb_nph
            k = Int(Tl[i])
            @NLconstraint(m, Vm[i]^2 + Vm[k]^2 - 2*Vm[i]*Vm[k]*cos(θ[i]-θ[k]) == Vll[i]^2 )
        end
        @constraint(m, 0.9 .<= Vm[Vm_idx] .<= 1.1 )

        ## Define bounds on active and reactive power
        Qg_idx = [4 8 12]; num_count= 1;
        if PV_en == 1
            @constraint(m, -Qmax .<= QgY .<= Qmax )
        end

        ## Calculate overall star-connected loads in system
        dumY_d = @NLexpression(m, [i=1:nl_y],  cos(Sy[i,8])*(Sy[i,3]*Vm[Int(Sy[i,7])]*cos(θ[Int(Sy[i,7])]) - Sy[i,4]*Vm[Int(Sy[i,7])]*sin(θ[Int(Sy[i,7])])) + sin(Sy[i,8])*(Sy[i,4]*Vm[Int(Sy[i,7])]*cos(θ[Int(Sy[i,7])]) + Sy[i,3]*Vm[Int(Sy[i,7])]*sin(θ[Int(Sy[i,7])]))  )
        dumY_q = @NLexpression(m, [i=1:nl_y], -sin(Sy[i,8])*(Sy[i,3]*Vm[Int(Sy[i,7])]*cos(θ[Int(Sy[i,7])]) - Sy[i,4]*Vm[Int(Sy[i,7])]*sin(θ[Int(Sy[i,7])])) + cos(Sy[i,8])*(Sy[i,4]*Vm[Int(Sy[i,7])]*cos(θ[Int(Sy[i,7])]) + Sy[i,3]*Vm[Int(Sy[i,7])]*sin(θ[Int(Sy[i,7])]))  )
        YloadP1 = @NLexpression(m, [i=1:nl_y], Sy[i,1]+dumY_d[i]+Sy[i,5]*Vm[Int(Sy[i,7])]^2 )
        YloadQ1 = @NLexpression(m, [i=1:nl_y], Sy[i,2]+dumY_q[i]+Sy[i,6]*Vm[Int(Sy[i,7])]^2 )
        YloadP2 = @NLexpression(m, [i=1:nlh], Sy[i+nl_y,1]+Sy[i+nl_y,3]*Vm[Int(Sy[i+nl_y,7])]+Sy[i+nl_y,5]*Vm[Int(Sy[i+nl_y,7])]^2 )
        YloadQ2 = @NLexpression(m, [i=1:nlh], Sy[i+nl_y,2]+Sy[i+nl_y,4]*Vm[Int(Sy[i+nl_y,7])]+Sy[i+nl_y,6]*Vm[Int(Sy[i+nl_y,7])]^2 )
        YloadP = [YloadP1; YloadP2]; YloadQ = [YloadQ1; YloadQ2];

        ## Calculate overall delta-connected loads in system
        dum_Vll_d = @NLexpression(m, [i=1:nld_nph], Vm[Int(Sd[i,7])]*cos(θ[Int(Sd[i,7])]) - Vm[Int(Tl[Int(Sd[i,7])])]*cos(θ[Int(Tl[Int(Sd[i,7])])]) )
        dum_Vll_q = @NLexpression(m, [i=1:nld_nph], Vm[Int(Sd[i,7])]*sin(θ[Int(Sd[i,7])]) - Vm[Int(Tl[Int(Sd[i,7])])]*sin(θ[Int(Tl[Int(Sd[i,7])])]) )
        dumD_d = @NLexpression(m, [i=1:nld_nph],  cos(Sd[i,8])*(Sd[i,3]*dum_Vll_d[i] - Sd[i,4]*dum_Vll_q[i]) + sin(Sd[i,8])*(Sd[i,4]*dum_Vll_d[i] + Sd[i,3]*dum_Vll_q[i])  )
        dumD_q = @NLexpression(m, [i=1:nld_nph], -sin(Sd[i,8])*(Sd[i,3]*dum_Vll_d[i] - Sd[i,4]*dum_Vll_q[i]) + cos(Sd[i,8])*(Sd[i,4]*dum_Vll_d[i] + Sd[i,3]*dum_Vll_q[i])  )
        DloadP = @NLexpression(m, [i=1:nld_nph], Sd[i,1]+dumD_d[i]/sqrt(3)+Sd[i,5]*Vll[Int(Sd[i,7])]^2/3 )
        DloadQ = @NLexpression(m, [i=1:nld_nph], Sd[i,2]+dumD_q[i]/sqrt(3)+Sd[i,6]*Vll[Int(Sd[i,7])]^2/3 )

        ## Convert delta to equivalent wye-connected loads
        @variable(m, Ill_d[i=1:nld_nph] ); @variable(m, Ill_q[i=1:nld_nph] )
        @variable(m, P_d2y[i=1:nldy_nph] ); @variable(m, Q_d2y[i=1:nldy_nph] )
        num_count= 1
        for i=1:nld
            j = Int(dload[i,1])
            for k=Int.(sum(bus_count[1:j-1])+1:sum(bus_count[1:j]))
                kk = Int(Tl[k])
                if k != kk
                    @NLconstraint(m, DloadP[num_count]==  (Vm[k]*cos(θ[k])-Vm[kk]*cos(θ[kk]))*Ill_d[num_count]+ (Vm[k]*sin(θ[k])-Vm[kk]*sin(θ[kk]))*Ill_q[num_count])
                    @NLconstraint(m, DloadQ[num_count]== -(Vm[k]*cos(θ[k])-Vm[kk]*cos(θ[kk]))*Ill_q[num_count]+ (Vm[k]*sin(θ[k])-Vm[kk]*sin(θ[kk]))*Ill_d[num_count])
                    num_count= num_count+ 1
                end
            end
        end
        Ip_d = @NLexpression(m, [i=1:nldy_nph], sum(Ill_d[j]*Tl2p[i,j] for j=1:nld_nph) )
        Ip_q = @NLexpression(m, [i=1:nldy_nph], sum(Ill_q[j]*Tl2p[i,j] for j=1:nld_nph) )
        num_count= 1
        for i=1:nld
            j = Int(dload[i,1])
            for k=Int.(sum(bus_count[1:j-1])+1:sum(bus_count[1:j]))
                @NLconstraint(m, P_d2y[num_count]==  Vm[k]*cos(θ[k])*Ip_d[num_count]+ Vm[k]*sin(θ[k])*Ip_q[num_count])
                @NLconstraint(m, Q_d2y[num_count]== -Vm[k]*cos(θ[k])*Ip_q[num_count]+ Vm[k]*sin(θ[k])*Ip_d[num_count])
                num_count= num_count+ 1
            end
        end

        ## Constraints on flow conservation (generation - demand = losses)
        Pij = @NLexpression(m, [i=1:nb_nph], sum(Vm[i]*Vm[j]*(G[i,j]*cos(θ[i]-θ[j])+B[i,j]*sin(θ[i]-θ[j])) for j=1:nb_nph) )
        Qij = @NLexpression(m, [i=1:nb_nph], sum(Vm[i]*Vm[j]*(G[i,j]*sin(θ[i]-θ[j])-B[i,j]*cos(θ[i]-θ[j])) for j=1:nb_nph) )
        for i in idx_noref
            @NLconstraint(m, sum(Cbg[i,k]*PgY[k] for k=1:ngy_inv) - sum(CblY[i,kk]*YloadP[kk] for kk=1:nly_nph) - sum(CblD[i,kkk]*P_d2y[kkk] for kkk=1:nldy_nph) == Pij[i] )
            @NLconstraint(m, sum(Cbg[i,k]*QgY[k] for k=1:ngy_inv) - sum(CblY[i,kk]*YloadQ[kk] for kk=1:nly_nph) - sum(CblD[i,kkk]*Q_d2y[kkk] for kkk=1:nldy_nph) == Qij[i] )
        end

        ## Define cost
        cost_polar(Vll);

        ## Solve and update solution
        @time status = optimize!(m); stat = string(termination_status(m))
        Vm0 = JuMP.value.(Vm); θ0 = JuMP.value.(θ); QgY0 = JuMP.value.(QgY); Vll0 = JuMP.value.(Vll);
        opf_err = length(Vm_idx) + length(QgY_idx)
        print("OPF_pol: (",opf_iter,") Violations = ",opf_err,", ",stat,"\n");
    end
    print_VU(Vm0,θ0,Vll0)
    return Vm0,θ0,QgY0,stat
end

#----------------------------------------------------------------------------------------#
