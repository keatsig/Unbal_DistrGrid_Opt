### FUNCTION TO CALCULATE OPTIMAL POWER FLOW (POLAR-CONVENTIONAL)
function TPOPF_pol(crit_node,Vm0,θ0,QgY0)

    ## Initialize model
    # start_point(1); Vm0 = s0[1:nb_nph]; θ0 = s0[1+nb_nph:2*nb_nph];
    # @time (Vm0,θ0,QgY0,stat) = FOT_opf_pol(crit_node,Vm0,θ0,QgY0,0)
    # (Vm0,θ0,QgY0,pf_iter) = FP_pf(Vm0,θ0,QgY0,50);
    # @time (Vm0,θ0,QgY0,stat) = FP_opf_rect(crit_node,Vm0,θ0,QgY0,0)
    # @time (Vm0,θ0,QgY0,stat) = BFS_opf_rect(crit_node,Vm0,θ0,QgY0,0)
    # @time (Vm0,θ0,QgY0,stat) = D3F_opf_rect(crit_node,Vm0,θ0,QgY0,0)
    Vll0 = zeros(nb_nph);
    for i=1:nb_nph
        k = Int(Tl[i])
        Vll0[i] =  sqrt(Vm0[i]^2 + Vm0[k]^2 - 2*Vm0[i]*Vm0[k]*cos(θ0[i]-θ0[k]) )
    end

    ## Generate load matrix
    ZIP_load()

    ## Calculate overall delta-connected loads in system
    DloadP0 = zeros(nld_nph); DloadQ0 = zeros(nld_nph); I_ll0 = zeros(nld_nph)*(0+im*0)
    for i=1:nld_nph
        j = Int(Tl[Int(Sd[i,7])]);
        dum1 = Vm0[Int(Sd[i,7])]*exp(im*θ0[Int(Sd[i,7])])-Vm0[j]*exp(im*θ0[j])
        dum2 = (Sd[i,3] + im*Sd[i,4])/(sqrt(3)*exp(im*Sd[i,8]))*dum1
        DloadP0[i] = Sd[i,1]+real(dum2)+Sd[i,5]*Vll0[Int(Sd[i,7])]^2/3
        DloadQ0[i] = Sd[i,2]+imag(dum2)+Sd[i,6]*Vll0[Int(Sd[i,7])]^2/3
    end
    num_count= 1;
    for i=1:nld
        j = Int(dload[i,1])
        for k=Int.(sum(bus_count[1:j-1])+1:sum(bus_count[1:j]))
            kk = Int(Tl[k])
            if k != kk
                I_ll0[num_count]= (DloadP0[num_count]-im*DloadQ0[num_count])/(Vm0[k]*exp(-im*θ0[k])-Vm0[kk]*exp(-im*θ0[kk]));
                num_count= num_count+ 1
            end
        end
    end
    Ill_d0 = real(I_ll0); Ill_q0 = imag(I_ll0);
    Ip0 = Tl2p*I_ll0;
    num_count= 1; P_d2y0 = zeros(nldy_nph); Q_d2y0 = zeros(nldy_nph);
    for i=1:nld
        j = Int(dload[i,1])
        for k=Int.(sum(bus_count[1:j-1])+1:sum(bus_count[1:j]))
            dum1 = Vm0[k]*exp(im*θ0[k])*conj(Ip0[num_count])
            P_d2y0[num_count]= real(dum1); Q_d2y0[num_count]= imag(dum1)
            num_count= num_count+ 1
        end
    end

    ## Define variables
    # global m = Model(optimizer_with_attributes(Ipopt.Optimizer,"warm_start_init_point"=>"yes","print_level"=>0, "tol"=>1e-7));
    global m = Model(optimizer_with_attributes(KNITRO.Optimizer,"honorbnds" => 1, "outlev" => 0, "presolve" => 0, "strat_warm_start" => 1, "algorithm" => 0));
    global Vm = @variable(m, [i=1:nb_nph], start = Vm0[i] )
    global θ = @variable(m,  [i=1:nb_nph], start = θ0[i] )
    global Vll = @variable(m, [i=1:nb_nph], start = Vll0[i] )
    global QgY = @variable(m, [i=1:ngy_inv], start = QgY0[i] )

    ## Define bounds on voltage magnitude and set slack bus angle to 0
    @constraint(m, Vm[idx_ref] .== 1 )
    @constraint(m, θ[idx_ref] .== [0;-2*pi/3;2*pi/3] )
    num_count= 1
    for i=1:nb
        if  bus[i,2] == 3
            num_count= num_count+ 3
        end
        for j = 3:2:7
            if bus[i,j] != 0  && bus[i,2] != 3
                @constraint(m, bus[i,j] <= Vm[num_count]<= bus[i,j+1] )
                num_count= num_count+ 1
            end
        end
    end
    for i=1:nb_nph
        k = Int(Tl[i])
        @NLconstraint(m, Vm[i]^2 + Vm[k]^2 - 2*Vm[i]*Vm[k]*cos(θ[i]-θ[k]) == Vll[i]^2 )
    end

    ## Define bounds on active and reactive power
    Pg_idx = [2 6 10];  Qg_idx = [4 8 12];
    PgY = spzeros(ngy_inv); num_count= 1;
    if PV_en == 1
        for i=1:ngy_inv
            j = Int(ygen1[i,1]); count1 = 1
            for k = 3:2:7
                if bus[j,k] != 0
                    PgY[num_count]= ygen1[i,Pg_idx[count1]]/baseMVA
                    @constraint(m, ygen1[i,Qg_idx[count1]]/baseMVA <= QgY[num_count]<= ygen1[i,Qg_idx[count1]+1]/baseMVA )
                    num_count= num_count+ 1;
                end
                count1 = count1 + 1;
            end
        end
    end

    ## Calculate overall star-connected loads in system
    dumY_d = @NLexpression(m, [i=1:nl_y],  cos(Sy[i,8])*(Sy[i,3]*Vm[Int(Sy[i,7])]*cos(θ[Int(Sy[i,7])]) - Sy[i,4]*Vm[Int(Sy[i,7])]*sin(θ[Int(Sy[i,7])])) + sin(Sy[i,8])*(Sy[i,4]*Vm[Int(Sy[i,7])]*cos(θ[Int(Sy[i,7])]) + Sy[i,3]*Vm[Int(Sy[i,7])]*sin(θ[Int(Sy[i,7])]))  )
    dumY_q = @NLexpression(m, [i=1:nl_y], -sin(Sy[i,8])*(Sy[i,3]*Vm[Int(Sy[i,7])]*cos(θ[Int(Sy[i,7])]) - Sy[i,4]*Vm[Int(Sy[i,7])]*sin(θ[Int(Sy[i,7])])) + cos(Sy[i,8])*(Sy[i,4]*Vm[Int(Sy[i,7])]*cos(θ[Int(Sy[i,7])]) + Sy[i,3]*Vm[Int(Sy[i,7])]*sin(θ[Int(Sy[i,7])]))  )
    # dumY_d = @NLexpression(m, [i=1:nl_y], Vm[Int(Sy[i,7])]*Sy[i,9]*cos(θ[Int(Sy[i,7])]+Sy[i,10]-Sy[i,8])  )
    # dumY_q = @NLexpression(m, [i=1:nl_y], Vm[Int(Sy[i,7])]*Sy[i,9]*sin(θ[Int(Sy[i,7])]+Sy[i,10]-Sy[i,8])  )
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
    # DloadP = @NLexpression(m, [i=1:nld_nph], Sd[i,1]+Sd[i,3]*Vm[Int(Sd[i,7])]+Sd[i,5]*Vm[Int(Sd[i,7])]^2 )
    # DloadQ = @NLexpression(m, [i=1:nld_nph], Sd[i,2]+Sd[i,4]*Vm[Int(Sd[i,7])]+Sd[i,6]*Vm[Int(Sd[i,7])]^2 )

    ## Convert delta to equivalent wye-connected loads
    @variable(m, Ill_d[i=1:nld_nph], start = Ill_d0[i] ); @variable(m, Ill_q[i=1:nld_nph], start = Ill_q0[i] )
    @variable(m, P_d2y[i=1:nldy_nph], start = P_d2y0[i] ); @variable(m, Q_d2y[i=1:nldy_nph], start = Q_d2y0[i] )
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

    ## Delta-connected sources
    # global Pgd = @variable(m, [i=1:ngd_nph], start = 0 )
    # global Qgd = @variable(m, [i=1:ngd_nph], start = 0 )
    # num_count= 1;
    # for i=1:ngd
    #     for k = 3:2:7
    #         @constraint(m, dgen[i,Pg_idx[num_count]]/baseMVA <= Pgd[num_count]<= ygen[i,Pg_idx[num_count]+1]/baseMVA )
    #         @constraint(m, dgen[i,Qg_idx[num_count]]/baseMVA <= Qgd[num_count]<= ygen[i,Qg_idx[num_count]+1]/baseMVA )
    #         num_count= num_count+ 1;
    #     end
    # end
    # @variable(m, Igen_ll_d[i=1:ngd_nph] );  @variable(m, Igen_ll_q[i=1:ngd_nph] )
    # @variable(m, Pg_d2y[i=1:ngdy_nph] );  @variable(m, Qg_d2y[i=1:ngdy_nph] )
    # num_count= 1
    # for i=1:ngd
    #     j = Int(dgen[i,1])
    #     for k=Int.(sum(bus_count[1:j-1])+1:sum(bus_count[1:j]))
    #         kk = Int(Tl[k])
    #         if k != kk
    #             @NLconstraint(m, Pgd[num_count]==  (Vm[k]*cos(θ[k])-Vm[kk]*cos(θ[kk]))*Igen_ll_d[num_count]+ (Vm[k]*sin(θ[k])-Vm[kk]*sin(θ[kk]))*Igen_ll_q[num_count])
    #             @NLconstraint(m, Qgd[num_count]== -(Vm[k]*cos(θ[k])-Vm[kk]*cos(θ[kk]))*Igen_ll_q[num_count]+ (Vm[k]*sin(θ[k])-Vm[kk]*sin(θ[kk]))*Igen_ll_d[num_count])
    #             num_count= num_count+ 1
    #         end
    #     end
    # end
    # Ipgen_d = @NLexpression(m, [i=1:ngdy_nph], sum(Igen_ll_d[j]*Tl2p[i,j] for j=1:ngd_nph) )
    # Ipgen_q = @NLexpression(m, [i=1:ngdy_nph], sum(Igen_ll_q[j]*Tl2p[i,j] for j=1:ngd_nph) )
    # num_count= 1
    # for i=1:ngd
    #     j = Int(dgen[i,1])
    #     for k=Int.(sum(bus_count[1:j-1])+1:sum(bus_count[1:j]))
    #         @NLconstraint(m, Pg_d2y[num_count]==  Vm[k]*cos(θ[k])*Ipgen_d[num_count]+ Vm[k]*sin(θ[k])*Ipgen_q[num_count])
    #         @NLconstraint(m, Qg_d2y[num_count]== -Vm[k]*cos(θ[k])*Ipgen_q[num_count]+ Vm[k]*sin(θ[k])*Ipgen_d[num_count])
    #         num_count= num_count+ 1
    #     end
    # end

    ## Constraints on flow conservation (generation - demand = losses)
    Pij = @NLexpression(m, [i=1:nb_nph], sum(Vm[i]*Vm[j]*(G[i,j]*cos(θ[i]-θ[j])+B[i,j]*sin(θ[i]-θ[j])) for j=1:nb_nph) )
    Qij = @NLexpression(m, [i=1:nb_nph], sum(Vm[i]*Vm[j]*(G[i,j]*sin(θ[i]-θ[j])-B[i,j]*cos(θ[i]-θ[j])) for j=1:nb_nph) )
    for i in idx_noref
        @NLconstraint(m, sum(Cbg[i,k]*PgY[k] for k=1:ngy_inv) - sum(CblY[i,kk]*YloadP[kk] for kk=1:nly_nph) - sum(CblD[i,kkk]*P_d2y[kkk] for kkk=1:nldy_nph) == Pij[i] )
        @NLconstraint(m, sum(Cbg[i,k]*QgY[k] for k=1:ngy_inv) - sum(CblY[i,kk]*YloadQ[kk] for kk=1:nly_nph) - sum(CblD[i,kkk]*Q_d2y[kkk] for kkk=1:nldy_nph) == Qij[i] )
    end

    ## Define cost
    cost_polar(Vll);

    ## Solve and print solution
    @time status = optimize!(m); println(termination_status(m));
    stat = string(termination_status(m))
    Vm0 = JuMP.value.(Vm); θ0 = JuMP.value.(θ); QgY0 = JuMP.value.(QgY); Vll0 = JuMP.value.(Vll);
    print_VU(Vm0,θ0,Vll0)
    return Vm0,θ0,QgY0,stat
end

#----------------------------------------------------------------------------------------#

### FUNCTION TO CALCULATE OPTIMAL POWER FLOW (RECTANGULAR-CONVENTIONAL)
function TPOPF_rect(crit_node,Vd0,Vq0,QgY0)

    ## Initialize model
    start_point(4); Vm0 = s0[1:nb_nph]; θ0 = s0[1+nb_nph:2*nb_nph];
    # @time (Vm0,θ0,QgY0,stat) = FOT_opf_pol(crit_node,Vm0,θ0,QgY0,0)
    # (Vm0,θ0,QgY0,pf_iter) = FP_pf(Vm0,θ0,QgY0,50);
    # @time (Vm0,θ0,QgY0,stat) = FP_opf_rect(crit_node,Vm0,θ0,QgY0,0)
    @time (Vm0,θ0,QgY0,stat) = BFS_opf_rect(crit_node,Vm0,θ0,QgY0,0)
    # @time (Vm0,θ0,QgY0,stat) = D3F_opf_rect(crit_node,Vm0,θ0,QgY0,0)
    Vd0 = Vm0.*cos.(θ0); Vq0 = Vm0.*sin.(θ0); Vll0 = zeros(nb_nph)
    for i=1:nb_nph
        k = Int(Tl[i])
        Vll0[i] =  sqrt(Vm0[i]^2 + Vm0[k]^2 - 2*Vm0[i]*Vm0[k]*cos(θ0[i]-θ0[k]) )
    end

    ## Calculate overall delta-connected loads in system
    DloadP0 = zeros(nld_nph); DloadQ0 = zeros(nld_nph); I_ll0 = zeros(nld_nph)*(0+im*0)
    for i=1:nld_nph
        j = Int(Tl[Int(Sd[i,7])]);
        dum1 = Vm0[Int(Sd[i,7])]*exp(im*θ0[Int(Sd[i,7])])-Vm0[j]*exp(im*θ0[j])
        dum2 = (Sd[i,3] + im*Sd[i,4])/(sqrt(3)*exp(im*Sd[i,8]))*dum1
        DloadP0[i] = Sd[i,1]+real(dum2)+Sd[i,5]*Vll0[Int(Sd[i,7])]^2/3
        DloadQ0[i] = Sd[i,2]+imag(dum2)+Sd[i,6]*Vll0[Int(Sd[i,7])]^2/3
    end
    num_count= 1;
    for i=1:nld
        j = Int(dload[i,1])
        for k=Int.(sum(bus_count[1:j-1])+1:sum(bus_count[1:j]))
            kk = Int(Tl[k])
            if k != kk
                I_ll0[num_count]= (DloadP0[num_count]-im*DloadQ0[num_count])/(Vm0[k]*exp(-im*θ0[k])-Vm0[kk]*exp(-im*θ0[kk]));
                num_count= num_count+ 1
            end
        end
    end
    Ill_d0 = real(I_ll0); Ill_q0 = imag(I_ll0);
    Ip0 = Tl2p*I_ll0;
    num_count= 1; P_d2y0 = zeros(nldy_nph); Q_d2y0 = zeros(nldy_nph);
    for i=1:nld
        j = Int(dload[i,1])
        for k=Int.(sum(bus_count[1:j-1])+1:sum(bus_count[1:j]))
            dum1 = Vm0[k]*exp(im*θ0[k])*conj(Ip0[num_count])
            P_d2y0[num_count]= real(dum1); Q_d2y0[num_count]= imag(dum1)
            num_count= num_count+ 1
        end
    end

    ## Define voltage magnitude and angle variables
    global m = Model(optimizer_with_attributes(Ipopt.Optimizer,"warm_start_init_point"=>"yes","print_level"=>5, "tol"=>1e-7));
    global Vd = @variable(m, [i=1:nb_nph], start = Vd0[i] )
    global Vq = @variable(m, [i=1:nb_nph], start = Vq0[i] )
    global Vll = @variable(m, [i=1:nb_nph], start = Vll0[i] )
    global QgY = @variable(m, [i=1:ngy_inv], start = QgY0[i] )

    ## Define bounds on voltage magnitude and set slack bus angle to 0
    global Vm = @NLexpression(m, [i=1:nb_nph], sqrt(Vd[i]^2+Vq[i]^2) )
    @constraint(m, Vd[idx_ref] .== [cos(0); cos(-2*pi/3); cos(2*pi/3)] )
    @constraint(m, Vq[idx_ref] .== [sin(0); sin(-2*pi/3); sin(2*pi/3)] )
    num_count = 1
    for i=1:nb
        if  bus[i,2] == 3
            num_count = num_count + 3
        end
        for j = 3:2:7
            if bus[i,j] != 0  && bus[i,2] != 3
                @constraint(m, bus[i,j]^2 <= Vd[num_count]^2+Vq[num_count]^2 <= bus[i,j+1]^2 )
                num_count = num_count + 1
            end
        end
    end
    global Vll = @variable(m, [i=1:nb_nph], start = Vll0[i] )
    for i=1:nb_nph
        k = Int(Tl[i])
        @NLconstraint(m, (Vd[i]-Vd[k])^2 + (Vq[i]-Vq[k])^2 == Vll[i]^2 )
    end

   ## Define bounds on active and reactive power
   Pg_idx = [2 6 10];  Qg_idx = [4 8 12];
   PgY = spzeros(ngy_inv); num_count = 1;
   if PV_en == 1
       for i=1:ngy_inv
           j = Int(ygen1[i,1]); count1 = 1
           for k = 3:2:7
               if bus[j,k] != 0
                   PgY[num_count] = ygen1[i,Pg_idx[count1]]/baseMVA
                   @constraint(m, ygen1[i,Qg_idx[count1]]/baseMVA <= QgY[num_count] <= ygen1[i,Qg_idx[count1]+1]/baseMVA )
                   num_count = num_count + 1;
               end
               count1 = count1 + 1;
           end
       end
   end

   ## Generate load matrix
   ZIP_load()

   ## Calculate overall star-connected loads in system
   dumY_d = @NLexpression(m, [i=1:nl_y],  cos(Sy[i,8])*(Sy[i,3]*Vd[Int(Sy[i,7])] - Sy[i,4]*Vq[Int(Sy[i,7])]) + sin(Sy[i,8])*(Sy[i,4]*Vd[Int(Sy[i,7])] + Sy[i,3]*Vq[Int(Sy[i,7])])  )
   dumY_q = @NLexpression(m, [i=1:nl_y], -sin(Sy[i,8])*(Sy[i,3]*Vd[Int(Sy[i,7])] - Sy[i,4]*Vq[Int(Sy[i,7])]) + cos(Sy[i,8])*(Sy[i,4]*Vd[Int(Sy[i,7])] + Sy[i,3]*Vq[Int(Sy[i,7])])  )
   YloadP1 = @NLexpression(m, [i=1:nl_y], Sy[i,1]+dumY_d[i]+Sy[i,5]*Vm[Int(Sy[i,7])]^2 )
   YloadQ1 = @NLexpression(m, [i=1:nl_y], Sy[i,2]+dumY_q[i]+Sy[i,6]*Vm[Int(Sy[i,7])]^2 )
   YloadP2 = @NLexpression(m, [i=1:nlh], Sy[i+nl_y,1]+Sy[i+nl_y,3]*Vm[Int(Sy[i+nl_y,7])]+Sy[i+nl_y,5]*Vm[Int(Sy[i+nl_y,7])]^2 )
   YloadQ2 = @NLexpression(m, [i=1:nlh], Sy[i+nl_y,2]+Sy[i+nl_y,4]*Vm[Int(Sy[i+nl_y,7])]+Sy[i+nl_y,6]*Vm[Int(Sy[i+nl_y,7])]^2 )
   YloadP = [YloadP1; YloadP2]; YloadQ = [YloadQ1; YloadQ2];

   ## Calculate overall delta-connected loads in system
   dum_Vll_d = @NLexpression(m, [i=1:nld_nph], Vd[Int(Sd[i,7])] - Vd[Int(Tl[Int(Sd[i,7])])] )
   dum_Vll_q = @NLexpression(m, [i=1:nld_nph], Vq[Int(Sd[i,7])] - Vq[Int(Tl[Int(Sd[i,7])])] )
   dumD_d = @NLexpression(m, [i=1:nld_nph],  cos(Sd[i,8])*(Sd[i,3]*dum_Vll_d[i] - Sd[i,4]*dum_Vll_q[i]) + sin(Sd[i,8])*(Sd[i,4]*dum_Vll_d[i] + Sd[i,3]*dum_Vll_q[i])  )
   dumD_q = @NLexpression(m, [i=1:nld_nph], -sin(Sd[i,8])*(Sd[i,3]*dum_Vll_d[i] - Sd[i,4]*dum_Vll_q[i]) + cos(Sd[i,8])*(Sd[i,4]*dum_Vll_d[i] + Sd[i,3]*dum_Vll_q[i])  )
   DloadP = @NLexpression(m, [i=1:nld_nph], Sd[i,1]+dumD_d[i]/sqrt(3)+Sd[i,5]*Vll[Int(Sd[i,7])]^2/3 )
   DloadQ = @NLexpression(m, [i=1:nld_nph], Sd[i,2]+dumD_q[i]/sqrt(3)+Sd[i,6]*Vll[Int(Sd[i,7])]^2/3 )

    ## Convert delta to equivalent wye-connected loads
    @variable(m, Ill_d[i=1:nld_nph], start = Ill_d0[i] ); @variable(m, Ill_q[i=1:nld_nph], start = Ill_q0[i] )
    global P_d2y = @variable(m, [i=1:nldy_nph], start = P_d2y0[i] );  @variable(m, Q_d2y[i=1:nldy_nph], start = Q_d2y0[i]  )
    num_count = 1
    for i=1:nld
        j = Int(dload[i,1])
        for k=Int.(sum(bus_count[1:j-1])+1:sum(bus_count[1:j]))
            kk = Int(Tl[k])
            if k != kk
                @NLconstraint(m, DloadP[num_count] ==  (Vd[k]-Vd[kk])*Ill_d[num_count] + (Vq[k]-Vq[kk])*Ill_q[num_count] )
                @NLconstraint(m, DloadQ[num_count] == -(Vd[k]-Vd[kk])*Ill_q[num_count] + (Vq[k]-Vq[kk])*Ill_d[num_count] )
                num_count = num_count + 1
            end
        end
    end
    Ip_d = @NLexpression(m, [i=1:nldy_nph], sum(Ill_d[j]*Tl2p[i,j] for j=1:nld_nph) )
    Ip_q = @NLexpression(m, [i=1:nldy_nph], sum(Ill_q[j]*Tl2p[i,j] for j=1:nld_nph) )
    num_count = 1
    for i=1:nld
        j = Int(dload[i,1])
        for k=Int.(sum(bus_count[1:j-1])+1:sum(bus_count[1:j]))
            @NLconstraint(m, P_d2y[num_count] ==  Vd[k]*Ip_d[num_count] + Vq[k]*Ip_q[num_count] )
            @NLconstraint(m, Q_d2y[num_count] == -Vd[k]*Ip_q[num_count] + Vq[k]*Ip_d[num_count] )
            num_count = num_count + 1
        end
    end

    ## Constraints on flow conservation (generation - demand = losses)
    Pij = @NLexpression(m, [i=1:nb_nph], sum(Vd[i]*( G[i,j]*Vd[j]-B[i,j]*Vq[j]) + Vq[i]*(B[i,j]*Vd[j]+G[i,j]*Vq[j]) for j=1:nb_nph) )
    Qij = @NLexpression(m, [i=1:nb_nph], sum(Vd[i]*(-B[i,j]*Vd[j]-G[i,j]*Vq[j]) + Vq[i]*(G[i,j]*Vd[j]-B[i,j]*Vq[j]) for j=1:nb_nph) )
    for i in idx_noref
        @NLconstraint(m, sum(Cbg[i,k]*PgY[k] for k=1:ngy_inv) - sum(CblY[i,kk]*YloadP[kk] for kk=1:nly_nph) - sum(CblD[i,kkk]*P_d2y[kkk] for kkk=1:nldy_nph) == Pij[i] )
        @NLconstraint(m, sum(Cbg[i,k]*QgY[k] for k=1:ngy_inv) - sum(CblY[i,kk]*YloadQ[kk] for kk=1:nly_nph) - sum(CblD[i,kkk]*Q_d2y[kkk] for kkk=1:nldy_nph) == Qij[i] )
    end

    ## Define cost
    cost_rect(Vll);

    ## Solve and print solution
    @time status = optimize!(m); println(termination_status(m));
    stat = string(termination_status(m))
    Vd0 = JuMP.value.(Vd); Vq0 = JuMP.value.(Vq); QgY0 = JuMP.value.(QgY);
    Vm0 = abs.(Vd0+im*Vq0); θ0 = angle.(Vd0+im*Vq0); Vll0 = JuMP.value.(Vll);
    print_VU(Vm0,θ0,Vll0)
    return Vm0,θ0,QgY0,stat

end

#----------------------------------------------------------------------------------------#
### FUNCTION TO CALCULATE LINEARIZED OPTIMAL POWER FLOW (POLAR-CONVENTIONAL)
function FOT_opf_pol(crit_node,Vm0,θ0,QgY0,dnc)

    ## Initialize values
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

    ## Start iterative linearized OPF
    max_iter = 1; opf_iter = 0; opf_err = 1; stat = "LOCALLY_SOLVED";
    YloadP = zeros(nly_nph); YloadQ = zeros(nly_nph);
    DloadP = zeros(nld_nph); DloadQ = zeros(nld_nph);
    I_ll = zeros(nld_nph)*(0+im*0); P_d2y = zeros(nldy_nph); Q_d2y = zeros(nldy_nph);
    Pl0 = zeros(nb_nph1); Ql0 = zeros(nb_nph1); Vll0 = zeros(nb_nph)
    Jp_Vm = spzeros(nb_nph1,nb_nph1); Jq_Vm = spzeros(nb_nph1,nb_nph1);
    Jp_θ = spzeros(nb_nph1,nb_nph1); Jq_θ = spzeros(nb_nph1,nb_nph1); global err_vec = zeros(max_iter)
    while opf_err>1e-4 && opf_iter<max_iter && (stat == "LOCALLY_SOLVED" || stat == "ALMOST_LOCALLY_SOLVED")

        ## Setup model and check for slow convergence
        opf_iter = opf_iter + 1;
        (Vm0,θ0,QgY0,pf_iter) = FP_pf(Vm0,θ0,QgY0,30); print("No of PF iterations: ", pf_iter,"\n")
        if (opf_iter == 11 && pf_iter >= 4)
            opf_iter = max_iter;
        elseif opf_iter>4 && ((err_vec[opf_iter-2:opf_iter-1] == err_vec[opf_iter-4:opf_iter-3]) || (err_vec[opf_iter-1] == err_vec[opf_iter-2] && err_vec[opf_iter-1] == err_vec[opf_iter-3]) || (err_vec[opf_iter-1] >= err_vec[opf_iter-2] >= err_vec[opf_iter-3]))
            opf_iter = max_iter;
        else
            for i=1:nb_nph
                k = Int(Tl[i])
                Vll0[i] =  sqrt(Vm0[i]^2 + Vm0[k]^2 - 2*Vm0[i]*Vm0[k]*cos(θ0[i]-θ0[k]) )
            end
            global m = Model(optimizer_with_attributes(Ipopt.Optimizer,"warm_start_init_point"=>"yes","print_level"=>0, "tol"=>1e-7));

            ## Define variables
            global Vm = @variable(m, [i=1:nb_nph], start = Vm0[i] )
            global θ = @variable(m,  [i=1:nb_nph], start = θ0[i] )
            global Vll = @variable(m, [i=1:nb_nph], start = Vll0[i] )
            global QgY = @variable(m, [i=1:ngy_inv], start = 0 )

            ## Define bounds on voltage magnitude and set slack bus angle to 0
            @constraint(m, Vm[idx_ref] .== 1 )
            @constraint(m, θ[idx_ref] .== [0;-2*pi/3;2*pi/3] )
            num_count= 1;
            for i=1:nb
                if  bus[i,2] == 3
                    num_count= num_count+ 3
                end
                for j = 3:2:7
                    if bus[i,j] != 0  && bus[i,2] != 3
                        @constraint(m, bus[i,j] <= Vm[num_count]<= bus[i,j+1] )
                        num_count= num_count+ 1
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
                # YloadP[i] = Sy[i,1]+Sy[i,3]*Vm0[Int(Sy[i,7])]+Sy[i,5]*Vm0[Int(Sy[i,7])]^2
                # YloadQ[i] = Sy[i,2]+Sy[i,4]*Vm0[Int(Sy[i,7])]+Sy[i,6]*Vm0[Int(Sy[i,7])]^2
            end
            for i=nl_y+1:nly_nph
                # dum1 = (Sy[i,3] + im*Sy[i,4])/(1*exp(im*Sy[i,8]))*(Vm0[Int(Sy[i,7])]*exp(im*θ0[Int(Sy[i,7])]))
                # YloadP[i] = Sy[i,1]+real(dum1)+Sy[i,5]*Vm0[Int(Sy[i,7])]^2
                # YloadQ[i] = Sy[i,2]+imag(dum1)+Sy[i,6]*Vm0[Int(Sy[i,7])]^2
                YloadP[i] = Sy[i,1]+Sy[i,3]*Vm0[Int(Sy[i,7])]+Sy[i,5]*Vm0[Int(Sy[i,7])]^2
                YloadQ[i] = Sy[i,2]+Sy[i,4]*Vm0[Int(Sy[i,7])]+Sy[i,6]*Vm0[Int(Sy[i,7])]^2
            end

            ## Calculate overall delta-connected loads in system
            for i=1:nld_nph
                j = Int(Tl[Int(Sd[i,7])]);
                dum1 = Vm0[Int(Sd[i,7])]*exp(im*θ0[Int(Sd[i,7])])-Vm0[j]*exp(im*θ0[j])
                dum2 = (Sd[i,3] + im*Sd[i,4])/(sqrt(3)*exp(im*Sd[i,8]))*dum1
                DloadP[i] = Sd[i,1]+real(dum2)+Sd[i,5]*Vll0[Int(Sd[i,7])]^2/3
                DloadQ[i] = Sd[i,2]+imag(dum2)+Sd[i,6]*Vll0[Int(Sd[i,7])]^2/3
                # DloadP[i] = Sd[i,1]+Sd[i,3]*Vm0[Int(Sd[i,7])]+Sd[i,5]*Vm0[Int(Sd[i,7])]^2
                # DloadQ[i] = Sd[i,2]+Sd[i,4]*Vm0[Int(Sd[i,7])]+Sd[i,6]*Vm0[Int(Sd[i,7])]^2
            end
            num_count= 1;
            for i=1:nld
                j = Int(dload[i,1])
                for k=Int.(sum(bus_count[1:j-1])+1:sum(bus_count[1:j]))
                    kk = Int(Tl[k])
                    if k != kk
                        I_ll[num_count]= (DloadP[num_count]-im*DloadQ[num_count])/(Vm0[k]*exp(-im*θ0[k])-Vm0[kk]*exp(-im*θ0[kk]));
                        num_count= num_count+ 1
                    end
                end
            end
            Ip = Tl2p*I_ll;
            num_count= 1
            for i=1:nld
                j = Int(dload[i,1])
                for k=Int.(sum(bus_count[1:j-1])+1:sum(bus_count[1:j]))
                    dum1 = Vm0[k]*exp(im*θ0[k])*conj(Ip[num_count])
                    P_d2y[num_count]= real(dum1); Q_d2y[num_count]= imag(dum1)
                    num_count= num_count+ 1
                end
            end

            ## Constraints on flow conservation (generation - demand = losses)
            count_i = 1;
            for i in idx_noref
                Pl0[count_i] = sum(Vm0[i]*Vm0[j]*(G[i,j]*cos(θ0[i]-θ0[j])+B[i,j]*sin(θ0[i]-θ0[j])) for j=1:nb_nph)
                Ql0[count_i] = sum(Vm0[i]*Vm0[j]*(G[i,j]*sin(θ0[i]-θ0[j])-B[i,j]*cos(θ0[i]-θ0[j])) for j=1:nb_nph)
                count_j = 1;
                for j in idx_noref
                    if i == j
                        Jp_Vm[count_i,count_i] =  Vm0[i]*Ym[i,i]*cos(Yθ[i,i]) + sum(Ym[i,k]*Vm0[k]*cos(θ0[i]-θ0[k]-Yθ[i,k]) for k=1:nb_nph);
                        Jq_Vm[count_i,count_i] = -Vm0[i]*Ym[i,i]*sin(Yθ[i,i]) + sum(Ym[i,k]*Vm0[k]*sin(θ0[i]-θ0[k]-Yθ[i,k]) for k=1:nb_nph);
                        if i == 1
                            Jp_θ[count_i,count_i] = -Vm0[i]*(sum(Ym[i,k]*Vm0[k]*sin(θ0[i]-θ0[k]-Yθ[i,k]) for k=2:nb_nph));
                            Jq_θ[count_i,count_i] =  Vm0[i]*(sum(Ym[i,k]*Vm0[k]*cos(θ0[i]-θ0[k]-Yθ[i,k]) for k=2:nb_nph));
                        elseif i == nb_nph
                            Jp_θ[count_i,count_i] = -Vm0[i]*(sum(Ym[i,k]*Vm0[k]*sin(θ0[i]-θ0[k]-Yθ[i,k]) for k=1:i-1));
                            Jq_θ[count_i,count_i] =  Vm0[i]*(sum(Ym[i,k]*Vm0[k]*cos(θ0[i]-θ0[k]-Yθ[i,k]) for k=1:i-1));
                        else
                            Jp_θ[count_i,count_i] = -Vm0[i]*(sum(Ym[i,k]*Vm0[k]*sin(θ0[i]-θ0[k]-Yθ[i,k]) for k=1:i-1) + sum(Ym[i,k]*Vm0[k]*sin(θ0[i]-θ0[k]-Yθ[i,k]) for k=i+1:nb_nph));
                            Jq_θ[count_i,count_i] =  Vm0[i]*(sum(Ym[i,k]*Vm0[k]*cos(θ0[i]-θ0[k]-Yθ[i,k]) for k=1:i-1) + sum(Ym[i,k]*Vm0[k]*cos(θ0[i]-θ0[k]-Yθ[i,k]) for k=i+1:nb_nph));
                        end
                    elseif i != j
                        Jp_θ[count_i,count_j]  =  Vm0[i]*Ym[i,j]*Vm0[j]*sin(θ0[i]-θ0[j]-Yθ[i,j]);
                        Jp_Vm[count_i,count_j] =  Vm0[i]*Ym[i,j]*cos(θ0[i]-θ0[j]-Yθ[i,j]);
                        Jq_θ[count_i,count_j]  = -Vm0[i]*Ym[i,j]*Vm0[j]*cos(θ0[i]-θ0[j]-Yθ[i,j]);
                        Jq_Vm[count_i,count_j] =  Vm0[i]*Ym[i,j]*sin(θ0[i]-θ0[j]-Yθ[i,j]);
                    end
                    count_j = count_j + 1;
                end
                count_i = count_i + 1;
            end
            ΔX = [Vm[idx_noref] - Vm0[idx_noref]; θ[idx_noref] - θ0[idx_noref]]
            J = [Jp_Vm Jp_θ; Jq_Vm Jq_θ];
            xP = Cbg_mod*PgY - CblY_mod*YloadP - CblD_mod*P_d2y - Pl0;
            xQ = Cbg_mod*QgY - CblY_mod*YloadQ - CblD_mod*Q_d2y - Ql0;
            @constraint(m, [xP; xQ] .== J*ΔX )
            # @constraint(m, Cbg_mod*PgY - CblY_mod*YloadP - CblD[idx_noref,:]*P_d2y .== Pl0 + [Jp_Vm Jp_θ]*ΔX )
            # @constraint(m, Cbg_mod*QgY - CblY_mod*YloadQ - CblD[idx_noref,:]*Q_d2y .== Ql0 + [Jq_Vm Jq_θ]*ΔX )

            ## Define cost
            cost_polar(Vll);

            ## Solve and update solution
            @time status = optimize!(m); stat = string(termination_status(m))
            V0 = Vm0.*exp.(im*θ0); Vm0 = JuMP.value.(Vm); θ0 = JuMP.value.(θ); QgY0 = JuMP.value.(QgY);
            V1 = Vm0.*exp.(im*θ0); opf_err = sum(abs.(V1-V0)); err_vec[opf_iter] = round(opf_err, RoundNearest, digits=4);
            print("FOT_OPF_pol: (",opf_iter,") Error = ",opf_err,", ",stat,"\n");
        end
    end

    ## Run power flow to convergence if OPF did not converge
    if opf_iter == max_iter
        println("\n****************OPF DID NOT CONVERGE!!*******************\n")
        # (Vm0,θ0,QgY0,pf_iter) = FP_pf(Vm0,θ0,QgY0,50); print("No of PF iterations: ", pf_iter,"\n");
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
    return Vm0,θ0,QgY0,stat
end

#----------------------------------------------------------------------------------------#

### FUNCTION TO CALCULATE OPTIMAL POWER FLOW (RECT-CURRENT & VOLTAGE FRAME)
function TPOPF_ivr(crit_node,Vm0,θ0,QgY0)

    ## Initialize model
    start_point(1); Vm0 = s0[1:nb_nph]; θ0 = s0[nb_nph+1:2*nb_nph];
    # (Vm0,θ0,QgY0,pf_iter) = FP_PF(Vm0,θ0,QgY0,50);
    Vll0 = zeros(nb_nph)
    for i=1:nb_nph
        k = Int(Tl[i])
        Vll0[i] =  sqrt(Vm0[i]^2 + Vm0[k]^2 - 2*Vm0[i]*Vm0[k]*cos(θ0[i]-θ0[k]) )
    end

    ## Define variables
    global m = Model(optimizer_with_attributes(Ipopt.Optimizer,"warm_start_init_point"=>"yes","print_level"=>0, "tol"=>1e-7));
    global Vm = @variable(m, [i=1:nb_nph], start = Vm0[i] )
    global θ = @variable(m,  [i=1:nb_nph], start = θ0[i] )
    global Vll = @variable(m, [i=1:nb_nph], start = Vll0[i] )
    global QgY = @variable(m, [i=1:ngy_inv], start = 0 )
    global Ibus_d = @variable(m,[i=1:nb_nph], start = 0 ); 
    global Ibus_q = @variable(m, [i=1:nb_nph], start = 0 );
    global Id_inj = @variable(m,[i=1:nb_nph], start = 0 ); 
    global Iq_inj = @variable(m, [i=1:nb_nph], start = 0 );
    global Pg_SS = @variable(m, [i=1:3], start = 0 )
    global Qg_SS = @variable(m, [i=1:3], start = 0 )

    ## Define bounds on voltage magnitude and set slack bus angle to 0
    @constraint(m, Vm[idx_ref] .== 1 )
    @constraint(m, θ[idx_ref] .== [0;-2*pi/3;2*pi/3] )
    @constraint(m, 0.9 .<= Vm[idx_noref] .<= 1.1 )
    for i=1:nb_nph
        k = Int(Tl[i])
        @NLconstraint(m, Vm[i]^2 + Vm[k]^2 - 2*Vm[i]*Vm[k]*cos(θ[i]-θ[k]) == Vll[i]^2 )
    end

    ## Define bounds on active and reactive power
    Pg_idx = [2 6 10];  Qg_idx = [4 8 12];
    PgY = spzeros(ngy_inv); num_count = 1;
    if PV_en == 1
        for i=1:ngy_inv
            j = Int(ygen1[i,1]); count1 = 1
            for k = 3:2:7
                if bus[j,k] != 0
                    PgY[num_count] = ygen1[i,Pg_idx[count1]]/baseMVA
                    @constraint(m, ygen1[i,Qg_idx[count1]]/baseMVA <= QgY[num_count] <= ygen1[i,Qg_idx[count1]+1]/baseMVA )
                    num_count = num_count + 1;
                end
                count1 = count1 + 1;
            end
        end
    end
  
    ## Generate load matrix
    ZIP_load()

    ## Calculate overall star-connected loads in system
    dumY_d = @NLexpression(m, [i=1:nl_y],  cos(Sy[i,8])*(Sy[i,3]*Vm[Int(Sy[i,7])]*cos(θ[Int(Sy[i,7])]) - Sy[i,4]*Vm[Int(Sy[i,7])]*sin(θ[Int(Sy[i,7])])) + sin(Sy[i,8])*(Sy[i,4]*Vm[Int(Sy[i,7])]*cos(θ[Int(Sy[i,7])]) + Sy[i,3]*Vm[Int(Sy[i,7])]*sin(θ[Int(Sy[i,7])]))  )
    dumY_q = @NLexpression(m, [i=1:nl_y], -sin(Sy[i,8])*(Sy[i,3]*Vm[Int(Sy[i,7])]*cos(θ[Int(Sy[i,7])]) - Sy[i,4]*Vm[Int(Sy[i,7])]*sin(θ[Int(Sy[i,7])])) + cos(Sy[i,8])*(Sy[i,4]*Vm[Int(Sy[i,7])]*cos(θ[Int(Sy[i,7])]) + Sy[i,3]*Vm[Int(Sy[i,7])]*sin(θ[Int(Sy[i,7])]))  )
    YloadP1 = @NLexpression(m, [i=1:nl_y], Sy[i,1]+dumY_d[i]+Sy[i,5]*Vm[Int(Sy[i,7])]^2 )
    YloadQ1 = @NLexpression(m, [i=1:nl_y], Sy[i,2]+dumY_q[i]+Sy[i,6]*Vm[Int(Sy[i,7])]^2 )
    YloadP2 = @NLexpression(m, [i=1:nlh], Sy[i+nl_y,1]+Sy[i+nl_y,3]*Vm[Int(Sy[i+nl_y,7])]+Sy[i+nl_y,5]*Vm[Int(Sy[i+nl_y,7])]^2 )
    YloadQ2 = @NLexpression(m, [i=1:nlh], Sy[i+nl_y,2]+Sy[i+nl_y,4]*Vm[Int(Sy[i+nl_y,7])]+Sy[i+nl_y,6]*Vm[Int(Sy[i+nl_y,7])]^2 )
    YloadP = [YloadP1; YloadP2]; YloadQ = [YloadQ1; YloadQ2];
    Ycap = @NLexpression(m, [i=1:nb_nph], Sh_cap[i]*Vm[i]^2 ) 

    ## Calculate overall delta-connected loads in system
    dum_Vll_d = @NLexpression(m, [i=1:nld_nph], Vm[Int(Sd[i,7])]*cos(θ[Int(Sd[i,7])]) - Vm[Int(Tl[Int(Sd[i,7])])]*cos(θ[Int(Tl[Int(Sd[i,7])])]) )
    dum_Vll_q = @NLexpression(m, [i=1:nld_nph], Vm[Int(Sd[i,7])]*sin(θ[Int(Sd[i,7])]) - Vm[Int(Tl[Int(Sd[i,7])])]*sin(θ[Int(Tl[Int(Sd[i,7])])]) )
    dumD_d = @NLexpression(m, [i=1:nld_nph],  cos(Sd[i,8])*(Sd[i,3]*dum_Vll_d[i] - Sd[i,4]*dum_Vll_q[i]) + sin(Sd[i,8])*(Sd[i,4]*dum_Vll_d[i] + Sd[i,3]*dum_Vll_q[i])  )
    dumD_q = @NLexpression(m, [i=1:nld_nph], -sin(Sd[i,8])*(Sd[i,3]*dum_Vll_d[i] - Sd[i,4]*dum_Vll_q[i]) + cos(Sd[i,8])*(Sd[i,4]*dum_Vll_d[i] + Sd[i,3]*dum_Vll_q[i])  )
    DloadP = @NLexpression(m, [i=1:nld_nph], Sd[i,1]+dumD_d[i]/sqrt(3)+Sd[i,5]*Vll[Int(Sd[i,7])]^2/3 )
    DloadQ = @NLexpression(m, [i=1:nld_nph], Sd[i,2]+dumD_q[i]/sqrt(3)+Sd[i,6]*Vll[Int(Sd[i,7])]^2/3 )
    
    ## Convert delta to equivalent wye-connected loads
    @variable(m, Ill_d[i=1:nld_nph], start = 0 ); @variable(m, Ill_q[i=1:nld_nph], start = 0 )
    @variable(m, P_d2y[i=1:nldy_nph], start = 0 ); @variable(m, Q_d2y[i=1:nldy_nph], start = 0 )
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

    ## KCL constraints
    @NLconstraint(m, [i in 1:nb_nph], Vm[i]*Id_inj[i] == ((sum(CblY[i,kk]*YloadP[kk] for kk=1:nly_nph) + sum(CblD[i,kkk]*P_d2y[kkk] for kkk=1:nldy_nph) - sum(Cbg[i,k]*PgY[k] for k=1:ngy_inv) - sum(Cbg_SS[i,k]*Pg_SS[k] for k=1:3))*cos(θ[i]) + (sum(CblY[i,kk]*YloadQ[kk] for kk=1:nly_nph) + Ycap[i] + sum(CblD[i,kkk]*Q_d2y[kkk] for kkk=1:nldy_nph) - sum(Cbg[i,k]*QgY[k] for k=1:ngy_inv))*sin(θ[i]) - sum(Cbg_SS[i,k]*Qg_SS[k] for k=1:3)) );
    @NLconstraint(m, [i in 1:nb_nph], Vm[i]*Iq_inj[i] == ((sum(CblY[i,kk]*YloadP[kk] for kk=1:nly_nph) + sum(CblD[i,kkk]*P_d2y[kkk] for kkk=1:nldy_nph) - sum(Cbg[i,k]*PgY[k] for k=1:ngy_inv) - sum(Cbg_SS[i,k]*Pg_SS[k] for k=1:3))*sin(θ[i]) - (sum(CblY[i,kk]*YloadQ[kk] for kk=1:nly_nph) + Ycap[i] + sum(CblD[i,kkk]*Q_d2y[kkk] for kkk=1:nldy_nph) - sum(Cbg[i,k]*QgY[k] for k=1:ngy_inv))*cos(θ[i]) - sum(Cbg_SS[i,k]*Qg_SS[k] for k=1:3)) );
    for i=1:nb
        ii = Int(sum(bus_count[1:i-1])+1):Int(sum(bus_count[1:i]))
        j = findall(x -> x==1, A_conn[:,i] ); neighbours = Array{Int}(undef, length(j)); 
        for (idx,val) in enumerate(j)
            neighbours[idx] = findfirst(x -> x==-1, A_conn[val,:] );
        end
        if isempty(neighbours)
            @constraint(m, Ibus_d[ii] .== Id_inj[ii] ) 
            @constraint(m, Ibus_q[ii] .== Iq_inj[ii] ) 
        else
            neighbours_1ph = Array{Int}(undef, Int(sum(bus_count[i] for i in neighbours))); num_count=1
            for (idx,val) in enumerate(neighbours)
                neighbours_1ph[num_count:num_count+Int(bus_count[val])-1] = Int(sum(bus_count[1:val-1])+1):Int(sum(bus_count[1:val]))
                num_count = Int(num_count+bus_count[val])
            end
            Isum_d = @expression(m, [k=1:length(ii)], sum(Ibus_d[kk] for kk in neighbours_1ph if bus_ϕ[kk] == bus_ϕ[ii[k]]) )
            Isum_q = @expression(m, [k=1:length(ii)], sum(Ibus_q[kk] for kk in neighbours_1ph if bus_ϕ[kk] == bus_ϕ[ii[k]]) )
            if i == Int(tf_branch[vr_idx[1],2])
                if tap_en == 0
                    if bus_ph[i] == "AN"
                        global vr_tap = 1 + 0.00625*tf_branch[vr_idx[1],6]
                    elseif bus_ph[i] == "BN"
                        global vr_tap = 1 + 0.00625*tf_branch[vr_idx[1],7]
                    elseif bus_ph[i] == "CN"
                        global vr_tap = 1 + 0.00625*tf_branch[vr_idx[1],8]
                    elseif bus_ph[i] == "ABN" || bus_ph[i] == "BAN"
                        global vr_tap = 1 .+ 0.00625*tf_branch[vr_idx[1],6:7]
                    elseif bus_ph[i] == "BCN" || bus_ph[i] == "CBN"
                        global vr_tap = 1 .+ 0.00625*tf_branch[vr_idx[1],7:8]
                    elseif bus_ph[i] == "ACN" || bus_ph[i] == "CAN"
                        global vr_tap = 1 .+ 0.00625*tf_branch[vr_idx[1],6:2:8]
                    else
                        global vr_tap = 1 .+ 0.00625*tf_branch[vr_idx[1],6:8]
                    end
                    @constraint(m, Ibus_d[ii] .== vr_tap.*(Id_inj[ii] + Isum_d) ) 
                    @constraint(m, Ibus_q[ii] .== vr_tap.*(Iq_inj[ii] + Isum_q) ) 
                else
                    if bus_ph[i] == "AN" || bus_ph[i] == "BN" || bus_ph[i] == "CN"
                        global vr_tap = @variable(m, start = 1 );
                    elseif bus_ph[i] == "ABN" || bus_ph[i] == "BAN" || bus_ph[i] == "BCN" || bus_ph[i] == "CBN" || bus_ph[i] == "ACN" || bus_ph[i] == "CAN"
                        global vr_tap = @variable(m, [i=1:2], start = 1 );
                    else
                        global vr_tap = @variable(m, [i=1:3], start = 1 );
                    end
                    @constraint(m, 0.9 .<= vr_tap .<= 1.1 )
                    for (idx,iii) in enumerate(ii)
                        @NLconstraint(m, Ibus_d[iii] == vr_tap[idx]*(Id_inj[iii] + Isum_d[idx]) ) 
                        @NLconstraint(m, Ibus_q[iii] == vr_tap[idx]*(Iq_inj[iii] + Isum_q[idx]) ) 
                    end
                end  
            else
                @constraint(m, Ibus_d[ii] .== Id_inj[ii] + Isum_d ) 
                @constraint(m, Ibus_q[ii] .== Iq_inj[ii] + Isum_q ) 
            end   
        end
    end

    ## KVL constraints
    count2 = 1;
    for i=1:nbr
        j = Int(branch[i,1]); jj = Int(sum(bus_count[1:j-1])+1):Int(sum(bus_count[1:j]))
        k = Int(branch[i,2]); kk = Int(sum(bus_count[1:k-1])+1):Int(sum(bus_count[1:k]))
        nz = findmin([length(jj) length(kk)])[1]
        jjjj = zeros(nz); kkkk = zeros(nz); num_count=1
        for jjj=jj
            for kkk=kk
                if bus_ϕ[jjj] == bus_ϕ[kkk]
                    jjjj[num_count] = jjj; kkkk[num_count] = kkk;
                    num_count=num_count+1
                end
            end
        end
        jj = Int.(jjjj); kk = Int.(kkkk); 
        Zbr = -pinv(Y[jj,kk]); Rbr = real(Zbr); Xbr = imag(Zbr);
        if tf_branch[i,3] == 8
            v_min= VR_config[1,3]*VR_config[1,2]/bus_Vnom[k]; 
            v_max= VR_config[1,4]*VR_config[1,2]/bus_Vnom[k];
             # @constraint(m, v_min.<= Vm[kk] .<= v_max)
            v_ref = 0.5*(VR_config[1,3]+VR_config[1,4]);            
            CT_s = 0.2; CT_p = VR_config[1,5]*CT_s;
            R_comp = VR_config[1,6]/CT_s; X_comp = VR_config[1,7]/CT_s;
            # Icomp_d = Ibus_d[kkk]*baseI[k]/VR_config[1,6]; Icomp_q = Ibus_q[kkk]*baseI[k]/VR_config[1,6];
            # Vcomp_d = Vm[kkk]*cos(θ[kkk])*bus_Vnom[k]/VR_config[1,2]; Vcomp_q = Vm[kkk]*sin(θ[kkk])*bus_Vnom[k]/VR_config[1,2];
            for (idx_r,jjj) in enumerate(jj)
                kkk = kk[idx_r]; 
                # Vreg_d = Vm[kkk]*cos(θ[kkk])*bus_Vnom[k]/VR_config[1,2] - (R_comp*(Ibus_d[kkk]*baseI[k]/VR_config[1,5]) - X_comp*(Ibus_q[kkk]*baseI[k]/VR_config[1,5]) );
                # Vreg_q = Vm[kkk]*sin(θ[kkk])*bus_Vnom[k]/VR_config[1,2] - (R_comp*(Ibus_q[kkk]*baseI[k]/VR_config[1,5]) + X_comp*(Ibus_d[kkk]*baseI[k]/VR_config[1,5]) );
                @NLconstraint(m, vr_tap[idx_r] == (v_ref-sqrt((Vm[kkk]*cos(θ[kkk])*bus_Vnom[k]/VR_config[1,2] - (R_comp*(Ibus_d[kkk]*baseI[k]/VR_config[1,5]) - X_comp*(Ibus_q[kkk]*baseI[k]/VR_config[1,5]) ))^2 + (Vm[kkk]*sin(θ[kkk])*bus_Vnom[k]/VR_config[1,2] - (R_comp*(Ibus_q[kkk]*baseI[k]/VR_config[1,5]) + X_comp*(Ibus_d[kkk]*baseI[k]/VR_config[1,5]) ))^2))/0.75*0.00625+1 );
                @NLconstraint(m, Vm[kkk]*cos(θ[kkk]) == vr_tap[idx_r]*Vm[jjj]*cos(θ[jjj]) - sum(Rbr[idx_r,idx_c]*Ibus_d[val] for (idx_c,val) in enumerate(kk)) + sum(Xbr[idx_r,idx_c]*Ibus_q[val] for (idx_c,val) in enumerate(kk)) )
                @NLconstraint(m, Vm[kkk]*sin(θ[kkk]) == vr_tap[idx_r]*Vm[jjj]*sin(θ[jjj]) - sum(Rbr[idx_r,idx_c]*Ibus_q[val] for (idx_c,val) in enumerate(kk)) - sum(Xbr[idx_r,idx_c]*Ibus_d[val] for (idx_c,val) in enumerate(kk)) )
            end
            count2 = count2 + 1;
        else
            for (idx_r,jjj) in enumerate(jj)
                kkk = kk[idx_r]; 
                @NLconstraint(m, Vm[kkk]*cos(θ[kkk]) == (Vm[jjj]*cos(θ[jjj]) - sum(Rbr[idx_r,idx_c]*Ibus_d[val] for (idx_c,val) in enumerate(kk)) + sum(Xbr[idx_r,idx_c]*Ibus_q[val] for (idx_c,val) in enumerate(kk))) )
                @NLconstraint(m, Vm[kkk]*sin(θ[kkk]) == (Vm[jjj]*sin(θ[jjj]) - sum(Rbr[idx_r,idx_c]*Ibus_q[val] for (idx_c,val) in enumerate(kk)) - sum(Xbr[idx_r,idx_c]*Ibus_d[val] for (idx_c,val) in enumerate(kk))) )
            end
         end
    end

    ## Define cost
    cost_polar(Vll);

    ## Solve and print solution
    @time status = optimize!(m); println(termination_status(m));
    sts = string(termination_status(m))
    Vm0 = JuMP.value.(Vm); θ0 = JuMP.value.(θ); QgY0 = JuMP.value.(QgY); Vll0 = JuMP.value.(Vll);
    global Ibus_d0 = JuMP.value.(Ibus_d); global Ibus_q0 = JuMP.value.(Ibus_q);
    print_VU(Vm0,θ0,Vll0)
    println("TOTAL TIME:")

    ## Adjust admittance matrix for new regulator taps
    if tap_en == 1
        j = Int(tf_branch[vr_idx[1],1]); jj = Int(sum(bus_count[1:j-1])+1):Int(sum(bus_count[1:j]))
        k = Int(tf_branch[vr_idx[1],2]); kk = Int(sum(bus_count[1:k-1])+1):Int(sum(bus_count[1:k]))
        global tap_ratio = (JuMP.value.(vr_tap) .- 1)./0.00625;
        tap_old = tf_branch[vr_idx[1],6:8]*0.00625 .+ 1;
        yt = 1/(tf_branch[vr_idx[1],4] + tf_branch[vr_idx[1],5]im);
        Y1_old = yt*Matrix{Float64}(I, 3, 3); tap_old = Diagonal(tap_old);
        tap_new = Diagonal(JuMP.value.(vr_tap));
        G[jj,jj] = G[jj,jj]+real(Y1_old).*tap_new^2-real(Y1_old).*tap_old^2;
        B[jj,jj] = B[jj,jj]+imag(Y1_old).*tap_new^2-imag(Y1_old).*tap_old^2;
        G[jj,kk] = real(-Y1_old).*tap_new;
        B[jj,kk] = imag(-Y1_old).*tap_new;
        G[kk,jj] = real(-Y1_old).*tap_new;
        B[kk,jj] = imag(-Y1_old).*tap_new;
    end
    global Y = G+im*B;
    return Vm0,θ0,QgY0,sts
end