### FUNCTION TO DEFINE COST (POLAR)
function cost_polar(Vll)

    ## Minimize nothing
    if obj_op == 1 || obj_op == 6
        @objective(m, Min, 0*M )

    ## Minimize VUF
    elseif obj_op == 2
        @NLexpression(m, Vn[i=1:n_VU], (sum(Vm[j]*cos(θ[j]+an[j]) for j=Int.(idx_bus_3ph1[i*3-2:i*3])))^2 + (sum(Vm[j]*sin(θ[j]+an[j]) for j=Int.(idx_bus_3ph1[i*3-2:i*3])))^2 )
        @NLexpression(m, Vp[i=1:n_VU], (sum(Vm[j]*cos(θ[j]+ap[j]) for j=Int.(idx_bus_3ph1[i*3-2:i*3])))^2 + (sum(Vm[j]*sin(θ[j]+ap[j]) for j=Int.(idx_bus_3ph1[i*3-2:i*3])))^2 )
        global cost_VUF = @NLexpression(m, sum(Vn[i]/Vp[i] for i=1:n_VU) )
        # global cost_VUF = @NLexpression(m, sum(Vn[i] for i=1:n_VU) )
        # @NLobjective(m, Min, M*cost_VUF )
        @NLobjective(m, Min, M*cost_VUF + Q_pen*sum(QgY[i]^2 for i=1:ngy_inv))

        ## Constraints on unbalance at all nodes
        # @NLexpression(m, Vn_all[i=1:nb_3ph], (sum(Vm[j]*cos(θ[j]+an[j]) for j=Int.(idx_bus_3ph2[i*3-2:i*3])))^2 + (sum(Vm[j]*sin(θ[j]+an[j]) for j=Int.(idx_bus_3ph2[i*3-2:i*3])))^2 )
        # @NLexpression(m, Vp_all[i=1:nb_3ph], (sum(Vm[j]*cos(θ[j]+ap[j]) for j=Int.(idx_bus_3ph2[i*3-2:i*3])))^2 + (sum(Vm[j]*sin(θ[j]+ap[j]) for j=Int.(idx_bus_3ph2[i*3-2:i*3])))^2 )
        # u = VUF_init[:,count_sch];
        # for i=1:nb_3ph
        #     if u[i] >= 5e-3
        #         @NLconstraint(m, Vn_all[i] <= Vp_all[i] * u[i]^2 )
        #     end
        # end

    ## Minimize power injections at substation
    elseif obj_op == 3
        Pij = @NLexpression(m, [i=1:3], sum(Vm[idx_ref[i]]*Vm[j]*(G[idx_ref[i],j]*cos(θ[idx_ref[i]]-θ[j])+B[idx_ref[i],j]*sin(θ[idx_ref[i]]-θ[j])) for j=1:nb_nph) )
        # Qij = @NLexpression(m, [i=1:3], sum(Vm[idx_ref[i]]*Vm[j]*(G[idx_ref[i],j]*sin(θ[idx_ref[i]]-θ[j])-B[idx_ref[i],j]*cos(θ[idx_ref[i]]-θ[j])) for j=1:nb_nph) )
        global cost_SS = @NLexpression(m, sum(Pij[i] for i=1:3) )
        @NLobjective(m, Min, cost_SS )

    ## Minimize PVUR
    elseif obj_op == 4
        @variable(m, zp1[i=1:n_VU*3] )
        @variable(m, zp2[i=1:n_VU] )
        @expression(m, Vavg[i=1:n_VU], sum(Vm[j] for j=Int.(idx_bus_3ph1[i*3-2:i*3]))/3 )
        @expression(m, Vdev[i=1:n_VU*3], Vm[Int(idx_bus_3ph1[i])]-Vavg[Int(ceil(i/3))] )
        for i=1:n_VU
            @constraint(m, zp1[i*3-2:i*3] .>=  Vdev[i*3-2:i*3] )
            @constraint(m, zp1[i*3-2:i*3] .>= -Vdev[i*3-2:i*3] )
            @constraint(m, zp2[i] * Vavg[i] .>= zp1[i*3-2:i*3] )
        end
        global cost_PVUR = @expression(m, sum(zp2[i] for i=1:n_VU) )
        @objective(m, Min, cost_PVUR*M )

    ## Minimize LVUR
    elseif obj_op == 5
        @variable(m, zl1[i=1:n_VU*3] )
        @variable(m, zl2[i=1:n_VU] )
        @NLexpression(m, Vlavg[i=1:n_VU], sum(Vll[j] for j=Int.(idx_bus_3ph1[i*3-2:i*3]))/3 )
        @NLexpression(m, Vldev[i=1:n_VU*3], Vll[Int(idx_bus_3ph1[i])]-Vlavg[Int(ceil(i/3))] )
        for i=1:n_VU*3
            @NLconstraint(m, zl1[i] >=  Vldev[i] )
            @NLconstraint(m, zl1[i] >= -Vldev[i] )
        end
        for i=1:n_VU
            @NLconstraint(m, zl2[i] * Vlavg[i] >= zl1[i*3-2] )
            @NLconstraint(m, zl2[i] * Vlavg[i] >= zl1[i*3-1] )
            @NLconstraint(m, zl2[i] * Vlavg[i] >= zl1[i*3] )
        end
        global cost_LVUR = @expression(m, sum(zl2[i] for i=1:n_VU) )
        @objective(m, Min, cost_LVUR*M )

    # Minimize losses with unbalance constraints
    elseif obj_op == 6

        ## Calculate losses
        Pij = @NLexpression(m, [i=1:nb_nph], sum(Vm[i]*Vm[j]*(G[i,j]*cos(θ[i]-θ[j])+B[i,j]*sin(θ[i]-θ[j])) for j=1:nb_nph) )
        global cost_loss = @NLexpression(m, sum(Pij[i] for i=1:nb_nph ) )
        @NLobjective(m, Min, cost_loss*M )

        ## VUF Constraints
        # @NLexpression(m, Vn_all[i=1:nb_3ph], (sum(Vm[j]*cos(θ[j]+an[j]) for j=Int.(idx_bus_3ph2[i*3-2:i*3])))^2 + (sum(Vm[j]*sin(θ[j]+an[j]) for j=Int.(idx_bus_3ph2[i*3-2:i*3])))^2 )
        # @NLexpression(m, Vp_all[i=1:nb_3ph], (sum(Vm[j]*cos(θ[j]+ap[j]) for j=Int.(idx_bus_3ph2[i*3-2:i*3])))^2 + (sum(Vm[j]*sin(θ[j]+ap[j]) for j=Int.(idx_bus_3ph2[i*3-2:i*3])))^2 )
        # @NLexpression(m, V0_all[i=1:nb_3ph], (sum(Vm[j]*cos(θ[j]) for j=Int.(idx_bus_3ph1[i*3-2:i*3])))^2 + (sum(Vm[j]*sin(θ[j]) for j=Int.(idx_bus_3ph1[i*3-2:i*3])))^2 )
        # VUF_lim = 0.02
        # num_count = 1
        # for i in idx_bus_3ph
            # @NLconstraint(m, Vn_all[num_count] <= Vp_all[num_count] * VUF_lim^2 )
            # @NLconstraint(m, V0_all[num_count] <= Vp_all[num_count] * VUF_lim^2 )
            # num_count = num_count + 1
        # end

        ## PVUR Constraints
        # PVUR_lim = 0.02
        # PVUR_mat = [PVUR_lim+2 PVUR_lim-1 PVUR_lim-1; PVUR_lim-1 PVUR_lim+2 PVUR_lim-1; PVUR_lim-1 PVUR_lim-1 PVUR_lim+2; PVUR_lim-2 PVUR_lim+1 PVUR_lim+1; PVUR_lim+1 PVUR_lim-2 PVUR_lim+1; PVUR_lim+1 PVUR_lim+1 PVUR_lim-2]
        # for i=1:nb_3ph
        #     j = Int.(idx_bus_3ph1[i*3-2:i*3])
        #     @constraint(m, PVUR_mat*Vm[j] .>= 0 )
        # end

        ## LVUR Constraints
        # @variable(m, Vll[i=1:nb_nph] >= 0, start = s0[i+2*nb_nph] )
        # for i=1:nb_nph
        #     k = Int(Tl[i])
        #     @NLconstraint(m, Vm[i]^2 + Vm[k]^2 - 2*Vm[i]*Vm[k]*cos(θ[i]-θ[k]) == Vll[i]^2 )
        # end
        # LVUR_lim = 0.03
        # LVUR_mat = [LVUR_lim+2 LVUR_lim-1 LVUR_lim-1; LVUR_lim-1 LVUR_lim+2 LVUR_lim-1; LVUR_lim-1 LVUR_lim-1 LVUR_lim+2; LVUR_lim-2 LVUR_lim+1 LVUR_lim+1; LVUR_lim+1 LVUR_lim-2 LVUR_lim+1; LVUR_lim+1 LVUR_lim+1 LVUR_lim-2]
        # for i=1:nb_3ph
        #     j = Int.(idx_bus_3ph1[i*3-2:i*3])
        #     @constraint(m, LVUR_mat*Vll[j] .>= 0 )
        # end

    ## Minimize Q-injections
    elseif obj_op == 7
        @objective(m, Min, sum(QgY[i] for i=1:ngy_inv) )

    ## Minimize linearized VUF
    elseif obj_op == 8
        @expression(m, VUF_lin[i=1:n_VU], 0 )
        for i=1:n_VU
            j=Int.(idx_bus_3ph1[i*3-2:i*3]); V_vec = [Vm[j];θ[j]];
            Vm_a = Vm0[j[1]]; Vm_b = Vm0[j[2]]; Vm_c = Vm0[j[3]];
            Va_a = θ0[j[1]]; Va_b = θ0[j[2]]; Va_c = θ0[j[3]];
            dum1 = [ (2*sin(Va_a)*(Vm_a*sin(Va_a) + Vm_b*sin(Va_b - (2*pi)/3) + Vm_c*sin(Va_c + (2*pi)/3)) + 2*cos(Va_a)*(Vm_a*cos(Va_a) + Vm_b*cos(Va_b - (2*pi)/3) + Vm_c*cos(Va_c + (2*pi)/3)))/((Vm_a*cos(Va_a) + Vm_b*cos(Va_b + (2*pi)/3) + Vm_c*cos(Va_c - (2*pi)/3))^2 + (Vm_a*sin(Va_a) + Vm_b*sin(Va_b + (2*pi)/3) + Vm_c*sin(Va_c - (2*pi)/3))^2) - ((2*sin(Va_a)*(Vm_a*sin(Va_a) + Vm_b*sin(Va_b + (2*pi)/3) + Vm_c*sin(Va_c - (2*pi)/3)) + 2*cos(Va_a)*(Vm_a*cos(Va_a) + Vm_b*cos(Va_b + (2*pi)/3) + Vm_c*cos(Va_c - (2*pi)/3)))*((Vm_a*cos(Va_a) + Vm_b*cos(Va_b - (2*pi)/3) + Vm_c*cos(Va_c + (2*pi)/3))^2 + (Vm_a*sin(Va_a) + Vm_b*sin(Va_b - (2*pi)/3) + Vm_c*sin(Va_c + (2*pi)/3))^2))/((Vm_a*cos(Va_a) + Vm_b*cos(Va_b + (2*pi)/3) + Vm_c*cos(Va_c - (2*pi)/3))^2 + (Vm_a*sin(Va_a) + Vm_b*sin(Va_b + (2*pi)/3) + Vm_c*sin(Va_c - (2*pi)/3))^2)^2, (2*sin(Va_b - (2*pi)/3)*(Vm_a*sin(Va_a) + Vm_b*sin(Va_b - (2*pi)/3) + Vm_c*sin(Va_c + (2*pi)/3)) + 2*cos(Va_b - (2*pi)/3)*(Vm_a*cos(Va_a) + Vm_b*cos(Va_b - (2*pi)/3) + Vm_c*cos(Va_c + (2*pi)/3)))/((Vm_a*cos(Va_a) + Vm_b*cos(Va_b + (2*pi)/3) + Vm_c*cos(Va_c - (2*pi)/3))^2 + (Vm_a*sin(Va_a) + Vm_b*sin(Va_b + (2*pi)/3) + Vm_c*sin(Va_c - (2*pi)/3))^2) - ((2*sin(Va_b + (2*pi)/3)*(Vm_a*sin(Va_a) + Vm_b*sin(Va_b + (2*pi)/3) + Vm_c*sin(Va_c - (2*pi)/3)) + 2*cos(Va_b + (2*pi)/3)*(Vm_a*cos(Va_a) + Vm_b*cos(Va_b + (2*pi)/3) + Vm_c*cos(Va_c - (2*pi)/3)))*((Vm_a*cos(Va_a) + Vm_b*cos(Va_b - (2*pi)/3) + Vm_c*cos(Va_c + (2*pi)/3))^2 + (Vm_a*sin(Va_a) + Vm_b*sin(Va_b - (2*pi)/3) + Vm_c*sin(Va_c + (2*pi)/3))^2))/((Vm_a*cos(Va_a) + Vm_b*cos(Va_b + (2*pi)/3) + Vm_c*cos(Va_c - (2*pi)/3))^2 + (Vm_a*sin(Va_a) + Vm_b*sin(Va_b + (2*pi)/3) + Vm_c*sin(Va_c - (2*pi)/3))^2)^2, (2*sin(Va_c + (2*pi)/3)*(Vm_a*sin(Va_a) + Vm_b*sin(Va_b - (2*pi)/3) + Vm_c*sin(Va_c + (2*pi)/3)) + 2*cos(Va_c + (2*pi)/3)*(Vm_a*cos(Va_a) + Vm_b*cos(Va_b - (2*pi)/3) + Vm_c*cos(Va_c + (2*pi)/3)))/((Vm_a*cos(Va_a) + Vm_b*cos(Va_b + (2*pi)/3) + Vm_c*cos(Va_c - (2*pi)/3))^2 + (Vm_a*sin(Va_a) + Vm_b*sin(Va_b + (2*pi)/3) + Vm_c*sin(Va_c - (2*pi)/3))^2) - ((2*sin(Va_c - (2*pi)/3)*(Vm_a*sin(Va_a) + Vm_b*sin(Va_b + (2*pi)/3) + Vm_c*sin(Va_c - (2*pi)/3)) + 2*cos(Va_c - (2*pi)/3)*(Vm_a*cos(Va_a) + Vm_b*cos(Va_b + (2*pi)/3) + Vm_c*cos(Va_c - (2*pi)/3)))*((Vm_a*cos(Va_a) + Vm_b*cos(Va_b - (2*pi)/3) + Vm_c*cos(Va_c + (2*pi)/3))^2 + (Vm_a*sin(Va_a) + Vm_b*sin(Va_b - (2*pi)/3) + Vm_c*sin(Va_c + (2*pi)/3))^2))/((Vm_a*cos(Va_a) + Vm_b*cos(Va_b + (2*pi)/3) + Vm_c*cos(Va_c - (2*pi)/3))^2 + (Vm_a*sin(Va_a) + Vm_b*sin(Va_b + (2*pi)/3) + Vm_c*sin(Va_c - (2*pi)/3))^2)^2, ((2*Vm_a*sin(Va_a)*(Vm_a*cos(Va_a) + Vm_b*cos(Va_b + (2*pi)/3) + Vm_c*cos(Va_c - (2*pi)/3)) - 2*Vm_a*cos(Va_a)*(Vm_a*sin(Va_a) + Vm_b*sin(Va_b + (2*pi)/3) + Vm_c*sin(Va_c - (2*pi)/3)))*((Vm_a*cos(Va_a) + Vm_b*cos(Va_b - (2*pi)/3) + Vm_c*cos(Va_c + (2*pi)/3))^2 + (Vm_a*sin(Va_a) + Vm_b*sin(Va_b - (2*pi)/3) + Vm_c*sin(Va_c + (2*pi)/3))^2))/((Vm_a*cos(Va_a) + Vm_b*cos(Va_b + (2*pi)/3) + Vm_c*cos(Va_c - (2*pi)/3))^2 + (Vm_a*sin(Va_a) + Vm_b*sin(Va_b + (2*pi)/3) + Vm_c*sin(Va_c - (2*pi)/3))^2)^2 - (2*Vm_a*sin(Va_a)*(Vm_a*cos(Va_a) + Vm_b*cos(Va_b - (2*pi)/3) + Vm_c*cos(Va_c + (2*pi)/3)) - 2*Vm_a*cos(Va_a)*(Vm_a*sin(Va_a) + Vm_b*sin(Va_b - (2*pi)/3) + Vm_c*sin(Va_c + (2*pi)/3)))/((Vm_a*cos(Va_a) + Vm_b*cos(Va_b + (2*pi)/3) + Vm_c*cos(Va_c - (2*pi)/3))^2 + (Vm_a*sin(Va_a) + Vm_b*sin(Va_b + (2*pi)/3) + Vm_c*sin(Va_c - (2*pi)/3))^2), ((2*Vm_b*sin(Va_b + (2*pi)/3)*(Vm_a*cos(Va_a) + Vm_b*cos(Va_b + (2*pi)/3) + Vm_c*cos(Va_c - (2*pi)/3)) - 2*Vm_b*cos(Va_b + (2*pi)/3)*(Vm_a*sin(Va_a) + Vm_b*sin(Va_b + (2*pi)/3) + Vm_c*sin(Va_c - (2*pi)/3)))*((Vm_a*cos(Va_a) + Vm_b*cos(Va_b - (2*pi)/3) + Vm_c*cos(Va_c + (2*pi)/3))^2 + (Vm_a*sin(Va_a) + Vm_b*sin(Va_b - (2*pi)/3) + Vm_c*sin(Va_c + (2*pi)/3))^2))/((Vm_a*cos(Va_a) + Vm_b*cos(Va_b + (2*pi)/3) + Vm_c*cos(Va_c - (2*pi)/3))^2 + (Vm_a*sin(Va_a) + Vm_b*sin(Va_b + (2*pi)/3) + Vm_c*sin(Va_c - (2*pi)/3))^2)^2 - (2*Vm_b*sin(Va_b - (2*pi)/3)*(Vm_a*cos(Va_a) + Vm_b*cos(Va_b - (2*pi)/3) + Vm_c*cos(Va_c + (2*pi)/3)) - 2*Vm_b*cos(Va_b - (2*pi)/3)*(Vm_a*sin(Va_a) + Vm_b*sin(Va_b - (2*pi)/3) + Vm_c*sin(Va_c + (2*pi)/3)))/((Vm_a*cos(Va_a) + Vm_b*cos(Va_b + (2*pi)/3) + Vm_c*cos(Va_c - (2*pi)/3))^2 + (Vm_a*sin(Va_a) + Vm_b*sin(Va_b + (2*pi)/3) + Vm_c*sin(Va_c - (2*pi)/3))^2), ((2*Vm_c*sin(Va_c - (2*pi)/3)*(Vm_a*cos(Va_a) + Vm_b*cos(Va_b + (2*pi)/3) + Vm_c*cos(Va_c - (2*pi)/3)) - 2*Vm_c*cos(Va_c - (2*pi)/3)*(Vm_a*sin(Va_a) + Vm_b*sin(Va_b + (2*pi)/3) + Vm_c*sin(Va_c - (2*pi)/3)))*((Vm_a*cos(Va_a) + Vm_b*cos(Va_b - (2*pi)/3) + Vm_c*cos(Va_c + (2*pi)/3))^2 + (Vm_a*sin(Va_a) + Vm_b*sin(Va_b - (2*pi)/3) + Vm_c*sin(Va_c + (2*pi)/3))^2))/((Vm_a*cos(Va_a) + Vm_b*cos(Va_b + (2*pi)/3) + Vm_c*cos(Va_c - (2*pi)/3))^2 + (Vm_a*sin(Va_a) + Vm_b*sin(Va_b + (2*pi)/3) + Vm_c*sin(Va_c - (2*pi)/3))^2)^2 - (2*Vm_c*sin(Va_c + (2*pi)/3)*(Vm_a*cos(Va_a) + Vm_b*cos(Va_b - (2*pi)/3) + Vm_c*cos(Va_c + (2*pi)/3)) - 2*Vm_c*cos(Va_c + (2*pi)/3)*(Vm_a*sin(Va_a) + Vm_b*sin(Va_b - (2*pi)/3) + Vm_c*sin(Va_c + (2*pi)/3)))/((Vm_a*cos(Va_a) + Vm_b*cos(Va_b + (2*pi)/3) + Vm_c*cos(Va_c - (2*pi)/3))^2 + (Vm_a*sin(Va_a) + Vm_b*sin(Va_b + (2*pi)/3) + Vm_c*sin(Va_c - (2*pi)/3))^2)]
            dum2 = [ 2*sin(Va_a)*(Vm_a*sin(Va_a) + Vm_b*sin(Va_b - (2*pi)/3) + Vm_c*sin(Va_c + (2*pi)/3)) + 2*cos(Va_a)*(Vm_a*cos(Va_a) + Vm_b*cos(Va_b - (2*pi)/3) + Vm_c*cos(Va_c + (2*pi)/3)), 2*sin(Va_b - (2*pi)/3)*(Vm_a*sin(Va_a) + Vm_b*sin(Va_b - (2*pi)/3) + Vm_c*sin(Va_c + (2*pi)/3)) + 2*cos(Va_b - (2*pi)/3)*(Vm_a*cos(Va_a) + Vm_b*cos(Va_b - (2*pi)/3) + Vm_c*cos(Va_c + (2*pi)/3)), 2*sin(Va_c + (2*pi)/3)*(Vm_a*sin(Va_a) + Vm_b*sin(Va_b - (2*pi)/3) + Vm_c*sin(Va_c + (2*pi)/3)) + 2*cos(Va_c + (2*pi)/3)*(Vm_a*cos(Va_a) + Vm_b*cos(Va_b - (2*pi)/3) + Vm_c*cos(Va_c + (2*pi)/3)), 2*Vm_a*cos(Va_a)*(Vm_a*sin(Va_a) + Vm_b*sin(Va_b - (2*pi)/3) + Vm_c*sin(Va_c + (2*pi)/3)) - 2*Vm_a*sin(Va_a)*(Vm_a*cos(Va_a) + Vm_b*cos(Va_b - (2*pi)/3) + Vm_c*cos(Va_c + (2*pi)/3)), 2*Vm_b*cos(Va_b - (2*pi)/3)*(Vm_a*sin(Va_a) + Vm_b*sin(Va_b - (2*pi)/3) + Vm_c*sin(Va_c + (2*pi)/3)) - 2*Vm_b*sin(Va_b - (2*pi)/3)*(Vm_a*cos(Va_a) + Vm_b*cos(Va_b - (2*pi)/3) + Vm_c*cos(Va_c + (2*pi)/3)), 2*Vm_c*cos(Va_c + (2*pi)/3)*(Vm_a*sin(Va_a) + Vm_b*sin(Va_b - (2*pi)/3) + Vm_c*sin(Va_c + (2*pi)/3)) - 2*Vm_c*sin(Va_c + (2*pi)/3)*(Vm_a*cos(Va_a) + Vm_b*cos(Va_b - (2*pi)/3) + Vm_c*cos(Va_c + (2*pi)/3))];
            VUF_lin[i] = dum1'*V_vec;
        end
        @objective(m, Min, M*sum(VUF_lin) )
    end
end
#----------------------------------------------------------------------------------------#

### FUNCTION TO DEFINE COST (RECTANGULAR)
function cost_rect(Vll)

    ## Minimize nothing
    if obj_op == 1 || obj_op == 6
        @objective(m, Min, 0*M )

    ## Minimize VUF
    elseif obj_op == 2
        @variable(m, Vn[i=1:n_VU], start = 0 );
        @variable(m, Vp[i=1:n_VU], start = 1 );
        for i=1:n_VU
            j=Int.(idx_bus_3ph1[i*3-2:i*3]);
            Vn_dq = a_rec[1:2,:]*[Vd[j]; Vq[j]]; Vp_dq = a_rec[3:4,:]*[Vd[j]; Vq[j]];
            @constraint(m, Vn[i] == Vn_dq[1]^2 + Vn_dq[2]^2)
            @constraint(m, Vp[i] == Vp_dq[1]^2 + Vp_dq[2]^2)
        end
        cost_VUF = @NLexpression(m, sum(Vn[i]/Vp[i] for i=1:n_VU) )
        @NLobjective(m, Min, M*cost_VUF + Q_pen*sum(QgY[i]^2 for i=1:ngy_inv) )

    ## Minimize power injections at substation
    elseif obj_op == 3
        Pij = @expression(m, [i=1:3], sum(Vd[idx_ref[i]]*( G[idx_ref[i],j]*Vd[j]-B[idx_ref[i],j]*Vq[j]) + Vq[idx_ref[i]]*(B[idx_ref[i],j]*Vd[j]+G[idx_ref[i],j]*Vq[j]) for j=1:nb_nph) )
        # Qij = @expression(m, [i=1:3], sum(Vd[idx_ref[i]]*(-B[idx_ref[i],j]*Vd[j]-G[idx_ref[i],j]*Vq[j]) + Vq[idx_ref[i]]*(G[idx_ref[i],j]*Vd[j]-B[idx_ref[i],j]*Vq[j]) for j=1:nb_nph) )
        global cost_SS = @expression(m, sum(Pij[i] for i=1:3) )
        @objective(m, Min, cost_SS )

    ## Minimize PVUR
    elseif obj_op == 4
        @variable(m, zp1[i=1:nb_3ph*3] )
        @variable(m, zp2[i=1:nb_3ph] )
        @expression(m, Vavg[i=1:nb_3ph], sum(Vm[j] for j=Int.(idx_bus_3ph1[i*3-2:i*3]))/3 )
        @expression(m, Vdev[i=1:nb_3ph*3], Vm[Int(idx_bus_3ph1[i])]-Vavg[Int(ceil(i/3))] )
         for i=1:nb_3ph
            @constraint(m, zp1[i*3-2:i*3] .>=  Vdev[i*3-2:i*3] )
            @constraint(m, zp1[i*3-2:i*3] .>= -Vdev[i*3-2:i*3] )
            @constraint(m, zp2[i] * Vavg[i] .>= zp1[i*3-2:i*3] )
        end
        global cost_PVUR = @expression(m, sum(zp2[i] for i=1:n_VU) )
        @objective(m, Min, cost_PVUR*M )

    ## Minimize LVUR
    elseif obj_op == 5
        @variable(m, zl1[i=1:nb_3ph*3] )
        @variable(m, zl2[i=1:nb_3ph] )
        @NLexpression(m, Vlavg[i=1:nb_3ph], sum(Vll[j] for j=Int.(idx_bus_3ph1[i*3-2:i*3]))/3 )
        @NLexpression(m, Vldev[i=1:nb_3ph*3], Vll[Int(idx_bus_3ph1[i])]-Vlavg[Int(ceil(i/3))] )
        for i=1:nb_3ph*3
            @NLconstraint(m, zl1[i] >=  Vldev[i] )
            @NLconstraint(m, zl1[i] >= -Vldev[i] )
        end
        for i=1:nb_3ph
            @NLconstraint(m, zl2[i] * Vlavg[i] >= zl1[i*3-2] )
            @NLconstraint(m, zl2[i] * Vlavg[i] >= zl1[i*3-1] )
            @NLconstraint(m, zl2[i] * Vlavg[i] >= zl1[i*3] )
        end
        global cost_LVUR = @expression(m, sum(zl2[i] for i=1:n_VU) )
        @objective(m, Min, cost_LVUR*M )

    # Minimize losses with unbalance constraints
    elseif obj_op == 6

        ## Calculate losses
        Pij = @expression(m, [i=1:nb_nph], sum(Vd[i]*( G[i,j]*Vd[j]-B[i,j]*Vq[j]) + Vq[i]*(B[i,j]*Vd[j]+G[i,j]*Vq[j]) for j=1:nb_nph) )
        global cost_loss = @expression(m, sum(Pij) )
        @objective(m, Min, cost_loss*M )

        # ## VUF Constraints
        # @NLexpression(m, Vn_all[i=1:nb_3ph], (sum(Vm[j]*cos(θ[j]+an[j]) for j=Int.(idx_bus_3ph2[i*3-2:i*3])))^2 + (sum(Vm[j]*sin(θ[j]+an[j]) for j=Int.(idx_bus_3ph2[i*3-2:i*3])))^2 )
        # @NLexpression(m, Vp_all[i=1:nb_3ph], (sum(Vm[j]*cos(θ[j]+ap[j]) for j=Int.(idx_bus_3ph2[i*3-2:i*3])))^2 + (sum(Vm[j]*sin(θ[j]+ap[j]) for j=Int.(idx_bus_3ph2[i*3-2:i*3])))^2 )
        # @NLexpression(m, V0_all[i=1:nb_3ph], (sum(Vm[j]*cos(θ[j]) for j=Int.(idx_bus_3ph1[i*3-2:i*3])))^2 + (sum(Vm[j]*sin(θ[j]) for j=Int.(idx_bus_3ph1[i*3-2:i*3])))^2 )
        # VUF_lim = 0.02
        # num_count = 1
        # for i in idx_bus_3ph
        #     @NLconstraint(m, Vn_all[num_count] <= Vp_all[num_count] * VUF_lim^2 )
        #     @NLconstraint(m, V0_all[num_count] <= Vp_all[num_count] * VUF_lim^2 )
        #     num_count = num_count + 1
        # end
        #
        # ## PVUR Constraints
        # PVUR_lim = 0.02
        # PVUR_mat = [PVUR_lim+2 PVUR_lim-1 PVUR_lim-1; PVUR_lim-1 PVUR_lim+2 PVUR_lim-1; PVUR_lim-1 PVUR_lim-1 PVUR_lim+2; PVUR_lim-2 PVUR_lim+1 PVUR_lim+1; PVUR_lim+1 PVUR_lim-2 PVUR_lim+1; PVUR_lim+1 PVUR_lim+1 PVUR_lim-2]
        # for i=1:nb_3ph
        #     j = Int.(idx_bus_3ph1[i*3-2:i*3])
        #     @constraint(m, PVUR_mat*Vm[j] .>= 0 )
        # end
        #
        # ## LVUR Constraints
        # @variable(m, Vll[i=1:nb_nph] >= 0, start = s0[i+2*nb_nph] )
        # for i=1:nb_nph
        #     k = Int(Tl[i])
        #     @NLconstraint(m, Vm[i]^2 + Vm[k]^2 - 2*Vm[i]*Vm[k]*cos(θ[i]-θ[k]) == Vll[i]^2 )
        # end
        # LVUR_lim = 0.03
        # LVUR_mat = [LVUR_lim+2 LVUR_lim-1 LVUR_lim-1; LVUR_lim-1 LVUR_lim+2 LVUR_lim-1; LVUR_lim-1 LVUR_lim-1 LVUR_lim+2; LVUR_lim-2 LVUR_lim+1 LVUR_lim+1; LVUR_lim+1 LVUR_lim-2 LVUR_lim+1; LVUR_lim+1 LVUR_lim+1 LVUR_lim-2]
        # for i=1:nb_3ph
        #     j = Int.(idx_bus_3ph1[i*3-2:i*3])
        #     @constraint(m, LVUR_mat*Vll[j] .>= 0 )
        # end

    ## Minimize Q-injections
    elseif obj_op == 7
        @objective(m, Min, sum(QgY[i] for i=1:ngy_inv) )

    ## Minimize linearized VUF
    elseif obj_op == 8
        @expression(m, VUF_lin[i=1:n_VU], 0 )
        for i=1:n_VU
            j=Int.(idx_bus_3ph1[i*3-2:i*3]); V_vec = [Vd[j];Vq[j]];
            Va_d = Vd0[j[1]]; Vb_d = Vd0[j[2]]; Vc_d = Vd0[j[3]];
            Va_q = Vq0[j[1]]; Vb_q = Vq0[j[2]]; Vc_q = Vq0[j[3]];
            dum1 = [ 2*Va_d - Vb_d - Vc_d + 3^(1/2)*Vb_q - 3^(1/2)*Vc_q, Vb_d/2 - Va_d + Vc_d/2 - (3^(1/2)*Vb_q)/2 + (3^(1/2)*Vc_q)/2 + 3^(1/2)*(Vb_q/2 - Va_q + Vc_q/2 + (3^(1/2)*Vb_d)/2 - (3^(1/2)*Vc_d)/2), Vb_d/2 - Va_d + Vc_d/2 - (3^(1/2)*Vb_q)/2 + (3^(1/2)*Vc_q)/2 - 3^(1/2)*(Vb_q/2 - Va_q + Vc_q/2 + (3^(1/2)*Vb_d)/2 - (3^(1/2)*Vc_d)/2), 2*Va_q - Vb_q - Vc_q - 3^(1/2)*Vb_d + 3^(1/2)*Vc_d, Vb_q/2 - Va_q + Vc_q/2 + (3^(1/2)*Vb_d)/2 - (3^(1/2)*Vc_d)/2 - 3^(1/2)*(Vb_d/2 - Va_d + Vc_d/2 - (3^(1/2)*Vb_q)/2 + (3^(1/2)*Vc_q)/2), Vb_q/2 - Va_q + Vc_q/2 + (3^(1/2)*Vb_d)/2 - (3^(1/2)*Vc_d)/2 + 3^(1/2)*(Vb_d/2 - Va_d + Vc_d/2 - (3^(1/2)*Vb_q)/2 + (3^(1/2)*Vc_q)/2)]
            VUF_lin[i] = dum1'*V_vec;
        end
        @objective(m, Min, 1e-10*sum(VUF_lin) )
    end
end
