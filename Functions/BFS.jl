### FUNCTION TO CALCULATE OPTIMAL POWER FLOW (BFS-rect)
function BFS_opf_rect(crit_node,Vm0,θ0,QgY0,dnc)

    ## Other inputs
    start_point(4); Vm0 = s0[1:nb_nph]; θ0 = s0[nb_nph+1:2*nb_nph];
    # (Vm0,θ0,QgY0,pf_iter) = FP_PF(Vm0,θ0,QgY0,50);
    Vd0 = Vm0.*cos.(θ0); Vq0 = Vm0.*sin.(θ0);

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
    Vll0 = zeros(nb_nph); opf_err = 1; opf_iter = 0;
    YloadP = zeros(nly_nph); YloadQ = zeros(nly_nph); stat = "LOCALLY_SOLVED"
    max_iter=30; DloadP = zeros(nld_nph); DloadQ = zeros(nld_nph); global err_vec = zeros(max_iter)
    I_ll = zeros(nld_nph)*(0+im*0); P_d2y = zeros(nldy_nph); Q_d2y = zeros(nldy_nph);
    while opf_err >= 1e-4 && opf_iter<max_iter && (stat == "LOCALLY_SOLVED" || stat == "ALMOST_LOCALLY_SOLVED")

        ## Setup model and check for slow convergence
        opf_iter = opf_iter + 1;
        # (Vm0,θ0,QgY0,pf_iter) = FP_pf(Vm0,θ0,QgY0,50); print("No of PF iterations: ", pf_iter,"\n");
        if opf_iter>4 && ((err_vec[opf_iter-2:opf_iter-1] == err_vec[opf_iter-4:opf_iter-3]) || (err_vec[opf_iter-1] == err_vec[opf_iter-2] && err_vec[opf_iter-1] == err_vec[opf_iter-3]) || (err_vec[opf_iter-1] >= err_vec[opf_iter-2] >= err_vec[opf_iter-3]))
            opf_iter = max_iter;
        elseif opf_iter>5 && err_vec[opf_iter-1] >10
            opf_iter = max_iter;
        else
            for i=1:nb_nph
                k = Int(Tl[i])
                Vll0[i] =  sqrt(Vm0[i]^2 + Vm0[k]^2 - 2*Vm0[i]*Vm0[k]*cos(θ0[i]-θ0[k]) )
            end
            global m = Model(optimizer_with_attributes(Ipopt.Optimizer,"warm_start_init_point"=>"yes","print_level"=>0, "tol"=>1e-7));

            ## Define variables
            global Vd = @variable(m, [i=1:nb_nph], start = Vd0[i] )
            global Vq = @variable(m, [i=1:nb_nph], start = Vq0[i] )
            global Vll = @variable(m, [i=1:nb_nph], start = Vll0[i] )
            global QgY = @variable(m, [i=1:ngy_inv], start = 0 )
            global Ibus_d = @variable(m,[i=1:nb_nph], start = 0 ); 
            global Ibus_q = @variable(m, [i=1:nb_nph], start = 0 );
            
            ## Define bounds on voltage magnitude and set slack bus angle to 0
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
            Ycap = im*Sh_cap.*Vm0.^2
            Yload = CblY*(YloadP + im*YloadQ) + Ycap;

            ## Calculate overall delta-connected loads in system
            for i=1:nld_nph
                j = Int(Tl[Int(Sd[i,7])]);
                dum1 = Vm0[Int(Sd[i,7])]*exp(im*θ0[Int(Sd[i,7])])-Vm0[j]*exp(im*θ0[j])
                dum2 = (Sd[i,3] + im*Sd[i,4])/(sqrt(3)*exp(im*Sd[i,8]))*dum1
                DloadP[i] = Sd[i,1]+real(dum2)+Sd[i,5]*Vll0[Int(Sd[i,7])]^2/3
                DloadQ[i] = Sd[i,2]+imag(dum2)+Sd[i,6]*Vll0[Int(Sd[i,7])]^2/3
            end
            num_count= 1;
            for i=1:nld
                j = Int(dload[i,1])
                for k=Int.(sum(bus_count[1:j-1])+1:sum(bus_count[1:j]))
                    kk = Int(Tl[k])
                    if k != kk
                        I_ll[num_count]= conj((DloadP[num_count]+im*DloadQ[num_count])/(Vm0[k]*exp(im*θ0[k])-Vm0[kk]*exp(im*θ0[kk])));
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
            Dload = CblD*(P_d2y + im*Q_d2y); YDload = Dload+Yload;

            ## Backward sweep
            V0 = Vm0.*exp.(im*θ0); S_inj = Dload+Yload; count3 = 1;
            Id_inj = ((real(YDload)-Cbg*PgY).*cos.(θ0) + (imag(YDload)-Cbg*QgY).*sin.(θ0))./Vm0;
            Iq_inj = ((real(YDload)-Cbg*PgY).*sin.(θ0) - (imag(YDload)-Cbg*QgY).*cos.(θ0))./Vm0;
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
                        @constraint(m, Ibus_d[ii] .== Id_inj[ii] + Isum_d ) 
                        @constraint(m, Ibus_q[ii] .== Iq_inj[ii] + Isum_q ) 
                    end   
                end
            end

            ## Forward sweep
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
                jj = Int.(jjjj); kk = Int.(kkkk); Zbr = -pinv(Y[jj,kk]);
                if tf_branch[i,3] == 8
                    @constraint(m, Vd[kk] .== vr_tap.*Vd[jj] - real(Zbr)*Ibus_d[kk] + imag(Zbr)*Ibus_q[kk] )
                    @constraint(m, Vq[kk] .== vr_tap.*Vq[jj] - real(Zbr)*Ibus_q[kk] - imag(Zbr)*Ibus_d[kk] )
                    count2 = count2 + 1;
                else
                    @constraint(m, Vd[kk] .== Vd[jj] - real(Zbr)*Ibus_d[kk] + imag(Zbr)*Ibus_q[kk] )
                    @constraint(m, Vq[kk] .== Vq[jj] - real(Zbr)*Ibus_q[kk] - imag(Zbr)*Ibus_d[kk] )
                end
                # Vm_jj = Vm1[jj].*exp.(im*θ1[jj]); ii = Int(line_prove[i,4])
                # if bus_ph[k] == "AN" || bus_ph[k] == "AS" || bus_ph[k] == "A"
                #     Vm_kk = A_mat[1,1,ii]*Vm_jj - B_mat[1,1,ii]*Ibr[count1:count1+length(kk)-1]
                # elseif bus_ph[k] == "BN" || bus_ph[k] == "BS" || bus_ph[k] == "B"
                #     Vm_kk = A_mat[2,2,ii]*Vm_jj - B_mat[2,2,ii]*Ibr[count1:count1+length(kk)-1]
                # elseif bus_ph[k] == "CN" || bus_ph[k] == "CS" || bus_ph[k] == "C"
                #     Vm_kk = A_mat[3,3,ii]*Vm_jj - B_mat[3,3,ii]*Ibr[count1:count1+length(kk)-1]
                # elseif bus_ph[k] == "ABN" || bus_ph[k] == "ABD" || bus_ph[k] == "BAN" || bus_ph[k] == "BAD" || bus_ph[k] == "BA" || bus_ph[k] == "AB"
                #     Vm_kk = A_mat[1:2,1:2,ii]*Vm_jj - B_mat[1:2,1:2,ii]*Ibr[count1:count1+length(kk)-1]
                # elseif bus_ph[k] == "BCN" || bus_ph[k] == "BCD" || bus_ph[k] == "CBN"  || bus_ph[k] == "CBD" || bus_ph[k] == "CB" || bus_ph[k] == "BC"
                #     Vm_kk = A_mat[2:3,2:3,ii]*Vm_jj - B_mat[2:3,2:3,ii]*Ibr[count1:count1+length(kk)-1]
                # elseif bus_ph[k] == "ACN" || bus_ph[k] == "ACD" || bus_ph[k] == "CAN" || bus_ph[k] == "CAD" || bus_ph[k] == "CA" || bus_ph[k] == "AC"
                #     Vm_kk = A_mat[1:2:3,1:2:3,ii]*Vm_jj - B_mat[1:2:3,1:2:3,ii]*Ibr[count1:count1+length(kk)-1]
                # else
                #     Vm_kk = A_mat[:,:,ii]*Vm_jj - B_mat[:,:,ii]*Ibr[count1:count1+length(kk)-1]
                # end
            end

            ## Define cost
            cost_rect(Vll);

            ## Solve and update solution
            @time status = optimize!(m); sts = string(termination_status(m))
            V0 = Vm0.*exp.(im*θ0); Vd0 = JuMP.value.(Vd); Vq0 = JuMP.value.(Vq); QgY0 = JuMP.value.(QgY);
            Vm0 = abs.(Vd0+im*Vq0); θ0 = angle.(Vd0+im*Vq0);
            V1 = Vm0.*exp.(im*θ0); opf_err = sum(abs.(V1-V0)); err_vec[opf_iter] = round(opf_err, RoundNearest, digits=4);
            print("BFS_OPF_rect: (",opf_iter,") Error = ",opf_err,"\n")
        end
    end

    ## Run power flow to convergence if OPF did not converge
    if opf_iter == max_iter
        println("\n****************OPF DID NOT CONVERGE!!*******************\n")
        (Vm0,θ0,QgY0,pf_iter) = FP_pf(Vm0,θ0,QgY0,50); print("No of PF iterations: ", pf_iter,"\n");
        if dnc == 0
            sts = "DID NOT CONVERGE";
        end
    end

    ## Calculate missing solutions
    Vph = Vm0.*exp.(im*θ0);
    Iph = Y*Vph; S_inj = Vph.*conj(Iph);
    P_inj = real(S_inj); global P_loss = sum(P_inj)*baseMVA*1e3
    for i=1:nb_nph
        k = Int(Tl[i])
        Vll0[i] =  sqrt(Vm0[i]^2 + Vm0[k]^2 - 2*Vm0[i]*Vm0[k]*cos(θ0[i]-θ0[k]) )
    end
    print_VU(Vm0,θ0,Vll0);
    println("TOTAL TIME:")
    return Vm0,θ0,QgY0,sts
end

#----------------------------------------------------------------------------------------#

### FUNCTION TO CALCULATE OPTIMAL POWER FLOW (BFS-rect)
function BFS_opf_rect1(crit_node,Vm0,θ0,QgY0,dnc)

    ## Other inputs
    start_point(4); Vm0 = s0[1:nb_nph]; θ0 = s0[nb_nph+1:2*nb_nph];
    # (Vm0,θ0,QgY0,pf_iter) = FP_PF(Vm0,θ0,QgY0,50);
    Vd0 = Vm0.*cos.(θ0); Vq0 = Vm0.*sin.(θ0);

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
    Vll0 = zeros(nb_nph); opf_err = 1; opf_iter = 0;
    YloadP = zeros(nly_nph); YloadQ = zeros(nly_nph); stat = "LOCALLY_SOLVED"
    max_iter=30; DloadP = zeros(nld_nph); DloadQ = zeros(nld_nph); global err_vec = zeros(max_iter)
    I_ll = zeros(nld_nph)*(0+im*0); P_d2y = zeros(nldy_nph); Q_d2y = zeros(nldy_nph);
    while opf_err >= 1e-4 && opf_iter<max_iter && (stat == "LOCALLY_SOLVED" || stat == "ALMOST_LOCALLY_SOLVED")

        ## Setup model and check for slow convergence
        opf_iter = opf_iter + 1;
        # (Vm0,θ0,QgY0,pf_iter) = FP_pf(Vm0,θ0,QgY0,50); print("No of PF iterations: ", pf_iter,"\n");
        if opf_iter>4 && ((err_vec[opf_iter-2:opf_iter-1] == err_vec[opf_iter-4:opf_iter-3]) || (err_vec[opf_iter-1] == err_vec[opf_iter-2] && err_vec[opf_iter-1] == err_vec[opf_iter-3]) || (err_vec[opf_iter-1] >= err_vec[opf_iter-2] >= err_vec[opf_iter-3]))
            opf_iter = max_iter;
        elseif opf_iter>5 && err_vec[opf_iter-1] >10
            opf_iter = max_iter;
        else
            for i=1:nb_nph
                k = Int(Tl[i])
                Vll0[i] =  sqrt(Vm0[i]^2 + Vm0[k]^2 - 2*Vm0[i]*Vm0[k]*cos(θ0[i]-θ0[k]) )
            end
            global m = Model(optimizer_with_attributes(Ipopt.Optimizer,"warm_start_init_point"=>"yes","print_level"=>0, "tol"=>1e-7));

            ## Define variables
            global Vd = @variable(m, [i=1:nb_nph], start = Vd0[i] )
            global Vq = @variable(m, [i=1:nb_nph], start = Vq0[i] )
            global Vll = @variable(m, [i=1:nb_nph], start = Vll0[i] )
            global QgY = @variable(m, [i=1:ngy_inv], start = 0 )     

            ## Define bounds on voltage magnitude and set slack bus angle to 0
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
            Ycap = im*Sh_cap.*Vm0.^2
            Yload = CblY*(YloadP + im*YloadQ) + Ycap;

            ## Calculate overall delta-connected loads in system
            for i=1:nld_nph
                j = Int(Tl[Int(Sd[i,7])]);
                dum1 = Vm0[Int(Sd[i,7])]*exp(im*θ0[Int(Sd[i,7])])-Vm0[j]*exp(im*θ0[j])
                dum2 = (Sd[i,3] + im*Sd[i,4])/(sqrt(3)*exp(im*Sd[i,8]))*dum1
                DloadP[i] = Sd[i,1]+real(dum2)+Sd[i,5]*Vll0[Int(Sd[i,7])]^2/3
                DloadQ[i] = Sd[i,2]+imag(dum2)+Sd[i,6]*Vll0[Int(Sd[i,7])]^2/3
            end
            num_count= 1;
            for i=1:nld
                j = Int(dload[i,1])
                for k=Int.(sum(bus_count[1:j-1])+1:sum(bus_count[1:j]))
                    kk = Int(Tl[k])
                    if k != kk
                        I_ll[num_count]= conj((DloadP[num_count]+im*DloadQ[num_count])/(Vm0[k]*exp(im*θ0[k])-Vm0[kk]*exp(im*θ0[kk])));
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
            Dload = CblD*(P_d2y + im*Q_d2y); YDload = Dload+Yload;

            ## Backward sweep
            V0 = Vm0.*exp.(im*θ0); count3 = 1;
            Sd_inj = real(YDload)-Cbg*PgY; Sq_inj = imag(YDload)-Cbg*QgY;
            Sbus_d = @expression(m,[i=1:nb_nph], 0*Vd[1]); Sbus_q = @expression(m, [i=1:nb_nph], 0*Vq[1]);
            for i=1:nbr
                j = Int(line_prove[nbr-i+1,2]); jj = Int(sum(bus_count[1:j-1])+1):Int(sum(bus_count[1:j]))
                k = Int(line_prove[nbr-i+1,1]); kk = Int(sum(bus_count[1:k-1])+1):Int(sum(bus_count[1:k]))
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
                Ybr = -Y[jj,kk]; Sbr = (V0[jj]-V0[kk]).*conj(Ybr*(V0[jj]-V0[kk]));
                Sbus_d[kk] = Sbus_d[kk] + Sd_inj[kk]; 
                Sbus_q[kk] = Sbus_q[kk] + Sq_inj[kk];
                Sbus_d[jj] = Sbus_d[kk] + Sbus_d[jj] + real(Sbr); 
                Sbus_q[jj] = Sbus_q[kk] + Sbus_q[jj] + imag(Sbr); 
            end
            Ibus_d = (Sbus_d.*cos.(θ0) + Sbus_q.*sin.(θ0))./Vm0;
            Ibus_q = (Sbus_d.*sin.(θ0) - Sbus_q.*cos.(θ0))./Vm0;

            ## Forward sweep
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
                jj = Int.(jjjj); kk = Int.(kkkk); Zbr = -pinv(Y[jj,kk]);
                if tf_branch[i,3] == 8
                    if bus_ph[k] == "AN"
                        vr_tap = 1 + 0.00625*tf_branch[vr_idx[count2],6]
                    elseif bus_ph[k] == "BN"
                        vr_tap = 1 + 0.00625*tf_branch[vr_idx[count2],7]
                    elseif bus_ph[k] == "CN"
                        vr_tap = 1 + 0.00625*tf_branch[vr_idx[count2],8]
                    elseif bus_ph[k] == "ABN" || bus_ph[k] == "BAN"
                        vr_tap = Diagonal(1 .+ 0.00625*tf_branch[vr_idx[count2],6:7])
                    elseif bus_ph[k] == "BCN" || bus_ph[k] == "CBN"
                        vr_tap = Diagonal(1 .+ 0.00625*tf_branch[vr_idx[count2],7:8])
                    elseif bus_ph[k] == "ACN" || bus_ph[k] == "CAN"
                        vr_tap = Diagonal(1 .+ 0.00625*tf_branch[vr_idx[count2],6:2:8])
                    else
                        vr_tap = Diagonal(1 .+ 0.00625*tf_branch[vr_idx[count2],6:8])
                    end
                    @constraint(m, Vd[kk] .== vr_tap*(Vd[jj] - real(Zbr)*Ibus_d[kk] + imag(Zbr)*Ibus_q[kk]) )
                    @constraint(m, Vq[kk] .== vr_tap*(Vq[jj] - real(Zbr)*Ibus_q[kk] - imag(Zbr)*Ibus_d[kk]) )
                    count2 = count2 + 1;
                else
                    @constraint(m, Vd[kk] .== Vd[jj] - real(Zbr)*Ibus_d[kk] + imag(Zbr)*Ibus_q[kk] )
                    @constraint(m, Vq[kk] .== Vq[jj] - real(Zbr)*Ibus_q[kk] - imag(Zbr)*Ibus_d[kk] )
                end
                # Vm_jj = Vm1[jj].*exp.(im*θ1[jj]); ii = Int(line_prove[i,4])
                # if bus_ph[k] == "AN" || bus_ph[k] == "AS" || bus_ph[k] == "A"
                #     Vm_kk = A_mat[1,1,ii]*Vm_jj - B_mat[1,1,ii]*Ibr[count1:count1+length(kk)-1]
                # elseif bus_ph[k] == "BN" || bus_ph[k] == "BS" || bus_ph[k] == "B"
                #     Vm_kk = A_mat[2,2,ii]*Vm_jj - B_mat[2,2,ii]*Ibr[count1:count1+length(kk)-1]
                # elseif bus_ph[k] == "CN" || bus_ph[k] == "CS" || bus_ph[k] == "C"
                #     Vm_kk = A_mat[3,3,ii]*Vm_jj - B_mat[3,3,ii]*Ibr[count1:count1+length(kk)-1]
                # elseif bus_ph[k] == "ABN" || bus_ph[k] == "ABD" || bus_ph[k] == "BAN" || bus_ph[k] == "BAD" || bus_ph[k] == "BA" || bus_ph[k] == "AB"
                #     Vm_kk = A_mat[1:2,1:2,ii]*Vm_jj - B_mat[1:2,1:2,ii]*Ibr[count1:count1+length(kk)-1]
                # elseif bus_ph[k] == "BCN" || bus_ph[k] == "BCD" || bus_ph[k] == "CBN"  || bus_ph[k] == "CBD" || bus_ph[k] == "CB" || bus_ph[k] == "BC"
                #     Vm_kk = A_mat[2:3,2:3,ii]*Vm_jj - B_mat[2:3,2:3,ii]*Ibr[count1:count1+length(kk)-1]
                # elseif bus_ph[k] == "ACN" || bus_ph[k] == "ACD" || bus_ph[k] == "CAN" || bus_ph[k] == "CAD" || bus_ph[k] == "CA" || bus_ph[k] == "AC"
                #     Vm_kk = A_mat[1:2:3,1:2:3,ii]*Vm_jj - B_mat[1:2:3,1:2:3,ii]*Ibr[count1:count1+length(kk)-1]
                # else
                #     Vm_kk = A_mat[:,:,ii]*Vm_jj - B_mat[:,:,ii]*Ibr[count1:count1+length(kk)-1]
                # end
            end

            ## Define cost
            cost_rect(Vll);

            ## Solve and update solution
            @time status = optimize!(m); sts = string(termination_status(m))
            V0 = Vm0.*exp.(im*θ0); Vd0 = JuMP.value.(Vd); Vq0 = JuMP.value.(Vq); QgY0 = JuMP.value.(QgY);
            Vm0 = abs.(Vd0+im*Vq0); θ0 = angle.(Vd0+im*Vq0);
            V1 = Vm0.*exp.(im*θ0); opf_err = sum(abs.(V1-V0)); err_vec[opf_iter] = round(opf_err, RoundNearest, digits=4);
            print("BFS_OPF_rect: (",opf_iter,") Error = ",opf_err,"\n")
        end
    end

    ## Run power flow to convergence if OPF did not converge
    if opf_iter == max_iter
        println("\n****************OPF DID NOT CONVERGE!!*******************\n")
        (Vm0,θ0,QgY0,pf_iter) = FP_pf(Vm0,θ0,QgY0,50); print("No of PF iterations: ", pf_iter,"\n");
        if dnc == 0
            sts = "DID NOT CONVERGE";
        end
    end

    ## Calculate missing solutions
    Vph = Vm0.*exp.(im*θ0);
    Iph = Y*Vph; S_inj = Vph.*conj(Iph);
    P_inj = real(S_inj); global P_loss = sum(P_inj)*baseMVA*1e3
    for i=1:nb_nph
        k = Int(Tl[i])
        Vll0[i] =  sqrt(Vm0[i]^2 + Vm0[k]^2 - 2*Vm0[i]*Vm0[k]*cos(θ0[i]-θ0[k]) )
    end
    print_VU(Vm0,θ0,Vll0);
    println("TOTAL TIME:")
    return Vm0,θ0,QgY0,sts
end

### FUNCTION TO CALCULATE OPTIMAL POWER FLOW (BFS-rect)
function BFS_opf_rect2(crit_node,Vm0,θ0,QgY0,dnc)

    ## Other inputs
    start_point(4); Vm0 = s0[1:nb_nph]; θ0 = s0[nb_nph+1:2*nb_nph];
    # (Vm0,θ0,QgY0,pf_iter) = FP_PF(Vm0,θ0,QgY0,50);
    Vd0 = Vm0.*cos.(θ0); Vq0 = Vm0.*sin.(θ0);

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
    Vll0 = zeros(nb_nph); opf_err = 1; opf_iter = 0;
    YloadP = zeros(nly_nph); YloadQ = zeros(nly_nph); stat = "LOCALLY_SOLVED"
    max_iter=30; DloadP = zeros(nld_nph); DloadQ = zeros(nld_nph); global err_vec = zeros(max_iter)
    I_ll = zeros(nld_nph)*(0+im*0); P_d2y = zeros(nldy_nph); Q_d2y = zeros(nldy_nph);
    while opf_err >= 1e-4 && opf_iter<max_iter && (stat == "LOCALLY_SOLVED" || stat == "ALMOST_LOCALLY_SOLVED")

        ## Setup model and check for slow convergence
        opf_iter = opf_iter + 1;
        # (Vm0,θ0,QgY0,pf_iter) = FP_pf(Vm0,θ0,QgY0,50); print("No of PF iterations: ", pf_iter,"\n");
        if opf_iter>4 && ((err_vec[opf_iter-2:opf_iter-1] == err_vec[opf_iter-4:opf_iter-3]) || (err_vec[opf_iter-1] == err_vec[opf_iter-2] && err_vec[opf_iter-1] == err_vec[opf_iter-3]) || (err_vec[opf_iter-1] >= err_vec[opf_iter-2] >= err_vec[opf_iter-3]))
            opf_iter = max_iter;
        elseif opf_iter>5 && err_vec[opf_iter-1] >10
            opf_iter = max_iter;
        else
            for i=1:nb_nph
                k = Int(Tl[i])
                Vll0[i] =  sqrt(Vm0[i]^2 + Vm0[k]^2 - 2*Vm0[i]*Vm0[k]*cos(θ0[i]-θ0[k]) )
            end
            global m = Model(optimizer_with_attributes(Ipopt.Optimizer,"warm_start_init_point"=>"yes","print_level"=>0, "tol"=>1e-7));

            ## Define variables
            global Vd = @variable(m, [i=1:nb_nph], start = Vd0[i] )
            global Vq = @variable(m, [i=1:nb_nph], start = Vq0[i] )
            global Vll = @variable(m, [i=1:nb_nph], start = Vll0[i] )
            global QgY = @variable(m, [i=1:ngy_inv], start = 0 )

            ## Define bounds on voltage magnitude and set slack bus angle to 0
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
            Ycap = Vm0.*exp.(im*θ0).*conj.(im*Sh_cap.*Vm0.*exp.(im*θ0)).*Vm0.^2
            Yload = CblY*(YloadP + im*YloadQ) + Ycap;

            ## Calculate overall delta-connected loads in system
            for i=1:nld_nph
                j = Int(Tl[Int(Sd[i,7])]);
                dum1 = Vm0[Int(Sd[i,7])]*exp(im*θ0[Int(Sd[i,7])])-Vm0[j]*exp(im*θ0[j])
                dum2 = (Sd[i,3] + im*Sd[i,4])/(sqrt(3)*exp(im*Sd[i,8]))*dum1
                DloadP[i] = Sd[i,1]+real(dum2)+Sd[i,5]*Vll0[Int(Sd[i,7])]^2/3
                DloadQ[i] = Sd[i,2]+imag(dum2)+Sd[i,6]*Vll0[Int(Sd[i,7])]^2/3
            end
            num_count= 1;
            for i=1:nld
                j = Int(dload[i,1])
                for k=Int.(sum(bus_count[1:j-1])+1:sum(bus_count[1:j]))
                    kk = Int(Tl[k])
                    if k != kk
                        I_ll[num_count]= conj((DloadP[num_count]+im*DloadQ[num_count])/(Vm0[k]*exp(im*θ0[k])-Vm0[kk]*exp(im*θ0[kk])));
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
            Dload = CblD*(P_d2y + im*Q_d2y); YDload = Dload+Yload;

            ## Backward sweep
            V0 = Vm0.*exp.(im*θ0); S_inj = Dload+Yload; count3 = 1;
            Id_inj = ((real(YDload)-Cbg*PgY).*cos.(θ0) + (imag(YDload)-Cbg*QgY).*sin.(θ0))./Vm0;
            Iq_inj = ((real(YDload)-Cbg*PgY).*sin.(θ0) - (imag(YDload)-Cbg*QgY).*cos.(θ0))./Vm0;
            Ibus_d = @expression(m,[i=1:nb_nph], 0*Vd[1]); Ibus_q = @expression(m, [i=1:nb_nph], 0*Vq[1]);
            for i=1:nbr
                j = Int(line_prove[nbr-i+1,2]); jj = Int(sum(bus_count[1:j-1])+1):Int(sum(bus_count[1:j]))
                k = Int(line_prove[nbr-i+1,1]); kk = Int(sum(bus_count[1:k-1])+1):Int(sum(bus_count[1:k]))
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
                jj = Int.(jjjj); kk = Int.(kkkk)
                if tf_branch[Int(line_prove[nbr-i+1,4]),3] == 8
                    if bus_ph[k] == "AN"
                        global vr_tap = 1 + 0.00625*tf_branch[vr_idx[count3],6]
                    elseif bus_ph[k] == "BN"
                        global vr_tap = 1 + 0.00625*tf_branch[vr_idx[count3],7]
                    elseif bus_ph[k] == "CN"
                        global vr_tap = 1 + 0.00625*tf_branch[vr_idx[count3],8]
                    elseif bus_ph[k] == "ABN" || bus_ph[k] == "BAN"
                        global vr_tap = Diagonal(1 .+ 0.00625*tf_branch[vr_idx[count3],6:7])
                    elseif bus_ph[k] == "BCN" || bus_ph[k] == "CBN"
                        global vr_tap = Diagonal(1 .+ 0.00625*tf_branch[vr_idx[count3],7:8])
                    elseif bus_ph[k] == "ACN" || bus_ph[k] == "CAN"
                        global vr_tap = Diagonal(1 .+ 0.00625*tf_branch[vr_idx[count3],6:2:8])
                    else
                        global vr_tap = Diagonal(1 .+ 0.00625*tf_branch[vr_idx[count3],6:8])
                    end
                        Ibus_d[kk] = vr_tap*(Ibus_d[kk] + Id_inj[kk]); 
                        Ibus_q[kk] = vr_tap*(Ibus_q[kk] + Iq_inj[kk]);
                        count3 = count3+1;
                else
                    Ibus_d[kk] = Ibus_d[kk] + Id_inj[kk]; 
                    Ibus_q[kk] = Ibus_q[kk] + Iq_inj[kk];
                end
                Ibus_d[jj] = Ibus_d[kk] + Ibus_d[jj]; 
                Ibus_q[jj] = Ibus_q[kk] + Ibus_q[jj];
            end

            ## Forward sweep
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
                jj = Int.(jjjj); kk = Int.(kkkk); Zbr = -pinv(Y[jj,kk]);
                if tf_branch[i,3] == 8
                    @constraint(m, Vd[kk] .== vr_tap*Vd[jj] - real(Zbr)*Ibus_d[kk] + imag(Zbr)*Ibus_q[kk] )
                    @constraint(m, Vq[kk] .== vr_tap*Vq[jj] - real(Zbr)*Ibus_q[kk] - imag(Zbr)*Ibus_d[kk] )
                    count2 = count2 + 1;
                else
                    @constraint(m, Vd[kk] .== Vd[jj] - real(Zbr)*Ibus_d[kk] + imag(Zbr)*Ibus_q[kk] )
                    @constraint(m, Vq[kk] .== Vq[jj] - real(Zbr)*Ibus_q[kk] - imag(Zbr)*Ibus_d[kk] )
                end
                # Vm_jj = Vm1[jj].*exp.(im*θ1[jj]); ii = Int(line_prove[i,4])
                # if bus_ph[k] == "AN" || bus_ph[k] == "AS" || bus_ph[k] == "A"
                #     Vm_kk = A_mat[1,1,ii]*Vm_jj - B_mat[1,1,ii]*Ibr[count1:count1+length(kk)-1]
                # elseif bus_ph[k] == "BN" || bus_ph[k] == "BS" || bus_ph[k] == "B"
                #     Vm_kk = A_mat[2,2,ii]*Vm_jj - B_mat[2,2,ii]*Ibr[count1:count1+length(kk)-1]
                # elseif bus_ph[k] == "CN" || bus_ph[k] == "CS" || bus_ph[k] == "C"
                #     Vm_kk = A_mat[3,3,ii]*Vm_jj - B_mat[3,3,ii]*Ibr[count1:count1+length(kk)-1]
                # elseif bus_ph[k] == "ABN" || bus_ph[k] == "ABD" || bus_ph[k] == "BAN" || bus_ph[k] == "BAD" || bus_ph[k] == "BA" || bus_ph[k] == "AB"
                #     Vm_kk = A_mat[1:2,1:2,ii]*Vm_jj - B_mat[1:2,1:2,ii]*Ibr[count1:count1+length(kk)-1]
                # elseif bus_ph[k] == "BCN" || bus_ph[k] == "BCD" || bus_ph[k] == "CBN"  || bus_ph[k] == "CBD" || bus_ph[k] == "CB" || bus_ph[k] == "BC"
                #     Vm_kk = A_mat[2:3,2:3,ii]*Vm_jj - B_mat[2:3,2:3,ii]*Ibr[count1:count1+length(kk)-1]
                # elseif bus_ph[k] == "ACN" || bus_ph[k] == "ACD" || bus_ph[k] == "CAN" || bus_ph[k] == "CAD" || bus_ph[k] == "CA" || bus_ph[k] == "AC"
                #     Vm_kk = A_mat[1:2:3,1:2:3,ii]*Vm_jj - B_mat[1:2:3,1:2:3,ii]*Ibr[count1:count1+length(kk)-1]
                # else
                #     Vm_kk = A_mat[:,:,ii]*Vm_jj - B_mat[:,:,ii]*Ibr[count1:count1+length(kk)-1]
                # end
            end

            ## Define cost
            cost_rect(Vll);

            ## Solve and update solution
            @time status = optimize!(m); sts = string(termination_status(m))
            V0 = Vm0.*exp.(im*θ0); Vd0 = JuMP.value.(Vd); Vq0 = JuMP.value.(Vq); QgY0 = JuMP.value.(QgY);
            Vm0 = abs.(Vd0+im*Vq0); θ0 = angle.(Vd0+im*Vq0);
            V1 = Vm0.*exp.(im*θ0); opf_err = sum(abs.(V1-V0)); err_vec[opf_iter] = round(opf_err, RoundNearest, digits=4);
            print("BFS_OPF_rect: (",opf_iter,") Error = ",opf_err,"\n")
        end
    end

    ## Run power flow to convergence if OPF did not converge
    if opf_iter == max_iter
        println("\n****************OPF DID NOT CONVERGE!!*******************\n")
        (Vm0,θ0,QgY0,pf_iter) = FP_pf(Vm0,θ0,QgY0,50); print("No of PF iterations: ", pf_iter,"\n");
        if dnc == 0
            sts = "DID NOT CONVERGE";
        end
    end

    ## Calculate missing solutions
    Vph = Vm0.*exp.(im*θ0);
    Iph = Y*Vph; S_inj = Vph.*conj(Iph);
    P_inj = real(S_inj); global P_loss = sum(P_inj)*baseMVA*1e3
    for i=1:nb_nph
        k = Int(Tl[i])
        Vll0[i] =  sqrt(Vm0[i]^2 + Vm0[k]^2 - 2*Vm0[i]*Vm0[k]*cos(θ0[i]-θ0[k]) )
    end
    print_VU(Vm0,θ0,Vll0);
    println("TOTAL TIME:")
    return Vm0,θ0,QgY0,sts
end
