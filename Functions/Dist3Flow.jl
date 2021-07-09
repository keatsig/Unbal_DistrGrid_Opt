### FUNCTION TO CALCULATE OPTIMAL POWER FLOW (DistFlow)
function D3F_opf_rect(crit_node,Vm0,θ0,QgY0,dnc)

    ## Other inputs
    # start_point(1); Vm0 = s0[1:nb_nph]; θ0 = s0[nb_nph+1:2*nb_nph];
    # (Vm0,θ0,QgY0,stat) = TP_OPF(crit_node)
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

    ## Start PF iterations
    Vll0 = zeros(nb_nph); opf_err = 1; opf_iter = 0;
    YloadP = zeros(nly_nph); YloadQ = zeros(nly_nph); stat = "LOCALLY_SOLVED"
    max_iter=1; DloadP = zeros(nld_nph); DloadQ = zeros(nld_nph); global err_vec = zeros(max_iter)
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
            V0 = Vm0.*exp.(im*θ0); y0 = Vm0.^2; γ=V0./conj(V0');
            for i=1:nb_nph
                k = Int(Tl[i])
                Vll0[i] =  sqrt(Vm0[i]^2 + Vm0[k]^2 - 2*Vm0[i]*Vm0[k]*cos(θ0[i]-θ0[k]) )
            end
            global m = Model(optimizer_with_attributes(Ipopt.Optimizer,"warm_start_init_point"=>"yes","print_level"=>0, "tol"=>1e-7));

            ## Define variables
            global Vd = @variable(m, [i=1:nb_nph], start = Vd0[i] )
            global Vq = @variable(m, [i=1:nb_nph], start = Vq0[i] )
            global Vll = @variable(m, [i=1:nb_nph], start = Vll0[i] )
            global y = @variable(m, [i=1:nb_nph], start = y0[i] )
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
            Ygen_d = Cbg*PgY; Ygen_q = Cbg*QgY;

            ## Calculate wye-connected loads
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

            ## Calculate delta-connected loads
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

            ## DistFlow
            V0 = Vm0.*exp.(im*θ0); count3 = 1; sd_inj = real(Dload+Yload)-Ygen_d; sq_inj = imag(Dload+Yload)-Ygen_q;
            Id_inj = ((real(YDload)-Cbg*PgY).*cos.(θ0) + (imag(YDload)-Cbg*QgY).*sin.(θ0))./Vm0;
            Iq_inj = ((real(YDload)-Cbg*PgY).*sin.(θ0) - (imag(YDload)-Cbg*QgY).*cos.(θ0))./Vm0;
            Ibus_d = @expression(m,[i=1:nb_nph], 0*Vd[1]); Ibus_q = @expression(m, [i=1:nb_nph], 0*Vq[1]);
            Sbus_d = @expression(m,[i=1:nb_nph], 0*Vd[1]); Sbus_q = @expression(m, [i=1:nb_nph], 0*Vq[1]);
            Ibus0 = zeros(nb_nph)*(0+im*0);
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
                Zbr = -pinv(Y[jj,kk]);
                if tf_branch[Int(line_prove[nbr-i+1,4]),3] == 8
                    if bus_ph[k] == "AN"
                        vr_tap = 1 + 0.00625*tf_branch[vr_idx[count3],6]
                    elseif bus_ph[k] == "BN"
                        vr_tap = 1 + 0.00625*tf_branch[vr_idx[count3],7]
                    elseif bus_ph[k] == "CN"
                        vr_tap = 1 + 0.00625*tf_branch[vr_idx[count3],8]
                    elseif bus_ph[k] == "ABN" || bus_ph[k] == "BAN"
                        vr_tap = Diagonal(1 .+ 0.00625*tf_branch[vr_idx[count3],6:7])
                    elseif bus_ph[k] == "BCN" || bus_ph[k] == "CBN"
                        vr_tap = Diagonal(1 .+ 0.00625*tf_branch[vr_idx[count3],7:8])
                    elseif bus_ph[k] == "ACN" || bus_ph[k] == "CAN"
                        vr_tap = Diagonal(1 .+ 0.00625*tf_branch[vr_idx[count3],6:2:8])
                    else
                        vr_tap = Diagonal(1 .+ 0.00625*tf_branch[vr_idx[count3],6:8])
                    end
                    Ibus_d[kk] = vr_tap*(Ibus_d[kk] + Id_inj[kk]); Ibus_q[kk] = vr_tap*(Ibus_q[kk] + Iq_inj[kk]);
                    count3 = count3+1;
                else
                    Ibus_d[kk] = Ibus_d[kk] + Id_inj[kk]; Ibus_q[kk] = Ibus_q[kk] + Iq_inj[kk];
                end
                Lbr = diag(Zbr*Ibus0[kk]*Ibus0[kk]');
                Ibus_d[jj] = Ibus_d[kk] + Ibus_d[jj]; Ibus_q[jj] = Ibus_q[kk] + Ibus_q[jj];
                Sbus_d[kk] = Sbus_d[kk] + sd_inj[kk]; Sbus_q[kk] = Sbus_q[kk] + sq_inj[kk];
                Sbus_d[jj] = Sbus_d[kk] + Sbus_d[jj] + real(Lbr); Sbus_q[jj] = Sbus_q[kk] + Sbus_q[jj] + imag(Lbr);
            end
            count2=1;
            for i=1:nbr
                j = Int(line_prove[i,2]); jj = Int(sum(bus_count[1:j-1])+1):Int(sum(bus_count[1:j]))
                k = Int(line_prove[i,1]); kk = Int(sum(bus_count[1:k-1])+1):Int(sum(bus_count[1:k]))
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
                jj = Int.(jjjj); kk = Int.(kkkk); γbr = γ[jj,kk];
                Zbr = -pinv(Y[jj,kk]); Hbr = real(diag(Zbr*Ibus0[kk]*Ibus0[kk]'*Zbr')); Zγ=conj(Zbr).*γbr;
                if tf_branch[Int(line_prove[i,4]),3] == 8
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
                    @constraint(m, y[kk] .== vr_tap.^2*y[jj]-Hbr-2*(real(Zγ)*Sbus_d[kk]-imag(Zγ)*Sbus_q[kk]) );
                    count2 = count2+1;
                else
                    @constraint(m, y[kk] .== y[jj]-Hbr-2*(real(Zγ)*Sbus_d[kk]-imag(Zγ)*Sbus_q[kk]) );
                end
            end

            ## Define cost
            cost_rect(Vll);

            ## Solve and update solution
            @time status = optimize!(m); sts = string(termination_status(m));
            QgY0 = JuMP.value.(QgY); Vm0 = sqrt.(JuMP.value.(y));
            (Vm0,θ0,QgY0,pf_iter) = FP_pf(Vm0,θ0,QgY0,50); print("No of PF iterations: ", pf_iter,"\n");
            V1 = Vm0.*exp.(im*θ0); opf_err = sum(abs.(V1-V0)); err_vec[opf_iter] = round(opf_err, RoundNearest, digits=4);
            print("D3F_OPF_rect: (",opf_iter,") Error = ",opf_err,"\n")
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
