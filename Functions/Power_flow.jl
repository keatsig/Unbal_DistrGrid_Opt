### FUNCTION TO CALCULATE POWER FLOW (FP)
function FP_pf(Vm0,θ0,QgY0,max_iter)

    ## Other inputs
    # start_point(4); Vm0 = s0[1:nb_nph]; θ0 = s0[nb_nph+1:2*nb_nph];
    # (Vm0,θ0,QgY0,stat) = TP_OPF(crit_node)

    ## Define active and reactive power of inverters
    Ygen = zeros(nb_nph1)
    if PV_en == 1
        Pg_idx = [2 6 10];  Qg_idx = [4 8 12];
        global PgY = zeros(ngy_inv); SgY = yinv_ratedS/baseMVA; num_count = 1;
        for i=1:ngy_inv
            j = Int(ygen1[i,1]); count1 = 1
            for k = 3:2:7
                if bus[j,k] != 0
                    PgY[num_count] = ygen1[i,Pg_idx[count1]]/baseMVA
                    if QgY0[num_count] > sqrt(SgY[num_count]^2 - PgY[num_count]^2)
                        QgY0[num_count] = sqrt(SgY[num_count]^2 - PgY[num_count]^2)
                    elseif QgY0[num_count] < -sqrt(SgY[num_count]^2 - PgY[num_count]^2)
                        QgY0[num_count] = -sqrt(SgY[num_count]^2 - PgY[num_count]^2)
                    end
                    num_count = num_count + 1;
                end
                count1 = count1 + 1;
            end
        end
        Ygen = Cbg_mod*(PgY + im*QgY0)
    end

    ## Generate load matrix
    ZIP_load()

    ## Start PF iterations
    V_new = zeros(nb_nph1); Vll0 = zeros(nb_nph); pf_err = 1; pf_iter = 0;
    V0 = Vm0[idx_ref].*exp.(im*θ0[idx_ref]); global w = -YLL_inv*YL0*V0;
    YloadP = zeros(nly_nph); YloadQ = zeros(nly_nph);
    DloadP = zeros(nld_nph); DloadQ = zeros(nld_nph);
    while pf_err >= 1e-4  && pf_iter<max_iter

        ## Initialize voltage
        V_old = Vm0[idx_noref].*exp.(im*θ0[idx_noref]);
        for i=1:nb_nph
            k = Int(Tl[i])
            Vll0[i] =  sqrt(Vm0[i]^2 + Vm0[k]^2 - 2*Vm0[i]*Vm0[k]*cos(θ0[i]-θ0[k]) )
        end
        V_cdY = diagm(0=>conj(V_old)); V_cdD = diagm(0=>H*conj(V_old));

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
        Yload = CblY_mod*(-YloadP - im*YloadQ);

        ## Calculate delta-connected loads
        for i=1:nld_nph
            j = Int(Tl[Int(Sd[i,7])]);
            dum1 = Vm0[Int(Sd[i,7])]*exp(im*θ0[Int(Sd[i,7])])-Vm0[j]*exp(im*θ0[j])
            dum2 = (Sd[i,3] + im*Sd[i,4])/(sqrt(3)*exp(im*Sd[i,8]))*dum1
            DloadP[i] = Sd[i,1]+real(dum2)+Sd[i,5]*Vll0[Int(Sd[i,7])]^2/3
            DloadQ[i] = Sd[i,2]+imag(dum2)+Sd[i,6]*Vll0[Int(Sd[i,7])]^2/3
        end
        Dload = CblD_FP*(-DloadP - im*DloadQ)

        ## Fixed-point linearization
        V_new = w+YLL_inv*(pinv(V_cdY)*conj(Yload+Ygen)+ H'*pinv(V_cdD)*conj(Dload))
        Vm_new= abs.(V_new); θ_new = angle.(V_new);
        Vm0 = [Vm_new[1:length(idx_noref1)]; Vm0[idx_ref]; Vm_new[length(idx_noref1)+1:end] ]
        θ0 = [θ_new[1:length(idx_noref1)]; θ0[idx_ref]; θ_new[length(idx_noref1)+1:end] ]

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
        #     start_point(4); Vm0 = s0[1:nb_nph]; θ0 = s0[nb_nph+1:2*nb_nph];
        # end

        ## Check error
        pf_err = sum(abs.(V_old-V_new))
        pf_iter = pf_iter + 1
        # print("PF- (",pf_iter,") Error = ",pf_err,"\n")
    end

    ## Calculate missing solutions
    Vph = Vm0.*exp.(im*θ0);
    Iph = (G+im*B)*Vph; S_inj = Vph.*conj(Iph);
    P_inj = real(S_inj); global P_loss = sum(P_inj)*baseMVA*1e3
    for i=1:nb_nph
        k = Int(Tl[i])
        Vll0[i] =  sqrt(Vm0[i]^2 + Vm0[k]^2 - 2*Vm0[i]*Vm0[k]*cos(θ0[i]-θ0[k]) )
    end
    print_VU(Vm0,θ0,Vll0);
    return Vm0,θ0,QgY0,pf_iter
end

#----------------------------------------------------------------------------------------#

### FUNCTION TO CALCULATE POWER FLOW (BFS)
function BFS_pf(Vm0,θ0,QgY0,max_iter)

    ## Other inputs
    # start_point(4); Vm0 = s0[1:nb_nph]; θ0 = s0[nb_nph+1:2*nb_nph];
    # (Vm0,θ0,QgY0,pf_iter) = FP_PF(Vm0,θ0,QgY0,50);

    ## Define active and reactive power of inverters
    Ygen = zeros(nb_nph)*(0+im*0)
    if PV_en == 1
        Pg_idx = [2 6 10];  Qg_idx = [4 8 12];
        global PgY = zeros(ngy_inv); SgY = yinv_ratedS/baseMVA; num_count = 1;
        for i=1:ngy_inv
            j = Int(ygen1[i,1]); count1 = 1
            for k = 3:2:7
                if bus[j,k] != 0
                    PgY[num_count] = ygen1[i,Pg_idx[count1]]/baseMVA
                    if QgY0[num_count] > sqrt(SgY[num_count]^2 - PgY[num_count]^2)
                        QgY0[num_count] = sqrt(SgY[num_count]^2 - PgY[num_count]^2)
                    elseif QgY0[num_count] < -sqrt(SgY[num_count]^2 - PgY[num_count]^2)
                        QgY0[num_count] = -sqrt(SgY[num_count]^2 - PgY[num_count]^2)
                    end
                    num_count = num_count + 1;
                end
                count1 = count1 + 1;
            end
        end
        Ygen = Cbg*(PgY + im*QgY0)
    end

    ## Generate load matrix
    ZIP_load()

    ## Start PF iterations
    Vll0 = zeros(nb_nph); pf_err = 1; pf_iter = 0; Ibr = zeros(nbr_nph)*(0+im*0);
    YloadP = zeros(nly_nph); YloadQ = zeros(nly_nph); Vm1 = zeros(nb_nph);
    DloadP = zeros(nld_nph); DloadQ = zeros(nld_nph); θ1 = zeros(nb_nph);
    I_ll = zeros(nld_nph)*(0+im*0); P_d2y = zeros(nldy_nph); Q_d2y = zeros(nldy_nph);
    while pf_err >= 1e-5 && pf_iter<max_iter

        ## Initialize voltage
        for i=1:nb_nph
            k = Int(Tl[i])
            Vll0[i] =  sqrt(Vm0[i]^2 + Vm0[k]^2 - 2*Vm0[i]*Vm0[k]*cos(θ0[i]-θ0[k]) )
        end

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
        Dload = CblD*(P_d2y + im*Q_d2y);

        ## Backward sweep
        V0 = Vm0.*exp.(im*θ0); S_inj = Dload+Yload-Ygen; count3 = 1;
        I_inj = conj(S_inj./V0); Ibus = zeros(nb_nph)*(0+im*0);
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
                Ibus[kk] = vr_tap*(Ibus[kk] + I_inj[kk]);
                count3 = count3+1;
            else
                Ibus[kk] = Ibus[kk] + I_inj[kk];
            end
            Ibus[jj] = Ibus[kk] + Ibus[jj];
        end

        ## Forward sweep
        Vm1[idx_ref] = Vm0[idx_ref]; θ1[idx_ref] = θ0[idx_ref]; count1 = 1;  count2=1
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
            jj = Int.(jjjj); kk = Int.(kkkk); Zbr = -pinv(Y[jj,kk]);
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
                Vm_kk = vr_tap*Vm1[jj].*exp.(im*θ1[jj]) - Zbr*Ibus[kk]
                count2 = count2+1;
            else
                Vm_kk = Vm1[jj].*exp.(im*θ1[jj]) - Zbr*Ibus[kk]
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
            Vm1[kk] = abs.(Vm_kk); θ1[kk] = angle.(Vm_kk);
            count1 = count1+length(kk);
        end

        ## Update solution
        V1 = Vm1.*exp.(im*θ1); pf_err = sum(abs.(V1-V0));
        Vm0 = Vm1; θ0 = θ1;
        for i=1:nb_nph
            k = Int(Tl[i])
            Vll0[i] =  sqrt(Vm0[i]^2 + Vm0[k]^2 - 2*Vm0[i]*Vm0[k]*cos(θ0[i]-θ0[k]) )
        end
        pf_iter = pf_iter + 1;
        # print("PF- (",pf_iter,") Error = ",pf_err,"\n")
    end

    ## Calculate missing solutions
    Vph = Vm0.*exp.(im*θ0);
    Iph = (G+im*B)*Vph; S_inj = Vph.*conj(Iph);
    P_inj = real(S_inj); global P_loss = sum(P_inj)*baseMVA*1e3
    for i=1:nb_nph
        k = Int(Tl[i])
        Vll0[i] =  sqrt(Vm0[i]^2 + Vm0[k]^2 - 2*Vm0[i]*Vm0[k]*cos(θ0[i]-θ0[k]) )
    end
    print_VU(Vm0,θ0,Vll0);
    return Vm0,θ0,QgY0,pf_iter
end

#----------------------------------------------------------------------------------------#

### FUNCTION TO CALCULATE POWER FLOW (CONVENTIONAL)
function NR_pf(Vm0,θ0,QgY0,max_iter)

    ## Initialize values
    # start_point(4); Vm0 = s0[1:nb_nph]; θ0 = s0[nb_nph+1:2*nb_nph];

    ## Define active and reactive power of inverters
    global PgY = zeros(ngy_inv);
    if PV_en == 1
        Pg_idx = [2 6 10];  Qg_idx = [4 8 12];
        global PgY = zeros(ngy_inv); SgY = yinv_ratedS/baseMVA; num_count = 1;
        for i=1:ngy_inv
            j = Int(ygen1[i,1]); count1 = 1
            for k = 3:2:7
                if bus[j,k] != 0
                    PgY[num_count] = ygen1[i,Pg_idx[count1]]/baseMVA
                    if QgY0[num_count] > sqrt(SgY[num_count]^2 - PgY[num_count]^2)
                        QgY0[num_count] = sqrt(SgY[num_count]^2 - PgY[num_count]^2)
                    elseif QgY0[num_count] < -sqrt(SgY[num_count]^2 - PgY[num_count]^2)
                        QgY0[num_count] = -sqrt(SgY[num_count]^2 - PgY[num_count]^2)
                    end
                    num_count = num_count + 1;
                end
                count1 = count1 + 1;
            end
        end
    end

    ## Generate load matrix
    ZIP_load()

    ## Start PF iterations
    pf_iter = 0; pf_err = 1;
    YloadP = zeros(nly_nph); YloadQ = zeros(nly_nph);
    DloadP = zeros(nld_nph); DloadQ = zeros(nld_nph);
    I_ll = zeros(nld_nph)*(0+im*0); P_d2y = zeros(nldy_nph); Q_d2y = zeros(nldy_nph);
    Pl0 = zeros(nb_nph1); Ql0 = zeros(nb_nph1); Vll0 = zeros(nb_nph)
    Jp_Vm = spzeros(nb_nph1,nb_nph1); Jq_Vm = spzeros(nb_nph1,nb_nph1);
    Jp_θ = spzeros(nb_nph1,nb_nph1); Jq_θ = spzeros(nb_nph1,nb_nph1);
    while pf_err>1e-4 && pf_iter<max_iter

        ## Initialize values
        pf_iter = pf_iter + 1; V0 = Vm0.*exp.(im*θ0);
        for i=1:nb_nph
            k = Int(Tl[i])
            Vll0[i] =  sqrt(Vm0[i]^2 + Vm0[k]^2 - 2*Vm0[i]*Vm0[k]*cos(θ0[i]-θ0[k]) )
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
        x0 = [Vm0[idx_noref]; θ0[idx_noref]];
        J_inv = inv(Matrix([Jp_Vm Jp_θ; Jq_Vm Jq_θ]));
        xP = Cbg_mod*PgY - CblY_mod*YloadP - CblD_mod*P_d2y - Pl0;
        xQ = Cbg_mod*QgY0 - CblY_mod*YloadQ - CblD_mod*Q_d2y - Ql0;
        x1 = x0 + J_inv*[xP; xQ];
        Vm_old = x1[1:nb_nph1]; θ_old = x1[1+nb_nph1:end];
        Vm0 = [Vm_old[1:length(idx_noref1)]; Vm0[idx_ref]; Vm_old[length(idx_noref1)+1:end] ]
        θ0 = [θ_old[1:length(idx_noref1)]; θ0[idx_ref]; θ_old[length(idx_noref1)+1:end] ]
        pf_err = abs(sum([xP; xQ]));
    end

    ## Calculate missing solutions
    for i=1:nb_nph
        k = Int(Tl[i])
        Vll0[i] =  sqrt(Vm0[i]^2 + Vm0[k]^2 - 2*Vm0[i]*Vm0[k]*cos(θ0[i]-θ0[k]) )
    end
    Vph = Vm0.*exp.(im*θ0);
    Iph = (G+im*B)*Vph; S_inj = Vph.*conj(Iph);
    P_inj = real(S_inj); global P_loss = sum(P_inj)*baseMVA*1e3
    print_VU(Vm0,θ0,Vll0);
    return Vm0,θ0,QgY0,pf_iter
end

#----------------------------------------------------------------------------------------#

### FUNCTION TO CALCULATE POWER FLOW (DistFlow)
function D3F_pf(Vm0,θ0,QgY0,max_iter)

    ## Other inputs
    start_point(1); Vm0 = s0[1:nb_nph]; θ0 = s0[nb_nph+1:2*nb_nph];
    # (Vm0,θ0,QgY0,stat) = TP_OPF(crit_node)

    ## Define active and reactive power of inverters
    Ygen = zeros(nb_nph)*(0+im*0)
    if ZIP_load == 1
        Pg_idx = [2 6 10];  Qg_idx = [4 8 12];
        global PgY = zeros(ngy_inv); SgY = yinv_ratedS/baseMVA; num_count = 1;
        for i=1:ngy_inv
            j = Int(ygen1[i,1]); count1 = 1
            for k = 3:2:7
                if bus[j,k] != 0
                    PgY[num_count] = ygen1[i,Pg_idx[count1]]/baseMVA
                    if QgY0[num_count] > sqrt(SgY[num_count]^2 - PgY[num_count]^2)
                        QgY0[num_count] = sqrt(SgY[num_count]^2 - PgY[num_count]^2)
                    elseif QgY0[num_count] < -sqrt(SgY[num_count]^2 - PgY[num_count]^2)
                        QgY0[num_count] = -sqrt(SgY[num_count]^2 - PgY[num_count]^2)
                    end
                    num_count = num_count + 1;
                end
                count1 = count1 + 1;
            end
        end
        Ygen = Cbg*(PgY + im*QgY0)
    end

    ## Generate load matrix
    ZIP_load()

    ## Start PF iterations
    Vll0 = zeros(nb_nph); pf_err = 1; pf_iter = 0;
    YloadP = zeros(nly_nph); YloadQ = zeros(nly_nph); Vm1 = zeros(nb_nph);
    DloadP = zeros(nld_nph); DloadQ = zeros(nld_nph); θ1 = zeros(nb_nph);
    I_ll = zeros(nld_nph)*(0+im*0); P_d2y = zeros(nldy_nph); Q_d2y = zeros(nldy_nph);
    while pf_err >= 1e-6 && pf_iter<max_iter

        ## Initialize voltage
        V0 = Vm0.*exp.(im*θ0); y0 = Vm0.^2; γ=V0./conj(V0');
        for i=1:nb_nph
            k = Int(Tl[i])
            Vll0[i] =  sqrt(Vm0[i]^2 + Vm0[k]^2 - 2*Vm0[i]*Vm0[k]*cos(θ0[i]-θ0[k]) )
        end

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
        Dload = CblD*(P_d2y + im*Q_d2y);

        ## DistFlow
        V0 = Vm0.*exp.(im*θ0); s_inj = Dload+Yload-Ygen; i_inj = conj(s_inj./V0);
        count3 = 1; Sbus = zeros(nb_nph)*(0+im*0); Ibus = zeros(nb_nph)*(0+im*0);
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
                Ibus[kk] = vr_tap*(Ibus[kk] + i_inj[kk]);
                count3 = count3+1;
            else
                Ibus[kk] = Ibus[kk] + i_inj[kk];
            end
            Lbr = diag(Zbr*Ibus[kk]*Ibus[kk]'); Ibus[jj] = Ibus[kk] + Ibus[jj];
            Sbus[kk] = Sbus[kk] + s_inj[kk]; Sbus[jj] = Sbus[kk] + Sbus[jj] + Lbr;
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
            Zbr = -pinv(Y[jj,kk]); Hbr = real(diag(Zbr*Ibus[kk]*Ibus[kk]'*Zbr'));
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
                y0[kk] = vr_tap.^2*y0[jj]-Hbr-2*real(conj(Zbr).*γbr*Sbus[kk]);
                count2 = count2+1;
            else
                y0[kk] = y0[jj]-Hbr-2*real(conj(Zbr).*γbr*Sbus[kk]);
            end
            θ0[kk] = θ0[jj] + angle.(y0[kk]+Zbr'*Sbus[kk])
        end

        ## Update solution
        Vm0 = sqrt.(y0); V1 = Vm0.*exp.(im*θ0);
        pf_err = sum(abs.(V1-V0)); V0 = V1;
        for i=1:nb_nph
            k = Int(Tl[i])
            Vll0[i] =  sqrt(Vm0[i]^2 + Vm0[k]^2 - 2*Vm0[i]*Vm0[k]*cos(θ0[i]-θ0[k]) )
        end
        pf_iter = pf_iter + 1;
        # print("PF- (",pf_iter,") Error = ",pf_err,"\n")
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
    return Vm0,θ0,QgY0,pf_iter
end

#----------------------------------------------------------------------------------------#
