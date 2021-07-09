### FUNCTION TO MAP GENERATORS AND LOADS TO BUS
function ZIP_load()
    ## Calculate overall star-connected loads in system
    global Sy=spzeros(nly_nph,8); num_count = 1; Py_idx = [2 3 4];
    for i=1:nly
        j = Int(yload[i,1]); count1 = 1;  count2 = 1
        for k = 3:2:7
            if bus[j,k] != 0
                Sy[num_count,1:6] = yload[i,Py_idx[count1]:3:Py_idx[count1]+15]/baseMVA
                Sy[num_count,7] = sum(bus_count[1:j-1]) + count2
                if bus_ϕ[Int(Sy[num_count,7])] == "B"
                    Sy[num_count,8] = -2*pi/3;
                elseif bus_ϕ[Int(Sy[num_count,7])] == "C"
                    Sy[num_count,8] = 2*pi/3;
                end
                num_count = num_count + 1; count2 = count2 + 1;
            end
        count1 = count1 + 1;
        end
    end

    ## Calculate overall delta-connected loads in system
    global Sd = spzeros(nld_nph,8); num_count = 1; Pd_idx = [2 3 4]; Pd_idx1 = 2:3:17;
    for i=1:nld
        j = Int(dload[i,1]); count1 = 1;  count2 = 1
        for k = 1:nld_count[j]
            Sd[num_count,1:6] = dload[i,Pd_idx[count1]:3:Pd_idx[count1]+15]/baseMVA
            if nld_count[j] == 1
                for kk=1:6
                    if !isempty(findall(idx->idx!=0,dload[i,Pd_idx1[kk]:Pd_idx1[kk]+2]))
                        dum1 = findall(idx->idx!=0,dload[i,Pd_idx1[kk]:Pd_idx1[kk]+2])[1]
                        Sd[num_count,kk] = dload[i,dum1+Pd_idx1[kk]-1]/baseMVA
                    end
                end
            end
            Sd[num_count,7] = sum(bus_count[1:j-1]) + count2
            if bus_ϕ[Int(Sd[num_count,7])] == "A"
                Sd[num_count,8] = pi/6;
            elseif bus_ϕ[Int(Sd[num_count,7])] == "B"
                Sd[num_count,8] = -pi/2;
            elseif bus_ϕ[Int(Sd[num_count,7])] == "C"
                Sd[num_count,8] = 5*pi/6;
            end
            num_count = num_count + 1; count2 = count2 + 1;
            count1 = count1 + 1;
        end
    end
end

#----------------------------------------------------------------------------------------#

### FUNCTION TO MAP GENERATORS AND LOADS TO BUS
function conn_bgl()

    ## Calculate A,B matrices
    z_file = open(string("data/",case_name,"/z_dump.xml")); z_dump = readlines(z_file)
    idx_A_f = findall(idx -> idx == "\t\t<A_matrix>", z_dump); idx_B_f = findall(idx -> idx == "\t\t<B_matrix>", z_dump)
    idx_c_f = findall(idx -> idx == "\t\t<c_matrix>", z_dump); idx_d_f = findall(idx -> idx == "\t\t<d_matrix>", z_dump)
    A_mat = zeros(3,3,nbr); B_mat = zeros(3,3,nbr)*(0+0*im);
    for i=1:nbr
        j = Int(branch[i,2]);
        if tf_branch[i,3] == 8
            A_mat[:,:,i] = Diagonal(1 .+ 0.00625*tf_branch[i,6:8])
            B_mat[:,:,i] = Matrix{Float64}(I, 3, 3)*0.00001/baseZ[j]
        end
        if tf_branch[i,3] == 1
            A_mat[:,:,i] = Matrix{Float64}(I, 3, 3)
            B_mat[:,:,i] = Matrix{Float64}(I, 3, 3)*(tf_branch[i,4]+im*tf_branch[i,5])
        end
        if tf_branch[i,3] == 0
            A_mat[:,:,i] = Matrix{Float64}(I, 3, 3)
            B_mat[:,:,i] = branch_matrix(z_dump,idx_B_f[i])/baseZ[j]
        end
    end

    ## Map generators to bus
    global ygen1= ygen[2:end,:];
    global Cbg = zeros(nb_nph,ngy_inv); count_c = 1
    for i=1:ngy_inv
        count_r = 1;
        for j=1:nb
            if bus[j,1] == ygen1[i,1]
                if bus_count[Int(bus[j,1])] == 3
                    Cbg[count_r:count_r+2,count_c:count_c+2] = Matrix{Float64}(I, 3, 3)
                    count_r = count_r+3; count_c = count_c+3;
                elseif bus_count[Int(bus[j,1])] == 2
                    Cbg[count_r:count_r+1,count_c:count_c+1] = Matrix{Float64}(I, 2, 2)
                    count_r = count_r+2; count_c = count_c+2
                elseif bus_count[Int(bus[j,1])] == 1
                    Cbg[count_r,count_c] = 1
                    count_r = count_r+1; count_c = count_c+1
                end
            else
                count_r = count_r+Int(bus_count[Int(bus[j,1])])
            end
        end
    end
    global Cbg_SS = zeros(nb_nph,3); count_c = 1
    for i=1:1
        count_r = 1;
        for j=1:nb
            if bus[j,1] == ygen[i,1]
                if bus_count[Int(bus[j,1])] == 3
                    Cbg_SS[count_r:count_r+2,count_c:count_c+2] = Matrix{Float64}(I, 3, 3)
                    count_r = count_r+3; count_c = count_c+3;
                elseif bus_count[Int(bus[j,1])] == 2
                    Cbg_SS[count_r:count_r+1,count_c:count_c+1] = Matrix{Float64}(I, 2, 2)
                    count_r = count_r+2; count_c = count_c+2
                elseif bus_count[Int(bus[j,1])] == 1
                    Cbg_SS[count_r,count_c] = 1
                    count_r = count_r+1; count_c = count_c+1
                end
            else
                count_r = count_r+Int(bus_count[Int(bus[j,1])])
            end
        end
    end
    global Cbg_mod = Cbg[idx_noref,:]

    ## Map wye-connected loads to bus
    global CblY = spzeros(nb_nph,nly_nph);
    if nly != 0
        count_c = 1;
        for i=1:nly
            count_r = 1;
            for j=1:nb
                if bus[j,1] == yload[i,1]
                    if bus_count[Int(bus[j,1])] == 3
                        CblY[count_r:count_r+2,count_c:count_c+2] = Matrix{Float64}(I, 3, 3)
                        count_r = count_r+3; count_c = count_c+3;
                    elseif bus_count[Int(bus[j,1])] == 2
                        CblY[count_r:count_r+1,count_c:count_c+1] = Matrix{Float64}(I, 2, 2)
                        count_r = count_r+2; count_c = count_c+2
                    elseif bus_count[Int(bus[j,1])] == 1
                        CblY[count_r,count_c] = 1
                        count_r = count_r+1; count_c = count_c+1
                    end
                else
                    count_r = count_r+Int(bus_count[Int(bus[j,1])])
                end
            end
        end
    end
    global CblY_mod = CblY[idx_noref,:]

    ## Map delta-connected loads to bus
    global CblD = spzeros(nb_nph,nldy_nph); global CblD_FP = spzeros(nbD_nph-3,nld_nph);
    global H = zeros(nbD_nph-3,nb_nph-3); global Tl2p = zeros(nldy_nph,nld_nph);
    if nld != 0
        count_c = 1;
        for i=1:nld
            count_r = 1;
            for j=1:nb
                if bus[j,1] == dload[i,1]
                    if bus_count[Int(bus[j,1])] == 3
                        CblD[count_r:count_r+2,count_c:count_c+2] = Matrix{Float64}(I, 3, 3)
                        count_r = count_r+3; count_c = count_c+3;
                    elseif bus_count[Int(bus[j,1])] == 2
                        CblD[count_r:count_r+1,count_c:count_c+1] = Matrix{Float64}(I, 2, 2)
                        count_r = count_r+2; count_c = count_c+2;
                    end
                else
                    count_r = count_r+Int(bus_count[Int(bus[j,1])])
                end
            end
        end

        count_c = 1;
        for i=1:nld
            count_r = 1;
            for j=1:nb
                if bus[j,1] == dload[i,1] && bus[j,2] !=3
                    if bus_count[Int(bus[j,1])] == 3
                       CblD_FP[count_r:count_r+2,count_c:count_c+2] = Matrix{Float64}(I, 3, 3)
                       count_r = count_r+3; count_c = count_c+3;
                    elseif bus_count[Int(bus[j,1])] == 2
                       CblD_FP[count_r,count_c] = 1
                       count_r = count_r+1; count_c = count_c+1;
                    end
                elseif bus[j,1] != dload[i,1] && bus[j,2] !=3
                    if bus_count[Int(bus[j,1])] == 3
                        count_r = count_r+3
                    elseif bus_count[Int(bus[j,1])] == 2
                        count_r = count_r+1
                    end
                end
            end
        end

        T1 = [1 -1 0; 0 1 -1; -1 0 1]
        count_r = 1; count_c = 1;
        for i=1:nb
            if bus_count[i] == 3 && bus[i,2] != 3
                H[count_r:count_r+2,count_c:count_c+2] = T1
                count_r = count_r+3; count_c = count_c+3;
            elseif bus_count[i] == 2 && (bus_ph[i] == "ABN" || bus_ph[i] == "ABD" || bus_ph[i] == "BAN" || bus_ph[i] == "BAD" || bus_ph[i] == "BA" || bus_ph[i] == "AB")
                H[count_r,count_c:count_c+1] = T1[1,1:2]
                count_r = count_r+1; count_c = count_c+2;
            elseif bus_count[i] == 2 && (bus_ph[i] == "BCN" || bus_ph[i] == "BCD" || bus_ph[i] == "CBN"  || bus_ph[i] == "CBD" || bus_ph[i] == "CB" || bus_ph[i] == "BC")
                H[count_r,count_c:count_c+1] = T1[2,2:3]
                count_r = count_r+1; count_c = count_c+2;
            elseif bus_count[i] == 2 && (bus_ph[i] == "ACN" || bus_ph[i] == "ACD" || bus_ph[i] == "CAN" || bus_ph[i] == "CAD" || bus_ph[i] == "CA" || bus_ph[i] == "AC")
                H[count_r,count_c:count_c+1] = T1[3,1:2:3]
                count_r = count_r+1; count_c = count_c+2;
            elseif bus_count[i] == 1
                count_c = count_c+1;
            end
        end

        T1 = [1 0 -1; -1 1 0; 0 -1 1]
        count_r = 1; count_c = 1
        for i=1:nld
            j = Int(dload[i,1]); count1 = 1;  count2 = 1
            k_ph = bus_ϕ[Int.(sum(bus_count[1:j-1])+1:sum(bus_count[1:j]))]
            if bus_count[j] == 3
                Tl2p[count_r:count_r+2,count_c:count_c+2] = T1
                count_r = count_r+3; count_c = count_c+3
            elseif bus_count[j] == 2 && (k_ph == ["A"; "B"] || k_ph == ["B"; "C"])
                Tl2p[count_r:count_r+1,count_c] = [1; -1]
                count_r = count_r+2; count_c = count_c+1;
            elseif bus_count[j] == 2 && k_ph == ["A"; "C"]
                Tl2p[count_r:count_r+1,count_c] = [-1; 1]
                count_r = count_r+2; count_c = count_c+1;
            end
        end
    end
    global CblD_mod = CblD[idx_noref,:]

    # global CbgD = spzeros(nb_nph,ngdy_nph);
    # if ngd != 0
    #     count_c = 1;
    #     for i=1:ngd
    #         count_r = 1;
    #         for j=1:nb
    #             if bus[j,1] == dgen[i,1]
    #                 if bus_count[Int(bus[j,1])] == 3
    #                     CbgD[count_r:count_r+2,count_c:count_c+2] = Matrix{Float64}(I, 3, 3)
    #                     count_r = count_r+3; count_c = count_c+3;
    #                 elseif bus_count[Int(bus[j,1])] == 2
    #                     CbgD[count_r:count_r+1,count_c:count_c+1] = Matrix{Float64}(I, 2, 2)
    #                     count_r = count_r+2; count_c = count_c+2;
    #                 end
    #             else
    #                 count_r = count_r+Int(bus_count[Int(bus[j,1])])
    #             end
    #         end
    #     end
    # end

    ## Map phase-phase connections to bus
    global Tl = spzeros(nb_nph); count = 1
    for i=1:nb
        if bus_count[Int(bus[i,1])] == 3
            Tl[count:count+2] = [count+1 count+2 count]
            count = count+3
        elseif bus_count[Int(bus[i,1])] == 2
            if bus_ϕ[count] == "A" && bus_ϕ[count+1] == "C"
                Tl[count:count+1] = [count count]
                count = count+2
            else
                Tl[count:count+1] = [count+1 count+1]
                count = count+2
            end
        elseif bus_count[Int(bus[i,1])] == 1
            Tl[count] = count
            count = count+1
        end
    end

    # Calculate symmetrical components
    global ap = zeros(nb_nph); global an = zeros(nb_nph)
    for i=1:nb_nph
        if bus_ϕ[i] == "B"
            ap[i] = 2*pi/3; an[i] = -2*pi/3
        end
        if bus_ϕ[i] == "C"
            ap[i] = -2*pi/3; an[i] = 2*pi/3
        end
    end
    global ad = cos(2*pi/3); global aq = sin(2*pi/3)
    a = cos(2*pi/3) + im*sin(2*pi/3);
    a_Vn_d = [1 real(a^2) real(a); 0 imag(a^2) imag(a)]; a_Vn_q = [0 imag(a) imag(a^2); 1 real(a) real(a^2)];
    a_Vp_d = [1 real(a) real(a^2); 0 imag(a) imag(a^2)]; a_Vp_q = [0 imag(a^2) imag(a); 1 real(a^2) real(a)];
    global a_rec = [a_Vn_d a_Vn_q; a_Vp_d a_Vp_q];
    global a_Vn_d1 = spzeros(2*n_VU,3*n_VU); global a_Vn_q1 = spzeros(2*n_VU,3*n_VU);
    global a_Vp_d1 = spzeros(2*n_VU,3*n_VU); global a_Vp_q1 = spzeros(2*n_VU,3*n_VU);
    for i=1:n_VU
        a_Vn_d1[i*2-1:2*i,i*3-2:3*i] = a_Vn_d; a_Vn_q1[i*2-1:2*i,i*3-2:3*i] = a_Vn_q;
        a_Vp_d1[i*2-1:2*i,i*3-2:3*i] = a_Vp_d; a_Vp_q1[i*2-1:2*i,i*3-2:3*i] = a_Vp_q;
    end
end

#----------------------------------------------------------------------------------------#

# FUNCTION TO CALCULATE INITIAL POINT
function start_point(init_op)

    if init_op == 1  # Flat start
        Vm0 = ones(nb_nph); Vll0 = sqrt(3)*ones(nb_nph); θ0 = zeros(nb_nph);
        for i=1:nb_nph
            if bus_ϕ[i] == "B"
                θ0[i] = -2*pi/3;
            end
            if bus_ϕ[i] == "C"
                θ0[i] = 2*pi/3;
            end
        end
        global s0 = [Vm0; θ0; Vll0; zeros(ngy_nph); zeros(ngy_nph)]

    elseif init_op == 2 # GridLabD-start
        Vm0 = zeros(nb_nph); θ0 = zeros(nb_nph); Vll0 = zeros(nb_nph);count=1
        for i=1:nb
            for j = 3:2:7
                if bus[i,j] != 0
                    Vm0[count] = bus[i,j+6]
                    θ0[count] = bus[i,j+7]*pi/180
                    count = count + 1
                end
            end
        end
        for i=1:nb_nph
            k = Int(Tl[i])
            Vll0[i] = sqrt.((Vm0[i].*cos.(θ0[i]) - Vm0[k].*cos.(θ0[k])).^2 +  (Vm0[i].*sin.(θ0[i]) - Vm0[k].*sin.(θ0[k])).^2)
        end
        Pg0 = zeros(ngy_nph); Qg0 = zeros(ngy_nph); count = 1
        for i=1:ngy
            j = Int(ygen[i,1])
            for k = 3:2:7
                if bus[j,k] != 0
                    Pg0[count] = ygen[i,k+11]
                    Qg0[count] = ygen[i,k+12]
                end
            end
        end
        global s0 = [Vm0; θ0; Vll0; Pg0; Qg0]  # [Vm θ Vll Pg Qg]

    elseif init_op == 3 # Previous soln-start
        Vm_init = readdlm(string("Vm_",case_name,"_soln.csv"), ',')
        Vm0 = Vm_init[:,1]; θ0 = Vm_init[:,2]
        Vll0 = zeros(nb_nph)
        for i=1:nb_nph
            k = Int(Tl[i])
            Vll0[i] = sqrt.((Vm0[i].*cos.(θ0[i]) - Vm0[k].*cos.(θ0[k])).^2 +  (Vm0[i].*sin.(θ0[i]) - Vm0[k].*sin.(θ0[k])).^2)
        end
        Pg0 = zeros(ngy_nph); Qg0 = readdlm(string("QgY0_",case_name,"_soln1.csv"), ',');
        global s0 = [Vm0; θ0; Vll0; Pg0; Qg0]

    elseif init_op == 4  # No-load start
            Vm_SS = ones(3); θ_SS = [0; -2*pi/3; 2*pi/3];
            V0 = Vm_SS.*exp.(im*θ_SS); global w = -YLL_inv*YL0*V0;
            Vm0 = zeros(nb_nph); θ0 = zeros(nb_nph); Vm0[idx_ref] = Vm_SS; θ0[idx_ref] = θ_SS;
            Vm0[idx_noref] = abs.(w); θ0[idx_noref] = angle.(w);
            Vll0 = zeros(nb_nph)
            for i=1:nb_nph
                k = Int(Tl[i])
                Vll0[i] =  sqrt(Vm0[i]^2 + Vm0[k]^2 - 2*Vm0[i]*Vm0[k]*cos(θ0[i]-θ0[k]) )
            end
            global s0 = [Vm0; θ0; Vll0; zeros(ngy_nph); zeros(ngy_nph)]
    end

end

#----------------------------------------------------------------------------------------#

# FUNCTION TO PRINT VUF DATA
function print_VU(V1,θ1,Vl)
    global soln_VUF = zeros(nb,3); Vpavg = zeros(nb); Vpdev = zeros(nb) ; Vlavg = zeros(nb); Vldev = zeros(nb)
    for i=1:nb
        Vn1 = 0; Vp1 = 0; V01 = 0
        if bus_count[i] == 3
            j=Int.(sum(bus_count[1:i-1])+1:sum(bus_count[1:i]))
            Vpavg[i] = sum(V1[j])/3
            Vpdev[i] = max(abs(V1[j[1]]-Vpavg[i]), abs(V1[j[2]]-Vpavg[i]), abs(V1[j[3]]-Vpavg[i]) )
            Vlavg[i] = sum(Vl[j])/3
            Vldev[i] = max(abs(Vl[j[1]]-Vlavg[i]), abs(Vl[j[2]]-Vlavg[i]), abs(Vl[j[3]]-Vlavg[i]) )
            Vn1 = (sqrt(sum(V1[j].*cos.(θ1[j]+an[j]) )^2 + sum(V1[j].*sin.(θ1[j]+an[j]) )^2 )/3)
            Vp1 = (sqrt(sum(V1[j].*cos.(θ1[j]+ap[j]) )^2 + sum(V1[j].*sin.(θ1[j]+ap[j]) )^2 )/3)
            V01 = (sqrt(sum(V1[j].*cos.(θ1[j]) )^2 + sum(V1[j].*sin.(θ1[j]) )^2 )/3)
            soln_VUF[i,1] = round(Vn1/Vp1*100, RoundNearest, digits=4)
            # soln_VUF[i,2] = round(V01/Vp1*100, RoundNearest, digits=4)
            soln_VUF[i,2] = round(Vldev[i]/Vlavg[i]*100, RoundNearest, digits=4)
            soln_VUF[i,3] = round(Vpdev[i]/Vpavg[i]*100, RoundNearest, digits=4)
        end
    end
end

#----------------------------------------------------------------------------------------#
