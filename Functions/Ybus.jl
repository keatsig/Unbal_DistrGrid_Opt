### FUNCTION TO MAKE YBUS MATRIX FROM BRANCH DATA

## Function to create ybus matrix for each transformer line
function ybus_submatrix(G1,B1,G2,B2,Y1,Y2,Y3,Y4,tap,ph_op)
    tap_pri = tap[1]; tap_sec = tap[2]
    if ph_op == 2 # AN
        Y1 = [Y1[1,1]]; Y2 = [Y2[1,1]]; Y3 = [Y3[1,1]]; Y4 = [Y4[1,1]];
    elseif ph_op == 3 # BN
        Y1 = [Y1[2,2]]; Y2 = [Y2[2,2]]; Y3 = [Y3[2,2]]; Y4 = [Y4[2,2]];
    elseif ph_op == 4 # CN
        Y1 = [Y1[3,3]]; Y2 = [Y2[3,3]]; Y3 = [Y3[3,3]]; Y4 = [Y4[3,3]];
    elseif ph_op == 5 # ABN
        Y1 = Y1[1:2,1:2]; Y2 = Y2[1:2,1:2]; Y3 = Y3[1:2,1:2]; Y4 = Y4[1:2,1:2];
    elseif ph_op == 6 # BCN
        Y1 = Y1[2:3,2:3]; Y2 = Y2[2:3,2:3]; Y3 = Y3[2:3,2:3]; Y4 = Y4[2:3,2:3];
    elseif ph_op == 7 # ACN
        Y1 = Y1[1:2:3,1:2:3]; Y2 = Y2[1:2:3,1:2:3]; Y3 = Y3[1:2:3,1:2:3]; Y4 = Y4[1:2:3,1:2:3];
    end
    G1_sub = spzeros(size(G1,1),size(G1,1)); B1_sub = spzeros(size(B1,1),size(B1,1))
    G1_sub = G1+real(Y1)/tap_pri^2;  B1_sub = B1+imag(Y1)/tap_pri^2;
    G2_sub = spzeros(size(G2,1),size(G2,1)); B2_sub = spzeros(size(B2,1),size(B2,1))
    G2_sub = G2+real(Y2)/tap_sec^2;  B2_sub = B2+imag(Y2)/tap_sec^2;
    G3_sub = spzeros(size(G1,1),size(G2,1)); B3_sub = spzeros(size(B1,1),size(B2,1))
    G4_sub = spzeros(size(G2,1),size(G1,1)); B4_sub = spzeros(size(B2,1),size(B1,1))
    G3_sub = real(Y3)/(tap_pri*tap_sec);  B3_sub = imag(Y3)/(tap_pri*tap_sec);
    G4_sub = real(Y4)/(tap_pri*tap_sec);  B4_sub = imag(Y4)/(tap_pri*tap_sec);
    return G1_sub,B1_sub,G2_sub,B2_sub,G3_sub,B3_sub,G4_sub,B4_sub
end

## Function to create overall system ybus matrix
function Ybus()
    global G = zeros(nb_nph,nb_nph); global B = zeros(nb_nph,nb_nph); global Bsh = zeros(nb_nph,nb_nph);
    for i = 1:nbr
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
        z = zeros(nz,nz)*(0+0*im); b = zeros(nz,nz)*(0+0*im)
        if ntf == 0 || tf_branch[i,3] == 0
            z11 = branch[i,3] + branch[i,4]im
            z12 = branch[i,5] + branch[i,6]im
            z13 = branch[i,7] + branch[i,8]im
            z22 = branch[i,9] + branch[i,10]im
            z23 = branch[i,11] + branch[i,12]im
            z33 = branch[i,13] + branch[i,14]im
            if abs(z11) != 0 && abs(z22) != 0 && abs(z33) != 0
                z = [z11 z12 z13; z12 z22 z23; z13 z23 z33]
                b = [branch[i,15] 0 0; 0 branch[i,16] 0; 0 0 branch[i,17]]/2
            elseif abs(z33) == 0 && abs(z11) != 0 && abs(z22) != 0
                z = [z11 z12; z12 z22]
                b = [branch[i,15] 0; 0 branch[i,16]]/2
            elseif abs(z22) == 0  && abs(z11) != 0 && abs(z33) != 0
                z = [z11 z13; z13 z33]
                b = [branch[i,15] 0; 0 branch[i,17]]/2
            elseif abs(z11) == 0 && abs(z33) != 0 && abs(z22) != 0
                z = [z22 z23; z23 z33]
                b = [branch[i,16] 0; 0 branch[i,17]]/2
            elseif abs(z11) == 0 && abs(z22) == 0 && abs(z33) != 0
                z = [z33]
                b = [branch[i,17]]/2
            elseif abs(z11) == 0 && abs(z33) == 0 && abs(z22) != 0
                z = [z22]
                b = [branch[i,16]]/2
            elseif abs(z33) == 0 && abs(z22) == 0 && abs(z11) != 0
                z = [z11]
                b = [branch[i,15]]/2
            end
            y = pinv(z);
            G[jj,jj] = G[jj,jj] + real(y)
            B[jj,jj] = B[jj,jj] + imag(y) + b
            G[kk,kk] = G[kk,kk] + real(y)
            B[kk,kk] = B[kk,kk] + imag(y) + b
            G[jj,kk] = -real(y); B[jj,kk] = -imag(y)
            G[kk,jj] = -real(y); B[kk,jj] = -imag(y)
            Bsh[jj,kk] = b; Bsh[kk,jj] = b
        else
            if ntf !=0
            op = tf_branch[i,3]
            zt = (tf_branch[i,4] + tf_branch[i,5]im)
            yt = 1/zt
            Y1 = yt*Matrix{Float64}(I, 3, 3)
            Y2 = yt*[2 -1 -1; -1 2 -1; -1 -1 2]/3
            Y3 = yt*[-1 1 0; 0 -1 1; 1 0 -1]/sqrt(3)
            if op == 1  # YNyn0
                (G[jj,jj],B[jj,jj],G[kk,kk],B[kk,kk],G[jj,kk],B[jj,kk],G[kk,jj],B[kk,jj]) = ybus_submatrix(G[jj,jj],B[jj,jj],G[kk,kk],B[kk,kk],Y1,Y1,-Y1,-Y1,tf_branch[i,6:7],tf_branch[i,8])

            elseif op == 2  # Yy0
                (G[jj,jj],B[jj,jj],G[kk,kk],B[kk,kk],G[jj,kk],B[jj,kk],G[kk,jj],B[kk,jj]) = ybus_submatrix(G[jj,jj],B[jj,jj],G[kk,kk],B[kk,kk],Y2,Y2,-Y2,-Y2,tf_branch[i,6:7],tf_branch[i,8])

            elseif op == 3  # YNd1
                (G[jj,jj],B[jj,jj],G[kk,kk],B[kk,kk],G[jj,kk],B[jj,kk],G[kk,jj],B[kk,jj]) = ybus_submatrix(G[jj,jj],B[jj,jj],G[kk,kk],B[kk,kk],Y1,Y2,Y3,conj(Y3'),tf_branch[i,6:7],tf_branch[i,8])

            elseif op == 4  # Yd1
                (G[jj,jj],B[jj,jj],G[kk,kk],B[kk,kk],G[jj,kk],B[jj,kk],G[kk,jj],B[kk,jj]) = ybus_submatrix(G[jj,jj],B[jj,jj],G[kk,kk],B[kk,kk],Y2,Y2,Y3,conj(Y3'),tf_branch[i,6:7],tf_branch[i,8])

            elseif op == 5  # Dyn1
                (G[jj,jj],B[jj,jj],G[kk,kk],B[kk,kk],G[jj,kk],B[jj,kk],G[kk,jj],B[kk,jj]) = ybus_submatrix(G[jj,jj],B[jj,jj],G[kk,kk],B[kk,kk],Y2,Y1,Y3,conj(Y3'),tf_branch[i,6:7],tf_branch[i,8])

            elseif op == 6  # Dyn11
                (G[jj,jj],B[jj,jj],G[kk,kk],B[kk,kk],G[jj,kk],B[jj,kk],G[kk,jj],B[kk,jj]) = ybus_submatrix(G[jj,jj],B[jj,jj],G[kk,kk],B[kk,kk],Y2,Y1,conj(Y3'),Y3,tf_branch[i,6:7],tf_branch[i,8])
            elseif op == 7  # Dd0
                (G[jj,jj],B[jj,jj],G[kk,kk],B[kk,kk],G[jj,kk],B[jj,kk],G[kk,jj],B[kk,jj]) = ybus_submatrix(G[jj,jj],B[jj,jj],G[kk,kk],B[kk,kk],Y2,Y2,-Y2,-Y2,tf_branch[i,6:7],tf_branch[i,8])
            elseif op == 8  # Regulator
                yt = 1/(tf_branch[i,4] + tf_branch[i,5]im);
                # yt = 1e5*(1+0*im)/(findmax(tap)[1])^2
                if bus_ph[k] == "AN"
                    tap = tf_branch[i,6]*0.00625 .+ 1;
                    Y1 = yt*Matrix{Float64}(I, 1, 1);
                elseif bus_ph[k] == "BN"
                    tap = tf_branch[i,7]*0.00625 .+ 1;
                    Y1 = yt*Matrix{Float64}(I, 1, 1);
                elseif bus_ph[k] == "CN"
                    tap = tf_branch[i,8]*0.00625 .+ 1;
                    Y1 = yt*Matrix{Float64}(I, 1, 1);
                elseif bus_ph[k] == "ABN" || bus_ph[k] == "BAN"
                    tap = [tf_branch[i,6], tf_branch[i,7]]*0.00625 .+ 1; tap = Diagonal(tap);
                    Y1 = yt*Matrix{Float64}(I, 2, 2);
                elseif bus_ph[k] == "BCN" || bus_ph[k] == "CBN"
                    tap = [tf_branch[i,7], tf_branch[i,8]]*0.00625 .+ 1; tap = Diagonal(tap);
                    Y1 = yt*Matrix{Float64}(I, 2, 2);
                elseif bus_ph[k] == "ACN" || bus_ph[k] == "CAN"
                    tap = [tf_branch[i,6], tf_branch[i,8]]*0.00625 .+ 1; tap = Diagonal(tap);
                    Y1 = yt*Matrix{Float64}(I, 2, 2);
                else
                    tap = [tf_branch[i,6], tf_branch[i,7], tf_branch[i,8]]*0.00625 .+ 1; tap = Diagonal(tap);
                    Y1 = yt*Matrix{Float64}(I, 3, 3);
                end
                G[jj,jj] = G[jj,jj]+real(Y1).*tap^2;
                B[jj,jj] = B[jj,jj]+imag(Y1).*tap^2;
                G[kk,kk] = G[kk,kk]+real(Y1);
                B[kk,kk] = B[kk,kk]+imag(Y1);
                G[jj,kk] = real(-Y1).*tap;
                B[jj,kk] = imag(-Y1).*tap;
                G[kk,jj] = real(-Y1).*tap;
                B[kk,jj] = imag(-Y1).*tap;
            end
            end
        end
    end

    ## Shunt capacitors
    global Sh_cap = zeros(nb_nph); num_count = 1
    for i=1:nbr
        if tf_branch[i,3] == 0
            j = Int(branch[i,1]); jj = Int(sum(bus_count[1:j-1])+1):Int(sum(bus_count[1:j]))
            k = Int(branch[i,2]); kk = Int(sum(bus_count[1:k-1])+1):Int(sum(bus_count[1:k]))
            dum1 = branch[i,15:17]; count1 = 0; ii=jj[1]; iii=kk[1]
            if length(findall(idx->idx=="A",[bus_ϕ[jj]; bus_ϕ[kk]])) == 2
                Sh_cap[ii+count1] = Sh_cap[ii+count1] + dum1[1]
                Sh_cap[iii+count1] = Sh_cap[iii+count1] + dum1[1]
                count1 = count1 + 1
            end
            if length(findall(idx->idx=="B",[bus_ϕ[jj]; bus_ϕ[kk]])) == 2
                Sh_cap[ii+count1] = Sh_cap[ii+count1] + dum1[2]
                Sh_cap[iii+count1] = Sh_cap[iii+count1] + dum1[1]
                count1 = count1 + 1
            end
            if length(findall(idx->idx=="C",[bus_ϕ[jj]; bus_ϕ[kk]])) == 2
                Sh_cap[ii+count1] = Sh_cap[ii+count1] + dum1[3]
                Sh_cap[iii+count1] = Sh_cap[iii+count1] + dum1[1]
            end
        end
    end

    ## Calculate submatrices of Ybus
    global Y = G+im*B; global Ym = abs.(Y); global Yθ = angle.(Y);
    global idx_ref = Int(sum(bus_count[1:bus_ref-1])+1):Int(sum(bus_count[1:bus_ref]))
    global idx_noref1 = 1:Int(sum(bus_count[1:bus_ref-1]))
    global idx_noref2 = Int(sum(bus_count[1:bus_ref])+1):nb_nph
    global idx_noref = [idx_noref1; idx_noref2]
    global nb_nph1 = nb_nph-3
    global Y00 = Y[idx_ref,idx_ref]
    global YLL = [Y[idx_noref1,idx_noref1] Y[idx_noref1,idx_noref2]; Y[idx_noref2,idx_noref1] Y[idx_noref2,idx_noref2] ];
    global YLL_inv = pinv(YLL); global Y0L = [Y[idx_ref,idx_noref1] Y[idx_ref,idx_noref2]]; global YL0 = Y0L';
end

#----------------------------------------------------------------------------------------#
