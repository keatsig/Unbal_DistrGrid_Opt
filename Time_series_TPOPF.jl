using JuMP, Ipopt, SparseArrays, LinearAlgebra, NamedArrays, DelimitedFiles, CSV, DataFrames, Statistics, PyPlot

# Include necessary functions
include("Functions/Ybus.jl");         # Function to calculate Ybus matrix (G,B)
include("Functions/Conventional.jl");   # Function to define constraints (polar)
include("Functions/Power_flow.jl");    # Function to calculate power flow
include("Functions/Cost.jl");          # Function to define cost (polar)
include("Functions/Misc.jl");                # Function to define start point,map input,print output
include("Functions/Parse_gld_allday.jl");      # Function to parse GridLabD system data
include("Functions/Fixed_point.jl");     # Function to implement fixed-point method (polar)
include("Functions/BFS.jl");             # Function to implement BFS method (polar)

## Main code
function TPOPF(case_name1,obj_op1)

    ## Load test case
    global case_name = case_name1; global obj_op = obj_op1;
    data_file = open(string("data/",case_name,"/data.xml")); data = readlines(data_file)
    z_file = open(string("data/",case_name,"/z_dump.xml")); z_dump = readlines(z_file)
    house_sch = CSV.read(string("data/Load_schedule_PS1.csv"),header=false,type=Float64,DataFrame)
    global house_sch = convert(Matrix,house_sch)/1e3;
    global n_house_row = size(house_sch,1); global n_house_col = size(house_sch,2)
    PV_sch = CSV.read(string("data/PV_schedule_PS1.csv"),header=false,type=Float64,DataFrame)
    global PV_sch = convert(Matrix,PV_sch)/1e3; global n_PV_col = size(PV_sch,2)
    # global VUF_init = readdlm(string("VUF_",case_name,"_init.csv"), ',')/100
    global QgY_init = readdlm(string("QgY0_",case_name,"_base.csv"), ',');
    # QgY_init = [zeros(size(QgY_init,1),120) QgY_init];

    ## Other inputs
    global PV_en = 1;   # 0-exclude, 1-include ZIP loads/PV inverters
    # crit_node = ["632"]
    crit_node = ["R1-12-47-1_node_359"]
    # crit_node = "all"
    global PV_rat = 0.020;
    global Q_pen = 1*1e-4; global ΔV_pen = 0*1e-2;

    ## Parse system data
    global n_sch = 1440;   # Time interval to run simulations
    global count_sch = 1;
    parse_gld(crit_node); global M = 1e3
    writedlm(string("House_",case_name,"_profile.csv"),house_prof , ',');
    writedlm(string("PV_",case_name,"_profile.csv"),PV_prof , ',');

    ## Create Y_bus and connection matrices
    Ybus(); conn_bgl();

    ## Initialize and start iterations
    start_point(4); Vm0 = s0[1:nb_nph]; θ0 = s0[nb_nph+1:2*nb_nph]; global QgY0 = zeros(ngy_inv);
    global Vm_soln = zeros(nb_nph,n_sch); global θ_soln = zeros(nb_nph,n_sch);
    global QgY0_soln = zeros(ngy_nph-3,n_sch); global Ploss_soln = zeros(n_sch);
    global soln_VUF = zeros(nb,3); global VU_soln = zeros(nb_3ph,3,n_sch); sts = "LOCALLY_SOLVED";
    while count_sch <= n_sch && (sts == "LOCALLY_SOLVED" || sts == "ALMOST_LOCALLY_SOLVED")

        ## Optimal power flow
        if count_sch%1 == 0
            # QgY0 = QgY_init[:,count_sch];
            # print("Time: ",hrs[count_sch],":",mins[count_sch],", ***OPTIMAL POWER-FLOW***\n")
            # @time (Vm0,θ0,QgY0,sts) = TPOPF_pol(crit_node,Vm0,θ0,QgY0)
            # @time (Vm0,θ0,QgY0,sts) = TPOPF_rect(crit_node,Vm0,θ0,QgY0)
            # @time (Vm0,θ0,QgY0,sts) = FOT_opf_pol(crit_node,Vm0,θ0,QgY0,0)
            # @time (Vm0,θ0,QgY0,sts) = FP_opf_rect(crit_node,Vm0,θ0,QgY0,0)
            # @time (Vm0,θ0,QgY0,sts) = FP_opf_pol(crit_node,Vm0,θ0,QgY0,0)
            # @time (Vm0,θ0,QgY0,sts) = BFS_opf_rect(crit_node,QgY0,0)
            if sts == "DID NOT CONVERGE"
                # start_point_pol(4); Vm0 = s0[1:nb_nph]; θ0 = s0[nb_nph+1:2*nb_nph]; QgY0 = zeros(ngy_nph-3);
                ΔV_pen = 0*1e-2; Q_pen = 1*1e-4;
                # @time (Vm0,θ0,QgY0,sts) = FOT_opf_pol(crit_node,Vm0,θ0,QgY0,1)
                # @time (Vm0,θ0,QgY0,sts) = FP_opf_rect(crit_node,Vm0,θ0,QgY0,1)
                # @time (Vm0,θ0,QgY0,sts) = FP_opf_pol(crit_node,Vm0,θ0,QgY0,1)
                # @time (Vm0,θ0,QgY0,sts) = BFS_opf_rect(crit_node,QgY0,1)
            end
            # (Vm0,θ0,QgY0,pf_iter) = NR_pf(Vm0,θ0,QgY0,50); print("NR iterations: ", pf_iter,"\n");
            (Vm0,θ0,QgY0,pf_iter) = FP_pf(Vm0,θ0,QgY0,50); print("FP iterations: ", pf_iter,"\n");

        ## Power flow
        else
            # QgY0 = QgY_init[:,count_sch];
            # (Vm0,θ0,QgY0,pf_iter) = NR_PF(Vm0,θ0,QgY0,50); print("NR iterations: ", pf_iter,"\n");
            # (Vm0,θ0,QgY0,pf_iter) = FP_PF(Vm0,θ0,QgY0,50);
            # print("Time: ",hrs[count_sch],":",mins[count_sch],", PF iterations:", pf_iter,"\n")
    end

        ## Store VU solution
        VU_soln[:,:,count_sch] = soln_VUF[idx_bus_3ph,:]; #Ploss_soln[count_sch] = P_loss;
        Vm_soln[:,count_sch] = Vm0; θ_soln[:,count_sch] = θ0; QgY0_soln[:,count_sch] = QgY0

        ## Update PV schedule
        yinv_powerP = PV_prof[:,count_sch+1];
        Qmax = sqrt.(yinv_ratedS.^2 .- yinv_powerP.^2); Qmin = -Qmax
        if obj_op == 1 || obj_op == 6
            Qmax = zeros(ngy_inv); Qmin = -Qmax;
        end
        for i=1:length(idx_ygen)
            if occursin("A",yinv_ph[i])
                ygen1[i,2] = yinv_powerP[i];   ygen1[i,3] = yinv_powerP[i];
                ygen1[i,4] = Qmin[i];          ygen1[i,5] = Qmax[i];
            elseif occursin("B",yinv_ph[i])
                ygen1[i,6] = yinv_powerP[i];   ygen1[i,7] = yinv_powerP[i];
                ygen1[i,8] = Qmin[i];          ygen1[i,9] = Qmax[i];
            else
                ygen1[i,10] = yinv_powerP[i];  ygen1[i,11] = yinv_powerP[i];
                ygen1[i,12] = Qmin[i];         ygen1[i,13] = Qmax[i];
            end
        end

        # Update load schedule
        # sch_Yfac = load_sch[r1_load+count_sch,r2_Yload];
        # sch_Dfac = load_sch[r1_load+count_sch,r2_Dload];
        # yload[1:nly-nlh,2:19] = yload[1:nly-nlh,2:19].*sch_Yfac
        # dload[:,2:19] = dload[:,2:19].*sch_Dfac
        Pbase = house_prof[:,count_sch+1];
        for i=1:nlh
            dum1 = [Pbase[i]*Pfrac[i] Pbase[i]*Ifrac[i] Pbase[i]*Zfrac[i]]
            dum2 = tan.(acos.([Ppf[i] Ipf[i] Zpf[i]])).*dum1
            if bus_ph[Int(hload[i,1])] == "AS" || bus_ph[Int(hload[i,1])] == "AN"
                yload[i+nly-nlh,2:6:14] = dum1; yload[i+nly-nlh,5:6:17] = dum2
            elseif bus_ph[Int(hload[i,1])] == "BS" || bus_ph[Int(hload[i,1])] == "BN"
                yload[i+nly-nlh,3:6:15] = dum1; yload[i+nly-nlh,6:6:18] = dum2
            elseif bus_ph[Int(hload[i,1])] == "CS" || bus_ph[Int(hload[i,1])] == "CN"
                yload[i+nly-nlh,4:6:16] = dum1; yload[i+nly-nlh,7:6:19] = dum2
            end
        end
        tot_load = (sum(sum(yload[:,i:i+2] for i=2:6:14)) + sum(sum(dload[:,i:i+2] for i=2:6:14)))
        PV_pen[count_sch+1] = sum(yinv_powerP)/tot_load*100;
        count_sch = count_sch + 1
    end


    ## Store other solutions in csv files
    writedlm(string("VUF_",case_name,"_soln.csv"),VU_soln[:,1,:] , ',');
    writedlm(string("Vm_",case_name,"_soln.csv"),Vm_soln, ',');
    writedlm(string("Va_",case_name,"_soln.csv"),θ_soln, ',');
    writedlm(string("QgY0_",case_name,"_soln.csv"),QgY0_soln , ',');
    writedlm(string("Ploss_",case_name,"_soln.csv"),Ploss_soln , ',');
    # writedlm(string("PVUR_",case_name,"_soln.csv"),VU_soln[:,3,:] , ',');
    # writedlm(string("LVUR_",case_name,"_soln.csv"),VU_soln[:,4,:] , ',');
end
