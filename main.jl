using JuMP, Ipopt, KNITRO, SparseArrays, LinearAlgebra, NamedArrays, DelimitedFiles, DataFrames, CSV, Statistics, PyPlot

# Include necessary functions
include("Functions/Ybus.jl");         # Function to calculate Ybus matrix (G,B)
include("Functions/Conventional.jl");   # Function to define constraints (polar)
include("Functions/Power_flow.jl");    # Function to calculate power flow
include("Functions/Cost.jl");          # Function to define cost (polar)
include("Functions/Misc.jl");                # Function to define start point,map input,print output
include("Functions/Parse_gld_1day.jl");      # Function to parse GridLabD system data
include("Functions/Fixed_point.jl");     # Function to implement fixed-point method
include("Functions/BFS.jl");             # Function to implement BFS method
include("Functions/Dist3Flow.jl");       # Function to implement Dist3Flow method
include("Functions/Active_constr.jl");   # Function to find active constraints

# Main code
function TPOPF(case_name1,obj_op1)

    ## Load test case
    global case_name = case_name1; global obj_op = obj_op1;
    global PV_en = 1;   # 0-exclude, 1-include ZIP loads/PV inverters
    if PV_en == 1
        house_sch = CSV.read(string("data/",case_name,"/House_",case_name,"_profile.csv"),header=false,type=Float64,DataFrame)
        global house_sch = convert(Matrix,house_sch);
        PV_sch = CSV.read(string("data/",case_name,"/PV_",case_name,"_profile.csv"),header=false,type=Float64,DataFrame)
        global PV_sch = convert(Matrix,PV_sch);
    end
    # global VUF_init = readdlm(string("VUF_",case_name,"_init.csv"), ',')/100
    # global QgY_init = readdlm(string("QgY0_",case_name,"_base.csv"), ',');
    # global QgY_init = readdlm("QgY0_R1-12.47-1_freq60_delay60.csv", ',');

    ## Other inputs
    global count_sch = 1;  # Time instant to run OPF
    # global crit_node = ["632"]
    # global crit_node = ["l_49"]
    # global crit_node = ["R1-12-47-1_node_359"]
    # global crit_node = ["R1-12-47-3_node_1"]
    global crit_node = "all"
    global PV_rat = 0.035;
    global Q_pen = 0*1e-5; global ΔV_pen = 0*1e-2;
    global tap_en = 0; # Enable taps as variables only for ivr formulation

    ## Parse system data
    parse_gld(); global M = 1e3

    ## Create Y_bus and connection matrices
    Ybus(); conn_bgl();

    # Initialize values
    start_point(4); global Vm0 = s0[1:nb_nph]; global θ0 = s0[nb_nph+1:2*nb_nph];
    global Vd0 = Vm0.*cos.(θ0); global Vq0 = Vm0.*sin.(θ0); global QgY0 = zeros(ngy_inv);
    stat = "LOCALLY_SOLVED";

    ## Power flow
    # QgY0 = QgY_init[:,count_sch];
    # (Vm0,θ0,QgY0,pf_iter) = NR_pf(Vm0,θ0,QgY0,50); print("No of NR iterations: ", pf_iter,"\n");
    # (Vm0,θ0,QgY0,pf_iter) = BFS_pf(Vm0,θ0,QgY0,100); print("No of BFS iterations: ", pf_iter,"\n");
    # (Vm0,θ0,QgY0,pf_iter) = FP_pf(Vm0,θ0,QgY0,50); print("No of FP iterations: ", pf_iter,"\n");
    # (Vm0,θ0,QgY0,pf_iter) = D3F_pf(Vm0,θ0,QgY0,50); print("No of DF iterations: ", pf_iter,"\n");

    ## Optimal power flow
    # (Vm0,θ0,QgY0,stat) = FOT_opf_pol(crit_node,Vm0,θ0,QgY0,0)
    # (Vm0,θ0,QgY0,stat) = FP_opf_rect(crit_node,Vm0,θ0,QgY0,0)
    # (Vm0,θ0,QgY0,stat) = BFS_opf_rect(crit_node,Vm0,θ0,QgY0,0)

    (Vm0,θ0,QgY0,stat) = TPOPF_pol(crit_node,Vm0,θ0,QgY0)
    # (Vm0,θ0,QgY0,stat) = TPOPF_rect(crit_node,Vm0,θ0,QgY0)
    # (Vm0,θ0,QgY0,stat) = TPOPF_ivr(crit_node,Vm0,θ0,QgY0)
   
    # (Vm0,θ0,QgY0,stat) = D3F_opf_rect(crit_node,Vm0,θ0,QgY0,0)
    # (Vm0,θ0,QgY0,stat) = TPOPF_act_constr(crit_node,Vm0,θ0,QgY0)

    if stat == "DID NOT CONVERGE"
        # start_point_pol(4); Vm0 = s0[1:nb_nph]; θ0 = s0[nb_nph+1:2*nb_nph]; QgY0 = zeros(ngy_nph-3);
        ΔV_pen = 0*1e-2; Q_pen = 1*1e-3;
        # @time (Vm0,θ0,QgY0,stat) = FOT_opf_pol(crit_node,Vm0,θ0,QgY0,1)
        # @time (Vm0,θ0,QgY0,stat) = FP_opf_rect(crit_node,Vm0,θ0,QgY0,1)
        # @time (Vm0,θ0,QgY0,stat) = FP_opf_pol(crit_node,Vm0,θ0,QgY0,1)
        # @time (Vm0,θ0,QgY0,stat) = BFS_opf_rect(crit_node,QgY0,1)
    end

    ## Print solution
    Vph = Vm0.*exp.(im*θ0); S_inj = Vph.*conj(Y*Vph);
    Vm1 = round.(Vm0,RoundNearest, digits=4); θ1 = round.(θ0*180/pi,RoundNearest, digits=4);
    P_inj = round.(real(S_inj),RoundNearest, digits=3);
    Q_inj = round.(imag(S_inj),RoundNearest, digits=3);
    global soln_VθPQ = NamedArray( [Vm1 θ1 P_inj Q_inj], ([i for i in bus_name], [:Vm_pu, :θ_deg, :Pinj_pu, :Qinj_pu]), ("Bus\n(no./ph)","Par"))
    # print("The substation power is ", sqrt(sum(P_inj[idx_ref].^2+Q_inj[idx_ref].^2)),"\n");
    # print("The substation active power is ", sum(P_inj[idx_ref]),"\n");
    # print("The network losses are ", sum(P_inj),"\n");
    # print("The total reactive power injection is ", sum(QgY0),"\n");
    print("------------------------------------------------------------\n")
    # @show soln_VθPQ
    print("------------------------------------------------------------\n")
    x = length(idx_bus_3ph)
    println("Total VUF   = ",round(sum(soln_VUF[:,1]),RoundNearest, digits=3),"%,  Average VUF   = ",round(sum(soln_VUF[:,1])/x,RoundNearest, digits=3),"%, Maximum VUF   = ",findmax(soln_VUF[:,1])[1],"%")
    println("Total LVUR  = ",round(sum(soln_VUF[:,2]),RoundNearest, digits=3),"%,  Average LVUR  = ",round(sum(soln_VUF[:,2])/x,RoundNearest, digits=3),"%, Maximum LVUR  = ",findmax(soln_VUF[:,2])[1],"%")
    println("Total PVUR  = ",round(sum(soln_VUF[:,3]),RoundNearest, digits=3),"%,  Average PVUR  = ",round(sum(soln_VUF[:,3])/x,RoundNearest, digits=3),"%, Maximum PVUR  = ",findmax(soln_VUF[:,3])[1],"%")
    # soln_VU = NamedArray( soln_VUF[idx_bus_3ph,:], ([i for i in bus_ID[idx_bus_3ph]], [:VUF_neg__pc, :LVUR_pc, :PVUR_pc]), ("Bus","Par"))
    # @show soln_VU

    solnV = zeros(nb,6); num_count=1
    for i=1:nb
        if occursin("A", bus_ph[i])
            solnV[i,1:2] = [Vm1[num_count] θ1[num_count]];
            num_count = num_count+1;
        end
        if occursin("B", bus_ph[i])
            solnV[i,3:4] = [Vm1[num_count] θ1[num_count]];
            num_count = num_count+1;
        end
        if occursin("C", bus_ph[i])
            solnV[i,5:6] = [Vm1[num_count] θ1[num_count]];
            num_count = num_count+1;
        end
    end

    ## Store other solutions in csv files
    # writedlm(string("QgY0_",case_name,"_soln1.csv"),QgY0 , ',');
    # writedlm(string("soln_",case_name,"_soln1.csv"),[bus_name soln_VθPQ] , ',');
    # writedlm("soln.csv",solnV , ',');
    # writedlm("idx_phA.csv",idx_phA , ','); writedlm("idx_phB.csv",idx_phB , ','); writedlm("idx_phC.csv",idx_phC , ',');


    ## Plot figures
    # figure(figsize=(5,5))
    # plot(Vm0,"-o"); grid("on");
    # xlabel("Three-phase nodes"); ylabel("Voltage [pu]")
end
