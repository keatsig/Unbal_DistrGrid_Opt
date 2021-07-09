# FUNCTION TO CHECK VOLTAGE UNBALANCE OF GRIDLABD SOLUTION
function VUcheck_glm(case_name)

data_file = open(string("data/",case_name,"/data.xml"))
data = readlines(data_file)
z_file = open(string("data/",case_name,"/z_dump.xml"))
z_dump = readlines(z_file);
    
## PARSE BUS DATA
idx_node_f = findall(idx -> idx == "\t\t\t<node>", data); idx_node_t = findall(idx -> idx == "\t\t\t</node>", data)
idx_load_f = findall(idx -> idx == "\t\t\t<load>", data); idx_load_t = findall(idx -> idx == "\t\t\t</load>", data)
idx_tnode_f = findall(idx -> idx == "\t\t\t<triplex_node>", data); idx_tnode_t = findall(idx -> idx == "\t\t\t</triplex_node>", data)
idx_meter_f = findall(idx -> idx == "\t\t\t<meter>", data); idx_meter_t = findall(idx -> idx == "\t\t\t</meter>", data)
idx_tmeter_f = findall(idx -> idx == "\t\t\t<triplex_meter>", data); idx_tmeter_t = findall(idx -> idx == "\t\t\t</triplex_meter>", data)
idx_bus_f = [idx_node_f; idx_load_f; idx_tnode_f]; idx_bus_t = [idx_node_t; idx_load_t; idx_tnode_t]

# Calculate total number of nodes
nb_n = length(idx_node_f); nb_l = length(idx_load_f); nb_tn = length(idx_tnode_f); 
nb_m = length(idx_meter_f); nb_tm = length(idx_tmeter_f);
global nb = nb_n+nb_l+nb_tn

# Create bus_ID matrix to map node-names to node-numbers
bus_ID_nn = get_data(data,nb_n,idx_node_f,idx_node_t,"<name>",3)
bus_ID_np = get_data(data,nb_n,idx_node_f,idx_node_t,"<parent>",3)
bus_ID_ln = get_data(data,nb_l,idx_load_f,idx_load_t,"<name>",3)
bus_ID_lp = get_data(data,nb_l,idx_load_f,idx_load_t,"<parent>",3)
bus_ID_tnn = get_data(data,nb_tn,idx_tnode_f,idx_tnode_t,"<name>",3)
bus_ID_tnp = get_data(data,nb_tn,idx_tnode_f,idx_tnode_t,"<parent>",3)
bus_ID_nam = [bus_ID_nn; bus_ID_ln; bus_ID_tnn]; bus_ID_par = [bus_ID_np; bus_ID_lp; bus_ID_tnp]; 
global bus_ID = bus_ID_nam;
bus_type = get_data(data,nb,idx_bus_f,idx_bus_t,"<bustype>",3)

# Impose minimum & maximum voltage constraints
bus_ph = get_data(data,nb,idx_bus_f,idx_bus_t,"<phases>",3)
global nb_nph = 0
for i=1:nb
    dum1 = length(split(split(bus_ph[i],r"S|D|N")[1],""))
    nb_nph = nb_nph + dum1
end
global bus_ϕ = repeat(["0"],nb_nph); global bus_count = zeros(nb)
count=1
for i=1:nb
    dum1 = split(split(bus_ph[i],r"S|D|N")[1],"")
    bus_count[i] = length(dum1)
    bus_ϕ[count:count+Int(bus_count[i])-1] = dum1
    count = Int(count+bus_count[i])
end
global bus_name = repeat(["0"],nb_nph); count = 1
for i=1:nb
    for j=Int.(sum(bus_count[1:i-1])+1:sum(bus_count[1:i]))
        bus_name[count] = string(count,":",bus_ID[i],"::",bus_ϕ[j])
        count = count + 1
    end
end
V_min = ones(nb,3)*0.9; V_max = ones(nb,3)*1.1;
for i=1:nb
    if bus_ph[i] == "AN" || bus_ph[i] == "AS" || bus_ph[i] == "A" 
        V_min[i,2:3] .= 0;   V_max[i,2:3] .= 0;   
    elseif bus_ph[i] == "BN" || bus_ph[i] == "BS" || bus_ph[i] == "B" 
        V_min[i,1:2:3] .= 0; V_max[i,1:2:3] .= 0; 
    elseif bus_ph[i] == "CN" || bus_ph[i] == "CS" || bus_ph[i] == "C" 
        V_min[i,1:2] .= 0;   V_max[i,1:2] .= 0;  
    elseif bus_ph[i] == "ABN" || bus_ph[i] == "ABD" || bus_ph[i] == "BAN" || bus_ph[i] == "BAD" || bus_ph[i] == "BA" || bus_ph[i] == "AB"
        V_min[i,3] = 0;      V_max[i,3] = 0;     
    elseif bus_ph[i] == "BCN" || bus_ph[i] == "BCD" || bus_ph[i] == "CBN"  || bus_ph[i] == "CBD" || bus_ph[i] == "CB" || bus_ph[i] == "BC" 
        V_min[i,1] = 0;      V_max[i,1] = 0;     
    elseif bus_ph[i] == "ACN" || bus_ph[i] == "ACD" || bus_ph[i] == "CAN" || bus_ph[i] == "CAD" || bus_ph[i] == "CA" || bus_ph[i] == "AC"
        V_min[i,2] = 0;      V_max[i,2] = 0;     
    end
end

# Read current powerflow solution from GridLabD (warm-start)
bus_ID_lc = get_data(data,nb_l,idx_load_f,idx_load_t,"<load_class>",3)
bus_Vnom = str_num(get_data(data,nb,idx_bus_f,idx_bus_t,"<nominal_voltage>",3))
for i=1:nb_l
    if bus_ID_lc[i] == "C"
       bus_Vnom[i+nb_n] = 480/sqrt(3)
    end
end
p2n = ["A" "B" "C"]; l2l = ["AB" "BC" "CA"];
global s0 = spzeros(nb,6);
for i=1:3
    s0[1:nb_n+nb_l,2*i-1:2*i] = warm_start(get_data(data,nb_n+nb_l,[idx_node_f; idx_load_f],[idx_node_t; idx_load_t],string("<voltage_",p2n[i],">"),3))
end
dum1 = warm_start(get_data(data,nb_tn,idx_tnode_f,idx_tnode_t,string("<voltage_1>"),3))
for i=1:nb_tn
    if dum1[i,2] >= -60 && dum1[i,2] <= 60
        s0[nb_n+nb_l+i,1:2] = dum1[i,:]
    elseif dum1[i,2] >= -180 && dum1[i,2] <= -60
        s0[nb_n+nb_l+i,3:4] = dum1[i,:]
    elseif dum1[i,2] >= 60 && dum1[i,2] <= 180
        s0[nb_n+nb_l+i,5:6] = dum1[i,:]
    end
end
s0[:,1:2:5] = s0[:,1:2:5]./ [bus_Vnom bus_Vnom bus_Vnom]

# Create bus matrix
bus = zeros(nb,14)
bus[:,1] = 1:nb 
for i=1:nb
    if bus_type[i] == "SWING"
        bus[i,2] = 3   
    elseif bus_type[i] == "PV"
        bus[i,2] = 2
    elseif bus_type[i] == "PQ"
        bus[i,2] = 1
    end
end
bus[:,3:2:7] = V_min
bus[:,4:2:8] = V_max
bus[:,9:14] = s0;
global idx_bus_3ph = findall(idx -> idx == 3,bus_count)
global nb_3ph = length(idx_bus_3ph); 
global idx_bus_3ph1 = spzeros(nb_3ph*3)
for i=1:nb_3ph
    j = idx_bus_3ph[i]
    idx_bus_3ph1[i*3-2:i*3] = sum(bus_count[1:j-1])+1:sum(bus_count[1:j])
end

global ap = zeros(nb_nph); global an = zeros(nb_nph)
for i=1:nb_nph
    if bus_ϕ[i] == "B"
        ap[i] = 2*pi/3; an[i] = -2*pi/3
    end
    if bus_ϕ[i] == "C"
        ap[i] = -2*pi/3; an[i] = 2*pi/3
    end
end
V1 = zeros(nb_nph); θ1 = zeros(nb_nph); count = 1
for i=1:nb
    if bus_count[i] == 3
        V1[count:count+2] = s0[i,1:2:5]
        θ1[count:count+2] = s0[i,2:2:6]*pi/180
        count = count+3
    elseif bus_count[i] == 2 && bus[i,3] != 0 && bus[i,5] != 0 && bus[i,7] == 0
        V1[count:count+1] = s0[i,1:2:3]
        θ1[count:count+1] = s0[i,2:2:4]*pi/180
        count = count+2
    elseif bus_count[i] == 2 && bus[i,5] != 0 && bus[i,7] != 0 && bus[i,3] == 0
        V1[count:count+1] = s0[i,3:2:5]
        θ1[count:count+1] = s0[i,4:2:6]*pi/180
        count = count+2
    elseif bus_count[i] == 2 && bus[i,3] != 0 && bus[i,7] != 0 && bus[i,5] == 0
        V1[count:count+1] = s0[i,1:4:5]
        θ1[count:count+1] = s0[i,2:4:6]*pi/180
        count = count+2        
    elseif bus_count[i] == 1 && bus[i,3] != 0 && bus[i,5] == 0 && bus[i,7] == 0
        V1[count] = s0[i,1]
        θ1[count] = s0[i,2]*pi/180
        count = count+1        
    elseif bus_count[i] == 1 && bus[i,3] == 0 && bus[i,5] != 0 && bus[i,7] == 0
        V1[count] = s0[i,3]
        θ1[count] = s0[i,4]*pi/180
        count = count+1          
    elseif bus_count[i] == 1 && bus[i,3] == 0 && bus[i,5] == 0 && bus[i,7] != 0
        V1[count] = s0[i,1]
        θ1[count] = s0[i,2]*pi/180
        count = count+1
    end
end
soln = zeros(nb); 
for i=1:nb
    Vn1 = 0; Vp1 = 0; 
    if bus_count[i] == 3
        j=Int.(sum(bus_count[1:i-1])+1:sum(bus_count[1:i]))    
        Vn1 = (sqrt(sum(V1[j].*cos.(θ1[j]+an[j]) )^2 + sum(V1[j].*sin.(θ1[j]+an[j]) )^2 )/3)
        Vp1 = (sqrt(sum(V1[j].*cos.(θ1[j]+ap[j]) )^2 + sum(V1[j].*sin.(θ1[j]+ap[j]) )^2 )/3)
        soln[i] = round(Vn1/Vp1*100, RoundNearest, digits=4)
    end
end
x=findall(idx->idx>=0.0001,soln);
println(sum(soln[x])," ",findmax(soln[x]))
plot(soln)
       
end
#----------------------------------------------------------------------------------------#

