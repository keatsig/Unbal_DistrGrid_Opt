# FUNCTION TO PARSE BUS DATA
function bus_parse(data::Dict{String,Any})

df_bus = data["bus"]
global nb = length(df_bus)
df_load = data["load"]
bus = zeros(nb,14)
for i=1:nb
   bus[i,1] = data["bus"][string(i)]["index"]
   bus[i,2] = data["bus"][string(i)]["bus_type"]
   bus[i,3] = data["bus"][string(i)]["vmin"][1]
   bus[i,4] = data["bus"][string(i)]["vmax"][1]
   bus[i,5] = data["bus"][string(i)]["vmin"][2]
   bus[i,6] = data["bus"][string(i)]["vmax"][2]
   bus[i,7] = data["bus"][string(i)]["vmin"][3]
   bus[i,8] = data["bus"][string(i)]["vmax"][3]
   bus[i,9] = data["bus"][string(i)]["vm"][1]
   bus[i,10] = data["bus"][string(i)]["va"][1]*180/pi
   bus[i,11] = data["bus"][string(i)]["vm"][2]
   bus[i,12] = data["bus"][string(i)]["va"][2]*180/pi
   bus[i,13] = data["bus"][string(i)]["vm"][3]
   bus[i,14] = data["bus"][string(i)]["va"][3]*180/pi
end
if isempty(findall(x->x==3,bus[:,2]))
   bus[1,2] = 3
end
   return bus
end

#----------------------------------------------------------------------------------------#

# FUNCTION TO PARSE GENERATOR DATA
function gen_parse(data::Dict{String,Any})

df_gen = data["gen"]
global ng = length(df_gen)
base = data["baseMVA"]
gen = zeros(ng,20)
for i=1:ng
        gen[i,1] = data["gen"][string(i)]["gen_bus"]
        gen[i,2] = data["gen"][string(i)]["pmin"][1]*base
        gen[i,3] = data["gen"][string(i)]["pmax"][1]*base
        gen[i,4] = data["gen"][string(i)]["qmin"][1]*base
        gen[i,5] = data["gen"][string(i)]["qmax"][1]*base
        gen[i,6] = data["gen"][string(i)]["pmin"][2]*base
        gen[i,7] = data["gen"][string(i)]["pmax"][2]*base
        gen[i,8] = data["gen"][string(i)]["qmin"][2]*base
        gen[i,9] = data["gen"][string(i)]["qmax"][2]*base
        gen[i,10] = data["gen"][string(i)]["pmin"][3]*base
        gen[i,11] = data["gen"][string(i)]["pmax"][3]*base
        gen[i,12] = data["gen"][string(i)]["qmin"][3]*base
        gen[i,13] = data["gen"][string(i)]["qmax"][3]*base
        gen[i,14] = data["gen"][string(i)]["pg"][1]*base
        gen[i,15] = data["gen"][string(i)]["qg"][1]*base
        gen[i,16] = data["gen"][string(i)]["pg"][2]*base
        gen[i,17] = data["gen"][string(i)]["qg"][2]*base
        gen[i,18] = data["gen"][string(i)]["pg"][3]*base
        gen[i,19] = data["gen"][string(i)]["qg"][3]*base
        gen[i,20] = data["gen"][string(i)]["gen_status"]
end
    return gen
end

#----------------------------------------------------------------------------------------#

# FUNCTION TO PARSE GENERATOR COST DATA
function gencost_parse(data::Dict{String,Any})

df_gencost = data["gen"]
global ngc = length(df_gencost)
base = data["baseMVA"]  
n=0
for i=1:ngc
    if data["gen"][string(i)]["ncost"] >= n
        n = data["gen"][string(i)]["ncost"]
    end
end
gencost = zeros(ngc,4+n)
for i=1:ngc
        gencost[i,1] = data["gen"][string(i)]["model"]
        gencost[i,2] = data["gen"][string(i)]["startup"]
        gencost[i,3] = data["gen"][string(i)]["shutdown"]
        gencost[i,4] = data["gen"][string(i)]["ncost"]
        for j=1:n
            gencost[i,j+4] = data["gen"][string(i)]["cost"][j]/base^(n-j)
        end        
end
    return gencost
end

#----------------------------------------------------------------------------------------#

# FUNCTION TO PARSE BRANCH DATA
function branch_parse(data::Dict{String,Any})

df_branch = data["branch"]
global nbr = length(df_branch)
base = data["baseMVA"]
branch = zeros(nbr,23)
for i=1:nbr
        branch[i,1] = data["branch"][string(i)]["f_bus"]
        branch[i,2] = data["branch"][string(i)]["t_bus"]
        branch[i,3] = data["branch"][string(i)]["br_r"][1,1]
        branch[i,4] = data["branch"][string(i)]["br_x"][1,1]
        branch[i,5] = data["branch"][string(i)]["br_r"][1,2]
        branch[i,6] = data["branch"][string(i)]["br_x"][1,2]
        branch[i,7] = data["branch"][string(i)]["br_r"][1,3]
        branch[i,8] = data["branch"][string(i)]["br_x"][1,3]
        branch[i,9] = data["branch"][string(i)]["br_r"][2,2]
        branch[i,10] = data["branch"][string(i)]["br_x"][2,2]
        branch[i,11] = data["branch"][string(i)]["br_r"][2,3]
        branch[i,12] = data["branch"][string(i)]["br_x"][2,3]
        branch[i,13] = data["branch"][string(i)]["br_r"][3,3]
        branch[i,14] = data["branch"][string(i)]["br_x"][3,3]
        branch[i,15] = data["branch"][string(i)]["b_fr"][1] + data["branch"][string(i)]["b_to"][1]
        branch[i,16] = data["branch"][string(i)]["b_fr"][2] + data["branch"][string(i)]["b_to"][2]
        branch[i,17] = data["branch"][string(i)]["b_fr"][3] + data["branch"][string(i)]["b_to"][3]
        branch[i,18] = data["branch"][string(i)]["rate_a"][1]*base
        branch[i,19] = data["branch"][string(i)]["rate_b"][1]*base
        branch[i,20] = data["branch"][string(i)]["rate_c"][1]*base
        branch[i,21] = data["branch"][string(i)]["angmin"][1]*180/pi
        branch[i,22] = data["branch"][string(i)]["angmax"][1]*180/pi
        branch[i,23] = data["branch"][string(i)]["br_status"]
end
    return branch
end

#----------------------------------------------------------------------------------------#