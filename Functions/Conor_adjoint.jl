using JuMP,Ipopt

Z=0.01+0.1*im
Y=inv(Z)
g=real(Y)
b=imag(Y)

myDict=Dict()

myDict["g"]=g
myDict["b"]=b
myDict["V_max"]=1.1
myDict["V_min"]=0.9
myDict["Pload"]=1
myDict["Qload"]=0
myDict["S_pv"]=0.5

# Initial Guess for NR powerflow
myDict["Var"]=Dict()
myDict["Var"]["V2"]=1
myDict["Var"]["T2"]=0


function Solve_PowerFlow(V1,P_PV,Q_PV)

    if V1<0.1
        println("V1 less that tol $(V1)")
        V1=0.9
    end
    V2= myDict["Var"]["V2"]
    T2=myDict["Var"]["T2"]

    Pload=myDict["Pload"]
    Qload=myDict["Qload"]
    g=myDict["g"]
    b=myDict["b"]

    for myiter=1:100

        f1= g*V2^2-V1*V2*(g*cos(T2)+b*sin(T2))+Pload-P_PV
        f2=-b*V2^2-V1*V2*(g*sin(T2)-b*cos(T2))+Qload-Q_PV

        f=[f1;f2]

        df1_dV2=  2*g*V2-V1*(g*cos(T2)+b*sin(T2))
        df2_dV2= -2*b*V2-V1*(g*sin(T2)-b*cos(T2))

        df1_dT2= -V1*V2*(-g*sin(T2)+b*cos(T2))
        df2_dT2= -V1*V2*( g*cos(T2)+b*sin(T2))

        J=[df1_dV2 df1_dT2;
           df2_dV2 df2_dT2;]

          dx = -inv(J)*f

          err = maximum(abs.(dx))

          if err <1e-12
              break
          end
          # println(err)
          # println(V2)
          # println(T2)

          V2=V2+dx[1]
          T2=T2+dx[2]

          if myiter==99
              println("No Convergence")
              # Make sure code throws an error!
              myDict["Var"]["V2"]=NaN
              myDict["Var"]["T2"]=NaN
              return
          end

    end

    myDict["Var"]["V2"]=V2
    myDict["Var"]["T2"]=T2


    df1_dV2=  2*g*V2-V1*(g*cos(T2)+b*sin(T2))
    df2_dV2= -2*b*V2-V1*(g*sin(T2)-b*cos(T2))

    df1_dT2= -V1*V2*(-g*sin(T2)+b*cos(T2))
    df2_dT2= -V1*V2*( g*cos(T2)+b*sin(T2))

    J=[df1_dV2 df1_dT2;
       df2_dV2 df2_dT2;]


    df1_dV1=-V2*(g*cos(T2)+b*sin(T2))
    df1_dP_PV=-1
    df1_dQ_PV=0

    df2_dV1=-V2*(g*sin(T2)-b*cos(T2))
    df2_dP_PV=0
    df2_dQ_PV=-1

    J2=[df1_dV1 df1_dP_PV df1_dQ_PV;
        df2_dV1 df2_dP_PV df2_dQ_PV;]

    dx_du=-inv((J))*J2

    return dx_du
    end

function Pinj_func(V1,P_PV,Q_PV)
    println("Calc Pinj at V1=$(round(V1,digits=2))     P_PV=$(round(P_PV,digits=2))     Q_PV=$(round(Q_PV,digits=2))")

    # Function to calulte Pinj

    # First solve the power flow for the control variables
    dx_du=Solve_PowerFlow(V1,P_PV,Q_PV)

    # Get the necessary values and calculate Pinj
    V2= myDict["Var"]["V2"]
    T2=myDict["Var"]["T2"]
    g=myDict["g"]
    b=myDict["b"]
    T1=0

    P12 = g*V1*V1-V1*V2*( g*cos(T1-T2)+b*sin(T1-T2))

    return P12
    end

function Pinj_func_diff(grad,V1,P_PV,Q_PV)
    println("Calc Pinj Grad at V1=$(round(V1,digits=2))     P_PV=$(round(P_PV,digits=2))     Q_PV=$(round(Q_PV,digits=2))")

     # same as function Pinj_func but calculate the sensitivity
    dx_du=Solve_PowerFlow(V1,P_PV,Q_PV)

    V2= myDict["Var"]["V2"]
    T2=myDict["Var"]["T2"]

    g=myDict["g"]
    b=myDict["b"]
    T1=0

    P12 = g*V1*V1-V1*V2*( g*cos(T1-T2)+b*sin(T1-T2))

    # Derivative wrt x and u
    dP_dV2=  -V1*( g*cos(T1-T2)+b*sin(T1-T2))
    dP_dT2=V1*V2*( -g*sin(T1-T2)+b*cos(T1-T2))
    dP_dV1=2*g*V1-V2*( g*cos(T1-T2)+b*sin(T1-T2))
    dp_dP_PV=0
    dp_dQ_PV=0

    dg_dx=[dP_dV2; dP_dT2]

    dg_du=[dP_dV1; dp_dP_PV; dp_dQ_PV]

    Q=dg_du+transpose(dx_du)*dg_dx

    grad[1]=Q[1]
    grad[2]=Q[2]
    grad[3]=Q[3]
    end

function Voltage2(V1,P_PV,Q_PV)
    println("Calc V2 at V1=$(round(V1,digits=2))     P_PV=$(round(P_PV,digits=2))     Q_PV=$(round(Q_PV,digits=2))")

    # Solve power flow to get value for V2
    dx_du=Solve_PowerFlow(V1,P_PV,Q_PV)
    V2=myDict["Var"]["V2"]
    return V2
    end

function diffVoltage2(g,V1,P_PV,Q_PV)
    println("Calc V2 Grad at V1=$(round(V1,digits=2))     P_PV=$(round(P_PV,digits=2))     Q_PV=$(round(Q_PV,digits=2))")

    # Calcualte sensitivy of V2 wrt u
    dx_du=Solve_PowerFlow(V1,P_PV,Q_PV)
    #println("Finished Sim")
    #println(dx_du)
    g[1]=dx_du[1,1]
    g[2]=dx_du[1,2]
    g[3]=dx_du[1,3]
    end

function Reduced_model()

    model = Model()
     # register functions and derivatives to model
    JuMP.register(model, :Voltage2, 3, Voltage2,diffVoltage2)
    JuMP.register(model, :Pinj_inner, 3, Pinj_func,Pinj_func_diff)


    @variable(model, myDict["V_min"]<=V1<=myDict["V_max"])
    @variable(model, P_PV)
    @variable(model, Q_PV)


    @NLconstraint(model, P_PV^2+Q_PV^2 <= myDict["S_pv"])
    @NLconstraint(model,myDict["V_min"]<= Voltage2(V1,P_PV,Q_PV))

    @NLobjective(model, Min, Pinj_inner(V1,P_PV,Q_PV))

     nlp_optimizer = with_optimizer(Ipopt.Optimizer,
     print_level=0,
     max_iter=100,
     tol=1e-3,
    # derivative_test="first-order",
     )
     JuMP.set_optimizer(model,nlp_optimizer)

     optimize!(model)

    println("\nReduced Space Results")
    println("-----------------------------")
    println("Obj is $(Pinj_func(JuMP.value(V1),JuMP.value(P_PV),JuMP.value(Q_PV)))")
    println("V1 is $(JuMP.value(V1))")
    println("P_PV is $(JuMP.value(P_PV))")
    println("Q_PV is $(JuMP.value(Q_PV))")
    println("V2 is $(myDict["Var"]["V2"])")
    end

# Reduced_model()




function test(x)
    s=0
    for i=1:2
        s = s+x^i
    end
    return s

end

function test_diff(x)
    ds=0
    for i=1:2
        ds = ds+i*x^(i-1)
    end
    return ds
end

function test_diff2(x)
    ds2=0
    for i=1:2
        ds2 = ds2+i*(i-1)*x^(i-2)
    end
    return ds2
end
