using MAT
using LinearAlgebra
using Printf
using Plots

# Load data
d = matread("p_acm.mat")

DJ = Int(d["DJ"])

d["Niter"] = 5000 #max iterations
d["err_tol"] = 1.0e-8 #tolerance
d["speed"] = 0.4 #update speed


#------------------------
function Solve_ACM(d)
    
    Niter = Int(d["Niter"])
    err_tol = d["err_tol"]
    T = Int(d["T"])
    tau = Int(d["tau"])
    speed = d["speed"]
    price = d["price"]
    DJ = Int(d["DJ"])
    nu = d["nu"]
    C = d["C"]
    theta = d["theta"]
    alpha = d["alpha"]
    psi = d["A"]
    bta = d["bta"]
    rho = d["rho"]
    
    L = d["L"]
    V = d["V"]
    L_star = copy(L)
    V_star = copy(V)
    
    err_big = true
    iter = 0
    
    diag_mask_T = repeat(Matrix{Bool}(I,DJ,DJ),1,1,T)
    
    #for scope 
    wage = zeros(DJ,T)    
    mij = zeros(DJ,DJ,T)
    cpi = ones(1,T)
    omega = zeros(DJ,T)
    err = [0.0, 0.0]
    
    while err_big && iter < Niter
        
        iter += 1
        
        # SOLUTION ALGORITHM BEGINS
        #if tau=0 then i_tau=[1, 2, 3, .. T-1, T]
        #if tau=1 then i_tau=[2, 3, .. T-1, T, T] thus last period repeats
        i_tau = vcat(collect((tau+1):(tau-1+T)),T) #index for convenience
        
        #Step 1: cpi and wages 
        cpi = exp.(sum(log.(price).*theta,dims=1))  #CPI      
        wage = (price./cpi).*psi.*alpha.*(L.^(rho.-1.0)).*(alpha.*(L.^rho).+(1.0.-alpha)).^((1.0./rho).-1.0) #Eq. 14
        
        #Step 2: moving probabilities
        mij_num = (reshape(bta.*V[:,i_tau],1,DJ,T).-C) ./ nu #Eq. 10, split into two
        mij = exp.(mij_num)./sum(exp.(mij_num),dims=2) #Eq. 10, , split into two
        
        #Step 3: option values
        omega = -nu.*log.(reshape(mij[diag_mask_T],DJ,T)) #Eq. 7
        
        #Step 4: implied values
        V_star[:,1:T] = wage.+bta.*V[:,i_tau].+omega #Eq.3
        
        #Step 5: implied labor allocations
        L_star[:,tau+1:T] = reshape(sum(reshape(L[:,1:T-tau], DJ, 1, T-tau).*mij[:,:,1:T-tau],dims=1),DJ,T-tau) #labor flow
        
        #Step 6: check error and update
        err1 = 100 * sum(abs.(L.-L_star))/T
        err2 = sum(abs.(V.-V_star))/T
        err = [err1, err2]
        
        L = (1-speed).*L.+speed.*L_star
        V = (1-speed).*V.+speed.*V_star
        
        if iter % 100 == 0
            err_big = (err1 > err_tol || err2 > err_tol)
            @printf("Iteration %d, Errors: %s\n", iter, err)
        end

        # SOLUTION ALGORITHM ENDS

    end #while
    
    return Dict(
        "V" => V,
        "L" => L,
        "wage" => wage,
        "mij" => mij,
        "err" => err,
        "iter" => iter,
        "omega" => omega,
        "cpi" => cpi
    )
end

#--------------------------------------

# STARTING STEADY STATE
d["T"] = 1 #horizon is irrelevant for ss
d["tau"] = 0 #tau=1 means steady state
d["price"] = d["price_ss"] #initial prices

# Make sure: Shapes are correct
d["L"] = ones(DJ, d["T"]) ./ DJ #guess, total labor=1
d["V"] = (1/(1-d["bta"])).*ones(DJ,d["T"]) #guess

println("Solving initial steady state")
rs = Solve_ACM(d)

# TRANSITION
d["T"] = 30 #horizon=30
d["tau"] = 1 #tau=1 means transition
d["price"] = repeat(d["price_ft"], 1, d["T"])

# Make sure: (1) L[:,1] is the initial labor alloc (2) Shapes are correct
d["V"] = repeat(rs["V"], 1, d["T"]) #copy steady state
d["L"] = repeat(rs["L"], 1, d["T"]) #copy steady state

println("Solving transition")
rt = Solve_ACM(d)

println("Figures")
Lt=[repeat(rs["L"],1,5) rt["L"]]
waget=[repeat(rs["wage"],1,5) rt["wage"]]

DJ = d["DJ"] 
T = d["T"] 
t = (1:(T+5)) .- 6

plot(t, Lt', 
     linewidth = 1.5,
     xlabel = "Time",
     ylabel = "Labor",
     title = "Labor allocation",
     legend = false,
     grid = true)

plot(t, waget', 
     linewidth = 1.5,
     xlabel = "Time",
     ylabel = "Wage",
     title = "Real wage",
     legend = false,
     grid = true)

