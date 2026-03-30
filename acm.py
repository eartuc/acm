import numpy as np
import scipy.io
import matplotlib.pyplot as plt

#Load data
d = scipy.io.loadmat("p_acm.mat", squeeze_me=True)

DJ = int(d["DJ"])

d["Niter"] = 5000
d["err_tol"] = 1.0e-8
d["speed"] = 0.4

#------------------------

def Solve_ACM(d):

    Niter = d["Niter"]
    err_tol = d["err_tol"]
    T = d["T"]
    tau = d["tau"]
    speed = d["speed"]
    price = d["price"]
    DJ = d["DJ"]
    nu = d["nu"]
    C = d["C"][:,:,np.newaxis]
    bta = d["bta"]
    theta = d["theta"][:,np.newaxis]
    alpha = d["alpha"][:,np.newaxis]
    psi = d["A"][:,np.newaxis]    
    rho = d["rho"][:,np.newaxis]
    
    L = d["L"]
    V = d["V"]
    L_star = np.copy(L)
    V_star = np.copy(V)
    
    err_big = True
    iter_count = 0

    while err_big and iter_count < Niter:

        iter_count += 1
        
        # SOLUTION ALGORITHM BEGINS
        i_tau = list(range(tau,tau-1+T))+[T-1] #index for convenience

        #Step 1: cpi and wages 
        cpi = np.exp(np.sum(np.log(price)*theta,axis=0,keepdims=True)) #CPI                
        wage = (price/cpi)*psi*alpha*(L**(rho-1))*(alpha*(L**rho)+(1-alpha))**((1/rho)-1) #Eq. 14

        #Step 2: moving probabilities
        mij_num = ((bta *V[:,i_tau])[np.newaxis,:,:]-C)/nu #Eq. 10, two parts       
        mij = np.exp(mij_num)/np.sum(np.exp(mij_num),axis=1,keepdims=True) #Eq. 10, two parts
        
        #Step 3: option values
        omega = -nu*np.log(np.diagonal(mij,axis1=0,axis2=1).T) #Eq. 7
        
        #Step 4: implied values
        V_star[:,0:T] = wage + bta*V[:,i_tau] + omega #Eq. 3
        
        #Step 5: implied labor allocations
        L_star[:,tau:T] = np.sum(L[:,0:T-tau][:,np.newaxis,:]*mij[:,:,0:T-tau],axis=0) #Labor flow
        
        #Step 6: check error and update
        err1 = 100*np.sum(np.abs(L-L_star))/T
        err2 = np.sum(np.abs(V-V_star))/T
        err = [err1, err2]
        
        L = (1-speed)*L+speed*L_star
        V = (1-speed)*V+speed*V_star
        
        if iter_count%100 == 0:
            err_big = (err1 > err_tol) or (err2 > err_tol)            
            print(f"Iteration {iter_count}, Errors: {err}")

        # SOLUTION ALGORITHM ENDS
    
    return {
        "V": V,
        "L": L,
        "wage": wage,
        "mij": mij,
        "err": err,
        "iter": iter_count,
        "cpi": cpi,
        "omega": omega
    }

#------------------------------------

# STARTING STEADY STATE
d["T"] = 1 #horizon is irrelevant for ss
d["tau"] = 0 #tau=0 means steady state

# Make sure: Shapes are correct
d["price"] = np.tile(d["price_ss"][:,np.newaxis],(1,d["T"]))
d["L"] = np.ones((DJ, d["T"])) / DJ #guess, total labor=1
d["V"] = (1/(1-d["bta"]))*np.ones((DJ,d["T"])) #guess

print("Solving initial steady state")
rs = Solve_ACM(d)

# TRANSITION
d["T"] = 30 #horizon=30
d["tau"] = 1 #tau=1 means transition

# Make sure: (1) L[:,0] is the initial labor alloc (2) Shapes are correct
d["price"] = np.tile(d["price_ft"][:,np.newaxis],(1,d["T"]))
d["V"] = np.tile(rs["V"],(1,d["T"])) #copy steady state
d["L"] = np.tile(rs["L"],(1,d["T"])) #copy ss

print("Solving transition")
rt = Solve_ACM(d)

print("Figures")
Lt = np.hstack( [ np.tile(rs["L"], (1, 5)), rt["L"]])
waget = np.hstack( [ np.tile(rs["wage"], (1, 5)), rt["wage"]])

DJ = d["DJ"]
T = d["T"]
t = np.arange(1, T + 6) - 6

plt.figure(figsize=(8, 5))
plt.plot(t, Lt.T, linewidth=1.5)
plt.xlabel('Time')
plt.ylabel('Labor')
plt.title('Labor allocation')
plt.grid(True)
plt.show()


plt.figure(figsize=(8, 5))
plt.plot(t, waget.T, linewidth=1.5)
plt.xlabel('Time')
plt.ylabel('Wage')
plt.title('Real wage')
plt.grid(True)
plt.show()


