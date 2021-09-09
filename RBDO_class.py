# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 08:39:11 2021

@author: lbrevaul
"""

import openturns as ot 
## Class to transform the Python function of evaluation of the Probabilistic inequality constraints using FORM (Pineq_eval_FORM) into a OpenTURNSPythonFunction with the appropriate input and output dimensions
class RBDO_PIneq_eval_FORM(ot.OpenTURNSPythonFunction):

    def __init__(self,RBDO_Problem):
        
        self.RBDO_Problem = RBDO_Problem
        self.dim_output = self.RBDO_Problem.NbPIneqCons+self.RBDO_Problem.NbDIneqCons  #dimension of the constraint vector of both probabilistic and deterministic constraints
        
        super(RBDO_PIneq_eval_FORM, self).__init__(self.RBDO_Problem.len_d+self.RBDO_Problem.len_p, self.dim_output)
                  
    def  _exec(self,XX):
        
        Y = self.RBDO_Problem.evalConsList_FORM(XX)

        return Y
## Class to transform the Python function of evaluation of the Probabilistic inequality constraints using invFORM (Pineq_eval_invFORM) into a OpenTURNSPythonFunction with the appropriate input and output dimensions
class RBDO_PIneq_eval_invFORM(ot.OpenTURNSPythonFunction):

    def __init__(self,RBDO_Problem):
        
        self.RBDO_Problem = RBDO_Problem
        self.dim_output = self.RBDO_Problem.NbPIneqCons+self.RBDO_Problem.NbDIneqCons  #dimension of the constraint vector of both probabilistic and deterministic constraints
        
        super(RBDO_PIneq_eval_invFORM, self).__init__(self.RBDO_Problem.len_d+self.RBDO_Problem.len_p, self.dim_output)
                  
    def  _exec(self,XX):
        
        Y = self.RBDO_Problem.evalConsList_invFORM(XX)

        return Y

## Class to transform the Python function of evaluation of the Probabilistic inequality constraints for SORA algorithm into a OpenTURNSPythonFunction with the appropriate input and output dimensions    
class RBDO_PIneq_eval_SORA(ot.OpenTURNSPythonFunction):

    def __init__(self,RBDO_Problem):
        
        self.RBDO_Problem = RBDO_Problem
        self.dim_output = self.RBDO_Problem.NbPIneqCons+self.RBDO_Problem.NbDIneqCons #dimension of the constraint vector of both probabilistic and deterministic constraints
        
        super(RBDO_PIneq_eval_SORA, self).__init__(self.RBDO_Problem.len_d+self.RBDO_Problem.len_p, self.dim_output)
                  
    def  _exec(self,XX):
        
        Y = self.RBDO_Problem.evalConsList_SORA(XX)

        return Y

## Class of RIA algorithm (Reliability Index Approach), takes as input:
# RBDO_Problem
# Solver
# InitialPoint
# Runs an optimization problem, and for each Probabilistic Inequality Constraint it runs a FORM method to estimate the probability of failure (and the reliability index)
# Careful, openTURNS optimization problem assumes inequality constraints >=0 as feasible   (opposite of classical RBDO problem which assumes constraint feasibility as <=0)
class RIA_Algorithm(object):

    def __init__(self,RBDO_Problem,Solver,InitialPoint):
        
        self.RBDO_Problem = RBDO_Problem
        self.RBDO_PIneq_eval_FORM = ot.Function(RBDO_PIneq_eval_FORM(self.RBDO_Problem))  #transformation into OpenTURNS function
        self.Solver = Solver
        self.InitialPoint = InitialPoint    

        
    def run(self):
        
        #Optimization problem definition
        algo = self.Solver 
        objective = self.RBDO_Problem.Objective     #objective function
        inequality_constraint = self.RBDO_PIneq_eval_FORM   #constraint functions
        bounds = self.RBDO_Problem.Bounds #research optimization bounds
        problem = ot.OptimizationProblem(objective)  #create an OpenTURNS OptimizationProblem
        problem.setMinimization(True)  #it is i a minimization problem
        problem.setInequalityConstraint(inequality_constraint) #set inequality constraints (both probabilistic and deterministic)
        problem.setBounds(bounds)  #set bounds
        
        algo.setProblem(problem)
        startingPoint = self.InitialPoint  #define the starting point of optimization problem
        algo.setStartingPoint(startingPoint) #set it into the algorithm
        algo.run() 
        
        self.result = algo.getResult()
        
        return self.result
    
    def get_optimum(self):
        ## Provide the optimal design point
        return self.result.getOptimalPoint()
        
    def get_foptimum(self):
        ## Provide the value of the objective function at the optimal design point
        return self.RBDO_Problem.Objective(self.result.getOptimalPoint())     
    
    def get_consoptimum(self):
        ## Provide the constraint values (Probabilistic and Deterministic)  at the optimal design point
        
        # Evaluation of the constraints (Probabilistic and Deterministic) on the optimal RIA design point 
        cons = -self.RBDO_PIneq_eval_FORM(self.result.getOptimalPoint()) ## the minus is to provide feasibility as standard RBDO problem (compared to openTURNS optimization problem)
        
        #Number of (Probabilistic and Deterministic) constraints
        NbPIneqCons = self.RBDO_Problem.NbPIneqCons
        NbDIneqCons = self.RBDO_Problem.NbDIneqCons
        
        #Value of the reliability index target for each Probabilistic inequality constraint
        beta_target_list = []
        for i in range(NbPIneqCons):
            beta_target_list.append(self.RBDO_Problem.PIneqCons_list[i].beta_target)
            
        #Computation of the probability of failure at the optimal design point for each Probabilistic inequality constraint
        dist_norm = ot.Normal()
        Pf_computed = []
        for i in range(NbPIneqCons):
            Pf_computed.append(dist_norm.computeCDF(cons[i]-beta_target_list[i]))
        
        #Computation of the Deterministic inequality constraints
        List_value_DIneqCons = []
        for i in range(NbDIneqCons):
            List_value_DIneqCons.append(-self.RBDO_Problem.DIneqCons_list[i](self.result.getOptimalPoint())[0])
        
        return   Pf_computed+List_value_DIneqCons

## Class of PMA algorithm (Probability Measure Approach), takes as input:
# RBDO_Problem
# Solver
# InitialPoint
# Runs an optimization problem, and for each Probabilistic Inequality Constraint it runs an inverse-FORM method to estimate the probability of failure (and the reliability index)
# Careful, openTURNS optimization problem assumes inequality constraints >=0 as feasible   (opposite of classical RBDO problem which assumes constraint feasibility as <=0)

class PMA_Algorithm(object):

    def __init__(self,RBDO_Problem,Solver,InitialPoint):
        self.RBDO_Problem = RBDO_Problem
        self.RBDO_PIneq_eval_invFORM = ot.Function(RBDO_PIneq_eval_invFORM(self.RBDO_Problem))  #transformation into OpenTURNS function
        self.Solver = Solver
        self.InitialPoint = InitialPoint
    
        
    def run(self):
        #idem RIA for the comments
        algo = self.Solver 
        objective = self.RBDO_Problem.Objective   #objective function
        inequality_constraint = self.RBDO_PIneq_eval_invFORM  #constraint functions
        bounds = self.RBDO_Problem.Bounds  #research optimization bounds
        problem = ot.OptimizationProblem(objective)  #create an OpenTURNS OptimizationProblem
        problem.setMinimization(True)  #it is i a minimization problem
        problem.setInequalityConstraint(inequality_constraint)  #set inequality constraints (both probabilistic and deterministic)
        problem.setBounds(bounds) #set bounds
        
        algo.setProblem(problem)
        # define starting point for optimization algorithm
        startingPoint = self.InitialPoint
        algo.setStartingPoint(startingPoint)
        algo.run()
        
        self.result = algo.getResult()
        
        return self.result
    
    def get_optimum(self):
        ## Provide the optimal design point        
        return self.result.getOptimalPoint()
        
    def get_foptimum(self):
        ## Provide the value of the objective function at the optimal design point        
        return self.RBDO_Problem.Objective(self.result.getOptimalPoint())     
    
    def get_consoptimum(self):
        ## Provide the constraint values (Probabilistic and Deterministic)  at the optimal design point
        
        # Evaluation of the constraints (Probabilistic and Deterministic) on the optimal PMA design point using invFORM        
        cons = -self.RBDO_PIneq_eval_invFORM(self.result.getOptimalPoint())
        NbPIneqCons = self.RBDO_Problem.NbPIneqCons
        NbDIneqCons = self.RBDO_Problem.NbDIneqCons

        #Value of the reliability index target for each Probabilistic inequality constraint
        beta_target_list = []
        for i in range(NbPIneqCons):
            beta_target_list.append(self.RBDO_Problem.PIneqCons_list[i].beta_target)
        
        #Computation of the probability of failure at the optimal design point for each Probabilistic inequality constraint
        dist_norm = ot.Normal()
        Pf_computed = []
        for i in range(NbPIneqCons):
            Pf_computed.append(dist_norm.computeCDF(cons[i]-beta_target_list[i]))
        
        #Computation of the Deterministic inequality constraints        
        List_value_DIneqCons = []
        for i in range(NbDIneqCons):
            List_value_DIneqCons.append(-self.RBDO_Problem.DIneqCons_list[i](self.result.getOptimalPoint())[0])

        return   Pf_computed+List_value_DIneqCons
        
## Class of SORA algorithm (Sequential Optimization and Reliability Analysis), takes as input:
# RBDO_Problem
# Solver
# InitialPoint
# It consists of a sequential loop of a deterministic optimization problem (with the uncertain variables fixed at the Most Probable Target Point), and the evaluation of the Probabilistic Inequality Constraints using inverse FORM for the design variables fixed at the optimal value of the deterministic optimization problem
# Careful, openTURNS optimization problem assumes inequality constraints >=0 as feasible   (opposite of classical RBDO problem which assumes constraint feasibility as <=0)
class SORA_Algorithm(object):

    def __init__(self,RBDO_Problem,Solver,InitialPoint):
        self.RBDO_Problem = RBDO_Problem
        self.Solver = Solver
        self.InitialPoint = InitialPoint
        self.RBDO_PIneq_eval_SORA = ot.Function(RBDO_PIneq_eval_SORA(self.RBDO_Problem))  #transformation into OpenTURNS function
    
    def MPTP_determination_list(self,XX):
        #Give the list of the Most Probable Target Point for each Probabilistic inequality constraint
        List_MPTP = []
        for i in range(self.RBDO_Problem.NbPIneqCons):
            List_MPTP.append(self.RBDO_Problem.PIneqCons_list[i].Eval_PIneqCons_invFORM(XX,self.RBDO_Problem.PIneqCons_list[i].solver_invFORM))
                            
        return List_MPTP
    
        
    def run(self):
        #Initialization of the sequential loop, definition of the shiftingVector and the MPTP values
        shiftingVector_List =[]
        z_MPTP_List = []
        for i in range(len(self.RBDO_Problem.PIneqCons_list)):
            self.RBDO_Problem.PIneqCons_list[i].Active_d_p(self.InitialPoint)
            self.RBDO_Problem.PIneqCons_list[i].len_X_Z(self.InitialPoint)
            len_X = self.RBDO_Problem.PIneqCons_list[i].len_X                
                
            shiftingVector_List.append([0.]*len_X)
            if self.RBDO_Problem.PIneqCons_list[i].dist_Z==None:
                z_MPTP_List.append(None)
            else:                
                z_MPTP_List.append(self.RBDO_Problem.PIneqCons_list[i].dist_Z.getMean())
        
        #Convergence criterion for SORA loop
        crit_eps_th = 0.001
        max_iter = 40
        crit_eps = 10.
        iter_ = 0
        optim_XX_previous = ot.Point(self.InitialPoint)

        while (crit_eps>=crit_eps_th) and (iter_<max_iter):           
            #Update of the shifting vector and the MPTP for Z variable
            for i in range(len(self.RBDO_Problem.PIneqCons_list)):
                self.RBDO_Problem.PIneqCons_list[i].shiftingVector = shiftingVector_List[i]
                self.RBDO_Problem.PIneqCons_list[i].z_MPTP =  z_MPTP_List[i]
            
            #Deterministic Optimization problem solving with X and Z fixed to the MPTP of the previous iteration and with the appropriate shifhting vector
            algo = ot.NLopt('LN_COBYLA')
            objective = self.RBDO_Problem.Objective    #objective function
            inequality_constraint = self.RBDO_PIneq_eval_SORA    #inequality constraints function (both probabilistic and deterministic)
            bounds = self.RBDO_Problem.Bounds 
            problem = ot.OptimizationProblem(objective)
            problem.setMinimization(True)
            problem.setInequalityConstraint(inequality_constraint)
            problem.setBounds(bounds)
            
            algo.setProblem(problem)
            startingPoint = self.InitialPoint
            algo.setStartingPoint(startingPoint)
            algo.run()
                
            self.result = algo.getResult()
            
            #Current deterministic optimum
            optimum_XX = self.result.getOptimalPoint()
            
            #Evaluation of the probabilistic inequality constraint with the design variable fixed at the current optimum
            #Calculation of the new MPTP and shifhting vector
            MPTP_xz_list = []
            z_MPTP_List = []
            shiftingVector_List = []
            for i in range(len(self.RBDO_Problem.PIneqCons_list)):
                MPTP_xz_courant = self.RBDO_Problem.PIneqCons_list[i].Eval_PIneqCons_invFORM(optimum_XX)[1]
                MPTP_xz_list.append(list(MPTP_xz_courant))
                Current_Index_Active_p = self.RBDO_Problem.PIneqCons_list[i].Index_Active_p
                self.RBDO_Problem.PIneqCons_list[i].len_X_Z(optimum_XX)
                len_X = self.RBDO_Problem.PIneqCons_list[i].len_X
                                
                x_MPTP = MPTP_xz_courant[:len_X]                
                z_MPTP = MPTP_xz_courant[len_X:]
                if z_MPTP.getDimension()==0:
                    z_MPTP_List.append(None)
                else:
                    z_MPTP_List.append(z_MPTP)
                shiftingVector_List.append(optimum_XX[Current_Index_Active_p]-x_MPTP)
        
            #Update of the stopping criterion
            crit_eps_ = optim_XX_previous-optimum_XX
            crit_eps = crit_eps_.norm()
            iter_ = iter_ +1
            optim_XX_previous = optimum_XX
        
        return self.result
    
    def get_optimum(self):
        ## Provide the optimal design point                
        return self.result.getOptimalPoint()
        
    def get_foptimum(self):
        ## Provide the value of the objective function at the optimal design point        
        return self.RBDO_Problem.Objective(self.result.getOptimalPoint())     
    
    def get_consoptimum(self):
        ## Provide the constraint values (Probabilistic and Deterministic)  at the optimal design point
        cons = [-a for a in self.RBDO_Problem.evalConsList_SORA(self.result.getOptimalPoint()) ]
        
        NbPIneqCons = self.RBDO_Problem.NbPIneqCons
        NbDIneqCons = self.RBDO_Problem.NbDIneqCons
        
        beta_target_list = []
        for i in range(NbPIneqCons):
            beta_target_list.append(self.RBDO_Problem.PIneqCons_list[i].beta_target)
        
        dist_norm = ot.Normal()
        Pf_computed = []
        for i in range(NbPIneqCons):
            Pf_computed.append(dist_norm.computeCDF(cons[i]-beta_target_list[i]))
        
        List_value_DIneqCons = []
        for i in range(NbDIneqCons):
            List_value_DIneqCons.append(cons[NbPIneqCons+i])
            
        return Pf_computed+List_value_DIneqCons
    
#Class to define the RBDO problem. It takes as input:
# - the objective function
# - the list of the probablistic inequality constraints PIneqCons_list
# - the list of the deterministic inequality constraints DIneqCons_list_
# - the bounds on the design variables Bounds
# - the lenght of the vector d: len_d
# - the lenght of the vecto p: len_p
class RBDO_Problem(object):

    def __init__(self,Objective,PIneqCons_list,DIneqCons_list_,Bounds,len_d,len_p):
        self.Objective = ot.Function(Objective)
        self.PIneqCons_list = PIneqCons_list
        self.DIneqCons_list = [ot.Function(DIneqCons_list_[i]) for i in range(len(DIneqCons_list_))]
        self.Bounds = Bounds
        self.len_d = len_d
        self.len_p = len_p
        self.NbPIneqCons = len(self.PIneqCons_list)
        self.NbDIneqCons = len(self.DIneqCons_list)
    
    #Evaluate the list of the constraints using FORM for the Probabilistic constraints
    def evalConsList_FORM(self,XX):
        
        List_value_PIneqCons = []
        for i in range(self.NbPIneqCons):
            List_value_PIneqCons.append(self.PIneqCons_list[i].Eval_PIneqCons_FORM(XX)[0])
        
        List_value_DIneqCons = []
        for i in range(self.NbDIneqCons):
            List_value_DIneqCons.append(self.DIneqCons_list[i](XX)[0])
            
        List_value_RBDO_cons = List_value_PIneqCons+List_value_DIneqCons
        
        return List_value_RBDO_cons

    #Evaluate the list of the constraints using inverse FORM for the Probabilistic constraints    
    def evalConsList_invFORM(self,XX):
        
        List_value_PIneqCons = []
        for i in range(self.NbPIneqCons):
            List_value_PIneqCons.append(self.PIneqCons_list[i].Eval_PIneqCons_invFORM(XX)[0])
        
        List_value_DIneqCons = []
        for i in range(self.NbDIneqCons):
            List_value_DIneqCons.append(self.DIneqCons_list[i](XX)[0])
            
        List_value_RBDO_cons = List_value_PIneqCons+List_value_DIneqCons
        
        return List_value_RBDO_cons

    #Evaluate the list of the constraints using inverse FORM for the Probabilistic constraints with the current design vector fixed
    def evalConsList_SORA(self,XX):
        
        List_value_PIneqCons = []
        for i in range(self.NbPIneqCons):
            
            if  self.PIneqCons_list[i].type_failure(1,2)==True:
                List_value_PIneqCons.append(self.PIneqCons_list[i].Eval_PIneqCons_SORA(XX)[0])
            else : 
                List_value_PIneqCons.append(-self.PIneqCons_list[i].Eval_PIneqCons_SORA(XX)[0])
        
        List_value_DIneqCons = []
        for i in range(self.NbDIneqCons):
            List_value_DIneqCons.append(self.DIneqCons_list[i](XX)[0])
            
        List_value_RBDO_cons = List_value_PIneqCons+List_value_DIneqCons
        
        return List_value_RBDO_cons

# Class to evaluate the Python constraint provided by the user and to convert to the formalism OpenTURNSPythonFunction 
class g_rbdopythonfunction(ot.OpenTURNSPythonFunction):

     def __init__(self,g_rbdo,Index_Active_d,XX,len_X,len_Z):
        
         self.g_rbdo = g_rbdo
         self.Index_Active_d = Index_Active_d
         self.XX = XX
         self.Active_d_value = [XX[i] for i in self.Index_Active_d]
         self.len_X = len_X
         self.len_Z = len_Z
         
         super(g_rbdopythonfunction, self).__init__(self.len_X+self.len_Z, 1)
                  
                
     def  _exec(self,W):

         Z = [W[i] for i in range(-self.len_Z,0)]
         X = [W[i] for i in range(0,len(W)-self.len_Z)]        
        
         Y = self.g_rbdo(self.Active_d_value,X,Z)

         return ot.Point(Y)
    
# Class to evaluate the Python constraint with SORA provided by the user and to convert to the formalism OpenTURNSPythonFunction with the correct shiftingVector 
class g_SORApythonfunction(ot.OpenTURNSPythonFunction):

     def __init__(self,g_rbdo,Index_Active_d,Index_Active_p,XX,shiftingVector,Z_MPTP):
        
         self.g_rbdo = g_rbdo
         self.Index_Active_d = Index_Active_d
         self.Index_Active_p = Index_Active_p
         self.shiftingVector = shiftingVector
         self.Z_MPTP = Z_MPTP
         
         super(g_SORApythonfunction, self).__init__(len(XX), 1)
                  
                
     def  _exec(self,XX):
         shiftingVector = self.shiftingVector
         Z_MPTP = self.Z_MPTP
         d = [XX[i] for i in self.Index_Active_d]
         p = [XX[i] for i in self.Index_Active_p]
         
         X = [p[i]-shiftingVector[i] for i in range(len(p))]
        
         Y = self.g_rbdo(d,X,Z_MPTP)

         return ot.Point(Y)

##Class to define the Probabilistic Constraints and to evaluate it using FORM or invFORM.
#It takes as input:
# - g: the constraint Python function g
# - Index_Active_d: the list of the index in the design vector XX=[d,p] for the active deterministic design variables d in g
# - Index_Active_p: the list of the index in the design vector XX=[d,p] for the active parameters p in g 
# - update_distX: Python function that takes as input p and provide a distribution for the random vector X with the updated hyperparameters p values
# - distZ: distribution of the random vector Z
# - Pf_target: the target maximal probablity for the probabilistic constraint
# - type_failure ot.Less() or ot.Greater()  depending if P[g<=0] or P[g>=0], by default it is ot.Less() meaning P[g<=0]
# - solver_invFORM: solver for the inverse FORM technique, only necessary for PMA and SORA algorithms. Choice of Inverse FORM Solver : ['SLSQP','AMV','CMV','HMV']
class PIneqCons(object):

    def __init__(self,g,Index_Active_d,Index_Active_p,update_distX,distZ,Pf_target,type_failure= ot.Less(),solver_invFORM = None,shiftingVector = None,z_MPTP = None):
        
        self.Pf_target = Pf_target
        self.update_distX = update_distX
        if distZ==[]:
            self.dist_Z = None
        else:
            self.dist_Z = distZ
        self.g = g
        self.dist_X = None
        self.Index_Active_d = Index_Active_d
        self.Index_Active_p = Index_Active_p        
        self.Active_d_value = None
        self.Active_p_value = None
        self.beta_target = list(-ot.Normal().computeQuantile(self.Pf_target))[0]
        self.solver_invFORM = solver_invFORM 
        self.shiftingVector = shiftingVector 
        self.z_MPTP = z_MPTP 
        self.len_X = None
        self.len_Z = None
        self.type_failure = type_failure

    #Update of the active d and p values for the constraint g considering the vector XX=[d,p]    
    def Active_d_p(self,XX):       
         ## define active d and p using  active_d and active_p  
         self.Active_d_value = [XX[i] for i in self.Index_Active_d]
         self.Active_p_value = [XX[i] for i in self.Index_Active_p]
        
    ## define the new distX using the function provided by user in update_distX                                     
    def Update_distX(self,p):
        self.dist_X = self.update_distX(p)

    ## define the len of the vector X and the vector Z
    def len_X_Z(self,XX):
        self.Active_d_p(XX)
        if self.update_distX==None:
            self.len_X = 0
        else:       
            self.Update_distX(self.Active_p_value)
            self.len_X = self.dist_X.getDimension()
        
        if self.dist_Z ==None:
            self.len_Z =0         
        else:
            self.len_Z = self.dist_Z.getDimension()
        
    ## evaluate the constraint using FORM and return -[-Beta_i + Beta_i_Target] to cope with openTURNS g()>=0 type of constraints
    ## it also returns the PhysicalSpaceDesignPoint based on FORM results
    def Eval_PIneqCons_FORM(self,XX):
       
        self.Active_d_p(XX)
        self.len_X_Z(XX)
        
        #Determination of the joint X,Z distribution for g
        if self.update_distX!=None:
            self.Update_distX(self.Active_p_value)
        
        if self.dist_Z ==None and self.update_distX!=None:
            jointdistribution_X_Z = self.dist_X 
        elif self.dist_Z !=None and self.update_distX!=None:
            jointdistribution_X_Z = ot.BlockIndependentDistribution([self.dist_X,self.dist_Z])
        elif self.dist_Z !=None  and self.update_distX==None :
            jointdistribution_X_Z = self.dist_Z
        else :
            raise ValueError('both dist X and dist Z are None, this is not a probabilistic constraint')
            
        
        #definition of the g_rbdofunction by fixing the design variables d
        g_rbdofunction = ot.Function(g_rbdopythonfunction(self.g,self.Index_Active_d,XX,\
                                              self.len_X,self.len_Z))
        # FORM    
        vect = ot.RandomVector(jointdistribution_X_Z)
        output = ot.CompositeRandomVector(g_rbdofunction, vect)
        event = ot.ThresholdEvent(output, self.type_failure, 0.0)
    
        solver = ot.AbdoRackwitz()
        algo = ot.FORM(solver, event, [1.]*jointdistribution_X_Z.getDimension())
        algo.run()
        result = algo.getResult()
        
        #Becareful if probability of failure is >0.5, openTURNS returns -beta
        if result.getEventProbability()<=0.5:
            y = result.getHasoferReliabilityIndex()
        else:
            y = -result.getHasoferReliabilityIndex()
    
        return [-(-y+self.beta_target), result.getPhysicalSpaceDesignPoint()]

    def Eval_PIneqCons_SORA(self,XX):
        ## evaluate the constraint using with shiftingVector and Z fixed at MPTP     
        g_SORAfunction = ot.Function(g_SORApythonfunction(self.g,self.Index_Active_d,\
                                                          self.Index_Active_p,XX,self.shiftingVector,self.z_MPTP))
        Y = g_SORAfunction(XX)
        
        return Y


    ## evaluate the constraint using FORM and return -Gi_minus to cope with openTURNS g()>=0 type of constraints
    ## it also returns the PhysicalSpaceDesignPoint based on FORM results
    def Eval_PIneqCons_invFORM(self,XX):
        self.Active_d_p(XX)
        
        # 4 types of algorithm to solve the inverse FORM --> optimization problem
        if self.solver_invFORM == 'SLSQP':
            return self.Eval_PIneqCons_invFORM_SLSQP(XX)    
        elif self.solver_invFORM == 'AMV':
            return self.Eval_PIneqCons_invFORM_AMV(XX)
        elif self.solver_invFORM == 'CMV':
            return self.Eval_PIneqCons_invFORM_CMV(XX)
        
        elif self.solver_invFORM == 'HMV':
            return self.Eval_PIneqCons_invFORM_HMV(XX)
                    
        else : raise ValueError('Error in invFORM solver, please provide an appropriate solver')
        
    #Evaluation of inverseFORM using SLSQP solver from NLOpt
    def Eval_PIneqCons_invFORM_SLSQP(self,XX):
        #Definition of the joint distribution of X, Z for the constraint
        #Determination of the joint X,Z distribution for g
        self.Active_d_p(XX)
        self.len_X_Z(XX)
        
        if self.update_distX!=None:
            self.Update_distX(self.Active_p_value)
        
        if self.dist_Z ==None and self.update_distX!=None:
            jointdistribution_X_Z = self.dist_X 
        elif self.dist_Z !=None and self.update_distX!=None:
            jointdistribution_X_Z = ot.BlockIndependentDistribution([self.dist_X,self.dist_Z])
        elif self.dist_Z !=None  and self.update_distX==None :
            jointdistribution_X_Z = self.dist_Z
        else :
            raise ValueError('both dist X and dist Z are None, this is not a probabilistic constraint')

        
        transformU_to_X = jointdistribution_X_Z.getInverseIsoProbabilisticTransformation() 
        
        
        g_rbdofunction = ot.Function(g_rbdopythonfunction(self.g,self.Index_Active_d,XX,\
                                              self.len_X,self.len_Z)) 
            
            
        if self.type_failure(1,2) == False:
            obj = lambda s: list(-g_rbdofunction(list(transformU_to_X(s))))
        else:
            obj = lambda s: list(g_rbdofunction(list(transformU_to_X(s))))
            
        cons = lambda s: [ot.Point(s).norm()-self.beta_target]
        
        dim = jointdistribution_X_Z.getDimension()   

        #Solving of the optimization problem corresponding to inverseFORM using SLSQP from NLopt
        algo = ot.NLopt('LD_SLSQP')
        objective = ot.PythonFunction(dim, 1, obj)
        equality_constraint = ot.PythonFunction(dim, 1, cons)
        bounds = ot.Interval([-7.]*dim, [7.]*dim)
        problem = ot.OptimizationProblem(objective)
        problem.setMinimization(True)
        problem.setEqualityConstraint(equality_constraint)
        problem.setBounds(bounds)
        algo.setProblem(problem)
        startingPoint = [2.0] * dim
        algo.setStartingPoint(startingPoint)
        algo.run()
        result = algo.getResult()
        iter_ = result.getEvaluationNumber()
        point_best = result.getOptimalPoint()
        
        
        Gi_minus = list(-objective(point_best))[0]
        
        return [-Gi_minus, transformU_to_X(point_best)]

    #Evaluation of inverseFORM using AMV algorithm
    def Eval_PIneqCons_invFORM_AMV(self,XX):
        #Determination of the joint X,Z distribution for g
        self.Active_d_p(XX)
        self.len_X_Z(XX)
        
        if self.update_distX!=None:
            self.Update_distX(self.Active_p_value)
        
        if self.dist_Z ==None and self.update_distX!=None:
            jointdistribution_X_Z = self.dist_X 
        elif self.dist_Z !=None and self.update_distX!=None:
            jointdistribution_X_Z = ot.BlockIndependentDistribution([self.dist_X,self.dist_Z])
        elif self.dist_Z !=None  and self.update_distX==None :
            jointdistribution_X_Z = self.dist_Z
        else :
            raise ValueError('both dist X and dist Z are None, this is not a probabilistic constraint')

        

        
        transformU_to_X = jointdistribution_X_Z.getInverseIsoProbabilisticTransformation() 
        
        
        g_rbdofunction = ot.Function(g_rbdopythonfunction(self.g,self.Index_Active_d,XX,\
                                              self.len_X,self.len_Z)) 
            
            
        if self.type_failure(1,2) == False:
            obj = lambda s: list(-g_rbdofunction(list(transformU_to_X(s))))
        else:
            obj = lambda s: list(g_rbdofunction(list(transformU_to_X(s))))
            
        cons = lambda s: [ot.Point(s).norm()-self.beta_target]
        
        dim = jointdistribution_X_Z.getDimension()        

        objective = ot.PythonFunction(dim, 1, obj)
        point_ini = ot.Point([0.0] * dim)
        point_old = point_ini
        crit_eps_u_th = 0.001
        crit_eps_u = 1.0
        iter_ = 0
        max_iter = 30
        while (crit_eps_u>=crit_eps_u_th) and (iter_<max_iter):
            
            grad = - objective.gradient(self.beta_target*point_old)
            grad_vec = ot.Point([grad[i,0] for i in range(grad.getNbRows())])
            grad_norm = grad_vec.norm()
            n = grad_vec/grad_norm
            point_new = self.beta_target*n
            
            eps_u = point_old-point_new
            crit_eps_u = eps_u.norm()
            point_old = point_new
            iter_ = iter_ +1 
    
            point_best = point_old
            
        Gi_minus = list(-objective(point_best))[0] 
        
        return [-Gi_minus, point_best]            
            
            
    #Evaluation of inverseFORM using CMV algorithm                     
    def Eval_PIneqCons_invFORM_CMV(self,XX):
        #Determination of the joint X,Z distribution for g
        self.Active_d_p(XX)
        self.len_X_Z(XX)
        
        if self.update_distX!=None:
            self.Update_distX(self.Active_p_value)
        
        if self.dist_Z ==None and self.update_distX!=None:
            jointdistribution_X_Z = self.dist_X 
        elif self.dist_Z !=None and self.update_distX!=None:
            jointdistribution_X_Z = ot.BlockIndependentDistribution([self.dist_X,self.dist_Z])
        elif self.dist_Z !=None  and self.update_distX==None :
            jointdistribution_X_Z = self.dist_Z
        else :
            raise ValueError('both dist X and dist Z are None, this is not a probabilistic constraint')

        

        
        transformU_to_X = jointdistribution_X_Z.getInverseIsoProbabilisticTransformation() 
        
        
        g_rbdofunction = ot.Function(g_rbdopythonfunction(self.g,self.Index_Active_d,XX,\
                                              self.len_X,self.len_Z)) 
            
            
            
        if self.type_failure(1,2) == False:
            obj = lambda s: list(-g_rbdofunction(list(transformU_to_X(s))))
        else:
            obj = lambda s: list(g_rbdofunction(list(transformU_to_X(s))))
        cons = lambda s: [ot.Point(s).norm()-self.beta_target]
        
        dim = jointdistribution_X_Z.getDimension()                 
                
        objective = ot.PythonFunction(dim, 1, obj)
        point_ini = ot.Point([0.0] * dim)
        point_0 = point_ini
        
        grad = - objective.gradient(self.beta_target*point_0)
        grad_vec = ot.Point([grad[i,0] for i in range(grad.getNbRows())])
        grad_norm = grad_vec.norm()
        n = grad_vec/grad_norm
        grad_0 = grad_vec
        point_1 = self.beta_target*n
        
        grad = - objective.gradient(self.beta_target*point_1)
        grad_vec = ot.Point([grad[i,0] for i in range(grad.getNbRows())])
        grad_norm = grad_vec.norm()
        n = grad_vec/grad_norm
        grad_1 = grad_vec
        point_2 = self.beta_target*n
        
        grad = - objective.gradient(self.beta_target*point_2)
        grad_vec = ot.Point([grad[i,0] for i in range(grad.getNbRows())])
        grad_2 = grad_vec
        
        point_list = [point_2,point_1,point_0]
        grad_list = [grad_2,grad_1,grad_0]
        
        crit_eps_u_th = 0.001
        crit_eps_u = 1.0
        iter_ = 0
        max_iter = 30    
        
        while (crit_eps_u>=crit_eps_u_th) and (iter_<max_iter):
            
            num = grad_list[0]/grad_list[0].norm()+grad_list[1]/grad_list[1].norm()+grad_list[2]/grad_list[2].norm()
            denom =  num.norm()
            
            point_new =self.beta_target*num/denom
            
            grad = - objective.gradient(point_new)
            grad_vec = ot.Point([grad[i,0] for i in range(grad.getNbRows())])
            
            grad_list = [grad_vec] + grad_list
            grad_list.pop()

            point_list = [point_new] + point_list
            point_list.pop()
            
            eps_u = point_list[0]-point_list[1]
            crit_eps_u = eps_u.norm()
            iter_ = iter_ +1 
    
            point_best = point_list[0]
            
        Gi_minus = list(-objective(point_best))[0]
        
        return [-Gi_minus, point_best] 

    #Evaluation of inverseFORM using HMV algorithm
    def Eval_PIneqCons_invFORM_HMV(self,XX):
        #Determination of the joint X,Z distribution for g
        self.Active_d_p(XX)
        self.len_X_Z(XX)
        
        if self.update_distX!=None:
            self.Update_distX(self.Active_p_value)
        
        if self.dist_Z ==None and self.update_distX!=None:
            jointdistribution_X_Z = self.dist_X 
        elif self.dist_Z !=None and self.update_distX!=None:
            jointdistribution_X_Z = ot.BlockIndependentDistribution([self.dist_X,self.dist_Z])
        elif self.dist_Z !=None  and self.update_distX==None :
            jointdistribution_X_Z = self.dist_Z
        else :
            raise ValueError('both dist X and dist Z are None, this is not a probabilistic constraint')

        

        
        transformU_to_X = jointdistribution_X_Z.getInverseIsoProbabilisticTransformation() 
        
        
        g_rbdofunction = ot.Function(g_rbdopythonfunction(self.g,self.Index_Active_d,XX,\
                                              self.len_X,self.len_Z)) 
            
            
            
        if self.type_failure(1,2) == False:
            obj = lambda s: list(-g_rbdofunction(list(transformU_to_X(s))))
        else:
            obj = lambda s: list(g_rbdofunction(list(transformU_to_X(s))))
        cons = lambda s: [ot.Point(s).norm()-self.beta_target]
        
        dim = jointdistribution_X_Z.getDimension()                 

                 
        objective = ot.PythonFunction(dim, 1, obj)
        point_ini = ot.Point([0.0] * dim)
        point_0 = point_ini
        
        grad = - objective.gradient(self.beta_target*point_0)
        grad_vec = ot.Point([grad[i,0] for i in range(grad.getNbRows())])
        grad_norm = grad_vec.norm()
        n = grad_vec/grad_norm
        grad_0 = grad_vec
        point_1 = self.beta_target*n
        
        grad = - objective.gradient(self.beta_target*point_1)
        grad_vec = ot.Point([grad[i,0] for i in range(grad.getNbRows())])
        grad_norm = grad_vec.norm()
        n = grad_vec/grad_norm
        grad_1 = grad_vec
        point_2 = self.beta_target*n
        
        grad = - objective.gradient(self.beta_target*point_2)
        grad_vec = ot.Point([grad[i,0] for i in range(grad.getNbRows())])
        grad_2 = grad_vec
        
        point_list = [point_2,point_1,point_0]
        grad_list = [grad_2,grad_1,grad_0]
        
        crit_eps_u_th = 0.001
        crit_eps_u = 1.0
        iter_ = 0
        max_iter = 30
        
        term1 = grad_list[0] - grad_list[1]
        term2 = grad_list[1] - grad_list[2]
        eta = term1.dot(term2)
        
        while (crit_eps_u>=crit_eps_u_th) and (iter_<max_iter):
            
            if eta>0:
                grad = - objective.gradient(self.beta_target*point_list[0])
                grad_vec = ot.Point([grad[i,0] for i in range(grad.getNbRows())])
                grad_norm = grad_vec.norm()
                n = grad_vec/grad_norm
                point_new =self.beta_target*n
                
                grad_list = [grad_vec] + grad_list
                grad_list.pop()

                point_list = [point_new] + point_list
                point_list.pop()
                
                term1 = grad_list[0] - grad_list[1]
                term2 = grad_list[1] - grad_list[2]
                eta = term1.dot(term2)
                
                eps_u = point_list[0]-point_list[1]
                crit_eps_u = eps_u.norm()
                iter_ = iter_ +1  
                point_best = point_list[0]
                
            if eta<=0:
                
                num = grad_list[0]/grad_list[0].norm()+grad_list[1]/grad_list[1].norm()+grad_list[2]/grad_list[2].norm()
                denom =  num.norm()
                
                point_new = self.beta_target*num/denom
                
                grad = - objective.gradient(point_new)
                grad_vec = ot.Point([grad[i,0] for i in range(grad.getNbRows())])
                
                grad_list = [grad_vec] + grad_list
                grad_list.pop()

                point_list = [point_new] + point_list
                point_list.pop()
                
                eps_u = point_list[0]-point_list[1]
                crit_eps_u = eps_u.norm()
                iter_ = iter_ +1 
        
                point_best = point_list[0]  
                    
        Gi_minus = list(-objective(point_best))[0]
        
        return [-Gi_minus, point_best]
        
# Class to define the deterministic constraint and to use the OpenTURNSPythonFunction formalism
# It takes as inputs:
# - g: the Python constraint
# - Index_Active_d: the list of the index in the design vector XX=[d,p] for the active deterministic design variables d in g
# - Index_Active_p: the list of the index in the design vector XX=[d,p] for the active parameters p in g 
# - len_XX: lenght of the vector XX
class DIneqCons(ot.OpenTURNSPythonFunction):

    def __init__(self,g,Index_Active_d,Index_Active_p,len_XX):
        
        self.DetermIneqCons = g
        self.Index_Active_d = Index_Active_d
        self.Index_Active_p = Index_Active_p 
        self.len_XX = len_XX
        
        super(DIneqCons, self).__init__(self.len_XX, 1)
                
            
    def Active_d_p(self,XX):       
        ## define active d and p using  active_d and active_p
        self.Active_d_value = [XX[i] for i in self.Index_Active_d]
        self.Active_p_value = [XX[i] for i in self.Index_Active_p]
         
        return self.Active_d_value,self.Active_p_value 
        
    def  _exec(self,XX):
        d_active, p_active = self.Active_d_p(XX)     
        
        Y = self.DetermIneqCons(d_active,p_active)
        
        return [-Y]
        
# Class to define the deterministic objective function and to use the OpenTURNSPythonFunction formalism
# It takes as inputs:        
# - obj: the Python objective function
# - Index_Active_d: the list of the index in the design vector XX=[d,p] for the active deterministic design variables d in g
# - Index_Active_p: the list of the index in the design vector XX=[d,p] for the active parameters p in g         
# - len_XX: lenght of the vector XX
class Objective(ot.OpenTURNSPythonFunction):

    def __init__(self,obj,Index_Active_d,Index_Active_p,len_XX):
        
        self.objective = obj
        self.Index_Active_d = Index_Active_d
        self.Index_Active_p = Index_Active_p
        self.len_XX   = len_XX      
        
        super(Objective, self).__init__(self.len_XX, 1)
                  
        
    def Active_d_p(self,XX):       
        ## define active d and p using  active_d and active_p
        self.Active_d_value = [XX[i] for i in self.Index_Active_d]
        self.Active_p_value = [XX[i] for i in self.Index_Active_p]
         
        return self.Active_d_value,self.Active_p_value
                
    def  _exec(self,XX):

        d_active, p_active = self.Active_d_p(XX)     
        
        Y = self.objective(d_active,p_active)

        return [Y]