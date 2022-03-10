import numpy as np
import qutip as qp
import scipy

class QuIPI_circuit():        
    def set_qumode_para(self, s, cut, max_cut):
        self.s = s # squeezing factor
        self.cut = cut # truction of Fock state for preparing the resource state
        self.max_cut = max_cut # Fock state space
    
    
    # Generate resource state |R,S>
    def gen_r_state(self):
        '''
        This function generate the qumode resource state in Fock space, 
        which is a squeezed position state but only have the positive
        amplitude. This function return a quantum object.
        '''
        
        def r_state_p(s, p):
            # resource state in momentum space
            ampl = np.sqrt(1/(np.sqrt(np.pi)*s))*np.exp(-p**2/(s**2*2)) # amplitude
            if p <= 0:
                # it only has positive position
                return 0
            else:
                return ampl
        
        def single_coeff(n, p):
            # single coefficient for transforming |p> to |n>
            den = (1j)**n*scipy.special.hermite(n)(p)*np.exp(-p**2/2)
            num = np.sqrt(2**n*np.math.factorial(n)*np.sqrt(np.pi))
            cnp = den/num
            return cnp
        
        def coeff_rstate(s, n):
            # integrate all the single coefficient with amplitude
            def integrand_real(p,s,n):
                return np.real(single_coeff(n,p))*r_state_p(s,p)
            def integrand_imag(p,s,n):
                return np.imag(single_coeff(n,p))*r_state_p(s,p)
            cn_real = scipy.integrate.quad(integrand_real,0,np.infty,args=(s,n))
            cn_imag = scipy.integrate.quad(integrand_imag,0,np.infty,args=(s,n))
            cn = cn_real[0] + 1j*cn_imag[0]
            return cn
        
        def gen_fock_vec(s,cut,max_cut):
            # generate the fock state vector
            c_vec = [coeff_rstate(s,n) for n in range(cut)]+[0 for i in range(max_cut-cut)]
            return c_vec
        
        c_vec = np.array(gen_fock_vec(self.s,self.cut,self.max_cut))
        r_state_fock = qp.Qobj(c_vec) # resource state in Fock space
        return r_state_fock
    
    
    # Ideal inverse evolution
    def ideal_inverse_evolution(self, H, b, K):
        '''
        Solve the ground state by iterative applying ideal
        inverse Hamiltonian on the initial state. It returns
        a list of energy and evolved state.
        
        Args:
            H (Qobj): Hamiltonian to be solved
            b (Qobj): initial state
            K (int): iteration step
        '''
        
        state_list = [b] # store the evolved state
        E_list = [qp.expect(H,b)] # store the Energy of the evolved state
        H_inv = H.inv()
        
        for i in range(K):
            b = (H_inv*b).unit()
            state_list.append(b)
            E = qp.expect(H, b)
            E_list.append(E)
        return E_list, state_list
    
    
    # Exact ground energy and state
    def exact_ground_E_and_state(self, H):
        E, states = H.eigenstates()
        E_g = E[0]
        g_state = states[0]
        return E_g, g_state
    
    
    # Perform the algorithm without considering Trotterdecomposition
    # The Hamiltonian evolution is directly applied
    def evolution(self, H, b, K):
        '''
        Solve the ground state by our algorithm, where the Hamiltonian
        is directly applied without considering a Trotter decomposition.
        It returns a list of Energy and state corresponding to each
        iteration.
        
        Args:
            H (Qobj): Hamiltonian to be solved
            b (Qobj): initial state
            K (int): iteration step
        '''
        
        # initlization
        N = len(H.dims[0]) # number of qubits
        b_den = qp.ket2dm(b) # density matrix of b
        r_state = self.gen_r_state() # resource state in Fock space
        r_den = qp.ket2dm(r_state) # density matrix of resource state
        E_list = [qp.expect(H,b)] # Energy list
        state_list = [b_den] # state list
        
        ide = qp.tensor([qp.qeye(2)]*N) # N qubits identity operator
        # It is a vacuum state being used to project the qumode.
        proj = qp.tensor(ide,qp.ket2dm(qp.basis(self.max_cut,0)))
        # Squeezing operator; squeezing parameter = log(squeezing factor)
        sque = qp.tensor(ide,qp.squeeze(self.max_cut,np.log(self.s)))
        P = qp.operators.momentum(self.max_cut) # momentum operator
        eHP = (-1j*qp.tensor(H, P)).expm() # exp(-i H P) to entangle qubit and qumode
        
        for i in range(K):
            b_r = qp.tensor(b_den, r_den)
            b_r = eHP*b_r*eHP.dag() # applying the unitary operator
            # project the qumode to the squeezed position state
            b_r = proj*sque.dag()*b_r*sque*proj
            # Partial trace the ancillary qumode and then normalization
            b_den = (b_r.ptrace([i for i in range(N)])).unit()
            state_list.append(b_den)
            E_list.append(qp.expect(H,b_den))
        
        return E_list, state_list
    
    
    # Evolution with First order Trotter decomposition
    def Trotter_evolution(self, H, h_list, b, K, S):
        '''
        Solve the ground state by our algorithm, where the Hamiltonian
        the first order Trotter decomposition is considered.
        It returns a list of Energy and state corresponding to each
        iteration.
        
        Args:
            H (Qobj): Hamiltonian to be solved
            h_list (list of Qobj): local operator list
            b (Qobj): initial state
            K (int): iteration step
            S (int): Trotter step
        '''
        
        # initlization
        N = len(H.dims[0]) # number of qubits
        b_den = qp.ket2dm(b) # density matrix of b
        r_state = self.gen_r_state() # resource state in Fock space
        r_den = qp.ket2dm(r_state) # density matrix of resource state
        E_list = [qp.expect(H,b)] # Energy list
        state_list = [b_den] # state list
        
        ide = qp.tensor([qp.qeye(2)]*N) # N qubits identity operator
        # It is a vacuum state being used to project the qumode.
        proj = qp.tensor(ide,qp.ket2dm(qp.basis(self.max_cut,0)))
        # Squeezing operator; squeezing parameter = log(squeezing factor)
        sque = qp.tensor(ide,qp.squeeze(self.max_cut,np.log(self.s)))
        P = qp.operators.momentum(self.max_cut) # momentum operator
        
        # Prepare local rotation e^{-i h P}
        ehP_list = []
        for i in range(len(h_list)):
            ehP = (-1j/S*qp.tensor(h_list[i], P)).expm()
            ehP_list.append(ehP)
            
        for i in range(K):
            b_r = qp.tensor(b_den, r_den)
            for s in range(S):
                for j in range(len(ehP_list)):
                    b_r = ehP_list[j]*b_r*ehP_list[j].dag()
            # project the qumode to the squeezed position state
            b_r = proj*sque.dag()*b_r*sque*proj
            # Partial trace the ancillary qumode and then normalization
            b_den = (b_r.ptrace([i for i in range(N)])).unit()
            state_list.append(b_den)
            E_list.append(qp.expect(H,b_den))
        
        return E_list, state_list
    
    
    # Define noise Channel
    # Both the input and output are quantum states
    # p is the probability of noise happened
    def bit_flip_channel(self, state, p):
        # It is a bit flip channel
        # Each qubit will flip with a probability of p
        N = len(state.dims[0]) # number of qubits
        # Identity operator with probability 1 - p
        E0 = qp.tensor(np.sqrt(1 - p)*qp.tensor([qp.qeye(2)]*N),qp.qeye(self.max_cut))
        Ele_list = [] # operator Element list
        Ele_list.append(E0)
        
        for i in range(N):
            ide = [qp.qeye(2)]*N # a list of identity operator
            ide[i] = qp.sigmax() # replace the i-th one to sigma X
            Ele_list.append(np.sqrt(p)*qp.tensor(qp.tensor(ide), qp.qeye(self.max_cut)))
        
        state = sum(Ele_list[i]*state*Ele_list[i].dag() for i in range(N+1))
        return state
    
    def phase_flip_channel(self, state, p):
        # It is a phase flip channel
        # Each qubit will flip with a probability of p
        N = len(state.dims[0]) # number of qubits
        # Identity operator with probability 1 - p
        E0 = qp.tensor(np.sqrt(1 - p)*qp.tensor([qp.qeye(2)]*N),qp.qeye(self.max_cut))
        Ele_list = [] # operator Element list
        Ele_list.append(E0)
        
        for i in range(N):
            ide = [qp.qeye(2)]*N # a list of identity operator
            ide[i] = qp.sigmaz() # replace the i-th one to sigma Z
            Ele_list.append(np.sqrt(p)*qp.tensor(qp.tensor(ide), qp.qeye(self.max_cut)))
        
        state = sum(Ele_list[i]*state*Ele_list[i].dag() for i in range(N+1))
        return state

    def loss_bosonic_channel(self, state, p):
        # It is a lossing bosonic noise on the qumode
        # The Fock state decay with a probablity p
        def one_m_p_ada(p, max_cut):
            # The operator (1-p)^{(a^{dag} a)/2} below Eq. 9
            mat = np.zeros([max_cut, max_cut])
            for i in range(max_cut):
                mat[i][i] = (1 - p)**(i/2)
            return qp.Qobj(mat)
        
        max_cut = self.max_cut
        N = len(state.dims[0]) # number of qubits
        Ele_list = [] # operator Element list
        ide = qp.tensor([qp.qeye(2)]*N)
        
        for k in range(max_cut):
            E_k = qp.tensor(ide, (1/np.sqrt(float(np.math.factorial(k)))*p**(k/2)
                   *one_m_p_ada(p, max_cut)*qp.destroy(max_cut)**k))
            Ele_list.append(E_k)
        
        state = sum(Ele_list[i]*state*Ele_list[i].dag() for i in range(max_cut))
        return state
    
    
    # Evolution with First order Trotter decomposition and under noise
    def Trotter_evolution_under_noise(self, H, h_list, b, K, S, p_n):
        '''
        Solve the ground state by our algorithm, where the Hamiltonian
        the first order Trotter decomposition is considered.
        The quantum noise is applied after each Trotter evolution.
        It returns a list of Energy and state corresponding to each
        iteration.
        
        Args:
            H (Qobj): Hamiltonian to be solved
            h_list (list of Qobj): local operator list
            b (Qobj): initial state
            K (int): iteration step
            S (int): Trotter step
            p_n (dict): probability of noise; a dictionary contains px, pz, pl
        '''
        # initlization
        N = len(H.dims[0]) # number of qubits
        b_den = qp.ket2dm(b) # density matrix of b
        r_state = self.gen_r_state() # resource state in Fock space
        r_den = qp.ket2dm(r_state) # density matrix of resource state
        E_list = [qp.expect(H,b)] # Energy list
        state_list = [b_den] # state list
        
        px = p_n['px'] # probablity of bit flip
        pz = p_n['pz'] # probablity of phase flip
        pl = p_n['pl'] # probablity of lossing bosonic
        
        ide = qp.tensor([qp.qeye(2)]*N) # N qubits identity operator
        # It is a vacuum state being used to project the qumode.
        proj = qp.tensor(ide,qp.ket2dm(qp.basis(self.max_cut,0)))
        # Squeezing operator; squeezing parameter = log(squeezing factor)
        sque = qp.tensor(ide,qp.squeeze(self.max_cut,np.log(self.s)))
        P = qp.operators.momentum(self.max_cut) # momentum operator
        
        # Prepare local rotation e^{-i h P}
        ehP_list = []
        for i in range(len(h_list)):
            ehP = (-1j/S*qp.tensor(h_list[i], P)).expm()
            ehP_list.append(ehP)
            
        for i in range(K):
            b_r = qp.tensor(b_den, r_den)
            for s in range(S):
                for j in range(len(ehP_list)):
                    b_r = ehP_list[j]*b_r*ehP_list[j].dag()
                    # if the probablity of noise > 0, it happens
                    if px != 0:
                        b_r = self.bit_flip_channel(b_r, px)
                    if pz != 0:
                        b_r = self.phase_flip_channel(b_r, pz)
                    if pl != 0:
                        b_r = self.loss_bosonic_channel(b_r, pl)
                        
            # project the qumode to the squeezed position state
            b_r = proj*sque.dag()*b_r*sque*proj
            # Partial trace the ancillary qumode and then normalization
            b_den = (b_r.ptrace([i for i in range(N)])).unit()
            state_list.append(b_den)
            E_list.append(qp.expect(H,b_den))
        
        return E_list, state_list