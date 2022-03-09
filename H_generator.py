import numpy as np
import qutip as qp

class H_generator():
    def build_Ising_model(self, N, a, J, shift_E):
        '''
        It constructs a transversed field Ising model,
        which presented in Eq. 10 in the article.
        
        Args:
            N: number of qubits.
            a: a vector with length N; tranverse field on each site.
            J: a matrix N*N; coupling strength for two sites.
            shift_E: shift-energy applied.
        '''
        
        sx_list = [] # i-th element is sigma x performing on i-th qubit
        sz_list = [] # i-th element is sigma z performing on i-th qubit
        
        for i in range(N):
            ope_list = [qp.qeye(2)]*N # a list of identity operator
            ope_list[i] = qp.sigmax() # replace the i-th element by sigma x
            sx_list.append(qp.tensor(ope_list))
            
            ope_list = [qp.qeye(2)]*N # a list of identity operator
            ope_list[i] = qp.sigmaz() # replace the i-th element by sigma z
            sz_list.append(qp.tensor(ope_list))
        
        # Transverse field term
        Hx = 0
        for i in range(N):
            Hx = Hx + a[i]*sx_list[i]
        
        # Interaction term
        Hz = 0
        for i in range(N):
            for j in range(i):
                Hz = Hz + J[i][j]*sz_list[i]*sz_list[j]
        
        # Total Hamiltonian
        H = Hx + Hz + shift_E
        return H
    
    
    def Ising_model_decomposition(self, N, a, J, shift_E):
        '''
        It gives the componenet of Ising model, which presented
        in Eq. 10 in the article. It can be used to apply the 
        Hamiltonian with Trotter decomposition.
        
        Args:
            N: number of qubits.
            a: a vector with length N; tranverse field on each site.
            J: a matrix N*N; coupling strength for two sites.
            shift_E: shift-energy applied.
        '''
        
        sx_list = [] # i-th element is sigma x performing on i-th qubit
        sz_list = [] # i-th element is sigma z performing on i-th qubit
        
        for i in range(N):
            ope_list = [qp.qeye(2)]*N # a list of identity operator
            ope_list[i] = qp.sigmax() # replace the i-th element by sigma x
            sx_list.append(qp.tensor(ope_list))
            
            ope_list = [qp.qeye(2)]*N # a list of identity operator
            ope_list[i] = qp.sigmaz() # replace the i-th element by sigma z
            sz_list.append(qp.tensor(ope_list))
        
        h_list = []
        # Transverse field term
        for i in range(N):
            h_list.append(a[i]*sx_list[i])
        
        # Interaction term
        for i in range(N):
            for j in range(i):
                h_list.append(J[i][j]*sz_list[i]*sz_list[j])
        
        # Energy of shift
        h_list.append(shift_E*qp.tensor([qp.qeye(2)]*N))
        return h_list
    
    
    def build_H2(self, ind, shift_E):
        '''
        It returns a Hydrogen molecular Hamiltonian and
        its corresponding bond distance.
        
        Args:
            ind: the index of data; the bond distance is from 0.2
            to 2.85 and the the interval is 0.05.
            shift_E: shift-energy applied.
        '''
        
        # Data from paper (arXiv:1512.06860v2) table 1: R, I, Z0, Z1, Z0Z1, X0X1, Y0Y1
        data = [
        [0.20,  2.8489, 0.5678, -1.4508, 0.6799, 0.0791, 0.0791],
        [0.25,  2.1868, 0.5449, -1.2870, 0.6719, 0.0798, 0.0798],
        [0.30,  1.7252, 0.5215, -1.1458, 0.6631, 0.0806, 0.0806],
        [0.35,  1.3827, 0.4982, -1.0226, 0.6537, 0.0815, 0.0815],
        [0.40,  1.1182, 0.4754, -0.9145, 0.6438, 0.0825, 0.0825],
        [0.45,  0.9083, 0.4534, -0.8194, 0.6336, 0.0835, 0.0835],
        [0.50,  0.7381, 0.4325, -0.7355, 0.6233, 0.0846, 0.0846],
        [0.55,  0.5979, 0.4125, -0.6612, 0.6129, 0.0858, 0.0858],
        [0.60,  0.4808, 0.3937, -0.5950, 0.6025, 0.0870, 0.0870],
        [0.65,  0.3819, 0.3760, -0.5358, 0.5921, 0.0883, 0.0883],
        [0.70,  0.2976, 0.3593, -0.4826, 0.5818, 0.0896, 0.0896],
        [0.75,  0.2252, 0.3435, -0.4347, 0.5716, 0.0910, 0.0910],
        [0.80,  0.1626, 0.3288, -0.3915, 0.5616, 0.0925, 0.0925],
        [0.85,  0.1083, 0.3149, -0.3523, 0.5518, 0.0939, 0.0939],
        [0.90,  0.0609, 0.3018, -0.3168, 0.5421, 0.0954, 0.0954],
        [0.95,  0.0193, 0.2895, -0.2845, 0.5327, 0.0970, 0.0970],
        [1.00, -0.0172, 0.2779, -0.2550, 0.5235, 0.0986, 0.0986],
        [1.05, -0.0493, 0.2669, -0.2282, 0.5146, 0.1002, 0.1002],
        [1.10, -0.0778, 0.2565, -0.2036, 0.5059, 0.1018, 0.1018],
        [1.15, -0.1029, 0.2467, -0.1810, 0.4974, 0.1034, 0.1034],
        [1.20, -0.1253, 0.2374, -0.1603, 0.4892, 0.1050, 0.1050],
        [1.25, -0.1452, 0.2286, -0.1413, 0.4812, 0.1067, 0.1067],
        [1.30, -0.1629, 0.2203, -0.1238, 0.4735, 0.1083, 0.1083],
        [1.35, -0.1786, 0.2123, -0.1077, 0.4660, 0.1100, 0.1100],
        [1.40, -0.1927, 0.2048, -0.0929, 0.4588, 0.1116, 0.1116],
        [1.45, -0.2053, 0.1976, -0.0792, 0.4518, 0.1133, 0.1133],
        [1.50, -0.2165, 0.1908, -0.0666, 0.4451, 0.1149, 0.1149],
        [1.55, -0.2265, 0.1843, -0.0549, 0.4386, 0.1165, 0.1165],
        [1.60, -0.2355, 0.1782, -0.0442, 0.4323, 0.1181, 0.1181],
        [1.65, -0.2436, 0.1723, -0.0342, 0.4262, 0.1196, 0.1196],
        [1.70, -0.2508, 0.1667, -0.0251, 0.4204, 0.1211, 0.1211],
        [1.75, -0.2573, 0.1615, -0.0166, 0.4148, 0.1226, 0.1226],
        [1.80, -0.2632, 0.1565, -0.0088, 0.4094, 0.1241, 0.1241],
        [1.85, -0.2684, 0.1517, -0.0015, 0.4042, 0.1256, 0.1256],
        [1.90, -0.2731, 0.1472,  0.0052, 0.3992, 0.1270, 0.1270],
        [1.95, -0.2774, 0.1430,  0.0114, 0.3944, 0.1284, 0.1284],
        [2.00, -0.2812, 0.1390,  0.0171, 0.3898, 0.1297, 0.1297],
        [2.05, -0.2847, 0.1352,  0.0223, 0.3853, 0.1310, 0.1310],
        [2.10, -0.2879, 0.1316,  0.0272, 0.3811, 0.1323, 0.1323],
        [2.15, -0.2908, 0.1282,  0.0317, 0.3769, 0.1335, 0.1335],
        [2.20, -0.2934, 0.1251,  0.0359, 0.3730, 0.1347, 0.1347],
        [2.25, -0.2958, 0.1221,  0.0397, 0.3692, 0.1359, 0.1359],
        [2.30, -0.2980, 0.1193,  0.0432, 0.3655, 0.1370, 0.1370],
        [2.35, -0.3000, 0.1167,  0.0465, 0.3620, 0.1381, 0.1381],
        [2.40, -0.3018, 0.1142,  0.0495, 0.3586, 0.1392, 0.1392],
        [2.45, -0.3035, 0.1119,  0.0523, 0.3553, 0.1402, 0.1402],
        [2.50, -0.3051, 0.1098,  0.0549, 0.3521, 0.1412, 0.1412],
        [2.55, -0.3066, 0.1078,  0.0572, 0.3491, 0.1422, 0.1422],
        [2.60, -0.3079, 0.1059,  0.0594, 0.3461, 0.1432, 0.1432],
        [2.65, -0.3092, 0.1042,  0.0614, 0.3433, 0.1441, 0.1441],
        [2.70, -0.3104, 0.1026,  0.0632, 0.3406, 0.1450, 0.1450],
        [2.75, -0.3115, 0.1011,  0.0649, 0.3379, 0.1458, 0.1458],
        [2.80, -0.3125, 0.0997,  0.0665, 0.3354, 0.1467, 0.1467],
        [2.85, -0.3135, 0.0984,  0.0679, 0.3329, 0.1475, 0.1475]]
        
        H = qp.tensor(qp.ket2dm(qp.Qobj(np.zeros(2))),qp.ket2dm(qp.Qobj(np.zeros(2))))
        H = H + data[ind][1]*qp.tensor(qp.qeye(2),qp.qeye(2))
        H = H + qp.tensor(data[ind][2]*qp.sigmaz(),qp.qeye(2))
        H = H + qp.tensor(qp.qeye(2),data[ind][3]*qp.sigmaz())
        H = H + data[ind][4]*qp.tensor(qp.sigmaz(),qp.sigmaz())
        H = H + data[ind][5]*qp.tensor(qp.sigmax(),qp.sigmax())
        H = H + data[ind][6]*qp.tensor(qp.sigmay(),qp.sigmay())
        H = H + shift_E*qp.tensor(qp.qeye(2),qp.qeye(2))
        
        return data[ind][0], H
    
    def build_Kitaev_ring(self, N, h, J, shift_E):
        '''
        It constructs a Kitaev ring, which presented in
        Eq. 12 in the article.
        
        Args:
            N: number of qubits.
            h: a vector with length N; tranverse field on each site.
            J: a matrix N*N; coupling strength for two sites.
            shift_E: shift-energy applied.
        '''
        
        sx_list=[] # i-th element is sigma x performing on i-th qubit
        sy_list=[] # i-th element is sigma y performing on i-th qubit
        sz_list=[] # i-th element is sigma z performing on i-th qubit
        
        for i in range(N):
            ope_list = [qp.qeye(2)]*N # a list of identity operator
            ope_list[i] = qp.sigmax() # replace the i-th element by sigma x
            sx_list.append(qp.tensor(ope_list))
            
            ope_list = [qp.qeye(2)]*N # a list of identity operator
            ope_list[i] = qp.sigmay() # replace the i-th element by sigma y
            sy_list.append(qp.tensor(ope_list))
            
            ope_list = [qp.qeye(2)]*N # a list of identity operator
            ope_list[i] = qp.sigmaz() # replace the i-th element by sigma z
            sz_list.append(qp.tensor(ope_list))
        
        # Hz term
        Hz = 0
        for i in range(N):
            Hz = Hz + sz_list[i]
        # Hxx term
        Hxx = 0
        for i in range(N-1):
            Hxx = Hxx + sx_list[i]*sx_list[i+1]
        # Hyz term
        Hyz = sy_list[0]
        for i in range(N-2):
            Hyz = Hyz*sz_list[i+1]
        Hyz = Hyz*sy_list[N-1]
        
        # Total Hamiltonian
        H = -h*Hz - J*Hxx - J*Hyz + shift_E
        return H