
import numpy as np
from scipy.linalg import solveh_banded, solve_banded
from scipy.sparse import diags
from math import log, ceil
from time import time
from datetime import timedelta

from biorad1_data import DataModel
from biorad1_output import OutputModel, OutputAppendError


class ComputingModel(DataModel):


   # *********** allocations ***********************
   
   def _states_allocation(self):
      i, n = self.isot_n, self.nodes_n

      # ********* flow ****************************
      self.f_h = np.zeros(n)
      self.f_theta_h = np.zeros(n)
      self.PR_f_theta_h = np.zeros(n)
      self.f_sat = np.zeros(n)
      self.f_K = np.zeros(n)
      self.f_cap = np.zeros(n)
      self.f_flux = np.zeros(n)
      self.PR_f_flux = np.zeros(n)

      # *********** transport ***********************
      self.r_c_water = np.zeros((i, n))
      self.PR_r_c_water = np.zeros((i, n))
      
      # *********** identity ******************************
      self.mat_ident = np.identity(n)


   # *******************************************************
   # ***************** flow ********************************

   def _gen_flow_iteration(self):

      # ****** constants definition **************
      nodes_n, Dx, Dt = self.nodes_n, self.nodes[1], self.Dt
      Dx2, Dt2 = 2 * Dx, 2 * Dt

      # ****** help calculation *******************
      theta_s_minus_theta_r = self.theta_s - self.theta_r
      abs_alpha_x_f_h = np.zeros(nodes_n)

      # ****** arrays allocation ***************** 
      mat_P_ab = np.zeros((2, nodes_n))
      mat_e, mat_d = mat_P_ab[0,1:], mat_P_ab[1]
      mat_f = np.zeros(nodes_n)

      # ***** Time 0 - initializations *******************
      i_bc_top, i_bc_bottom = 0, 0
      self.f_h[:] = self.flow_init_heads

      # ****** iterations *******************
      while True:

         abs_alpha_x_f_h[:] = np.absolute(self.alpha * self.f_h)

         # ******* water retention curve ***********
         self.f_theta_h[:] = np.where(self.f_h < 0,
              self.theta_r + theta_s_minus_theta_r /  \
              (1 + abs_alpha_x_f_h ** self.n) ** self.m,
              self.theta_s )
        
         # ****** saturation ****************
         self.f_sat[:] = (self.f_theta_h - self.theta_r) / theta_s_minus_theta_r

         # ****** conductivity ***************************
         self.f_K[:] = self.Ks * self.f_sat**self.l_pc * (1 - (1 - self.f_sat**(1/self.m))**self.m)**2

         # ****** capacity ******************************
         self.f_cap[:] = np.where(self.f_h < 0, \
                theta_s_minus_theta_r *  \
                ( self.alpha * self.m * self.n ) *  \
                ( abs_alpha_x_f_h ** (self.n-1) ) /  \
                ( (1 + abs_alpha_x_f_h ** self.n) ** (self.m+1) )  \
             , 0 )
         
         # ************** stop ***********
         yield

         # ***** calculation of iteration solving **********
         mat_d[1:-1] = Dx * self.f_cap[1:-1] / Dt + (self.f_K[:-2] + 2 * self.f_K[1:-1] + self.f_K[2:]) / Dx2
         mat_f[1:-1] = Dx * self.f_cap[1:-1] * self.f_h[1:-1] / Dt - Dx * (self.f_theta_h[1:-1] - self.PR_f_theta_h[1:-1]) / Dt +  \
                       (self.f_K[2:] - self.f_K[:-2]) / 2 - self.flow_src_flux[1:-1]
         mat_e[1:-1] = - (self.f_K[1:-2] + self.f_K[2:-1]) / Dx2

         # ********* top boundary conditions *********
         while self.t >= self.flow_bc_top[i_bc_top+1][0]: i_bc_top += 1
         dir_h, neum_flux = self.flow_bc_top[i_bc_top][1:]
         if dir_h is not None:
            mat_d[-1], mat_e[-1], mat_f[-1] = 1, 0, dir_h
            mat_f[-1] = mat_f[-1]
            mat_f[-2] += dir_h * (self.f_K[-2] + self.f_K[-1]) / Dx2
         else:
            mat_d[-1] = Dx * self.f_cap[-1] / Dt2 + (self.f_K[-1] + self.f_K[-2]) / Dx2
            mat_e[-1] = - (self.f_K[-1] + self.f_K[-2]) / Dx2
            mat_f[-1] = Dx * self.f_cap[-1] * self.f_h[-1] / Dt2 - Dx * (self.f_theta_h[-1] - self.PR_f_theta_h[-1]) / Dt2  \
                        - (self.f_K[-1] + self.f_K[-2]) / 2 - neum_flux

         # ********* bottom boundary conditions *******
         while self.t >= self.flow_bc_bottom[i_bc_bottom+1][0]: i_bc_bottom += 1
         dir_h, neum_flux = self.flow_bc_bottom[i_bc_bottom][1:]
         if dir_h is not None:
            mat_d[0], mat_e[0], mat_f[0] = 1, 0, dir_h
            mat_f[0] = mat_f[0]
            mat_f[1] += dir_h * (self.f_K[0] + self.f_K[1]) / Dx2
         else:
            mat_d[0] = (self.f_K[0] + self.f_K[1]) / Dx2
            mat_e[0] = - (self.f_K[0] + self.f_K[1]) / Dx2
            mat_f[0] = (self.f_K[0] + self.f_K[1]) / 2 + neum_flux

         # *********** calculate_f_h *****************
         self.f_h[:] = solveh_banded(mat_P_ab, self.mat_ident, overwrite_ab=True, check_finite=False) @ mat_f


   def _gen_flow_states_after_last_iteration(self):

      # ****** constants definition **************
      Dx, Dt = self.nodes[1], self.Dt

      # ******* flux of init time on highest node *******
      self.f_flux[-1] = - (self.f_K[-1] + self.f_K[-2]) / 2 * ((self.f_h[-1] - self.f_h[-2]) / Dx + 1)

      while True:

         # ******** flux without highest node ***************
         self.f_flux[1:-1] =  \
                  - (self.f_K[2:] + self.f_K[1:-1]) / 4 * ((self.f_h[2:] - self.f_h[1:-1]) / Dx + 1)  \
                  - (self.f_K[1:-1] + self.f_K[:-2]) / 4 * ((self.f_h[1:-1] - self.f_h[:-2]) / Dx + 1)

         self.f_flux[0] = - (self.f_K[1] + self.f_K[0]) / 2 * ((self.f_h[1] - self.f_h[0]) / Dx + 1)

         # **** stop and get CR ***********
         yield np.max(np.absolute(self.f_flux) * Dt / self.f_theta_h / Dx)
         
         # **** flux on highest node for t > 0 **************
         self.f_flux[-1] = - (self.f_K[-1] + self.f_K[-2]) / 2 * ((self.f_h[-1] - self.f_h[-2]) / Dx + 1)  \
               - Dx / 2 * (self.f_theta_h[-1] - self.PR_f_theta_h[-1]) / Dt



   # ********************************************************
   # ************ transport **********************************


   def _gen_trans_step(self):

      # ********** constants definitions ****************
      isot_n, nodes_n, Dt, Dx = self.isot_n, self.nodes_n, self.Dt, self.nodes[1]
      eps = self.trans_numerical_scheme
      Dx12Dt, epsB, epsE, epsF, PepsB, PepsE, PepsF = Dx/12/Dt, eps/6, eps/Dx/2, eps*Dx/12, (1-eps)/6, (1-eps)/Dx/2, (1-eps)*Dx/12

      # *********** help calculation ********************
      theta_s_pow2 = self.theta_s**2
      density_x_dist_coef = self.density * self.isot_dist_coef

      # *********** arrays allocations *******************
      diff_coef = np.zeros((isot_n, nodes_n))
      B, E, F, G, Rk, PR_B, PR_E, PR_F, PR_G, PR_Rk = [np.zeros((isot_n, nodes_n)) for _ in range(10)]

      mat_P_ab = np.zeros((isot_n, 3, nodes_n))
      mat_P_u = mat_P_ab[:,0,1:]
      mat_P_d = mat_P_ab[:,1,:]
      mat_P_l = mat_P_ab[:,2,:-1]

      mat_T_d = np.zeros((isot_n, nodes_n))
      mat_T_u, mat_T_l = np.zeros((isot_n, nodes_n - 1)), np.zeros((isot_n, nodes_n - 1))
      mat_T_l_d_u_diags = [[mat_T_l[i], mat_T_d[i], mat_T_u[i]] for i in range(isot_n)]

      mat_R, mat_s = np.zeros(nodes_n), np.zeros(nodes_n)

      mat_inv_P = np.zeros((nodes_n, nodes_n))
    
      hm_top, hm_bottom = np.zeros((isot_n,1)), np.zeros((isot_n,1))
      hm_sources = np.zeros((isot_n, nodes_n - 2))

      # **** matrix of isotope conversion after Dt ********
      isot_conversion_rate = np.where(self.isot_half_life > -1, log(2) * Dt / self.isot_half_life, 0)
      sum_conversion_rate = - np.sum(isot_conversion_rate, axis=1).reshape((isot_n,1))
      
      # ****** TIME 0 - boundary conditions  **********
      i_bc_top, i_bc_bottom = 0, 0

      # ****** TIME 0 - initial conditions  *****************
      self.r_c_water[:] = self.trans_init_c
      self.r_c_water[:,:self.trans_sat_c_nodes] = self.trans_bc_bottom[:,[0]]

      # ****** TIME 0 - diff coefficient *******************
      if self.trans_tortuosity: 
            diff_coef[:] = (self.trans_dispersivity * np.absolute(self.f_flux) + self.f_theta_h * self.isot_diff_coef * self.f_theta_h**2.3333333 / theta_s_pow2) / self.f_theta_h
      else: diff_coef[:] = (self.trans_dispersivity * np.absolute(self.f_flux) + self.f_theta_h * self.isot_diff_coef) / self.f_theta_h

      # ****** TIME 0 -  matrix B, E, F, G ****************
      E[:] = self.f_theta_h * diff_coef
      B[:] = self.f_flux
      F[:] = - sum_conversion_rate * (self.f_theta_h + density_x_dist_coef)
      Rk[:] = 1 + density_x_dist_coef / self.f_theta_h
      for i in range(isot_n): G[i] = np.sum(isot_conversion_rate[:i,[i]] * (self.f_theta_h + density_x_dist_coef[:i]) * self.r_c_water[:i], axis=0)

      while True:

         # ******* result of previous step ******************
         B, PR_B = PR_B, B
         E, PR_E = PR_E, E
         F, PR_F = PR_F, F
         G, PR_G = PR_G, G
         Rk, PR_Rk = PR_Rk, Rk

         # ******** stop and get Pe *******************
         yield (np.absolute(self.f_flux) * Dx / (self.f_theta_h * diff_coef)).max()

         # ******* diff coefficient *******************
         if self.trans_tortuosity: 
               diff_coef[:] = (self.trans_dispersivity * np.absolute(self.f_flux) + self.f_theta_h * self.isot_diff_coef * self.f_theta_h**2.3333333 / theta_s_pow2) / self.f_theta_h
         else: diff_coef[:] = (self.trans_dispersivity * np.absolute(self.f_flux) + self.f_theta_h * self.isot_diff_coef) / self.f_theta_h
         
         # ******** matrix B, E, F, G **********
         E[:] = self.f_theta_h * diff_coef
         B[:] = self.f_flux
         F[:] = - sum_conversion_rate * (self.f_theta_h + density_x_dist_coef)
         Rk[:] = 1 + density_x_dist_coef / self.f_theta_h
         
         # ******* equation 8.26 ********************
         mat_P_d[:,0] = Dx12Dt * (3 * self.f_theta_h[0] * Rk[:,0] + self.f_theta_h[1] * Rk[:,1])  \
            + epsE  * (E[:,0] + E[:,1])  \
            + epsB  * (2 * B[:,0] + B[:,1])  \
            + epsF  * (3 * F[:,0] + F[:,1])
         mat_T_d[:,0] = Dx12Dt * (3 * self.PR_f_theta_h[0] * PR_Rk[:,0] + self.PR_f_theta_h[1] * PR_Rk[:,1])  \
            - PepsE * (PR_E[:,0] + PR_E[:,1])  \
            - PepsB * (2 * PR_B[:,0] + PR_B[:,1])  \
            - PepsF * (3 * PR_F[:,0] + PR_F[:,1])

         # ******* equation 8.27 *********************
         mat_P_l[:] = Dx12Dt * (self.f_theta_h[:-1] * Rk[:,:-1] + self.f_theta_h[1:] * Rk[:,1:])  \
            - epsE  * (E[:,:-1] + E[:,1:])  \
            - epsB  * (2 * B[:,:-1] + B[:,1:])  \
            + epsF * (F[:,:-1] + F[:,1:])
         mat_T_l[:] = Dx12Dt * (self.PR_f_theta_h[:-1] * PR_Rk[:,:-1] + self.PR_f_theta_h[1:] * PR_Rk[:,1:])  \
            + PepsE * (PR_E[:,:-1] + PR_E[:,1:])  \
            + PepsB * (2 * PR_B[:,:-1] + PR_B[:,1:])  \
            - PepsF * (PR_F[:,:-1] + PR_F[:,1:])
        
         # ******* equation 8.28 ***********************
         mat_P_d[:,1:-1] = Dx12Dt * (self.f_theta_h[:-2] * Rk[:,:-2] + 6 * self.f_theta_h[1:-1] * Rk[:,1:-1] + self.f_theta_h[2:] * Rk[:,2:])  \
            + epsE * (E[:,:-2] + 2 * E[:,1:-1] + E[:,2:])  \
            + epsB * (B[:,2:] - B[:,:-2])  \
            + epsF * (F[:,:-2] + 6 * F[:,1:-1] + F[:,2:])
         mat_T_d[:,1:-1] = Dx12Dt * (self.PR_f_theta_h[:-2] * PR_Rk[:,:-2] + 6 * self.PR_f_theta_h[1:-1] * PR_Rk[:,1:-1] + self.PR_f_theta_h[2:] * PR_Rk[:,2:])  \
            - PepsE * (PR_E[:,:-2] + 2 * PR_E[:,1:-1] + PR_E[:,2:])  \
            - PepsB * (PR_B[:,2:] - PR_B[:,:-2])  \
            - PepsF * (PR_F[:,:-2] + 6 * PR_F[:,1:-1] + PR_F[:,2:])
 
         # ******* equation 8.29 *********************
         mat_P_u[:] = Dx12Dt * (self.f_theta_h[:-1] * Rk[:,:-1] + self.f_theta_h[1:] * Rk[:,1:])  \
            - epsE * (E[:,:-1] + E[:,1:])  \
            + epsB * (2 * B[:,1:] + B[:,:-1])  \
            + epsF * (F[:,:-1] + F[:,1:])
         mat_T_u[:] = Dx12Dt * (self.PR_f_theta_h[:-1] * PR_Rk[:,:-1] + self.PR_f_theta_h[1:] * PR_Rk[:,1:])  \
            + PepsE * (PR_E[:,:-1] + PR_E[:,1:])  \
            - PepsB * (2 * PR_B[:,1:] + PR_B[:,:-1])  \
            - PepsF * (PR_F[:,:-1] + PR_F[:,1:])

         # ******* equation 8.30 ********************
         mat_P_d[:,-1] = Dx12Dt * (self.f_theta_h[-2] * Rk[:,-2] + 3 * self.f_theta_h[-1] * Rk[:,-1])  \
            + epsE  * (E[:,-2] + E[:,-1])  \
            - epsB  * (B[:,-2] + 2 * B[:,-1])  \
            + epsF  * (F[:,-2] + 3 * F[:,-1])
         mat_T_d[:,-1] = Dx12Dt * (self.PR_f_theta_h[-2] * PR_Rk[:,-2] + 3 * self.PR_f_theta_h[-1] * PR_Rk[:,-1])  \
            - PepsE * (PR_E[:,-2] + PR_E[:,-1])  \
            + PepsB * (PR_B[:,-2] + 2 * PR_B[:,-1])  \
            - PepsF * (PR_F[:,-2] + 3 * PR_F[:,-1])

         # ****** set bc time index **********************
         while self.t >= self.trans_bc_top_times[i_bc_top + 1]: i_bc_top += 1
         while self.t >= self.trans_bc_bottom_times[i_bc_bottom+1]: i_bc_bottom += 1

         # ******* tansport over top ****************************
         if self.f_flux[-1] < 0:            
            hm_top[:] = - self.f_flux[-1] * self.trans_bc_top[:,[i_bc_top]]
         else:
            hm_top[:] = - self.f_flux[-1] * self.PR_r_c_water[:,[-1]]
        
         # ******* transport over bottom *************************
         if self.f_flux[0] < 0:
            hm_bottom[:] = eps * self.f_flux[0] * self.PR_r_c_water[:,[0]]  \
                         + (1-eps) * self.PR_f_flux[0] * self.PR_r_c_water[:,[0]]
         else: 
            hm_bottom[:] = self.trans_bc_bottom[:,[i_bc_bottom]]  \
                         * (eps * self.f_flux[0] + (1-eps) * self.PR_f_flux[0])
         
         # ****** transport over sources *************************
         hm_sources[:] = np.where(self.flow_src_flux[1:-1] > 0, - self.flow_src_flux[1:-1] * Dx * self.PR_r_c_water[:,1:-1], 0)

         for i in range(isot_n):

            # *********** matrix G ***********************
            G[i] = np.sum(isot_conversion_rate[:i,[i]] * (self.f_theta_h + density_x_dist_coef[:i]) * self.r_c_water[:i], axis=0)
          
            # ******* vector R - equation 8.25 **************
            mat_s[:] = (eps * G[i] + (1-eps) * PR_G[i]) / 6
            mat_R[0] = Dx * (2 * mat_s[0] + mat_s[1]) + hm_bottom[i]
            mat_R[1:-1] = Dx * (mat_s[:-2] + 4 * mat_s[1:-1] + mat_s[2:]) + hm_sources[i]
            mat_R[-1] = Dx * (mat_s[-2] + 2 * mat_s[-1]) + hm_top[i]

            # ******* new concentrations ********************
            mat_inv_P[:] = solve_banded((1,1), mat_P_ab[i], self.mat_ident, overwrite_ab=True, check_finite=False)
            self.r_c_water[i] = mat_inv_P  \
               @ ( diags(mat_T_l_d_u_diags[i], (-1,0,1))  \
                     @ self.PR_r_c_water[i] + mat_R )

            # ****** saturation concentration **************
            self.r_c_water[i,:self.trans_sat_c_nodes] = self.trans_bc_bottom[i,i_bc_bottom]



   # *********** actual states to prior **************

   def _actual_states_to_prior_states(self):
      self.PR_f_theta_h[:] = self.f_theta_h
      self.PR_f_flux[:] = self.f_flux
      self.PR_r_c_water[:] = self.r_c_water



   # *************************************************
   # ************** main computation *****************


   def compute(self):

      print ("Starting of simulation ...")

      # ************ initialize ************
      self._states_allocation()
      gen_flow_iter = self._gen_flow_iteration()
      gen_flow_step = self._gen_flow_states_after_last_iteration()
      gen_tran_step = self._gen_trans_step()
      try:
         self.output_model = OutputModel(self)
      except OutputAppendError as e:
         print(e)
         return False


      # ******** time 0 ********************
      print("Computing of time 0 ...")
      self.t = 0
      computing_start_time = time()
      progress = 0
      next(gen_flow_iter)
      Cr_max, Cr_max_t = next(gen_flow_step), 0
      Pe_max, Pe_max_t = next(gen_tran_step), 0
      self._actual_states_to_prior_states()
      self.output_model.write_state()
      
      # ******** cycle over steps **************
      print("Computing next time steps ...")
      while self.t < self.sim_time:
         
         # ****** new time ************
         self.t += self.Dt

         # ****** flow step ***********
         for i in range(self.flow_iter_count): next(gen_flow_iter)
         r = next(gen_flow_step)   
         if r > Cr_max: Cr_max, Cr_max_t = r, self.t
         if r > 100: print("Warning time {:g}: Cr = {:g}".format(self.t, r))

         # ****** transport step ***********
         r = next(gen_tran_step)
         if r > Pe_max: Pe_max, Pe_max_t = r, self.t

         # ******** write results ********** 
         if self.t >= self.output_model.out_t:
            self.output_model.write_state()
       
         # ******* progress ********************
         if 100 * self.t // self.sim_time > progress:
            progress = 100 * self.t // self.sim_time
            print("{:g}% calculated \r".format(progress), end="")

         # **** actual states as prior states for next step ***
         self._actual_states_to_prior_states()

      self.output_model.finish()
      
      ct = time() - computing_start_time
      print("Computation was finished. OK.")
      print("Computing time: " + str(timedelta(seconds=ct)) + ".")
      print("Maximum Cr in computed time {:g}: {:g}".format(Cr_max_t / self.u_time_coef, Cr_max))
      print("Maximum Pe in computed time {:g}: {:g}".format(Pe_max_t / self.u_time_coef, Pe_max))
      
      return True

