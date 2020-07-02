
import numpy as np
from math import ceil
from sys import float_info

from biorad1_vangenuchten import VanGenuchten
from biorad1_base import BaseModel


class DataModel(BaseModel):

   @property
   def nodes_n(self):
      return len(self.nodes)

   @property
   def isot_n(self):
      return len(self.isot_names)

   
   # ************* units ****************

   def _load_units(self, main_key):                 
      # ********* length ***************
      length_units = {'mm':0.001, 'cm':0.01, 'dm':0.1, 'm':1}
      if not self.readdictvalue(length_units, main_key + ['length']): return False
      self.u_length, self.u_length_coef = self.lrvalue
      # ********* mass *****************
      mass_units = {'ng':1e-12, 'ug':1e-9, 'mg':1e-6, 'g':1e-3, 'kg':1}
      if not self.readdictvalue(mass_units, main_key + ['mass']): return False
      self.u_mass, self.u_mass_coef = self.lrvalue
      # ********* time ******************
      time_units = {'s':1/3600/24, 'min':1/60/24, 'h':1/24, 'day':1, 'year':365}
      if not self.readdictvalue(time_units, main_key + ['time']): return False
      self.u_time, self.u_time_coef = self.lrvalue
      print("Units loaded. OK.")
      return True


   # *********** horizons *****************************

   def _model_horizons(self, horizons_list_key, F_INIT=False):
      if not self.readlist(horizons_list_key, min_items_count = 2 if F_INIT else 1): return False
      n = len(self.lrvalue)
      H = []

      # ************ load bottoms ***********************
      i_top_h = None
      for i in range(n):
         if F_INIT:
            if self.key_exists_and_not_None(horizons_list_key + [i, 'top_head']): 
               i_top_h = i
               n -= 1
               continue
         if not self.readfloat(0, 1e30, horizons_list_key + [i, "bottom"]): return False
         H.append([i, self.lrvalue * self.u_length_coef, None, None, None])

      # ******** sorting by bottoms ********************
      def bottom_sort(b): return b[1]
      H.sort(key = bottom_sort)

      # ******** validate sorted list ******************
      if H[0][1] > 0:
         self.error("Bottom with zero value is required.", key=horizons_list_key, value="")
         return False
      if not F_INIT:
         if H[-1][1] > self.nodes[-1]:
            self.error("Bottom must be less or equal than model height.", key=horizons_list_key + [H[-1][0]], value=H[-1][1])
            return False
      else:
         if H[-1][1] >= self.nodes[-1]:
            self.error("Bottom must be less than model height.", key=horizons_list_key + [H[-1][0]], value=H[-1][1])
            return False      
      for i in range(n-1):
         if H[i][1] == H[i+1][1]:
            self.error("Bottom value is duplicate with another bottom.")
            return False

      # ***************** top **************
      for i in range(n-1): H[i][2] = H[i+1][1]
      H[-1][2] = self.nodes[-1]

      # ***************** indexes **********
      for h in H: h[3] = np.where(self.nodes >= h[1])[0][0]
      for i in range(n-1): H[i][4] = H[i+1][3]
      H[-1][4] = len(self.nodes)

      # ***************** ignore **********
      H_2 = []
      for i,h in enumerate(H):
         if h[3] < h[4]: H_2.append(h)
         else: self.warning(str(i) + " item will be ignored.")

      if F_INIT: 
         if i_top_h is None:
            self.error('Item with "top_head" was not found.', key=main_key, value="")
            return False
         H_2.append(i_top_h)

      # *** return [i, bottom, top, bottom_index, top_index+1] ****
      # *** top head index at the end, if flow initial condition are loaded
      return H_2
      

   # *********** mesh and materials ***************************

   def _load_mesh(self, main_key, vgfilename):
      # ************* geometry ***********************
      if not self.readfloat(1e-30, 1e30, main_key + ["height"]): return False
      h = self.lrvalue * self.u_length_coef
      if not self.readfloat(1e-30, 1e30, main_key + ["element_height"]): return False
      dz = self.lrvalue * self.u_length_coef
      if h / dz < 3:
         self.error("Element height is too big relative model height.")
         return False
      self.nodes = np.arange(0, h + dz, dz)
      if self.nodes[-1] > h:
         self.warning("The height of model was changed from {:g} {:s} to {:g} {:s}.".format(h / self.u_length_coef, self.u_length, self.nodes[-1] / self.u_length_coef, self.u_length))
      # ******** horizons **********************************
      self.theta_r, self.theta_s, self.alpha, self.Ks, self.n, self.density = [np.zeros(self.nodes_n) for _ in range(6)]
      self.horizon_id = np.zeros(self.nodes_n, dtype='int')
      hrz_key = main_key + ["horizons"]
      horizons = self._model_horizons(hrz_key)
      if not horizons: return False
      self.horizons_names = [None] * len(horizons)
      vg = None
      for i, b, t, ib, it in horizons:
         hik = hrz_key + [i]
         # ********* nodes params *****************
         self.horizon_id[ib:it] = i
         if not self.readstr(hik+["parameters_mode"]): return False
         params_mode = self.lrvalue
         if params_mode == "van_genuchten":
            if not self.readfloat(1, 1e5, hik+["density_kg_m3"]): return False
            self.density[ib:it] = self.lrvalue
            if not self.readfloat(0,1,hik+["theta_r"]): return False
            self.theta_r[ib:it] = self.lrvalue
            if not self.readfloat(0,1,hik+["theta_s"]): return False
            self.theta_s[ib:it] = self.lrvalue
            if not self.readfloat(0,1e10,hik+["alpha"]): return False
            self.alpha[ib:it] = self.lrvalue / self.u_length_coef
            if not self.readfloat(0,1e10,hik+["n"]): return False
            self.n[ib:it] = self.lrvalue
            if not self.readfloat(0,1e30,hik+["Ks"]): return False
            self.Ks[ib:it] = self.lrvalue * self.u_length_coef / self.u_time_coef
         elif params_mode == "granular_structure":
            if vg is None: 
               vg = VanGenuchten()
               if not vg.load_file(vgfilename): return False
            if not self.readfloat(1,1e5,hik+["density_kg_m3"]): return False
            self.density[ib:it] = density = self.lrvalue
            if not self.readfloat(0,100,hik+["sand"]): return False
            sand = self.lrvalue
            if not self.readfloat(0,100,hik+["silt"]): return False
            silt = self.lrvalue
            if not self.readfloat(0,100,hik+["clay"]): return False
            clay = self.lrvalue
            v = vg.get_params(density, sand, silt, clay)
            if v is None: return False
            self.theta_r[ib:it], self.theta_s[ib:it], self.alpha[ib:it], self.n[ib:it], self.Ks[ib:it] = v
         elif params_mode == "material":
            if vg is None: 
               vg = VanGenuchten()
               if not vg.load_file(vgfilename): return False
            if not self.readstr(hik+["material"]): return False
            material = self.lrvalue
            v = vg.get_params_by_meterial(material)
            if v is None: return False
            self.density[ib:it] = v[0]
            self.theta_r[ib:it], self.theta_s[ib:it], self.alpha[ib:it], self.n[ib:it], self.Ks[ib:it] = v[3:]
         else:
            self.error("Invalid parameters mode.")
            return False
         # ********* horizon name *****************
         if self.key_exists_and_not_None(hik + ["name"]):
            if not self.readstr(hik + ["name"]): return False
            self.horizons_names[i] = self.lrvalue
         else:
            self.horizons_names[i] = material if params_mode == "material" else "horizon {:g} - {:g} {:s}".format(b / self.u_length_coef, t / self.u_length_coef, self.u_length)
      # ******** next parameters ****************
      self.m = 1 - 1 / self.n
      self.l_pc = np.full((self.nodes_n), 0.5)
      print("Mesh loaded. OK.")
      return True

   
   # ************* flow initial conditions ************************

   def _load_flow_init_heads(self, main_key): 
      self.flow_init_heads = np.zeros(self.nodes_n)
      h_inits = self._model_horizons(main_key, F_INIT=True)
      if not h_inits: return False
      if not self.readfloat(-1e35, 1e35, main_key + [h_inits[-1], "top_head"]): return False
      t_h = self.lrvalue * self.u_length_coef
      h_inits.reverse()
      for i, b, t, ib, it in h_inits[1:]:
         if not self.readfloat(-1e35,1e35,main_key + [i, "head"]): return False
         b_h = self.lrvalue * self.u_length_coef
         # ****** head for nodes  **************************
         k = (t_h - b_h) / (t - b)
         for i in range(ib, it):            
            self.flow_init_heads[i] = k * (self.nodes[i] - b) + b_h
         t_h = b_h
      print("Flow initial conditions loaded. OK.")
      return True


   # ************* flow boundary conditions **************

   def _load_flow_bc_top(self, main_key):
      if not self.readlist(main_key): return False
      bc_count = len(self.lrvalue)
      self.flow_bc_top = []
      for i in range(bc_count):
         if not self.readfloat(0,1e50, main_key + [i, "time"]): return False
         t = self.lrvalue * self.u_time_coef
         if i == 0:
            if t != 0:
               self.error("First time of boundary condition must be zero.")
               return False
         elif t <= old_t:
            self.error("Time of top boudary condition must be more than previovs time.")
            return False
         if not self.readstr(main_key + [i, "type"]): return False
         if self.lrvalue == "dirichlet":
            if not self.readfloat(-1e35, 1e35, main_key + [i, "head"]): return False
            self.flow_bc_top.append((t, self.lrvalue * self.u_length_coef, None))
         elif self.lrvalue == "neumann":
            if not self.readfloat(-1e35, 1e35, main_key + [i, "flux"]): return False
            self.flow_bc_top.append((t, None, self.lrvalue * self.u_length_coef / self.u_time_coef))
         else:
            self.error("Value must be 'dirichlet' or 'neumann'.")
            return False
         old_t = t
      self.flow_bc_top.append((float_info.max, -1, -1))
      print("Flow top boundary conditions loaded. OK.")
      return True


   def _load_flow_bc_bottom(self, main_key):
      if not self.readlist(main_key): return False
      bc_count = len(self.lrvalue)
      self.flow_bc_bottom = []
      for i in range(bc_count):
         if not self.readfloat(0,1e50, main_key + [i, "time"]): return False
         t = self.lrvalue * self.u_time_coef
         if i == 0:
            if t != 0:
               self.error("First time of boundary condition must be zero.")
               return False
         elif t <= old_t:
            self.error("Time of bottom boudary condition must be more than previovs time.")
            return False
         if not self.readstr(main_key + [i, "type"]): return False
         if self.lrvalue == "dirichlet":
            if not self.readfloat(-1e35, 1e35, main_key + [i, "head"]): return False
            self.flow_bc_bottom.append((t, self.lrvalue * self.u_length_coef, None))
         elif self.lrvalue == "neumann":
            if not self.readfloat(-1e35, 1e35, main_key + [i, "flux"]): return False
            self.flow_bc_bottom.append((t, None, self.lrvalue * self.u_length_coef / self.u_time_coef))
         else:
            self.error("Value must be 'dirichlet' or 'neumann'.")
            return False
         old_t = t
      self.flow_bc_bottom.append((float_info.max, -1, -1))
      print("Flow bottom boundary conditions loaded. OK.")
      return True


   #  ********** flow sources ***************************

   def _load_flow_sources(self, main_key):
      self.flow_src_flux = np.zeros(self.nodes_n)
      sources = self._model_horizons(main_key)
      if not sources: return False
      for i, b, t, ib, it in sources:
         if not self.readfloat(-1e35, 1e35, main_key + [i, "flux_of_height_unit"]): return False
         self.flow_src_flux[ib:it] = self.lrvalue / self.u_time_coef
      print("Flow sources loaded. OK.")
      return True


   # ************ isotopes ************* ***********

   def _load_isotopes(self, main_key):
      if not self.readlist(main_key): return False
      n = len(self.lrvalue)
      self.isot_names = []
      self.isot_diff_coef = np.zeros((n,1))
      self.isot_dist_coef = np.zeros((n,1))
      for i in range(n):
         if not self.readstr(main_key + [i, 'name']): return False
         if self.isot_names.count(self.lrvalue) > 0:
            self.error("Isotope with this name is already used.")
            return False
         self.isot_names.append(self.lrvalue)
         if not self.readfloat(0,1e10, main_key + [i, 'dist_coef_m3_kg']): return False
         self.isot_dist_coef[i,0] = self.lrvalue
         if not self.readfloat(0,1e8,main_key + [i, 'diff_coef_m2_s']): return False 
         self.isot_diff_coef[i,0] = self.lrvalue * 3600 * 24
      print("Isotopes Loaded. OK.")
      return True


   def _readstr_isotope_index(self, isotope_key):         
      if not self.readstr(isotope_key): return False
      try:
         i = self.isot_names.index(self.lrvalue)
      except ValueError:
         self.error("Isotope was not found in isotopes list.")
         return False
      self.lrvalue = i
      return True


   def _load_isotopes_half_life(self, main_key):
      self.isot_half_life = np.full((self.isot_n,self.isot_n), -1, dtype="float")
      if not self.key_exists_and_not_None(main_key):
        self.warning("No isotope half life was set.")
        return True
      if not self.readlist(main_key): return False
      hl_count = len(self.lrvalue)
      for i in range(hl_count):
         if not self._readstr_isotope_index(main_key+[i,'isotope']): return False
         from_i = self.lrvalue        
         if not self._readstr_isotope_index(main_key+[i,'new_isotope']): return False
         to_i = self.lrvalue
         if not (from_i < to_i):
            self.error("Isotope can emerge from another isotope with lower ordering in isotpes list.", value=self.isot_names[to_i])
            return False
         if self.isot_half_life[from_i,to_i] > -1:
            self.error('Half life from isotope "'+self.isot_names[from_i]+'" to isotope "'+self.isot_names[to_i]+'" is already set.', value=self.isot_names[to_i])
            return False
         if not self.readfloat(0,1e60,main_key+[i,'half_life']): return False
         if self.lrvalue > 0:
           self.isot_half_life[from_i,to_i] = self.lrvalue * self.u_time_coef
      print("Isotopes half life loaded. OK.")
      return True


   # ******* transport global parameters ***********

   def _load_trans_global_params(self, main_key):
      # *********** tortuosity *****************
      if not self.readbool(main_key + ["tortuosity"]): return False
      self.trans_tortuosity = self.lrvalue
      # *********** dispersivity ***************
      if not self.readfloat(0, 1e30, main_key + ["dispersivity"]): return False
      self.trans_dispersivity = self.lrvalue * self.u_length_coef      
      # *********** numerical scheme ***********
      num_scheme = {"explicit":0, "crank_nicolson":0.5, "implicit":1}
      if not self.readdictvalue(num_scheme, main_key + ["numerical_scheme"]): return False
      self.trans_numerical_scheme = self.lrvalue[1]
      # ***** saturated zone concentration *****
      sat_c_key = main_key + ['saturated_zone_concentration']
      if not self.readbool(sat_c_key + ['apply']): return False
      if self.lrvalue:
         if not self.readfloat(0, 1e30, sat_c_key + ["height"]): return False
         z = self.lrvalue * self.u_length_coef
         if z > self.nodes[-2]:
            self.error("Value is too big relative model height.")
            return False
         self.trans_sat_c_nodes = np.where(self.nodes <= z)[0][-1] + 1
      else:
         self.trans_sat_c_nodes = 0
      print("Transport global parameters loaded. OK.")
      return True


   # ******* transport initial condition *********************

   def _load_trans_init_c(self, main_key):
      self.trans_init_c = np.zeros((self.isot_n,self.nodes_n))
      if not self.readlist(main_key): return False
      isot_init_count = len(self.lrvalue)
      for k in range(isot_init_count):
         if not self._readstr_isotope_index(main_key + [k, 'isotope']): return False
         isot_index = self.lrvalue
         isot_init_key = main_key + [k, 'concentration_in_water']
         isot_inits = self._model_horizons(isot_init_key)
         if not isot_inits: return False
         for i, b, t, ib, it in isot_inits:
            if not self.readfloat(0,1e30,isot_init_key + [i, 'c']): return False
            self.trans_init_c[isot_index,ib:it] = self.lrvalue * self.u_mass_coef / self.u_length_coef
      if self.trans_sat_c_nodes > -1:
         self.warning("Initial conditions to {:g} m will be ignored because saturated zone concentration is applied.".format(self.nodes[self.trans_sat_c_nodes]))
      print("Transport initial conditions loaded. OK.")
      return True


   # ********* transport boundary conditions **********

   def _load_trans_bc_top(self, main_key):
      if not self.readlist(main_key): return False
      isot_list = self.lrvalue
      times = {}
      for k in range(len(isot_list)):
         if not self._readstr_isotope_index(main_key + [k,'isotope']): return False
         isot_index = self.lrvalue
         isot_bc_key = main_key + [k,'time_function']
         if not self.readlist(isot_bc_key): return False
         isot_bc_list = self.lrvalue    
         for i in range(len(isot_bc_list)):
            if not self.readfloat(0,1e50,isot_bc_key+[i,"time"]): return False
            t = self.lrvalue * self.u_time_coef
            if i == 0:
               if t != 0:
                  self.error("First time of boundary condition must be zero.")
                  return False
            elif t <= old_t:
               self.error("Time of boudary condition must be more than previovs time.")
               return False               
            if not self.readfloat(0, 1e30, isot_bc_key+[i,"c_flux"]): return False
            if times.get(t, None) is None: times[t] = [-1 if t > 0 else 0] * self.isot_n
            times[t][isot_index] = self.lrvalue * self.u_mass_coef / self.u_length_coef
            old_t = t
      # ******** creating array *********************
      self.trans_bc_top = np.zeros((self.isot_n, len(times)))
      self.trans_bc_top_times = list(times.keys())
      self.trans_bc_top_times.sort()
      for ti, t in enumerate(self.trans_bc_top_times):    
         m = times[t]
         for k in range(self.isot_n):
            self.trans_bc_top[k,ti] = m[k] if m[k] != -1 else self.trans_bc_top[k,ti-1]
      self.trans_bc_top_times.append(float_info.max)
      print("Transport top boundary condition loaded. OK.")
      return True


   def _load_trans_bc_bottom(self, main_key):
      if not self.readlist(main_key): return False
      isot_list = self.lrvalue
      times = {}
      for k in range(len(isot_list)):
         if not self._readstr_isotope_index(main_key + [k,'isotope']): return False
         isot_index = self.lrvalue
         isot_bc_key = main_key + [k,'time_function']
         if not self.readlist(isot_bc_key): return False
         isot_bc_list = self.lrvalue    
         for i in range(len(isot_bc_list)):
            if not self.readfloat(0,1e50,isot_bc_key+[i,"time"]): return False
            t = self.lrvalue * self.u_time_coef
            if i == 0:
               if t != 0:
                  self.error("First time of boundary condition must be zero.")
                  return False
            elif t <= old_t:
               self.error("Time of boudary condition must be more than previovs time.")
               return False               
            if not self.readfloat(0, 1e30, isot_bc_key+[i,"c_flux"]): return False
            if times.get(t,None) is None: times[t] = [-1 if t > 0 else 0] * self.isot_n
            times[t][isot_index] = self.lrvalue * self.u_mass_coef / self.u_length_coef
            old_t = t
      # ****** creating boundary array **************
      self.trans_bc_bottom = np.zeros((self.isot_n, len(times)))
      self.trans_bc_bottom_times = list(times.keys())
      self.trans_bc_bottom_times.sort()
      for ti, t in enumerate(self.trans_bc_bottom_times):    
         m = times[t]
         for k in range(self.isot_n):
            self.trans_bc_bottom[k,ti] = m[k] if m[k] != -1 else self.trans_bc_bottom[k,ti-1]
      self.trans_bc_bottom_times.append(float_info.max)
      print("Transport bottom boundary condition loaded. OK.")
      return True


   # ********** simulation parameters ********************

   def _load_simulation_params(self, main_key):
      if not self.readfloat(1e-50,1e50,main_key + ['simulation_time']): return False
      self.sim_time = self.lrvalue * self.u_time_coef
      if not self.readint(1, 128, main_key + ['flow_iteration_count']): return False
      self.flow_iter_count = self.lrvalue
      if not self.readfloat(1e-50, 1e50, main_key + ['output_step_time']): return False
      self.output_step_time = self.lrvalue * self.u_time_coef
      if self.output_step_time > self.sim_time:
         print("Output step time must be less or equal than simulation time.")
         return False
      if not self.readfloat(1e-50,1e50,main_key + ['Dt']): return False
      self.Dt = self.lrvalue * self.u_time_coef
      if self.Dt > self.output_step_time:
         print("Computation step must be less or equal than output step time.")
         return False
      print("Simulation parameters loaded. OK.")
      return True


   # ************* outputs *********************************

   def _load_output_params(self, main_key):
      entities = {
         'nodes':['c_water', 'c_rock', 'test'],
         'elements':['V', 'c_water', 'c_rock', 'm_water', 'm_rock'] , 
         'summary':['V', 'm_water', 'm_rock']
         }
      isotope_required = {'c_water', 'c_rock', 'm_water' , 'm_rock'}
      file_formats = {'gmesh_v2_ASCII', 'txt', 'screen'}
      # ******** outputs list **************************
      if not self.readlist(main_key): return False
      out_list_count = len(self.lrvalue)
      self.output_params = []
      for k in range(out_list_count):
         # ********** entity ***********************
         if not self.readdictvalue(entities, main_key + [k, "entity"]): return False
         entity, q_list = self.lrvalue
         # ********** quantity **********************
         if not self.readstr_selection(q_list, main_key + [k, "physical_quantity"]): return False
         quantity = self.lrvalue
         # ******** isotopes list *******************
         isotopes = []
         if quantity in isotope_required:
            if not self.key_exists_and_not_None(main_key + [k, 'isotopes']):
               for i,m in enumerate(self.isot_names):
                  isotopes.append((i, m))
            else:
               if not self.readlist(main_key + [k, 'isotopes']): return False
               n = len(self.lrvalue)
               for i in range(n):
                  if not self._readstr_isotope_index(main_key + [k, 'isotopes', i]): return False
                  isotopes.append((self.lrvalue, self.isot_names[self.lrvalue]))
         # *********** file format **********************
         if not self.readstr_selection(file_formats, main_key + [k, 'file_format']): return False
         file_format = self.lrvalue
         # *********** file name ************************
         if not self.readstr(main_key + [k, 'file_name']): return False
         file_name = self.lrvalue
         # *********** creating output generator *******
         self.output_params.append((entity, quantity, isotopes, file_format, file_name))
      print("Output parameters loaded. OK.")
      return True


   # ****************************************************
   # *********** loading data model *********************
  
   def load_data_model(self):      
      print("Creating data model ...")
      if not self._load_units(['units']): return False
      if not self._load_mesh(['mesh'], 'vangenuchten.csv'): return False
      if not self._load_flow_init_heads(['flow','initial_conditions']): return False
      if not self._load_flow_bc_top(['flow','top_boundary_conditions']): return False
      if not self._load_flow_bc_bottom(['flow','bottom_boundary_conditions']): return False
      if not self._load_flow_sources(['flow','sources']): return False      
      if not self._load_isotopes(['transport','isotopes']): return False
      if not self._load_isotopes_half_life(['transport','isotopes_half_life']): return False
      if not self._load_trans_global_params(['transport']): return False
      if not self._load_trans_init_c(['transport','initial_conditions']): return False
      if not self._load_trans_bc_top(['transport','top_boundary_conditions']): return False
      if not self._load_trans_bc_bottom(['transport','bottom_boundary_conditions']): return False
      if not self._load_simulation_params(['simulation_parameters']): return False
      if not self._load_output_params(['outputs']): return False      
      print ("Data model was created. OK.")
      return True




