
import numpy as np

class OutputAppendError(Exception): pass

class OutputModel:

   def __init__(self, computing_model):
      self.computing_model = computing_model
      self.generators = []
      print("Creating output generators ...")
      for i, o in enumerate(self.computing_model.output_params):
         if not self.append_output(*o): raise OutputAppendError("Error of loading %d. output parameters." %(i))
         print("{:g}. Output generator for file {:s} was created. OK.".format(i+1, o[4]))
      # ********* initial values ************
      self.out_t_index = 0
      self.out_t = 0
      # *** alocate array for interpolated states *****
      i, n = self.computing_model.isot_n, self.computing_model.nodes_n
      self.out_f_theta_h = np.zeros(n)
      self.out_r_c_water = np.zeros((i,n))

      
   # ******* gmesh support methods *************************

   def _gmesh_v2_ASCII_mesh(self, f):
      c = self.computing_model
      f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
      f.write("$PhysicalNames\n"+str(len(c.horizons_names))+"\n")
      for i, name in enumerate(c.horizons_names):
         print(1, i, '"' + name + '"', file=f)
      f.write("$EndPhysicalNames\n")
      f.write("$Nodes\n"+str(c.nodes_n)+"\n")
      for i in range(c.nodes_n):
         print (i, 0, 0, "{:g}".format(c.nodes[i] / c.u_length_coef), file=f)
      f.write("$EndNodes\n")
      f.write("$Elements\n"+str(c.nodes_n - 1)+"\n")
      for i in range(c.nodes_n - 1):
         print (i, 1, 2, c.horizon_id[i], 0, i, i + 1, file=f)
      f.write("$EndElements\n")

      
   def _gmesh_v2_ASCII_NodeData(self, f, data_name, data_time, data_time_index, data_number_format, data):
       print("$NodeData", 1, data_name, 1, data_time, 3, data_time_index, 1, len(data), sep='\n', file=f)
       for i, d in enumerate(data): print(i, data_number_format.format(d), file=f)
       f.write("$EndNodeData\n")


   def _gmesh_v2_ASCII_ElementData(self, f, data_name, data_time, data_time_index, data_number_format, data):
       print("$ElementData", 1, data_name, 1, data_time, 3, data_time_index, 1, len(data), sep='\n', file=f)
       for i, d in enumerate(data): print(i, data_number_format.format(d), file=f)
       f.write("$EndElementData\n")
      
       

   # ******* generators control *********************

   def write_state(self):
      cp = self.computing_model
      T = cp.t - self.out_t

      # *** states interpolation ********
      self.out_f_theta_h = (cp.f_theta_h - cp.PR_f_theta_h) / cp.Dt * T + cp.PR_f_theta_h
      self.out_r_c_water[:] = (cp.r_c_water - cp.PR_r_c_water) / cp.Dt * T + cp.PR_r_c_water

      # *** writing of outputs for actual time ***
      for g in self.generators: next(g)
      print("{:g}. {:s} was written.".format(self.out_t / cp.u_time_coef, cp.u_time))

      # ******* set future output time *********
      self.out_t += cp.output_step_time
      self.out_t_index += 1


   def finish(self):
      for g in self.generators: g.close()
      print("Output generators were finished. OK.")



   # ******* generators instances creations ***************

   def append_output(self, entity, quantity, isotopes, file_format, file_name):
      g = self.generators.append
      p = (entity, quantity, file_format)
      if   p == ('nodes', 'c_water', 'gmesh_v2_ASCII'): g(self._out_nodes_c_water_gmesh_v2_ASCII(isotopes, file_name))
      elif p == ('nodes', 'c_rock',  'gmesh_v2_ASCII'): g(self._out_nodes_c_rock_gmesh_v2_ASCII (isotopes, file_name))

      elif p == ('nodes', 'test',    'screen'):         g(self._out_nodes_test_screen())


      elif p == ('elements', 'V',    'gmesh_v2_ASCII'): g(self._out_elements_V_gmesh_v2_ASCII(file_name))

      elif p == ('summary', 'V', 'txt'): g(self._out_summary_V_txt (file_name))
      elif p == ('summary', 'm_water', 'txt'): g(self._out_summary_m_water_txt (isotopes, file_name))
      elif p == ('summary', 'm_rock', 'txt'): g(self._out_summary_m_rock_txt (isotopes, file_name))


      else:
         print("Invalid combination of output: " + str(p) + ".")
         return False
      return True



   # *******************************************************
   # ********** nodes generators ***************************

   def _out_nodes_c_water_gmesh_v2_ASCII(self, isotopes, file_name):
      cp = self.computing_model
      f = open(file_name, "w")
      data = np.zeros((cp.isot_n, cp.nodes_n))
      self._gmesh_v2_ASCII_mesh(f)
      while True:
         data[:] = self.out_r_c_water / cp.u_mass_coef * cp.u_length_coef
         for i, name in isotopes:
            self._gmesh_v2_ASCII_NodeData(f, name, self.out_t / cp.u_time_coef, self.out_t_index, "{:.4e}", data[i])
         f.flush()      
         try: yield
         except GeneratorExit:
            f.close(); return


   def _out_nodes_c_rock_gmesh_v2_ASCII(self, isotopes, file_name):
      cp = self.computing_model
      data = np.zeros((cp.isot_n, cp.nodes_n))
      f = open(file_name, "w")
      self._gmesh_v2_ASCII_mesh(f)
      while True:
         data[:] = self.out_r_c_water * cp.isot_dist_coef
         for i, name in isotopes:
            self._gmesh_v2_ASCII_NodeData(f, name, self.out_t / cp.u_time_coef, self.out_t_index, "{:.4e}", data[i])
         f.flush()      
         try: yield
         except GeneratorExit:
            f.close(); return





   # *****************************************************
   # ************** elements generators ******************

   def _out_elements_V_gmesh_v2_ASCII(self, file_name):
      cp = self.computing_model
      data = np.zeros(cp.nodes_n - 1)
      Dx2u = cp.nodes[1] / 2 / cp.u_length_coef
      f = open(file_name, "w")
      self._gmesh_v2_ASCII_mesh(f)
      while True:
         data[:] = Dx2u * (self.out_f_theta_h[:-1] + self.out_f_theta_h[1:])
         self._gmesh_v2_ASCII_ElementData(f, "V_water", self.out_t / cp.u_time_coef, self.out_t_index, "{:.4e}", data)
         f.flush()      
         try: yield
         except GeneratorExit:
            f.close(); return




   # **************************************************
   # **************** summary generators **************

   def _out_summary_V_txt(self, file_name):
      cp = self.computing_model
      f = open(file_name, "w")
      print("time\tV", file=f)
      Dx2 = cp.nodes[1] / 2
      while True:
         V = np.sum(self.out_f_theta_h[:-1] + self.out_f_theta_h[1:]) * Dx2 / cp.u_length_coef
         print("{:g}\t{:.4e}".format(self.out_t / cp.u_time_coef, V), file=f)
         f.flush()      
         try: yield
         except GeneratorExit:
            f.close(); return


   def _out_summary_m_water_txt(self, isotopes, file_name):
      cp = self.computing_model
      f = open(file_name, "w")
      print("time\t" + "\t".join([n for i, n in isotopes]), file=f)
      i_isot = [i for i, n in isotopes]
      m_nodes = np.zeros((cp.isot_n, cp.nodes_n))
      m_sum = np.zeros(cp.isot_n)
      Dx = cp.nodes[1]
      while True:
         m_nodes[:] = self.out_r_c_water[i_isot] * self.out_f_theta_h * Dx
         m_sum[:] = (0.5 * m_nodes[i_isot,0] + np.sum(m_nodes[i_isot,1:-1],axis=1) + 0.5 * m_nodes[i_isot,-1]) / cp.u_mass_coef
         print("{:g}".format(self.out_t / cp.u_time_coef) + (len(i_isot)*"\t{:.4e}").format(*m_sum), file=f)
         f.flush()      
         try: yield
         except GeneratorExit:
            f.close(); return


   def _out_summary_m_rock_txt(self, isotopes, file_name):
      cp = self.computing_model
      f = open(file_name, "w")
      print("time\t" + "\t".join([n for i, n in isotopes]), file=f)
      i_isot = [i for i, n in isotopes]
      m_nodes = np.zeros((cp.isot_n, cp.nodes_n))
      m_sum = np.zeros(cp.isot_n)
      Dx = cp.nodes[1]
      while True:
         m_nodes[:] = self.out_r_c_water[i_isot] * cp.isot_dist_coef[i_isot] * cp.density * Dx
         m_sum[:] = (0.5 * m_nodes[i_isot,0] + np.sum(m_nodes[i_isot,1:-1],axis=1) + 0.5 * m_nodes[i_isot,-1]) / cp.u_mass_coef
         print("{:g}".format(self.out_t / cp.u_time_coef) + (len(i_isot)*"\t{:.4e}").format(*m_sum), file=f)
         f.flush()      
         try: yield
         except GeneratorExit:
            f.close(); return







   # ******************************************************
   # *********** test generators *************************
   # ***** POUZE PRO TESTOVACI UCELY !!! ****************

   def _out_nodes_test_screen(self):
      cp = self.computing_model
      

      while True:
         if cp.t == 1000 * cp.Dt:
            print(cp.t)
            for i in range(cp.nodes_n-1, -1, -1):
               print(i,
                     "{:.2f}".format(cp.nodes[i]),
                     "{:.8e}".format(cp.f_h[i]),
                     "{:.8e}".format(cp.f_theta_h[i]),
                     "{:.8e}".format(cp.f_K[i]),
                     "{:.8e}".format(cp.f_flux[i]),
                     "{:.8e}".format(cp.r_c_water[0,i]),
                     "{:.8e}".format(cp.r_c_water[10,i]),

                     sep='\t')

        
                  
         yield
         






