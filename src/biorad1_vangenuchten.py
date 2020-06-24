
import numpy as np
from pathlib import Path


base_path = Path(__file__).parent

class VanGenuchten:

   def __init__(self):
      self.materials = {}

      
   # ***** nacteni VG parametru ze souboru ***********

   def load_file(self, vgfile):
      try:
         f = open((base_path / "../database/vangenuchten.csv").resolve(),"r")
      except:
         print('Error with opening file "'+vgfile+'".')
         return
      with f:
         print('Loading file "'+vgfile+'" ...')
         # ******** head of file *******************
         f.readline()
         # ******** initial states *****************
         R = []
         hustota, pisek, prach = set(), set(), set()         
         i = 2
         # **** reading of row and grouping of sand, silt and clay
         s = f.readline()
         while s != "":
            s = s.strip()
            if s != "" and not s.startswith("%"):
               col_str = [c.strip() for c in s.split(",")]
               if len(col_str) != 9:
                  print("Error with column count on the line %d." %(i))
                  return False
               try:
                  col_float = list(map(float,col_str[1:]))
               except:
                  print("Error with reading float number on the line %d." %(i))
                  return False
               # *** converting alpha 1/cm --> 1/m ****
               col_float[5] = col_float[5] * 100
               # *** converting Ks cm/day --> m/day *****
               col_float[7] = col_float[7] / 100
               R.append(col_float)
               if col_str[0] != "":
                  self.materials[col_str[0]] = col_float
               hustota.add(col_float[0])
               pisek.add(col_float[1])
               prach.add(col_float[2])
            s = f.readline()
            i += 1
      self.hustota = list(hustota)
      self.hustota.sort()
      self.pisek = list(pisek)
      self.pisek.sort()
      self.prach = list(prach)
      self.prach.sort()
      # ********* nacteni VG hodnot ******************
      self.vg = np.full((len(self.hustota)+1,len(self.pisek)+1,len(self.prach)+1,5), -1, dtype=np.float64)
      for col in R:
         self.vg[self.hustota.index(col[0]),
                 self.pisek.index(col[1]),
                 self.prach.index(col[2]),
                 :] = col[3:]
      print("File was loaded.")
      return True


   # ********** vypocet VG parametru *******************

   def _uprava_VG_vstupu(self, VG_hustota, VG_pisek, VG_prach, VG_jil):
      if VG_pisek < 0 or VG_pisek > 100:
         print("Sand rate is out of interval <0,100>.")
         return
      if VG_prach < 0 or VG_prach > 100:
         print("Silt rate is out of interval <0,100>.")
         return
      if VG_jil < 0 or VG_jil > 100: 
         print("Clay rate is out of interval <0,100>.")
         return
      n = VG_pisek + VG_prach + VG_jil 
      if n == 100:
         prach, pisek = VG_prach, VG_pisek
      else:
         prach, pisek = VG_prach*100/n, VG_pisek*100/n
         print("The rate of sand and silt was changed.")
      if VG_hustota < 500:
         hustota = 500
         print("Low density was changed to 500kg/m3.")
      elif VG_hustota > 2000:
         hustota = 2000
         print("High density was changed to 2000kg/m3.")
      else:
         hustota = VG_hustota
      return hustota, pisek, prach


   def _interpolace_VG_params(self, VG_hustota, VG_pisek, VG_prach):
      # ****** nalezeni indexu 0 ************************
      h0 = len(self.hustota)-1
      while self.hustota[h0] > VG_hustota: h0 -= 1
      s0 = len(self.pisek)-1
      while self.pisek[s0] > VG_pisek: s0 -= 1
      p0 = len(self.prach)-1
      while self.prach[p0] > VG_prach: p0 -= 1
      # ********** kontrola platnosti ******************
      x = -1
      if self.vg[h0,s0,p0,0] == x \
         or (VG_hustota > self.hustota[h0] and self.vg[h0+1,s0,p0,0]==x) \
         or (VG_pisek > self.pisek[s0] and self.vg[h0,s0+1,p0,0]==x) \
         or (VG_prach > self.prach[p0] and self.vg[h0,s0,p0+1,0]==x):
            print("Vstupní data pro interpolaci VG parametrů nejsou v platném definičním oboru.")
            return
      # ********** provedeni interpolace ***************
      v = []
      for i in range(5):
         v0 = self.vg[h0,s0,p0,i]
         v.append(v0 + \
            (VG_hustota - self.hustota[h0]) * \
               (self.vg[h0+1,s0,p0,i] - v0) / (self.hustota[h0+1] - self.hustota[h0]) + \
            (VG_pisek - self.pisek[s0]) * \
               (self.vg[h0,s0+1,p0,i] - v0) / (self.pisek[s0+1] - self.pisek[s0]) + \
            (VG_prach - self.prach[p0]) * \
               (self.vg[h0,s0,p0+1,i] - v0) / (self.prach[p0+1] - self.prach[p0])
         )
      return v


   def get_params(self, VG_hustota, VG_pisek, VG_prach, VG_jil):
      """
      Return list of parameters in following order:
      [theta_r, theta_s, alpha[1/m-1], n, Ks[m/day]
      """
      v = self._uprava_VG_vstupu(VG_hustota, VG_pisek, VG_prach, VG_jil)
      if v is None: return
      return self._interpolace_VG_params(v[0], v[1], v[2])
      

   def get_params_by_meterial(self, material_name):
      """
      Return list of parameters in following order:
      [density, sand, silt, theta_r, theta_s, alpha[1/m-1], n, Ks[m/day]]
      """
      v = self.materials.get(material_name, None)
      if v is not None: return v
      else:
         print("Material '"+material_name+"' was not found.")
         return

      
