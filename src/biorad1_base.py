
from yaml import safe_load


class BaseModel:

           
   def keystr(self, key):
      return " --> ".join((str(k) if type(k) != int else str(k+1)+". item" for k in key))


   def key_exists_and_not_None(self, key):
      value = self.input_data
      for p in key:
         try:
            value = value[p]
         except KeyError:
            return False
      return value is not None


   # ********** informations printing ***********
      
   def warning(self, *text):
      print("Warning:", *text)


   def error(self, text, key=None, value=None):
      print("Error:")
      if key is None: print("Key: " + self.keystr(self.lrkey))
      else: print("Key: " + self.keystr(key))
      if value is None: print("Value: " + str(self.lrvalue))
      else: print("Value: " + str(value))
      print(text)


   # ******* reading any type ****************

   def readvalue(self, key, allow_None=False):
      value = self.input_data
      for p in key:
         try:
            value = value[p]
         except KeyError:
            print('Error:')
            print('Key "'+self.keystr(key)+'" was not found.')
            return False
         self.lrkey = key
         self.lrvalue = value
         if not allow_None and self.lrvalue is None:
            self.error("Value must be set.")
            return False
      return True


   # ********* reading specific types ***************

   def readfloat(self, x_min, x_max, key, open_interval=False):
      if not self.readvalue(key): return False
      try:
         self.lrvalue = float(self.lrvalue)
      except:     
         self.error("Value must be valid float number.")
         return False
      if open_interval:
         if self.lrvalue <= x_min or self.lrvalue >= x_max:
            self.error("Value must be from interval ("+str(x_min)+";"+str(x_max)+").")
            return False
      else:
         if self.lrvalue < x_min or self.lrvalue > x_max:
            self.error("Value must be from interval <"+str(x_min)+";"+str(x_max)+">.")
            return False
      return True


   def readint(self, x_min, x_max, key, open_interval=False):
      if not self.readvalue(key): return False
      try:
         self.lrvalue = int(self.lrvalue)
      except:     
         self.error("Value must be valid integer number.")
         return False
      if open_interval:
         if self.lrvalue <= x_min or self.lrvalue >= x_max:
            self.error("Value must be from interval ("+str(x_min)+";"+str(x_max)+").")
            return False
      else:
         if self.lrvalue < x_min or self.lrvalue > x_max:
            self.error("Value must be from interval <"+str(x_min)+";"+str(x_max)+">.")
            return False
      return True


   def readstr(self, key, allow_None_or_empty=False):
      if not self.readvalue(key, allow_None=allow_None_or_empty): return False
      if self.lrvalue is not None:
         try: self.lrvalue = str(self.lrvalue)
         except ValueError:
            self.error("The value must be valid string.")
            return False
         if self.lrvalue.strip() == "" and allow_None_or_empty:
            self.error("String value is empty.")
            return False
      return True


   def readbool(self, key):
      if not self.readvalue(key): return False
      e = False
      if type(self.lrvalue) == str:
         s = self.lrvalue.lower()
         if s in {'yes', 'true'}: self.lrvalue = True
         elif s in {'no', 'false'}: self.lrvalue = False
         else : e = True
      elif type(self.lrvalue) == bool:
         pass
      else:
         e = True
      if e:
         self.error('Boolean value "yes", "no", "true" or "false" is required.')
         return False
      return True


   def readlist(self, key, min_items_count=1):
      if not self.readvalue(key): return False
      if type(self.lrvalue) != list:
         self.error("Value must be sequence.")
         return False
      if len(self.lrvalue) < min_items_count:
         self.error("The sequence must have minimal %d items" %(min_items_count))
         return False
      return True


   def readdictvalue(self, dict, key):
      if not self.readstr(key): return False
      k = self.lrvalue
      v = dict.get(k, None)
      if v is None:
         self.error("Value can contain: " + str(list(dict.keys()))[1:-1] + ".")
         return False
      self.lrvalue = (k, v)
      return True


   def readstr_selection(self, ls, key):
      if not self.readstr(key): return False
      if self.lrvalue not in ls:
         self.error("Value can contain: " + str(ls)[1:-1] + ".")
         return False
      return True


   # ***************************************************
   # ******** load input_data *************************

   def load_yaml(self, yaml_file_name):
      print('Loading file "' + yaml_file_name + '" ...')
      try:
         f = open(yaml_file_name, "r", encoding="utf8")
         self.input_data = safe_load(f)
      except Exception as e:
         print(e)
         print('Error of "' + yaml_file_name + '" loading.')
         return False
      print('File "' + yaml_file_name + '" was loaded. OK.')
      return True

