# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 15:08:49 2019

@author: jiri.landa
"""

hybhen=u"\u002D" #spojovnik (zakladni pomlcky, minus)

#biorad keywords begin
#consumption
wat_ing_path="ing_water" #identificator of water ingestion pathway
soil_ing_path="ing_soil" #identificator of water ingestion pathway
root_ing_path="root_vegetables" #identificator of root vegetable ingestion pathway
leafy_ing_path="leafy_vegetables" #identificator of leafy vegetable ingestion pathway
potatoes_ing_path="potatoes" #identificator of potatoes ingestion pathway
mushrooms_ing_path="mushrooms" #identificator of mushroom ingestion pathway
fish_ing_path="fish" #identificator of mushroom ingestion pathway
#plants parameters
plant_env="plant_field" #oznaceni rostlin na poli
plant_hum_ingestion="ing_plant_prod_w_resp" #ingesce rostliny clovekem
#animals
pasture="feed" #identification of pasture plant product for animal feed
cattle="cattle" #cattle animal
#index identificators
conf_nuclides="nuclides" #identificator of nuclides in config file
conf_msh_el="msh_el" #identificator of mash elements in config file
conf_acts_type="activities" #identificator of activities type in config file
conf_baskets="baskets" #identificator of basket type in config file
#index of activities
groundwater="water_ground"
surfacewater="water_surface"
soil="soil_field"
wood="soil_forest"
total="total"
#biorad keywords end
endline="\n"
#sum of iput_data_error

list_data_errors=[]

import yaml
import numpy as np
import sqlite3
import re

def get_nucl_name(nuclname,hyphen='',total=total):
    if nuclname==total:
        return nuclname
    l=1
    if 'a'<=nuclname[l]<='z':
        l=2
    q=l
    if nuclname[q]<'1' or '9' <nuclname[q]:
        q=q+1
    return nuclname[0:l]+hyphen+nuclname[q:len(nuclname)]

def get_ch_el_name(nuclname):
    l=1
    if 'a'<=nuclname[l]<='z':
        l=2
    return nuclname[0:l]

def get_atomic_weight_from_nucl_name(nuclide,prt=False):
    ns=re.findall("[-+]?\d*\.\d+|\d+", nuclide)
    #weight=0.050 # i.e. 50 g/mol = 0.050 kg/mol
    #weight=0.056 # i.e. 56 g/mol = 0.056 kg/mol  Fe56
    weight=0.120 # i.e. 120 g/mol = 0.120 kg/mol Sn120
    if len(ns)>0 and float(ns[0])>0:
        weight=float(ns[0])/1000
    return weight

def get_atomic_weight(prot_file,db_cur,nuclide_index,errors,prt=False):
    aw=np.zeros((len(nuclide_index)),dtype=float)
    for rn in nuclide_index:
        rn_ndb=True
        sql_ask="select * from nuclide where symbol='"+rn+"'"
        sql_ask_e="\""+sql_ask+"\""
        db_cur.execute(sql_ask)
        if len(db_cur.fetchall())>0:
            sql_ask="select * from w_nuclide_parameters where nuclide='"+rn+"' and param_identificator='mass'"
            if prt:
                print(sql_ask)
            db_cur.execute(sql_ask)
            rec=db_cur.fetchall()
            if len(rec)>0:
                rn_ndb=False
                col_index=get_db_ask_column_names_index(db_cur,sql_ask,prt)
                if prt:
                    print(col_index)
                    print(rec)
                    print(rec[0])
                    print(rec[0][col_index["value"]])
                    print(aw[nuclide_index[rn]])
                aw[nuclide_index[rn]]=rec[0][col_index["value"]]*rec[0][col_index["conversion"]]
            else:
                sql_ask_e="\""+sql_ask+"\""
                
        if rn_ndb:
            if prt:
                print(type(rn))
                print(type(nuclide_index[rn]))
            aw[nuclide_index[rn]]=get_atomic_weight_from_nucl_name(rn,prt)
            sql_ask_e+=" is NULL"+endline+"the weight of the nuclide set to "+str(aw[nuclide_index[rn]])+" [kg.mol-1]"
            errors=write_input_data_error(prot_file,sql_ask_e,errors,"nucl_mass",prt)
    return aw,errors

def get_half_life(prot_file,db_cur,nuclide_index,a2s,errors,prt=False):
    hl=np.zeros((len(nuclide_index)),dtype=float)
    for rn in nuclide_index:
        rn_ndb=True
        sql_ask="select * from nuclide where symbol='"+rn+"'"
        sql_ask_e="\""+sql_ask+"\""
        db_cur.execute(sql_ask)
        if len(db_cur.fetchall())>0:
            sql_ask="select * from w_nuclide_parameters where nuclide='"+rn+"' and param_identificator='half_life'"
            if prt:
                print(sql_ask)
            db_cur.execute(sql_ask)
            rec=db_cur.fetchall()
            if len(rec)>0:
                rn_ndb=False
                col_index=get_db_ask_column_names_index(db_cur,sql_ask,prt)
                if prt:
                    print(col_index)
                    print(rec)
                    print(rec[0])
                    print(rec[0][col_index["value"]])
                    print(hl[nuclide_index[rn]])
                hl[nuclide_index[rn]]=rec[0][col_index["value"]]*rec[0][col_index["conversion"]]
            else:
                sql_ask_e="\""+sql_ask+"\""
                
        if rn_ndb:
            hlinf=1e20 # [a]
            if prt:
                print(type(rn))
                print(type(nuclide_index[rn]))
            hl[nuclide_index[rn]]=hlinf*a2s
            sql_ask_e+=" is NULL"+endline+"the half life of the nuclide set to "+str(hl[nuclide_index[rn]])+" [s] (i.e. "+str(hlinf)+" [a])"
            errors=write_input_data_error(prot_file,sql_ask_e,errors,"nucl_mass",prt)
    return hl,errors

def set_specific_activity(t12,M,Na,prt=False):
    sa=np.zeros((len(t12)),dtype=float)
    numerator=Na*np.log(2)
    for i in range(len(t12)):
        sa[i]=numerator/(t12[i]*M[i])
    return sa

def set_activity_from_concentration(activities,spec_activity,prt=False):
    act=activities*0
    dimensions=act.shape
    
    for i in range(dimensions[0]):
        for j in range(dimensions[1]):
            for k in range(dimensions[2]):
                for l in range(dimensions[3]):
                    act[i][j][k][l]=activities[i][j][k][l]*spec_activity[k]
    return act


def read_index(data,index_ident,prt=False):
    index={}
    i=-1        
    for d in data[index_ident]:
        i+=1
        index[d]=i
        if prt:
            print(d)
    return index

def read_activities_conf(data,prt=False):
    conf=False
    if data["activity_type"]=="all":
        conf=True
    return conf

def read_activities_conf2(data,prt=False):
    return data["activity_type"]

def read_cont_data_type(data,prt=False):
    key="contamination_data_type"
    data_type="concentration"
    if key in data:
        if data[key]=="act" or data[key]=="activity" or data[key]=="activities":
            data_type="activity"
    return data_type


def get_env_index(br_db_cur,prt=False):
    index_rec=[None,'','']
    index={}
    
    sql_ask="SELECT * FROM s_env"
    col_index=get_db_ask_column_names_index(br_db_cur,sql_ask,prt)
    br_db_cur.execute(sql_ask)
    rows = br_db_cur.fetchall()

    i=-1
    for row in rows:
        i+=1
        if prt:
            print(row,type(row))
        index_rec[0]=i
        index_rec[1]=row[col_index["name_cz"]]
        index_rec[2]=row[col_index["name_en"]]
        
        index[row[col_index["sub"]]]=index_rec.copy()
    return index


def get_env_type(br_db_cur,env_index,prt=False):
    index={}
    
    sql_ask="SELECT * FROM s_env"
    col_index=get_db_ask_column_names_index(br_db_cur,sql_ask,prt)

    for env in env_index:
        sql_ask="SELECT * FROM s_env where sub='"+str(env)+"'"
        br_db_cur.execute(sql_ask)
        rec = br_db_cur.fetchall()

        if prt:
            print(rec,type(rec))
        index[rec[0][col_index["sub"]]]=rec[0][col_index["super"]]

    return index


def get_env_param_index(br_db_cur,prt=False):
    index_rec=[None,'','']
    index={}
    env_params=[
    "density",
    "porosity",
    "humidity",
    "dust_level",
    "deposit",
    "dilution"
    ]
    
    sql_ask="SELECT * FROM Parameter"
    col_index=get_db_ask_column_names_index(br_db_cur,sql_ask,prt)
    if prt:
        print(col_index)

    i=-1
    for param in env_params:
        i+=1
        sql_ask="SELECT * FROM Parameter where identificator ='"+param+"'"
        if prt:
            print(sql_ask)
        br_db_cur.execute(sql_ask)
        rec = br_db_cur.fetchall()

        if prt:
            print("rec[0]:",rec[0])
            print("rec[0][0]:",rec[0][0])
            print(rec,type(rec))
            print(col_index["name_cz"])
            print(rec[0][col_index["name_cz"]])
            
        index_rec[0]=i
        index_rec[1]=rec[0][col_index["name_cz"]]
        index_rec[2]=rec[0][col_index["name_en"]]
        
        index[rec[0][col_index["identificator"]]]=index_rec.copy()

    return index

def get_env_params(prot_file,br_db_cur,env_index,env_par_index,prt=False):
    prt_function_string="get environmental parameters"
    prot_file.write(prt_function_string+" begin"+endline)
    env_params=np.zeros((len(env_index),len(env_par_index)),dtype=float)

    sql_ask="SELECT * FROM w_parameters"
    col_index=get_db_ask_column_names_index(br_db_cur,sql_ask,prt)
    if prt:
        print(col_index)
    for env in env_index:
        for par in env_par_index:
            sql_ask="SELECT * FROM w_parameters where identificator='"+str(env)+"' and param_identificator='"+str(par)+"'"
            br_db_cur.execute(sql_ask)
            rec= br_db_cur.fetchall()
            if prt:
                print("printx:",env_index[env][0],env_par_index[par][0])
                print("printy:",col_index["value"],col_index["conversion"])
            env_params[env_index[env][0]][env_par_index[par][0]]=rec[0][col_index["value"]]*rec[0][col_index["conversion"]]
            #env_params[env_index[env][0]][env_par_index[par][0]]=1
    prot_file.write(prt_function_string+" end"+endline+endline)
    return env_params

def get_Kd(prot_file,br_db_cur,env_index,nuclide_index,env_type,errors,prt=False):
    prt_function_string="get distribution coefficients (Kd)"
    prot_file.write(prt_function_string+" begin"+endline)
    Kd=np.zeros((len(env_index),len(nuclide_index)),dtype=float)
    sql_ask="SELECT * FROM w_tf"
    col_index=get_db_ask_column_names_index(br_db_cur,sql_ask,prt)
    for env in env_index:
        if env_type[env]=="soil":
            for rn in nuclide_index:
                element=get_ch_el_name(rn)
                sql_ask="SELECT * FROM w_tf where element='"+element+"' and identificator='"+str(env)+"' and literature='IAEA'"
                br_db_cur.execute(sql_ask)
                rec= br_db_cur.fetchall()
                if len(rec)>0:
                    Kd[env_index[env][0]][nuclide_index[rn]]=rec[0][col_index["tf"]]*rec[0][col_index["conversion"]]
                else:
                    errors=write_input_data_error(prot_file,sql_ask,errors,"db",prt)
                    Kd[env_index[env][0]][nuclide_index[rn]]=0

    prot_file.write(prt_function_string+" end"+endline+endline)
    return Kd, errors

def nucl_rep_name(nuclide_index,hypten,prt=False):
    index_new={}
    for n in nuclide_index:
        if prt:
            print(n,nuclide_index[n])
        index_new[get_nucl_name(str(n),hypten)] = nuclide_index[n]
    if prt:
        print("ni:\n",nuclide_index)
    return index_new   

def read_activity(prot_file,activity_file_name,element_index,activity_index,nuclide_index,errors,prt=False):
    prt_function_string="read activities"
    prot_file.write(prt_function_string+" begin"+endline)
    with open(activity_file_name) as bioconf:
        data = yaml.load(bioconf, Loader=yaml.FullLoader)
    if prt:
        print(data["times"])
        print(len(data["times"]))
        print(data.values())
    
    times=[]
    for time in data["times"]:
        times.append(time)

    activities=np.zeros((len(element_index),len(activity_index),len(nuclide_index),len(times)),dtype=float)


    for element in element_index:
        if prt:
             print(element)
        for activity in activity_index:
            for nuclide in nuclide_index:
                for i in range(len(times)):
                    try:
                        activities[element_index[element]][activity_index[activity][0]][nuclide_index[nuclide]][i]=\
                        data["elements"][str(element)][str(activity)][str(nuclide)][i]
                    except KeyError:
                        errorstring="element: "+str(element)+", activity type: "+str(activity)+", nuclide: "+str(nuclide)+", time: "+str(times[i])
                        errors=write_input_data_error(prot_file,errorstring,errors,"act",prt)
                        activities[element_index[element]][activity_index[activity][0]][nuclide_index[nuclide]][i]=0
                    if prt:
                        print(element_index[element],activity_index[str(activity)],nuclide_index[str(nuclide)],i)
                        print(element,i,element_index[element],activity_index[activity],nuclide_index[nuclide])

    if prt:
        print("data elements:",len(data["elements"]))
        print(data["elements"])
    
        for element in data["elements"]:
            #print("element:",element,"elementype:",type(element))
            pass
    
        print("activities:")
        print(activities)
        print("activities end")
        
    bioconf.close()
    prot_file.write(prt_function_string+" end"+endline+endline)
    return times, activities, errors

def read_activity_from_gw(prot_file,activity_file_name,element_index,activity_index,nuclide_index,
                          env_params_index,env_params,Kd,errors,prt=False):
    prt_function_string="read activities from groundwater"
    prot_file.write(prt_function_string+" begin"+endline)
    with open(activity_file_name) as bioconf:
        data = yaml.load(bioconf, Loader=yaml.FullLoader)
    if prt:
        print(data["times"])
        print(len(data["times"]))
        print(data.values())
    
    times=[]
    for time in data["times"]:
        times.append(time)

    activities=np.zeros((len(element_index),len(activity_index),len(nuclide_index),len(times)),dtype=float)


    for element in element_index:
        if prt:
             print(element)
        activity=groundwater
        if prt:
            print("activity:",activity,str(activity))
        for nuclide in nuclide_index:
            for i in range(len(times)):
                if prt:
                    print(data)
                    print("elements")
                    print(str(element))
                    print(str(activity))
                    print(str(nuclide))
                    print(i)
                    print(data["elements"])
                    print(data["elements"][str(element)])
                    print(data["elements"][str(element)][str(activity)])
                    try:
                        print(data["elements"][str(element)][str(activity)][str(nuclide)][i])
                    except KeyError:
                        print("error")
                    #print(error[0])
                try:
                    activities[element_index[element]][activity_index[activity][0]][nuclide_index[nuclide]][i]=\
                    data["elements"][str(element)][str(activity)][str(nuclide)][i]
                except KeyError:
                    errorstring="element: "+str(element)+", activity type: "+str(activity)+", nuclide: "+str(nuclide)+", time: "+str(times[i])
                    errors=write_input_data_error(prot_file,errorstring,errors,"act",prt)
                    activities[element_index[element]][activity_index[activity][0]][nuclide_index[nuclide]][i]=0
                    if prt:
                        print(errorstring)
                if prt:
                    #print("d2x",d2)
                    #print(activity[element_index[str(element)]][activity_index[str(activity)]][nuclide_index[str(nuclide)]][i])
                    print(element_index[element],activity_index[str(activity)],nuclide_index[str(nuclide)],i)
                    print(element,i,element_index[element],activity_index[activity],nuclide_index[nuclide])
        
        for activity in activity_index:
            if activity!=groundwater:
                for nuclide in nuclide_index:
                    for i in range(len(times)):
                        if activity==surfacewater:
                            activities[element_index[element]][activity_index[activity][0]][nuclide_index[nuclide]][i]=\
                            activities[element_index[element]][activity_index[groundwater][0]][nuclide_index[nuclide]][i]/\
                            env_params[activity_index[activity][0]][env_params_index["dilution"][0]]
                        else:
                            activities[element_index[element]][activity_index[activity][0]][nuclide_index[nuclide]][i]=\
                            activities[element_index[element]][activity_index[groundwater][0]][nuclide_index[nuclide]][i]/\
                            env_params[activity_index[activity][0]][env_params_index["dilution"][0]]*\
                            (
                            Kd[activity_index[activity][0]][nuclide_index[nuclide]]*env_params[activity_index[activity][0]][env_params_index["density"][0]]
                            +env_params[activity_index[activity][0]][env_params_index["porosity"][0]]*env_params[activity_index[activity][0]][env_params_index["humidity"][0]]
                            )/\
                            env_params[activity_index[activity][0]][env_params_index["density"][0]]
                        if prt:
                            #print("d2x",d2)
                            #print(activity[element_index[str(element)]][activity_index[str(activity)]][nuclide_index[str(nuclide)]][i])
                            print(element_index[element],activity_index[str(activity)],nuclide_index[str(nuclide)],i)
                            print(element,i,element_index[element],activity_index[activity],nuclide_index[nuclide])

    if prt:
        print("data elements:",len(data["elements"]))
        print(data["elements"])
    
        for element in data["elements"]:
            #print("element:",element,"elementype:",type(element))
            #for act in element["groundwater_activity"]:
            pass
    
        print("activities:")
        print(activities)
        print("activities end")
        
    bioconf.close()
    prot_file.write(prt_function_string+" end"+endline+endline)
    return times, activities, errors


def read_activity_from_water(prot_file,activity_file_name,element_index,activity_index,nuclide_index,
                          env_params_index,env_params,Kd,act_conf,errors,prt=False):
    prt_function_string="read activities from groundwater"
    prot_file.write(prt_function_string+" begin"+endline)
    with open(activity_file_name) as bioconf:
        data = yaml.load(bioconf, Loader=yaml.FullLoader)
    if prt:
        print(data["times"])
        print(len(data["times"]))
        print(data.values())
    
    times=[]
    for time in data["times"]:
        times.append(time)

    activities=np.zeros((len(element_index),len(activity_index),len(nuclide_index),len(times)),dtype=float)


    for element in element_index:
        if prt:
             print(element)
             
        #read groundwater activity     
        activity=groundwater
        if prt:
            print("activity:",activity,str(activity))
        for nuclide in nuclide_index:
            for i in range(len(times)):
                if prt:
                    print(data)
                    print("elements")
                    print(str(element))
                    print(str(activity))
                    print(str(nuclide))
                    print(i)
                    print(data["elements"])
                    print(data["elements"][str(element)])
                    print(data["elements"][str(element)][str(activity)])
                    try:
                        print(data["elements"][str(element)][str(activity)][str(nuclide)][i])
                    except KeyError:
                        print("error")
                    #print(error[0])
                try:
                    activities[element_index[element]][activity_index[activity][0]][nuclide_index[nuclide]][i]=\
                    data["elements"][str(element)][str(activity)][str(nuclide)][i]
                except KeyError:
                    errorstring="element: "+str(element)+", activity type: "+str(activity)+", nuclide: "+str(nuclide)+", time: "+str(times[i])
                    errors=write_input_data_error(prot_file,errorstring,errors,"act",prt)
                    activities[element_index[element]][activity_index[activity][0]][nuclide_index[nuclide]][i]=0
                    if prt:
                        print(errorstring)
                if prt:
                    #print("d2x",d2)
                    #print(activity[element_index[str(element)]][activity_index[str(activity)]][nuclide_index[str(nuclide)]][i])
                    print(element_index[element],activity_index[str(activity)],nuclide_index[str(nuclide)],i)
                    print(element,i,element_index[element],activity_index[activity],nuclide_index[nuclide])
        
        #get surfacewater activity
        activity=surfacewater
        if act_conf=="water":
            for nuclide in nuclide_index:
                for i in range(len(times)):
                    try:
                        activities[element_index[element]][activity_index[activity][0]][nuclide_index[nuclide]][i]=\
                        data["elements"][str(element)][str(activity)][str(nuclide)][i]
                    except KeyError:
                        errorstring="element: "+str(element)+", activity type: "+str(activity)+", nuclide: "+str(nuclide)+", time: "+str(times[i])
                        errors=write_input_data_error(prot_file,errorstring,errors,"act",prt)
                        activities[element_index[element]][activity_index[activity][0]][nuclide_index[nuclide]][i]=0
                        if prt:
                            print(errorstring)
        else:
            for nuclide in nuclide_index:
                for i in range(len(times)):
                    activities[element_index[element]][activity_index[activity][0]][nuclide_index[nuclide]][i]=\
                    activities[element_index[element]][activity_index[groundwater][0]][nuclide_index[nuclide]][i]/\
                    env_params[activity_index[activity][0]][env_params_index["dilution"][0]]
        
        #get soils activities
        #activity=surfacewater
        for activity in activity_index:
            if activity!=groundwater and activity!=surfacewater:
                for nuclide in nuclide_index:
                    for i in range(len(times)):
                        activities[element_index[element]][activity_index[activity][0]][nuclide_index[nuclide]][i]=\
                        activities[element_index[element]][activity_index[surfacewater][0]][nuclide_index[nuclide]][i]*\
                        (
                        Kd[activity_index[activity][0]][nuclide_index[nuclide]]*env_params[activity_index[activity][0]][env_params_index["density"][0]]
                        +env_params[activity_index[activity][0]][env_params_index["porosity"][0]]*env_params[activity_index[activity][0]][env_params_index["humidity"][0]]
                        )/\
                        env_params[activity_index[activity][0]][env_params_index["density"][0]]
                        if prt:
                            #print("d2x",d2)
                            #print(activity[element_index[str(element)]][activity_index[str(activity)]][nuclide_index[str(nuclide)]][i])
                            print(element_index[element],activity_index[str(activity)],nuclide_index[str(nuclide)],i)
                            print(element,i,element_index[element],activity_index[activity],nuclide_index[nuclide])

    if prt:
        print("data elements:",len(data["elements"]))
        print(data["elements"])
    
        for element in data["elements"]:
            #print("element:",element,"elementype:",type(element))
            #for act in element["groundwater_activity"]:
            pass
    
        print("activities:")
        print(activities)
        print("activities end")
        
    bioconf.close()
    prot_file.write(prt_function_string+" end"+endline+endline)
    return times, activities, errors



def write_input_data_error(prot_file,sql_ask,errors,error_type,prt=False):
    errors+=1
    if prt:
        print("error no.",errors)
    errorstring="data error no. "+str(errors)+":"+endline
    if error_type=="db":
        errorstring+="the QSL ask: \""+sql_ask+"\" is NULL"+endline
        errorstring+="the parameter value set to  0."+endline
    if error_type=="act":
        errorstring+="acitity/concentration for next items not exist in activity input file:"+endline
        errorstring+="the record: \""+sql_ask+"\" is NULL"+endline
        errorstring+="the parameter value set to  0."+endline
    if error_type=="bas_cons":
        errorstring+="the QSL ask: \""+sql_ask+"\" is NULL"+endline
        errorstring+="the parameter value set to same as for 'basic' basket"+endline
    if error_type=="basket":
        errorstring+="the basket: '"+sql_ask+"' is not in the database"+endline
        errorstring+="the basket name in all languages set same as basket identificator"+endline
    if error_type=="nucl_mass":
        errorstring+="the QSL ask: "+sql_ask+endline
    prot_file.write(errorstring)
    list_data_errors.append(errorstring+endline)
    return errors
        

def get_db_ask_column_names_index(br_db_cur,sql_ask,prt=False):
    col_index = {}

    br_db_cur.execute(sql_ask)
    col_names=br_db_cur.description
    if prt:
        print(col_names)

    i=-1
    for row in col_names:
        if prt:
            print(row[0])
        i+=1
        col_index[str(row[0])]=i

    return col_index


def get_db_pathways3(br_db_cur,path_type,prt=False):
    if prt:
        print("get_db_pathways3:")
    br_path_index={}
    pathway_names=[0,'','']
    
    br_db_cur.execute("SELECT identificator FROM Pathways where "+path_type+"=1")
    rows = br_db_cur.fetchall()
    i=-1
    for row in rows:
        i+=1
        path=row[0]
        pathway_names[0]=i
        br_db_cur.execute("SELECT name_cz FROM Pathways where identificator='"+str(path)+"'")
        pathway_names[1]=br_db_cur.fetchall()[0][0]
        br_db_cur.execute("SELECT name_en FROM Pathways where identificator='"+str(path)+"'")
        pathway_names[2]=br_db_cur.fetchall()[0][0]
        br_path_index[path]=pathway_names.copy()
        if prt:
            print(row)
            print(row[0],i)
            print(pathway_names)
    
    if prt:
        print("len len:",len(br_path_index),len(pathway_names))
        print(pathway_names)
        print(br_path_index)
        print("get_db_pathways3 end")
    return br_path_index


def get_baskets_consumptions(prot_file,br_db_cur,basket_index,path_index,errors,prt=False):
    prt_function_string="read consumptions by baskets"
    prot_file.write(prt_function_string+" begin"+endline)
    bas_comps=np.zeros((len(basket_index),len(path_index)),dtype=float)
    sql_ask="SELECT * FROM w_baskets"
    col_index=get_db_ask_column_names_index(br_db_cur,sql_ask,prt)
    for b in basket_index:
        for p in path_index:
            sql_ask="SELECT * FROM w_baskets where basket_identificator='"+str(b)+"' and prod_identificator='"+str(p)+"'"
            if prt:
                print(sql_ask)
                print(basket_index[b])
                print(path_index[p][0])
            br_db_cur.execute(sql_ask)
            #rec=br_db_cur.fetchall()[0]
            rec=br_db_cur.fetchall()
            if len(rec)>0:
                rec=rec[0]
                bas_comps[basket_index[b]][path_index[p][0]]=rec[col_index["value"]]*rec[col_index["conversion"]]
            else:
                #errors=write_input_data_error(prot_file,sql_ask,errors,"db",prt)
                errors=write_input_data_error(prot_file,sql_ask,errors,"bas_cons",prt)
                #bas_comps[basket_index[b]][path_index[p][0]]=0
                sql_ask="SELECT * FROM w_baskets where prod_identificator = '"+str(p)+"' AND basket_identificator = 'basic'"
                br_db_cur.execute(sql_ask)
                rec=br_db_cur.fetchall()
                bas_comps[basket_index[b]][path_index[p][0]]=rec[0][col_index["value"]]*rec[0][col_index["conversion"]]

    prot_file.write(prt_function_string+" end"+endline+endline)
    return bas_comps, errors


def get_basket_index(prot_file,br_db_cur,errors,prt=False):
    if prt:
        print("get_basket_index:")
    br_basket_index={}
    basket_names=[0,'','']
    
    br_db_cur.execute("SELECT identificator FROM Basket_type")
    rows = br_db_cur.fetchall()
    i=-1
    for row in rows:
        i+=1
        path=row[0]
        basket_names[0]=i
        br_db_cur.execute("SELECT name_cz FROM Basket_type where identificator='"+str(path)+"'")
        basket_names[1]=br_db_cur.fetchall()[0][0]
        br_db_cur.execute("SELECT name_en FROM Basket_type where identificator='"+str(path)+"'")
        basket_names[2]=br_db_cur.fetchall()[0][0]
        br_basket_index[path]=basket_names.copy()
        if prt:
            print(row)
            print(row[0],i)
            print(basket_names)
    
    if prt:
        print("len len:",len(br_basket_index),len(basket_names))
        print(basket_names)
        print(br_basket_index)
        print("get_db_pathways3 end")
    return br_basket_index, errors


def get_a2s(br_db_cur):
    br_db_cur.execute("SELECT conversion FROM Unit where unit='a'")
    return(br_db_cur.fetchall()[0][0])

def get_Na(br_db_cur,prt=False):
    sql_ask="SELECT * FROM w_constants where identificator='Na'"
    col_index=get_db_ask_column_names_index(br_db_cur,sql_ask,prt)
    br_db_cur.execute(sql_ask)
    rec=br_db_cur.fetchall()[0]
    return rec[col_index["value"]]*rec[col_index["conversion"]]

def get_consumption(prot_file,br_db_cur,br_basket_index,path_name,errors,prt=False):
    prt_function_string="read consumptions for '"+path_name+" by baskets"
    prot_file.write(prt_function_string+" begin"+endline)
    consumpt_ing=np.zeros(len(br_basket_index),dtype=float)
    sql_ask="SELECT * FROM w_baskets"
    col_index=get_db_ask_column_names_index(br_db_cur,sql_ask,prt)
    
    if prt:
        print("basket index",br_basket_index)
    i=-1
    for basket in br_basket_index:
        i+=1
        sql_ask="SELECT * FROM w_baskets where prod_identificator = '"+path_name+"' AND basket_identificator = '"+str(basket)+"'" 
        if prt:
            print(sql_ask)
        br_db_cur.execute(sql_ask)
        rec=br_db_cur.fetchall()
        if len(rec)>0:
            consumpt_ing[i]=rec[0][col_index["value"]]*rec[0][col_index["conversion"]]
        else:
            errors=write_input_data_error(prot_file,sql_ask,errors,"bas_cons",prt)
            sql_ask="SELECT * FROM w_baskets where prod_identificator = '"+path_name+"' AND basket_identificator = 'basic'"
            br_db_cur.execute(sql_ask)
            rec=br_db_cur.fetchall()
            consumpt_ing[i]=rec[0][col_index["value"]]*rec[0][col_index["conversion"]]
    if prt:
        print("consumption ",path_name,":")
        print(consumpt_ing)
        print(path_name,"end")
    prot_file.write(prt_function_string+" end"+endline)
    return consumpt_ing, errors

def get_plant_consumption(prot_file,br_db_cur,br_basket_index,plant_index,errors,prt=False):
    prt_function_string="get plants consumptions"
    prot_file.write(prt_function_string+" begin"+endline)
    consumpt_ing=np.zeros((len(br_basket_index),len(plant_index)),dtype=float)
    for plant in plant_index:
        if prt:
            print(plant,plant[2])
        if plant_index[plant][2]:
            #consumpt_ing[plant_index[plant][0]]=get_consumption(br_db_cur,br_basket_index,str(plant),prt)
            consumpt, errors=get_consumption(prot_file,br_db_cur,br_basket_index,str(plant),errors,prt)
            for i in range(len(consumpt)):
                consumpt_ing[i][plant_index[plant][0]]=consumpt[i]
    prot_file.write(prt_function_string+" end"+endline+endline)
    return consumpt_ing, errors

def get_animal_prod_consumption(prot_file,br_db_cur,br_basket_index,animal_prod_index,errors,prt=False):
    prt_function_string="get animal products consumptions"
    prot_file.write(prt_function_string+" begin"+endline)
    consumpt_ing=np.zeros((len(br_basket_index),len(animal_prod_index)),dtype=float)
    for ap in animal_prod_index:
        consumpt, errors=get_consumption(prot_file,br_db_cur,br_basket_index,str(ap),errors,prt)
        for i in range(len(consumpt)):
                consumpt_ing[i][animal_prod_index[ap][0]]=consumpt[i]
    prot_file.write(prt_function_string+" end"+endline+endline)
    return consumpt_ing, errors


def get_water_prod_consumption(prot_file,br_db_cur,br_basket_index,wat_prod_index,errors,prt=False):
    prt_function_string="get water products consumptions"
    prot_file.write(prt_function_string+" begin"+endline)
    consumpt_ing=np.zeros((len(br_basket_index),len(wat_prod_index)),dtype=float)
    for wp in wat_prod_index:
        #consumpt_ing[plant_index[plant][0]]=get_consumption(br_db_cur,br_basket_index,str(plant),prt)
        consumpt,errors=get_consumption(prot_file,br_db_cur,br_basket_index,str(wp),errors,prt)
        for i in range(len(consumpt)):
            consumpt_ing[i][wat_prod_index[wp]]=consumpt[i]
    prot_file.write(prt_function_string+" end"+endline+endline)
    return consumpt_ing, errors


def get_hing(prot_file,br_db_cur,br_nuclide_index,errors,prt=False):
    prt_function_string="get dose conversion factors for ingestions"
    prot_file.write(prt_function_string+" begin"+endline)
    sql_ask="SELECT * FROM w_hing_max"
    col_index=get_db_ask_column_names_index(br_db_cur,sql_ask,prt)
    if prt:
        print("br_nuclide_index=",br_nuclide_index)
    hing=np.zeros(len(br_nuclide_index),dtype=float)
    for nuclide in br_nuclide_index:
        if prt:
            print("nuclide=",nuclide)
        
        sql_ask="SELECT * FROM w_hing_max where nuclide='"+nuclide+"'"
        br_db_cur.execute(sql_ask)
        rec=br_db_cur.fetchall()
        if len(rec)>0:
            rec=rec[0]
            hing[br_nuclide_index[nuclide]]=rec[col_index["dcf_hing_adult"]]*rec[col_index["conversion"]]
        else:
            errors=write_input_data_error(prot_file,sql_ask,errors,"db",prt)
            hing[br_nuclide_index[nuclide]]=0
            
    if prt:
        print(hing)
    prot_file.write(prt_function_string+" end"+endline+endline)
    return hing, errors


def get_hinh(prot_file,br_db_cur,br_nuclide_index,errors,prt=False):
    prt_function_string="get dose conversion factors for inhalation"
    prot_file.write(prt_function_string+" begin"+endline)
    sql_ask="SELECT * FROM w_hinh_max"
    col_index=get_db_ask_column_names_index(br_db_cur,sql_ask,prt)
    if prt:
        print("br_nuclide_index=",br_nuclide_index)
    hinh=np.zeros(len(br_nuclide_index),dtype=float)
    for nuclide in br_nuclide_index:
        if prt:
            print("nuclide=",nuclide)
        
        sql_ask="SELECT * FROM w_hinh_max where nuclide='"+nuclide+"'"
        br_db_cur.execute(sql_ask)
        rec=br_db_cur.fetchall()
        if len(rec)>0:
            rec=rec[0]
            hinh[br_nuclide_index[nuclide]]=rec[col_index["dcf_hinh_max_adult"]]*rec[col_index["conversion"]]
        else:
            errors=write_input_data_error(prot_file,sql_ask,errors,"db",prt)
            hinh[br_nuclide_index[nuclide]]=0
    if prt:
        print(hinh)
    return hinh, errors

def get_hext(prot_file,br_db_cur,ext_ident,br_nuclide_index,errors,prt=False):
    prot_file.write("get external irradiation dose coefficients for "+ext_ident+" begin"+endline)
    if prt:
        print("br_nuclide_index=",br_nuclide_index)
    hext=np.zeros(len(br_nuclide_index),dtype=float)
    sql_ask="SELECT * FROM w_hext"
    col_index=get_db_ask_column_names_index(br_db_cur,sql_ask,prt)
    br_db_cur.execute(sql_ask)

    for nuclide in br_nuclide_index:
        sql_ask="SELECT * FROM w_hext where identificator='"+ext_ident+"' and nuclide='"+nuclide+"'"
        if prt:
            print("nuclide=",nuclide)
            print(sql_ask)
        br_db_cur.execute(sql_ask)
        rec=br_db_cur.fetchall()
        if len(rec)>0:
            rec=rec[0]
            hext[br_nuclide_index[nuclide]]=rec[col_index["effective"]]*rec[col_index["conversion"]]
        else:
            errors=write_input_data_error(prot_file,sql_ask,errors,"db",prt)
            hext[br_nuclide_index[nuclide]]=0
            
    if prt:
        print(hext)
    prot_file.write("get external irradiation dose coefficients for "+ext_ident+" end"+endline+endline)
    return hext, errors


def get_plants_index(prot_file,br_db_cur,prt=False):

    index_rec=[None,'',False]
    index = {}
    
    sql_ask="SELECT * FROM s_plant_products"
    col_index=get_db_ask_column_names_index(br_db_cur,sql_ask,prt)
    br_db_cur.execute(sql_ask)
    rows = br_db_cur.fetchall()

    i=-1
    for row in rows:
        i+=1
        if prt:
            print(row,type(row))
        index_rec[0]=i
        if row[col_index["type_plant"]]==plant_env:
            env=soil
        else:
            env=wood
        index_rec[1]=env
        
        if row[col_index["super"]]==plant_hum_ingestion:
            index_rec[2]=True
            
        index[row[col_index["sub"]]]=index_rec.copy()

    return index

def get_resuspension_factor(prot_file,br_db_cur,prt=False):
    lfi="leafy_vegetables"

    sql_ask="SELECT * FROM w_parameters"
    col_index=get_db_ask_column_names_index(br_db_cur,sql_ask,prt)


    sql_ask="SELECT * FROM w_parameters where identificator='soil_field' and param_identificator='deposit'"
    br_db_cur.execute(sql_ask)
    rec = br_db_cur.fetchall()
    if prt:
        print("col:",col_index)
        print("val i:",col_index["value"])
        print("rec:",rec)
        print("rec val:",rec[0][col_index["value"]])
        print("conv:",rec[0][col_index["conversion"]])
    rf=rec[0][col_index["value"]]*rec[0][col_index["conversion"]]

    sql_ask="SELECT * FROM w_parameters where identificator='"+lfi+"' and param_identificator='vegetation'"
    br_db_cur.execute(sql_ask)
    rec = br_db_cur.fetchall()
    rf*=rec[0][col_index["value"]]*rec[0][col_index["conversion"]]
    
    sql_ask="SELECT * FROM w_parameters where identificator='"+lfi+"' and param_identificator='yeld'"
    br_db_cur.execute(sql_ask)
    rec = br_db_cur.fetchall()
    rf/=(rec[0][col_index["value"]]*rec[0][col_index["conversion"]])

    if prt:
        print("resuspension_factor:",rf)
    return rf

def get_animal_product_index(prot_file,br_db_cur,prt=False):

    index_rec=[None,'']
    index = {}
    
    sql_ask="SELECT * FROM s_animal_products"
    col_index=get_db_ask_column_names_index(br_db_cur,sql_ask,prt)
    br_db_cur.execute(sql_ask)
    rows = br_db_cur.fetchall()

    i=-1
    for row in rows:
        i+=1
        if prt:
            print("xxx",row,type(row))
        index_rec[0]=i
        index_rec[1]=row[col_index["animal"]]
        
        index[row[col_index["sub"]]]=index_rec.copy()

    return index

def get_animal_index(prot_file,br_db_cur,prt=False):
    index = {}
    
    sql_ask="SELECT animal FROM s_animals"
    br_db_cur.execute(sql_ask)
    rows = br_db_cur.fetchall()

    i=-1
    for row in rows:
        i+=1
        if prt:
            print("row ai:",row,type(row))
            print("row ai:",row[0],type(row[0]))
        #index[row]=i
        index[row[0]]=i

    return index

def get_animal_consumption_index(prot_file,br_db_cur,prt=False):
    index = {}
    
    sql_ask="select identificator from parameter where identificator like '%consump%'"
    br_db_cur.execute(sql_ask)
    rows = br_db_cur.fetchall()

    i=-1
    for row in rows:
        i+=1
        if prt:
            print(row,type(row))
        #index[row]=i
        index[row[0]]=i

    return index

def get_animal_consumption(prot_file,br_db_cur,animal_index,animal_consump_index,prt=False):
    a_consumption=np.zeros((len(animal_index),len(animal_consump_index)),dtype=float)
    
    for ai in animal_index:
        for aci in animal_consump_index:
            #sql_ask="select value, conversion from w_parameters where identificator='"+str(ai[0])+"' and param_identificator='"+str(aci[0])+"'"
            sql_ask="select value, conversion from w_parameters where identificator='"+ai+"' and param_identificator='"+aci+"'"
            if prt:
                print(sql_ask)
            br_db_cur.execute(sql_ask)
            row = br_db_cur.fetchall()[0]
            if prt:
                print(row)
            a_consumption[animal_index[ai]][animal_consump_index[aci]]=row[0]*row[1]
    
    return a_consumption

def get_animal_consumption_acitivies(prot_file,element_index,animal_index,nuclide_index,times,
                                     animal_consump_index,animal_consump,activity_index,activities,plant_index,plant_tff,env_params_index,env_params,prt=False):
    aa_consumption=np.zeros((len(element_index),len(animal_index),len(nuclide_index),len(times)),dtype=float)
    for e in range(len(element_index)):    
        for ai in animal_index:        
            for ni in nuclide_index:
                for t in range(len(times)):
                    if prt:
                        print("xxx0",ni,nuclide_index[ni],activity_index[groundwater][0])
                        print("xxx1:",activity_index)
                        #print("xxx2:",activity_index["soil_forest"])
                        print("xxx3:",animal_index)
                        print("xxx4:",ai,animal_index[ai])
                        #print("xxx5:",animal_index["pig"])
                        print("xxx6:",animal_consump_index)
                        #print("xxx7:",animal_consump_index["consump_water"])
                        #print("xxx8:",animal_consump[animal_index[ai]][animal_consump_index["consump_water"]])
                    #water:
                    aa_consumption[e][animal_index[ai]][nuclide_index[ni]][t]+=\
                    activities[e][activity_index[groundwater][0]][nuclide_index[ni]][t]*\
                    animal_consump[animal_index[ai]][animal_consump_index["consump_water"]]
                    #soil
                    aa_consumption[e][animal_index[ai]][nuclide_index[ni]][t]+=\
                    activities[e][activity_index[soil][0]][nuclide_index[ni]][t]*\
                    animal_consump[animal_index[ai]][animal_consump_index["consump_soil"]]
                    #pasture:
                    aa_consumption[e][animal_index[ai]][nuclide_index[ni]][t]+=\
                    activities[e][activity_index[soil][0]][nuclide_index[ni]][t]*plant_tff[plant_index["feed"][0]][nuclide_index[ni]]*\
                    animal_consump[animal_index[ai]][animal_consump_index["consump_feed"]]
                     #air:
                    aa_consumption[e][animal_index[ai]][nuclide_index[ni]][t]+=\
                    activities[e][activity_index[soil][0]][nuclide_index[ni]][t]*env_params[activity_index[soil][0]][env_params_index["dust_level"][0]]*\
                    animal_consump[animal_index[ai]][animal_consump_index["consump_air"]]
    
    return aa_consumption


def get_wat_prod_index(prot_file,br_db_cur,prt=False):

    index = {}
    
    sql_ask="SELECT * FROM s_water_products"
    col_index=get_db_ask_column_names_index(br_db_cur,sql_ask,prt)
    br_db_cur.execute(sql_ask)
    rows = br_db_cur.fetchall()

    i=-1
    for row in rows:
        i+=1
        index[row[col_index["sub"]]]=i

    return index

def get_dry_matter(prot_file,br_db_cur,plant_id,prt=False):
    prot_file.write("get "+plant_id+" dry matter begin"+endline)
    br_db_cur.execute("SELECT value FROM w_parameters where identificator='"+plant_id+"'")
    dm=br_db_cur.fetchall()[0][0]
    br_db_cur.execute("SELECT conversion FROM w_parameters where identificator='"+plant_id+"'")
    dm*=br_db_cur.fetchall()[0][0]
    prot_file.write("get "+plant_id+" dry matter begin"+endline)
    return dm


def get_plant_tf_fresh(prot_file,br_db_cur,plant_index,br_nuclide_index,errors,prt=False):
    prot_file.write("get plants transfer factors begin"+endline)
    ptf=np.zeros((len(plant_index),len(br_nuclide_index)),dtype=float)
    sql_ask="SELECT * FROM w_tf"
    col_index=get_db_ask_column_names_index(br_db_cur,sql_ask)
    if prt:
        print("col_index['tf']",col_index["tf"])        
    
    for plant in plant_index:
        dry_matter=get_dry_matter(prot_file,br_db_cur,plant,prt)
        prot_file.write("dry matter "+plant+" "+str(dry_matter)+endline)
        for nuclide in br_nuclide_index:
            elem=get_ch_el_name(nuclide)
            #br_db_cur.execute("SELECT * FROM w_tf where element='"+elem+"' and identificator='"+plant+"'")
            sql_ask="SELECT * FROM w_tf where element='"+elem+"' and identificator='"+plant+"'"
            br_db_cur.execute(sql_ask)
            row = br_db_cur.fetchall()
            if len(row)>0:
                value=row[0][col_index['tf']]
                conversion=row[0][col_index["conversion"]]
                dm=row[0][col_index["dry_matter"]]
                if  (dm==None or dm==0):
                    dm=dry_matter
            else:
                errors=write_input_data_error(prot_file,sql_ask,errors,"db",prt)
                value=0
                conversion=0
                dm=0
            ptf[plant_index[plant][0],br_nuclide_index[nuclide]]=value*conversion*dm

    if prt:
        print(ptf)
    prot_file.write("get plants transfer factors end"+endline)
    return ptf, errors

def get_wat_prod_tf(prot_file,br_db_cur,wat_prod_index,br_nuclide_index,errors,prt=False):
    prot_file.write("get water products transfer factors begin"+endline)
    wptf=np.zeros((len(wat_prod_index),len(br_nuclide_index)),dtype=float)
    sql_ask="SELECT * FROM w_tf"
    col_index=get_db_ask_column_names_index(br_db_cur,sql_ask)
    if prt:
        print("col_index['tf']",col_index["tf"])        
    
    for wp in wat_prod_index:
        for nuclide in br_nuclide_index:
            elem=get_ch_el_name(nuclide)
            sql_ask="SELECT * FROM w_tf where element='"+elem+"' and identificator='"+wp+"'"
            #br_db_cur.execute("SELECT * FROM w_tf where element='"+elem+"' and identificator='"+wp+"'")
            br_db_cur.execute(sql_ask)
            row = br_db_cur.fetchall()
            if len(row)>0:
                value=row[0][col_index['tf']]
                conversion=row[0][col_index["conversion"]]
            else:
                errors=write_input_data_error(prot_file,sql_ask,errors,"db",prt)
                value=0
                conversion=0
            wptf[wat_prod_index[wp],br_nuclide_index[nuclide]]=value*conversion

    if prt:
        print(wptf)
    prot_file.write("get water products transfer factors end"+endline)
    return wptf, errors

def get_animal_prod_tf(prot_file,br_db_cur,animal_prod_index,br_nuclide_index,errors,prt=False):
    prot_file.write("get water products transfer factors begin"+endline)
    aptf=np.zeros((len(animal_prod_index),len(br_nuclide_index)),dtype=float)
    sql_ask="SELECT * FROM w_tf"
    col_index=get_db_ask_column_names_index(br_db_cur,sql_ask)
    if prt:
        print("col_index['tf']",col_index["tf"])        
    
    for ap in animal_prod_index:
        for nuclide in br_nuclide_index:
            elem=get_ch_el_name(nuclide)
            if prt:
                print(ap)
            sql_ask="SELECT * FROM w_tf where element='"+elem+"' and identificator='"+str(ap)+"'"
            #br_db_cur.execute("SELECT * FROM w_tf where element='"+elem+"' and identificator='"+str(ap)+"'")
            br_db_cur.execute(sql_ask)
            row = br_db_cur.fetchall()
            if len(row):
                value=row[0][col_index['tf']]
                conversion=row[0][col_index["conversion"]]
            else:
                errors=write_input_data_error(prot_file,sql_ask,errors,"db",prt)
                value=0
                conversion=0
            aptf[animal_prod_index[ap][0],br_nuclide_index[nuclide]]=value*conversion

    if prt:
        print(aptf)
    prot_file.write("get water products transfer factors end"+endline)
    return aptf, errors


def collect_basket_index(prot_file,basket_index,basket_index_db,errors,prt=False):
    prt_function_string="get baskets index"
    prot_file.write(prt_function_string+" begin"+endline)
    if prt:
        print("collect_basket_index:")
    index={}
    index_names=[0,'','']

    i=-1        
    for bi in basket_index:
        i+=1
        #index_names[0]=bi[0]
        index_names[0]=basket_index[str(bi)]
        if bi in basket_index_db:
            index_names[1]=basket_index_db[str(bi)][1]
            index_names[2]=basket_index_db[str(bi)][2]
        else:
            errors=write_input_data_error(prot_file,str(bi),errors,"basket",prt)
            index_names[1]=index_names[2]=str(bi)
            
        index[str(bi)]=index_names.copy()
    if prt:
        print(index)
        print("collect_basket_index end")
    prot_file.write(prt_function_string+" end"+endline+endline)
    return index,errors


def doses_calc3(prot_file,msh_el_index,basket_index,nuclide_index,path_index,times,activities,activity_index,
               env_params_index,env_params,
               consumption_index,consumption_all,
               hing,hinh,hext_soil,hext_air,hext_water,wat_ing,soil_ing,
               plant_index,plant_tff,plant_ing,lf_resp_factor,
               animal_prod_index,animal_prod_ing,animal_index,animal_prod_tf,animal_consumption_activity,
               wat_prod_index,wat_prod_tf,wat_prod_ing,
               a2s=1,prt=False):
    doses=np.zeros((len(msh_el_index),len(basket_index),len(nuclide_index),len(path_index),len(times)), dtype=float)
    
    for el in msh_el_index:
        i=msh_el_index[el]
        for bask in basket_index:
            j=basket_index[bask]
            for rn in nuclide_index:
                k=nuclide_index[rn]
                for m in range(len(times)):
                    doses[i,j,k,path_index[wat_ing_path][0],m]=activities[i][activity_index[groundwater][0]][k][m]*hing[k]*wat_ing[j]
                    doses[i,j,k,path_index[soil_ing_path][0],m]=activities[i][activity_index[soil][0]][k][m]*hing[k]*soil_ing[j]
                    
                    #ingesce rostlinnych produktu 
                    for pi in plant_index:
                        l=plant_index[pi][0]
                        if plant_index[pi][2]:
                            doses[i,j,k,path_index[pi][0],m]=activities[i][activity_index[plant_index[pi][1]][0]][k][m]*hing[k]*plant_ing[j][l]*plant_tff[l][k]
                    
                    #ingesce z resuspenze prachu na listove zelenine:
                    pi="resuspension"
                    ei="soil_field"
                    l=plant_index["leafy_vegetables"][0]
                    doses[i,j,k,path_index[pi][0],m]=activities[i][activity_index[ei][0]][k][m]*hing[k]*plant_ing[j][l]*lf_resp_factor

                    if prt:
                        print(i,j,k,m)
                        print("activity:",activities[i][activity_index[groundwater][0]][k][m])
                        print("hing:",hing[k])
                        print("wat_ing:",wat_ing[j])
                        print(activities[i][activity_index[groundwater][0]][k][m]*hing[k]*wat_ing[j])
                        print(activities[i][activity_index[soil][0]][k][m]*hing[k]*soil_ing[j])

                    #ingesce zivocisnych produktu:
                    for ap in animal_prod_index:
                        l=animal_prod_index[ap][0]
                        animal=animal_prod_index[ap][1]
                        ai=animal_index[animal]
                        if prt:
                            print("animal:",l,animal,ai)
                            print(animal_prod_index)
                            print(animal_prod_index[ap])
                        doses[i,j,k,path_index[ap][0],m]=animal_prod_tf[l][k]*animal_consumption_activity[i][ai][k][m]*hing[k]*animal_prod_ing[j][l]
                    
                    #ingesce produktu vodniho hospodarstvi:  
                    for wp in wat_prod_index:
                        l=wat_prod_index[wp]
                        doses[i,j,k,path_index[wp][0],m]=activities[i][activity_index[surfacewater][0]][k][m]*hing[k]*wat_prod_ing[j][l]*wat_prod_tf[l][k]
                        
                    #inhalace
                    #inhalace na poli
                    pi="inh_field"
                    ot="sr_field"
                    doses[i,j,k,path_index[pi][0],m]=activities[i][activity_index[soil][0]][k][m]*\
                    env_params[activity_index[soil][0]][env_params_index["dust_level"][0]]*\
                    consumption_all[j][consumption_index["air"][0]]*consumption_all[j][consumption_index[ot][0]]*hinh[k]
                    #inhalace v lese
                    pi="inh_forest"
                    ot="sr_forest"
                    doses[i,j,k,path_index[pi][0],m]=activities[i][activity_index[wood][0]][k][m]*\
                    env_params[activity_index[wood][0]][env_params_index["dust_level"][0]]*\
                    consumption_all[j][consumption_index["air"][0]]*consumption_all[j][consumption_index[ot][0]]*hinh[k]
                    
                    #zevni ozareni
                    #zevni ozareni pole
                    pi="exp_field"
                    ot="sr_field"
                    env=soil
                    #zevni ozareni pole z pudy
                    doses[i,j,k,path_index[pi][0],m]=activities[i][activity_index[env][0]][k][m]/\
                    env_params[activity_index[env][0]][env_params_index["density"][0]]*\
                    consumption_all[j][consumption_index[ot][0]]*hext_soil[k]
                    #zevni ozareni pole ze vzduchu
                    doses[i,j,k,path_index[pi][0],m]+=activities[i][activity_index[env][0]][k][m]*\
                    env_params[activity_index[env][0]][env_params_index["dust_level"][0]]*\
                    consumption_all[j][consumption_index[ot][0]]*hext_air[k]
                    #zevni ozareni les
                    pi="exp_forest"
                    ot="sr_forest"
                    env=wood
                    #zevni ozareni les z pudy
                    doses[i,j,k,path_index[pi][0],m]=activities[i][activity_index[env][0]][k][m]/\
                    env_params[activity_index[env][0]][env_params_index["density"][0]]*\
                    consumption_all[j][consumption_index[ot][0]]*hext_soil[k]
                    #zevni ozareni les ze vzduchu
                    doses[i,j,k,path_index[pi][0],m]+=activities[i][activity_index[env][0]][k][m]*\
                    env_params[activity_index[env][0]][env_params_index["dust_level"][0]]*\
                    consumption_all[j][consumption_index[ot][0]]*hext_air[k]
                    #zevni ozareni voda
                    pi="exp_water"
                    ot="sr_water"
                    env=surfacewater
                    doses[i,j,k,path_index[pi][0],m]=activities[i][activity_index[env][0]][k][m]*\
                    consumption_all[j][consumption_index[ot][0]]*hext_water[k]

    
    doses*=a2s
    return doses


def get_sum_cross_matrix(prot_file,br_db_cur,path_index,path_index_sum,prt=False):
    prot_file.write("get sum cross matrix begin"+endline)
    prot_file.write(str(path_index)+endline)
    prot_file.write(str(path_index_sum)+endline)
    matrix=np.zeros((len(path_index_sum),len(path_index)), dtype=float)
    for psum in path_index_sum:
        prot_file.write(str(psum)+endline)
        if prt:
            print("path=",psum)
            print("path=",path_index_sum[psum][0])
        br_db_cur.execute("SELECT P2.identificator FROM Pathways_structure PS JOIN Pathways P1 ON PS.Pathways_1_id = P1.id JOIN Pathways P2 ON PS.Pathways_2_id = P2.id WHERE P1.identificator = '"+psum+"' and P2.calculated = 1")
        calc_pathways = br_db_cur.fetchall()
        for pcalc in calc_pathways:
            if prt:
                print("pcalc=",pcalc)
                print("pcalc[0]=",pcalc[0])
            prot_file.write(str(pcalc)+endline)
            matrix[path_index_sum[psum][0]][path_index[pcalc[0]][0]]=1
            if prt:
                print(path_index_sum[psum][1],path_index[pcalc[0]][1])
                continue
    prot_file.write("get sum cross matrix end"+endline)
    return matrix


def get_doses_all_sum(prot_file,doses,sum_matrix,nelems,nbaskets,nnuclides,path_index,ntimes,path_index_sum,prt=False):
    prot_file.write("doses all + sum"+endline)
    path_index_all=path_index.copy()
    npaths=len(path_index)
    for pi in path_index_sum:
        path_index_all[pi]=path_index_sum[pi]
        path_index_all[pi][0]+=npaths
        continue
    doses_all_sum=np.zeros((nelems,nbaskets,nnuclides,len(path_index_all),ntimes), dtype=float)
    for i in range(nelems):
        for j in range(nbaskets):
            for k in range(nnuclides-1):
                for l in range(len(path_index)):
                    for m in range(ntimes):
                        doses_all_sum[i][j][k][l][m]=doses[i][j][k][l][m]
                        doses_all_sum[i][j][nnuclides-1][l][m]+=doses[i][j][k][l][m]
                        for n in range(npaths,len(path_index_all)):
                            doses_all_sum[i][j][k][n][m]+=doses[i][j][k][l][m]*sum_matrix[n-npaths][l]
                            doses_all_sum[i][j][nnuclides-1][n][m]+=doses[i][j][k][l][m]*sum_matrix[n-npaths][l]
                            if prt:
                                #print(n)
                                continue
                            continue
                        continue
    
    return doses_all_sum,path_index_all
    

def biorad2_v05(config_file_name,prot_file,prt=False):
    biorad2_print_string="computation of doses, version 2.05"
    prot_file.write(biorad2_print_string+" begin:"+endline+endline)
    errors=0
    if prt:
        print("print test")

    #cteni kofigurace a nazvu poli
    bioconf=open(config_file_name,encoding='utf8')
    bioconfdata = yaml.load(bioconf, Loader=yaml.FullLoader)
    if prt:
        print(type(bioconf),type(bioconfdata))
    bioconf.close()
    
    activity_file_name=str(bioconfdata["activity_file_name"][0])
    database_file_name=str(bioconfdata["database_file_name"][0])
    results_file_name=str(bioconfdata["results_file_name"][0])
    if prt:
        print(activity_file_name,'\n',database_file_name,'\n',results_file_name,sep='')
    
    br_nuclide_index=read_index(bioconfdata,conf_nuclides,prt)
    br_nuclide_index=nucl_rep_name(br_nuclide_index,'',prt)
    br_nuclide_index_tot=br_nuclide_index.copy()
    br_nuclide_index_tot[total]=len(br_nuclide_index)
    if prt:
        print("nuclides with total")
        print(br_nuclide_index_tot)
    
    br_msh_el_index=read_index(bioconfdata,conf_msh_el,prt)

    #cteni databaze
    br_db=sqlite3.connect(database_file_name)
    br_db_cur=br_db.cursor()
    
    br_activity_index=get_env_index(br_db_cur,prt)
    if prt:
        print(br_activity_index)
        
    br_env_type=get_env_type(br_db_cur,br_activity_index)
    if prt:
        print("br_env_type:",br_env_type)
    #prt=False
    
    br_env_params_index=get_env_param_index(br_db_cur,prt)
    if prt:
        print("env. params. index:\n",br_env_params_index)
    
    br_env_params=get_env_params(prot_file,br_db_cur,br_activity_index,br_env_params_index,prt)
    if prt:
        print("env. params.:\n",br_env_params)
        for env in br_activity_index:
            e=br_activity_index[env][0]
            print(env)
            for param in br_env_params_index:
                p=br_env_params_index[param][0]
                print(param,br_env_params[e][p])
                
    br_Kd, errors=get_Kd(prot_file,br_db_cur,br_activity_index,br_nuclide_index,br_env_type,errors,prt)
    if prt:
        print("Kd:\n",br_Kd)
        for rn in br_nuclide_index:
            print(rn,br_Kd[br_activity_index[wood][0]][br_nuclide_index[rn]])
        
    br_basket_index=read_index(bioconfdata,conf_baskets,prt)
    
    #cteni aktivit
    #act_conf=read_activities_conf(bioconfdata,prt)
    act_conf=read_activities_conf2(bioconfdata,prt)
    if prt:
        print("activities all:",act_conf)
        #print(xxx)

    #if act_conf:
    if act_conf=="all":
        times,activities,errors=read_activity(prot_file,activity_file_name,br_msh_el_index,br_activity_index,br_nuclide_index,errors,prt)
    else:
        #times,activities,errors=read_activity_from_gw(prot_file,activity_file_name,br_msh_el_index,br_activity_index,br_nuclide_index,
        #                                       br_env_params_index,br_env_params,br_Kd,errors,prt)
        times,activities,errors=read_activity_from_water(prot_file,activity_file_name,br_msh_el_index,br_activity_index,br_nuclide_index,
                                               br_env_params_index,br_env_params,br_Kd,act_conf,errors,prt)

    a2s=get_a2s(br_db_cur)
    if prt:
        print("a2s:",a2s)

    cont_data_type_conf=read_cont_data_type(bioconfdata,prt)
    if prt:
        print(cont_data_type_conf)
        #print(xxx)
        
    if cont_data_type_conf=="concentration":
        if prt:
            print("concentration recalculate")
            print(cont_data_type_conf)
            #printf(xxxx)
        Na=get_Na(br_db_cur,prt)
        if prt:
            print("Na:",Na)
        for nuclide in br_nuclide_index:
            nw=get_atomic_weight_from_nucl_name(nuclide,prt)
            print(nuclide,nw)
        #print(error)
        br_nuclide_weight,errors=get_atomic_weight(prot_file,br_db_cur,br_nuclide_index,errors,prt)
        br_nuclide_half_life,errors=get_half_life(prot_file,br_db_cur,br_nuclide_index,a2s,errors,prt)
        br_nuclide_spec_activity=set_specific_activity(br_nuclide_half_life,br_nuclide_weight,Na,prt)
        activities=set_activity_from_concentration(activities,br_nuclide_spec_activity,prt)


    br_path_index_3=get_db_pathways3(br_db_cur,"calculated",prt)
    
    br_consumption_index=get_db_pathways3(br_db_cur,"in_basket",prt)
    br_consumption_all,errors=get_baskets_consumptions(prot_file,br_db_cur,br_basket_index,br_consumption_index,errors,prt)
    if prt:
        print("br_consumption_index:\n",br_consumption_index)
        for c in br_consumption_index:
            print(c,br_consumption_all[0][br_consumption_index[c][0]])
    
    br_path_index_sum=get_db_pathways3(br_db_cur,"sum",prt)
    if prt:
        print("br_path_index_sum:\n",br_path_index_sum)

    hing, errors=get_hing(prot_file,br_db_cur,br_nuclide_index,errors,prt)
    hinh, errors=get_hinh(prot_file,br_db_cur,br_nuclide_index,errors,prt)

    hext_soil,errors=get_hext(prot_file,br_db_cur,"exp_soil_inf",br_nuclide_index,errors,prt)
    hext_air,errors=get_hext(prot_file,br_db_cur,"exp_air",br_nuclide_index,errors,prt)
    hext_water,errors=get_hext(prot_file,br_db_cur,"exp_water",br_nuclide_index,errors,prt)
    
    if prt:
        for rn in br_nuclide_index:
            n=br_nuclide_index[rn]
            print(rn,hext_water[n])

    wat_ing, errors=get_consumption(prot_file,br_db_cur,br_basket_index,wat_ing_path,errors,prt)
    soil_ing, errors=get_consumption(prot_file,br_db_cur,br_basket_index,soil_ing_path,errors,prt)
    
    plant_index=get_plants_index(prot_file,br_db_cur,prt)
    if prt:
        print("plant_index:\n",plant_index)
        print("plant_index['mushroom']:\n",plant_index["mushrooms"])
    plant_tff, errors=get_plant_tf_fresh(prot_file,br_db_cur,plant_index,br_nuclide_index,errors,prt)
    plant_ing, errors=get_plant_consumption(prot_file,br_db_cur,br_basket_index,plant_index,errors,prt)

    if prt:
        for p in plant_index:
            pi=plant_index[p][0]
            print(p)
            for rn in br_nuclide_index:
                n=br_nuclide_index[rn]
                print(rn,plant_tff[pi][n])

    
    lf_resp_factor=get_resuspension_factor(prot_file,br_db_cur,prt)
    
    wat_prod_index=get_wat_prod_index(prot_file,br_db_cur,prt)
    wat_prod_tf, errors=get_wat_prod_tf(prot_file,br_db_cur,wat_prod_index,br_nuclide_index,errors,prt)
    wat_prod_ing, errors =get_water_prod_consumption(prot_file,br_db_cur,br_basket_index,wat_prod_index,errors,prt)

    if prt:
        for wp in wat_prod_index:
            wpi=wat_prod_index[wp]
            print(wp)
            for rn in br_nuclide_index:
                n=br_nuclide_index[rn]
                print(rn,wat_prod_tf[wpi][n])
    
    animal_prod_index=get_animal_product_index(prot_file,br_db_cur,prt)
    if prt:
        print("\nanimal_prod_index",animal_prod_index)
    animal_index=get_animal_index(prot_file,br_db_cur,prt)
    if prt:
        print("animal_index",animal_index)
    animal_consump_index=get_animal_consumption_index(prot_file,br_db_cur,prt)
    if prt:
        print("animal_consump_index",animal_consump_index)
    animal_consumption=get_animal_consumption(prot_file,br_db_cur,animal_index,animal_consump_index,prt)
    if prt:
        print("animal_consumption_index",animal_consump_index)
        for a in animal_index:
            ai=animal_index[a]
            print(a)
            for ac in animal_consump_index:
                aci=animal_consump_index[ac]
                print(animal_consumption[ai][aci])

    animal_consumption_activity=get_animal_consumption_acitivies(prot_file,br_msh_el_index,
                                                                 animal_index,br_nuclide_index,times,animal_consump_index,animal_consumption,
                                                                 br_activity_index,activities,plant_index,plant_tff,br_env_params_index,br_env_params,prt)
    if prt:
        print("animal_consumption_activity",animal_consumption_activity)
    animal_prod_ing,errors=get_animal_prod_consumption(prot_file,br_db_cur,br_basket_index,animal_prod_index,errors,prt)
    if prt:
        print("animal product:",animal_prod_index)
        print("water product",wat_prod_index)

    animal_prod_tf,errors=get_animal_prod_tf(prot_file,br_db_cur,animal_prod_index,br_nuclide_index,errors,prt)
    if prt:
        for ap in animal_prod_index:
            api=animal_prod_index[ap][0]
            print(ap)
            for rn in br_nuclide_index:
                n=br_nuclide_index[rn]
                print(rn,animal_prod_tf[api][n])

    if prt:
        for a in animal_index:
            ai=animal_index[a]
            print(a)
            for rn in br_nuclide_index:
                n=br_nuclide_index[rn]
                print(rn,animal_consumption_activity[0][ai][n][0])

    br_basket_index_db,errors=get_basket_index(prot_file,br_db_cur,errors,prt)
    sum_cross_matrix=get_sum_cross_matrix(prot_file,br_db_cur,br_path_index_3,br_path_index_sum,prt)    
    br_db.close()
    
    br_basket_index_col,errors=collect_basket_index(prot_file,br_basket_index,br_basket_index_db,errors,prt)
    if prt:
        print("br_basket_index:")
        print(br_basket_index)
        print("br_basket_index_col:")
        print(br_basket_index_col)

    prot_file.write("doses calculations"+endline)

    doses=doses_calc3(prot_file,br_msh_el_index,br_basket_index,br_nuclide_index,br_path_index_3,times,activities,br_activity_index,
               br_env_params_index,br_env_params,
               br_consumption_index,br_consumption_all,
               hing,hinh,hext_soil,hext_air,hext_water,wat_ing,soil_ing,
               plant_index,plant_tff,plant_ing,lf_resp_factor,
               animal_prod_index,animal_prod_ing,animal_index,animal_prod_tf,animal_consumption_activity,
               wat_prod_index,wat_prod_tf,wat_prod_ing,
               a2s,prt)

    doses_all_sum,path_index_all=get_doses_all_sum(prot_file,doses,sum_cross_matrix,len(br_msh_el_index),len(br_basket_index),
                                               len(br_nuclide_index_tot),br_path_index_3,len(times),br_path_index_sum,prt=True)
    if prt:
        print("path all sum: ",len(path_index_all))
        print(path_index_all)
        print("path all sum end")

        for p in path_index_all:
            pi=path_index_all[p][0]
            print(p)
            for rn in br_nuclide_index_tot:
                n=br_nuclide_index_tot[rn]
                print(rn,doses_all_sum[0][0][n][pi][0])

    np.savez(results_file_name,doses=doses_all_sum,times=times,elem_index=br_msh_el_index,basket_index=br_basket_index_col,nuclide_index=br_nuclide_index_tot,path_index=path_index_all,env_index=br_activity_index,activities=activities)

    prot_file.write(endline+"total input data errors number is: "+str(errors)+endline)
    if errors>0:
        prot_file.write("list of input data erors:"+endline+endline)
    else:
        prot_file.write(endline)
    for i in range(len(list_data_errors)):
        prot_file.write(list_data_errors[i])
        
    prot_file.write(biorad2_print_string+" end:"+endline+endline)
    return errors

def read_results2(results_file_name):
    doses=np.load(results_file_name)["doses"]
    times=np.load(results_file_name)["times"]
    br_msh_el_index=np.load(results_file_name,allow_pickle='TRUE')["elem_index"][()]
    br_basket_index=np.load(results_file_name,allow_pickle='TRUE')["basket_index"][()]
    br_nuclide_index=np.load(results_file_name,allow_pickle='TRUE')["nuclide_index"][()]
    br_path_index=np.load(results_file_name,allow_pickle='TRUE')["path_index"][()]
    env_index=np.load(results_file_name,allow_pickle='TRUE')["env_index"][()]
    activities=np.load(results_file_name,allow_pickle='TRUE')["activities"][()]
    return doses, times, br_msh_el_index, br_basket_index, br_nuclide_index, br_path_index, env_index, activities


def get_baskets(database_name,prt=False):

    index_names=[None,'','']
    index = {}

    br_db=sqlite3.connect(database_name)
    br_db_cur=br_db.cursor()
    #sql_ask="SELECT * FROM Basket_type"
    sql_ask="SELECT * FROM s_used_baskets"
    col_index=get_db_ask_column_names_index(br_db_cur,sql_ask)
    br_db_cur.execute(sql_ask)
    rows = br_db_cur.fetchall()

    for row in rows:
        if prt:
            print(row,type(row))
        index_names[1]=row[col_index["name_cz"]]
        index_names[2]=row[col_index["name_en"]]
        index[row[col_index["identificator"]]]=index_names.copy()
    br_db.close()
    return index

def get_environments(database_name,prt=False):
    index_names=[None,'','']
    index = {}
    br_db=sqlite3.connect(database_name)
    br_db_cur=br_db.cursor()
    
    eindex=get_env_index(br_db_cur)
    for ei in eindex:
        index_names[1]=eindex[ei][1]
        index_names[2]=eindex[ei][2]
        index[ei]=index_names.copy()
    br_db.close()
    return index
