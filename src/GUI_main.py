import tkinter as tk
from tkinter.filedialog import askopenfilename
from tkinter import ttk
import yaml
import sqlite3 as sq
import biorad2_library_v0072 as brl
import numpy as np
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
import matplotlib.pyplot as plt
from biorad1_main import execute_Biorad1
from process_UZ_results import *
from pathlib import Path


base_path = Path(__file__).parent

global scenario_var


class MyDumper(yaml.Dumper):
    def increase_indent(self, flow=False, indentless=False):
        return super(MyDumper, self).increase_indent(flow, False)


class ScrollableFrame(ttk.Frame):
    def __init__(self, container, *args, **kwargs):
        super().__init__(container, *args, **kwargs)
        canvas = tk.Canvas(self, width=600, height=400)
        scrollbar_1 = ttk.Scrollbar(self, orient="vertical", command=canvas.yview)
        scrollbar_2 = ttk.Scrollbar(self, orient="horizontal", command=canvas.xview)
        self.scrollable_frame = ttk.Frame(canvas)

        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(
                scrollregion=canvas.bbox("all")
            )
        )

        canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")

        canvas.configure(yscrollcommand=scrollbar_1.set)
        canvas.configure(xscrollcommand=scrollbar_2.set)

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar_1.pack(side="right", fill="y")
        scrollbar_2.pack(side="bottom", fill="x")


def what_up():
    about_w = tk.Toplevel()
    about_w.title('About BIORAD')
    about_label_text = tk.StringVar()
    about_label_text.set('BioRad version 1.0.4\n\nDatabase version 22\n\nCopyright (C) 2020 Technical University of '
                         'Liberec. All rights reserved.\n\nThis program is '
                         'free software; you can redistribute it and/or modify it under the terms of the \nGNU General '
                         'Public License version 3 as published by the Free Software Foundation. \n'
                         '(http://www.gnu.org/licenses/gpl-3.0.en.html)\n\nThis program is distributed in the hope that'
                         ' it will be useful, but WITHOUT ANY WARRANTY; \nwithout even the implied warranty of '
                         'MERCHANTABILITY or FITNESS FOR A PARTICULAR \nPURPOSE. See the GNU General '
                         'Public License for more details.\n\nCreated with a support from the Technology Agency of the '
                         'Czech Republic \nwithin the EPSILON program through the project TH03030274 \nSoftware for '
                         'evaluation of radionuclide transport on the geosphere/biosphere '
                         'interface and its impact on man.\n\n')
    tk.Label(about_w, textvariable=about_label_text).pack()

    def quit_this():
        about_w.destroy()

    ttk.Button(about_w, text='OK', command=quit_this).pack()


def show_me():
    class Visual(tk.Tk):
        def __init__(self):
            super(Visual, self).__init__()
            self.title('Biosphere module results')
            w, h = self.winfo_screenwidth(), self.winfo_screenheight()
            self.geometry("%dx%d+0+0" % (w, h))

    def select_all(varr):
        for p in varr:
            p.set(True)

    def deselect_all(varr):
        for p in varr:
            p.set(False)

    def quit_this():
        visual.quit()
        visual.destroy()

    visual = Visual()
    file = (base_path / "../outputs/Biosphere_module_results.npz").resolve()
    visual.labelFrame_3 = ttk.LabelFrame(visual, text='Results have been computed for the following basket:')
    visual.labelFrame_3.grid(column=0, row=0, padx=20, pady=20, sticky=tk.W)
    basket = np.load(file, allow_pickle='TRUE')["basket_index"][()]
    basket_text = tk.StringVar(visual.labelFrame_3, '...')
    tk.Label(visual.labelFrame_3, textvariable=basket_text).pack(anchor=tk.CENTER)
    for b in basket:
        basket_text.set(b)

    visual.labelFrame_1 = ttk.LabelFrame(visual, text='Select paths to show:')
    visual.labelFrame_1.grid(column=0, row=1, padx=20, pady=20, sticky=tk.W)
    paths = np.load(file, allow_pickle='TRUE')["path_index"][()]
    path_var = []
    path_name = []
    index = -1
    for item in paths.values():
        index += 1
        path_var.append(tk.BooleanVar(visual, True))
        path_name.append(item[2])
        c = ttk.Checkbutton(visual.labelFrame_1, text=item[2], variable=path_var[index])
        c.pack(anchor=tk.W)
    ttk.Button(visual.labelFrame_1, text='Select all', command=lambda: select_all(path_var)).pack(side=tk.LEFT)
    ttk.Button(visual.labelFrame_1, text='Deselect all', command=lambda: deselect_all(path_var)).pack(side=tk.LEFT)

    visual.labelFrame_2 = ttk.LabelFrame(visual, text='Select nuclides to show:')
    visual.labelFrame_2.grid(column=2, row=1, padx=20, pady=20, sticky=tk.NW)
    nuclides = np.load(file, allow_pickle='TRUE')["nuclide_index"][()]
    nucl_var = []
    nucl_name = []
    index = -1
    for item in nuclides:
        index += 1
        nucl_var.append(tk.BooleanVar(visual, True))
        nucl_name.append(item)
        c = ttk.Checkbutton(visual.labelFrame_2, text=item, variable=nucl_var[index])
        c.pack(anchor=tk.W)
    ttk.Button(visual.labelFrame_2, text='Select all', command=lambda: select_all(nucl_var)).pack(side=tk.LEFT)
    ttk.Button(visual.labelFrame_2, text='Deselect all', command=lambda: deselect_all(nucl_var)).pack(side=tk.LEFT)

    visual.labelFrame_5 = ttk.LabelFrame(visual, text='Plot area:')
    visual.labelFrame_5.grid(column=1, row=1, rowspan=3, padx=20, pady=20, sticky=tk.W)

    times = np.load(file)["times"]
    doses = np.load(file)["doses"]

    fig = plt.Figure(figsize=(8, 6), dpi=100)
    canvas = FigureCanvasTkAgg(fig, master=visual.labelFrame_5)
    canvas.draw()
    toolbar = NavigationToolbar2Tk(canvas, visual.labelFrame_5)
    toolbar.update()

    def on_key_press(event):
        print("you pressed {}".format(event.key))
        key_press_handler(event, canvas, toolbar)

    canvas.mpl_connect("key_press_event", on_key_press)
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    def plot_selected():
        fig.clf()
        lines = []
        labels = []
        doses_csv = []
        for i in range(0, len(nucl_name)):
            if nucl_var[i].get():
                print(nucl_name[i])
                for j in range(0, len(path_name)):
                    if path_var[j].get():
                        print(nucl_name[i])
                        print(path_name[j])
                        line = fig.add_subplot(1, 1, 1, xlabel='Time [' + time_unit + ']', ylabel='Dose [Sv/y]',
                                               yscale='log').plot(times, doses[0][0][i][j][:])[0]
                        doses_csv.append(doses[0][0][i][j][:])
                        lines.append(line)
                        labels.append(nucl_name[i] + ' - ' + path_name[j])
        fig.legend(lines, labels, loc="center right", borderaxespad=0.1, fontsize=5)
        fig.subplots_adjust(right=0.7)
        canvas.draw()

        file_path = (base_path / "../outputs/BioRad_results_export.csv").resolve()
        with open(file_path, 'w') as csv_dump:
            csv_dump.write('Time;')
            for i in range(0, len(labels)):
                if i == len(labels)-1:
                    csv_dump.write(labels[i] + '\n')
                else:
                    csv_dump.write(labels[i] + '; ')
            for i in range(0, len(times)):
                csv_dump.write(str(times[i]) + '; ')
                for j in range(0, len(doses_csv)):
                    if j == len(doses_csv)-1:
                        csv_dump.write(str(doses_csv[j][i]) + '\n')
                    else:
                        csv_dump.write(str(doses_csv[j][i]) + '; ')

    visual.labelFrame_4 = ttk.LabelFrame(visual)
    visual.labelFrame_4.grid(column=1, row=0, padx=20, pady=20, sticky=tk.W)
    ttk.Button(visual.labelFrame_4, text='Plot and save as CSV', command=plot_selected).pack(side=tk.LEFT)
    ttk.Button(visual.labelFrame_4, text='Quit', command=quit_this).pack(side=tk.LEFT)
    visual.mainloop()


def open_transport_file():
    label_text.set('Reading transport file.\n This may take a while if the file is large...')
    filename = askopenfilename(filetypes=[("GMSH files", "*.msh")])
    print(filename)
    galerkin = True
    with open(filename) as trans_file:
        line = trans_file.readline()
        while line != '$EndElements\n':
            line = trans_file.readline()
        line = trans_file.readline()
        if line == '$ElementNodeData\n':
            print('Transport results from DG method.')
        else:
            galerkin = False
            print('Transport results from OS method.')
    fields = []
    with open(filename) as trans_file:
        for line in trans_file:
            if line.rstrip() == '$ElementNodeData' or line.rstrip() == '$ElementData':
                trans_file.readline()
                field = trans_file.readline().rstrip()[1:-6]
                if field in fields:
                    break
                else:
                    fields.append(field)
    top = tk.Toplevel()
    top.title('Tracer selection.')
    top_label_text = tk.StringVar()
    top_label_text.set('The following is available in transport output.\nPlease choose what you are interested in.')
    top_label = tk.Label(top, textvariable=top_label_text).pack()
    checkboxlist = []
    varlist = []
    selected = []
    element_ID = []

    def resolve_checkbox():
        result = []
        element_ID.append(int(e.get()))
        for j in range(0, len(checkboxlist)):
            result.append(varlist[j].get())
        for item in fields:
            if result[fields.index(item)] == 1:
                selected.append(item)
        top.quit()
        top.destroy()

    for i in range(0, len(fields)):
        varlist.append(tk.BooleanVar())
        c = tk.Checkbutton(top, text=fields[i], variable=varlist[i], onvalue=True, offvalue=False).pack()
        checkboxlist.append(c)
    tk.Label(top, text='Enter ID of element you want to evaluate.').pack()
    e = tk.Entry(top)
    e.insert(0, '267')
    e.pack()
    tk.Label(top, text='\n').pack()
    button = tk.Button(top, text='Done.', width=25, command=resolve_checkbox).pack()
    top.mainloop()
    if len(selected) == 0:
        label_text.set('You have selected nothing. Please, try again.')
    else:
        file_path = (base_path / "../outputs/temp.txt").resolve()
        temp_file = open(file_path, 'w')
        for item in selected:
            temp_file.write(item + '\n')
        temp_file.close()

    with open(filename) as trans_file:
        file_path = (base_path / "../inputs/saturated_concentrations.csv").resolve()
        with open(file_path, 'w') as output_file:
            output_file.write(str(element_ID) + '\n')
            output_file.write('time;')
            for k in range(0, len(selected)):
                if k != len(selected) - 1:
                    output_file.write(selected[k] + ';')
                else:
                    output_file.write(selected[k] + '\n')
            times = []
            concentrations = []
            conc_t = []
            for line in trans_file:
                if line == '$ElementNodeData\n' or line == '$ElementData\n':
                    trans_file.readline()
                    field = trans_file.readline().rstrip()[1:-6]
                    if field in selected:
                        trans_file.readline()
                        time = float(trans_file.readline().rstrip())
                        if len(times) == 0:
                            times.append(time)
                        if time != times[-1]:
                            times.append(time)
                        trans_file.readline()
                        trans_file.readline()
                        trans_file.readline()
                        nr = int(trans_file.readline().rstrip())
                        for i in range(0, nr):
                            line = trans_file.readline().rstrip().split(' ')
                            if int(line[0]) == element_ID[0]:
                                if not galerkin:
                                    concentration = float(line[1])
                                    if concentration < 0:
                                        concentration = 0
                                    conc_t.append(concentration)
                                    if len(conc_t) == len(selected):
                                        concentrations.append(conc_t)
                                        output_file.write(str(times[-1]) + ';')
                                        for k in range(0, len(conc_t)):
                                            if k != len(conc_t) - 1:
                                                output_file.write(str(conc_t[k]) + ';')
                                            else:
                                                output_file.write(str(conc_t[k]) + '\n')
                                        conc_t = []
                                else:
                                    val_nr = int(line[1])
                                    value = 0
                                    for j in range(0, val_nr):
                                        value = value + float(line[-(j + 1)])
                                    value = value / val_nr
                                    if value < 0:
                                        value = 0
                                    conc_t.append(value)
                                    if len(conc_t) == len(selected):
                                        concentrations.append(conc_t)
                                        output_file.write(str(times[-1]) + ';')
                                        for k in range(0, len(conc_t)):
                                            if k != len(conc_t) - 1:
                                                output_file.write(str(conc_t[k]) + ';')
                                            else:
                                                output_file.write(str(conc_t[k]) + '\n')
                                        conc_t = []
    label_text.set('Transport file read.')


def call_biorad_1():
    b1 = tk.Toplevel()
    b1.title('Unsaturated zone module user options')
    ww, hh = b1.winfo_screenwidth(), b1.winfo_screenheight()
    b1.geometry("%dx%d+0+0" % (ww, hh))

    labelFrame_1 = ttk.LabelFrame(b1, text='Select units:')
    labelFrame_1.grid(column=0, row=0, padx=20, pady=20, sticky=tk.W)
    ttk.Label(labelFrame_1, text='Length: ').grid(column=0, row=0)
    length_c = tk.StringVar(b1, value="m")
    ttk.Combobox(labelFrame_1, values=["mm", "cm", "dm", "m"], state='readonly', textvariable=length_c)\
        .grid(column=1, row=0)
    mass_c = tk.StringVar(b1, value="kg")
    ttk.Label(labelFrame_1, text='Mass: ').grid(column=0, row=1)
    ttk.Combobox(labelFrame_1, values=["ng", "ug", "mg", "g", "kg"], state='readonly', textvariable=mass_c)\
        .grid(column=1, row=1)
    time_c = tk.StringVar(b1, value="s")
    ttk.Label(labelFrame_1, text='Time: ').grid(column=0, row=2)
    ttk.Combobox(labelFrame_1, values=["s", "h", "day", "year"], state='readonly', textvariable=time_c)\
        .grid(column=1, row=2)

    labelFrame_2 = ttk.LabelFrame(b1, text='Simulation time parameters:')
    labelFrame_2.grid(column=0, row=1, padx=20, pady=20, sticky=tk.W)

    period_e = tk.StringVar(b1, value='0')
    ttk.Label(labelFrame_2, text='Simulation period: ').grid(column=0, row=0)
    ttk.Entry(labelFrame_2, textvariable=period_e).grid(column=1, row=0)
    step_e = tk.StringVar(b1, value='0')
    ttk.Label(labelFrame_2, text='Simulation time step: ').grid(column=0, row=1)
    ttk.Entry(labelFrame_2, textvariable=step_e).grid(column=1, row=1)
    iter_e = tk.StringVar(b1, value='0')
    ttk.Label(labelFrame_2, text='Flow iterations per step: ').grid(column=0, row=2)
    ttk.Entry(labelFrame_2, textvariable=iter_e).grid(column=1, row=2)
    out_e = tk.StringVar(b1, value='0')
    ttk.Label(labelFrame_2, text='Output time step: ').grid(column=0, row=3)
    ttk.Entry(labelFrame_2, textvariable=out_e).grid(column=1, row=3)

    def geometry(event):
        global trash
        for item in trash:
            item.destroy()
        trash = []
        global heights
        heights = []
        group_count = int(group_nr.get())
        labelFrame_4 = ttk.LabelFrame(b1, text='Geometry layers (bottom to top): ')
        labelFrame_4.grid(column=0, row=3, padx=20, pady=20, sticky=tk.W, columnspan=2)
        trash.append(labelFrame_4)
        for i in range(0, group_count):
            t = ttk.Label(labelFrame_4, text='Layer number '+str(i+1)+': ')
            t.grid(column=2*i, row=0)
            trash.append(t)
            h = ttk.Label(labelFrame_4, text='Enter layer height: ')
            h.grid(column=2*i, row=1)
            trash.append(h)
            height = tk.StringVar(b1, value='0')
            heights.append(height)
            m = ttk.Entry(labelFrame_4, textvariable=heights[i])
            m.grid(column=2*i+1, row=1)
            trash.append(m)
            mode_label = ttk.Label(labelFrame_4, text='UZ parameters from: ')
            mode_label.grid(column=2*i, row=2)
            trash.append(mode_label)

        def parameters(eff, ID, text):
            global trash_1, trash_2, trash_3
            global sand_perc_1, sand_perc_2, sand_perc_3, silt_perc_1, silt_perc_2, silt_perc_3, clay_perc_1, \
                clay_perc_2, clay_perc_3, dens_1, dens_2, dens_3, soil_sel_1, soil_sel_2, soil_sel_3, theta_r_1, \
                theta_r_2, theta_r_3, theta_s_1, theta_s_2, theta_s_3, alpha_1, alpha_2, alpha_3, en_1, en_2, en_3, \
                ks_1, ks_2, ks_3, v_dens_1, v_dens_2, v_dens_3

            if ID == 1:
                for item in trash_1:
                    item.destroy()
                if text == 'Grain structure':
                    label_a = ttk.Label(labelFrame_4, text='sand [%]: ')
                    label_a.grid(column=0, row=3)
                    sand_perc_1 = tk.StringVar(b1)
                    value_a = ttk.Entry(labelFrame_4, textvariable=sand_perc_1)
                    value_a.grid(column=1, row=3)
                    label_b = ttk.Label(labelFrame_4, text='silt [%]: ')
                    label_b.grid(column=0, row=4)
                    silt_perc_1 = tk.StringVar(b1)
                    value_b = ttk.Entry(labelFrame_4, textvariable=silt_perc_1)
                    value_b.grid(column=1, row=4)
                    label_c = ttk.Label(labelFrame_4, text='clay [%]: ')
                    label_c.grid(column=0, row=5)
                    clay_perc_1 = tk.StringVar(b1)
                    value_c = ttk.Entry(labelFrame_4, textvariable=clay_perc_1)
                    value_c.grid(column=1, row=5)
                    label_d = ttk.Label(labelFrame_4, text='density [kg/m3]: ')
                    label_d.grid(column=0, row=6)
                    dens_1 = tk.StringVar(b1)
                    value_d = ttk.Entry(labelFrame_4, textvariable=dens_1)
                    value_d.grid(column=1, row=6)
                    trash_1.append(label_a)
                    trash_1.append(label_b)
                    trash_1.append(label_c)
                    trash_1.append(label_d)
                    trash_1.append(value_a)
                    trash_1.append(value_b)
                    trash_1.append(value_c)
                    trash_1.append(value_d)
                if text == 'Soil type':
                    label_a = ttk.Label(labelFrame_4, text='Select soil type: ')
                    label_a.grid(column=0, row=3)
                    soil_sel_1 = tk.StringVar(b1)
                    sel = ttk.Combobox(labelFrame_4, state='readonly', textvariable=soil_sel_1,
                                       values=['Sand', 'LoamySand', 'SandyLoam', 'Loam', 'Silt', 'SiltLoam',
                                               'SandyClayLoam', 'ClayLoam', 'SiltyClayLoam', 'SandyClay',
                                               'SiltyClay', 'Clay'])
                    sel.grid(column=1, row=3)
                    trash_1.append(sel)
                    trash_1.append(label_a)
                if text == 'van Genuchten':
                    label_a = ttk.Label(labelFrame_4, text='Theta_r: ')
                    label_a.grid(column=0, row=3)
                    theta_r_1 = tk.StringVar(b1)
                    value_a = ttk.Entry(labelFrame_4, textvariable=theta_r_1)
                    value_a.grid(column=1, row=3)
                    label_b = ttk.Label(labelFrame_4, text='Theta_s: ')
                    label_b.grid(column=0, row=4)
                    theta_s_1 = tk.StringVar(b1)
                    value_b = ttk.Entry(labelFrame_4, textvariable=theta_s_1)
                    value_b.grid(column=1, row=4)
                    label_c = ttk.Label(labelFrame_4, text='Alpha: ')
                    label_c.grid(column=0, row=5)
                    alpha_1 = tk.StringVar(b1)
                    value_c = ttk.Entry(labelFrame_4, textvariable=alpha_1)
                    value_c.grid(column=1, row=5)
                    label_d = ttk.Label(labelFrame_4, text='n: ')
                    label_d.grid(column=0, row=6)
                    en_1 = tk.StringVar(b1)
                    value_d = ttk.Entry(labelFrame_4, textvariable=en_1)
                    value_d.grid(column=1, row=6)
                    label_e = ttk.Label(labelFrame_4, text='Ks: ')
                    label_e.grid(column=0, row=7)
                    ks_1 = tk.StringVar(b1)
                    value_e = ttk.Entry(labelFrame_4, textvariable=ks_1)
                    value_e.grid(column=1, row=7)
                    label_f = ttk.Label(labelFrame_4, text='Density [kg/m3]: ')
                    label_f.grid(column=0, row=8)
                    v_dens_1 = tk.StringVar(b1)
                    value_f = ttk.Entry(labelFrame_4, textvariable=v_dens_1)
                    value_f.grid(column=1, row=8)
                    trash_1.append(label_a)
                    trash_1.append(label_b)
                    trash_1.append(label_c)
                    trash_1.append(label_d)
                    trash_1.append(label_e)
                    trash_1.append(label_f)
                    trash_1.append(value_a)
                    trash_1.append(value_b)
                    trash_1.append(value_c)
                    trash_1.append(value_d)
                    trash_1.append(value_e)
                    trash_1.append(value_f)

                for item in trash_1:
                    trash.append(item)

            if ID == 2:
                for item in trash_2:
                    item.destroy()
                if text == 'Grain structure':
                    label_a = ttk.Label(labelFrame_4, text='sand [%]: ')
                    label_a.grid(column=2, row=3)
                    sand_perc_2 = tk.StringVar(b1)
                    value_a = ttk.Entry(labelFrame_4, textvariable=sand_perc_2)
                    value_a.grid(column=3, row=3)
                    label_b = ttk.Label(labelFrame_4, text='silt [%]: ')
                    label_b.grid(column=2, row=4)
                    silt_perc_2 = tk.StringVar(b1)
                    value_b = ttk.Entry(labelFrame_4, textvariable=silt_perc_2)
                    value_b.grid(column=3, row=4)
                    label_c = ttk.Label(labelFrame_4, text='clay [%]: ')
                    label_c.grid(column=2, row=5)
                    clay_perc_2 = tk.StringVar(b1)
                    value_c = ttk.Entry(labelFrame_4, textvariable=clay_perc_2)
                    value_c.grid(column=3, row=5)
                    label_d = ttk.Label(labelFrame_4, text='density [kg/m3]: ')
                    label_d.grid(column=2, row=6)
                    dens_2 = tk.StringVar(b1)
                    value_d = ttk.Entry(labelFrame_4, textvariable=dens_2)
                    value_d.grid(column=3, row=6)
                    trash_2.append(label_a)
                    trash_2.append(label_b)
                    trash_2.append(label_c)
                    trash_2.append(label_d)
                    trash_2.append(value_a)
                    trash_2.append(value_b)
                    trash_2.append(value_c)
                    trash_2.append(value_d)
                if text == 'Soil type':
                    label_a = ttk.Label(labelFrame_4, text='Select soil type: ')
                    label_a.grid(column=2, row=3)
                    soil_sel_2 = tk.StringVar(b1)
                    sel = ttk.Combobox(labelFrame_4, state='readonly', textvariable=soil_sel_2,
                                       values=['Sand', 'LoamySand', 'SandyLoam', 'Loam', 'Silt', 'SiltLoam',
                                               'SandyClayLoam', 'ClayLoam', 'SiltyClayLoam', 'SandyClay',
                                               'SiltyClay', 'Clay'])
                    sel.grid(column=3, row=3)
                    trash_2.append(sel)
                    trash_2.append(label_a)
                if text == 'van Genuchten':
                    label_a = ttk.Label(labelFrame_4, text='Theta_r: ')
                    label_a.grid(column=2, row=3)
                    theta_r_2 = tk.StringVar(b1)
                    value_a = ttk.Entry(labelFrame_4, textvariable=theta_r_2)
                    value_a.grid(column=3, row=3)
                    label_b = ttk.Label(labelFrame_4, text='Theta_s: ')
                    label_b.grid(column=2, row=4)
                    theta_s_2 = tk.StringVar(b1)
                    value_b = ttk.Entry(labelFrame_4, textvariable=theta_s_2)
                    value_b.grid(column=3, row=4)
                    label_c = ttk.Label(labelFrame_4, text='Alpha: ')
                    label_c.grid(column=2, row=5)
                    alpha_2 = tk.StringVar(b1)
                    value_c = ttk.Entry(labelFrame_4, textvariable=alpha_2)
                    value_c.grid(column=3, row=5)
                    label_d = ttk.Label(labelFrame_4, text='n: ')
                    label_d.grid(column=2, row=6)
                    en_2 = tk.StringVar(b1)
                    value_d = ttk.Entry(labelFrame_4, textvariable=en_2)
                    value_d.grid(column=3, row=6)
                    label_e = ttk.Label(labelFrame_4, text='Ks: ')
                    label_e.grid(column=2, row=7)
                    ks_2 = tk.StringVar(b1)
                    value_e = ttk.Entry(labelFrame_4, textvariable=ks_2)
                    value_e.grid(column=3, row=7)
                    label_f = ttk.Label(labelFrame_4, text='Density [kg/m3]: ')
                    label_f.grid(column=2, row=8)
                    v_dens_2 = tk.StringVar(b1)
                    value_f = ttk.Entry(labelFrame_4, textvariable=v_dens_2)
                    value_f.grid(column=3, row=8)
                    trash_2.append(label_a)
                    trash_2.append(label_b)
                    trash_2.append(label_c)
                    trash_2.append(label_d)
                    trash_2.append(label_e)
                    trash_2.append(label_f)
                    trash_2.append(value_a)
                    trash_2.append(value_b)
                    trash_2.append(value_c)
                    trash_2.append(value_d)
                    trash_2.append(value_e)
                    trash_2.append(value_f)

                for item in trash_2:
                    trash.append(item)

            if ID == 3:
                for item in trash_3:
                    item.destroy()
                if text == 'Grain structure':
                    label_a = ttk.Label(labelFrame_4, text='sand [%]: ')
                    label_a.grid(column=4, row=3)
                    sand_perc_3 = tk.StringVar(b1)
                    value_a = ttk.Entry(labelFrame_4, textvariable=sand_perc_3)
                    value_a.grid(column=5, row=3)
                    label_b = ttk.Label(labelFrame_4, text='silt [%]: ')
                    label_b.grid(column=4, row=4)
                    silt_perc_3 = tk.StringVar(b1)
                    value_b = ttk.Entry(labelFrame_4, textvariable=silt_perc_3)
                    value_b.grid(column=5, row=4)
                    label_c = ttk.Label(labelFrame_4, text='clay [%]: ')
                    label_c.grid(column=4, row=5)
                    clay_perc_3 = tk.StringVar(b1)
                    value_c = ttk.Entry(labelFrame_4, textvariable=clay_perc_3)
                    value_c.grid(column=5, row=5)
                    label_d = ttk.Label(labelFrame_4, text='density [kg/m3]: ')
                    label_d.grid(column=4, row=6)
                    dens_3 = tk.StringVar(b1)
                    value_d = ttk.Entry(labelFrame_4, textvariable=dens_3)
                    value_d.grid(column=5, row=6)
                    trash_3.append(label_a)
                    trash_3.append(label_b)
                    trash_3.append(label_c)
                    trash_3.append(label_d)
                    trash_3.append(value_a)
                    trash_3.append(value_b)
                    trash_3.append(value_c)
                    trash_3.append(value_d)
                if text == 'Soil type':
                    label_a = ttk.Label(labelFrame_4, text='Select soil type: ')
                    label_a.grid(column=4, row=3)
                    soil_sel_3 = tk.StringVar(b1)
                    sel = ttk.Combobox(labelFrame_4, state='readonly', textvariable=soil_sel_3,
                                       values=['Sand', 'LoamySand', 'SandyLoam', 'Loam', 'Silt', 'SiltLoam',
                                               'SandyClayLoam', 'ClayLoam', 'SiltyClayLoam', 'SandyClay',
                                               'SiltyClay', 'Clay'])
                    sel.grid(column=5, row=3)
                    trash_3.append(sel)
                    trash_3.append(label_a)
                if text == 'van Genuchten':
                    label_a = ttk.Label(labelFrame_4, text='Theta_r: ')
                    label_a.grid(column=4, row=3)
                    theta_r_3 = tk.StringVar(b1)
                    value_a = ttk.Entry(labelFrame_4, textvariable=theta_r_3)
                    value_a.grid(column=5, row=3)
                    label_b = ttk.Label(labelFrame_4, text='Theta_s: ')
                    label_b.grid(column=4, row=4)
                    theta_s_3 = tk.StringVar(b1)
                    value_b = ttk.Entry(labelFrame_4, textvariable=theta_s_3)
                    value_b.grid(column=5, row=4)
                    label_c = ttk.Label(labelFrame_4, text='Alpha: ')
                    label_c.grid(column=4, row=5)
                    alpha_3 = tk.StringVar(b1)
                    value_c = ttk.Entry(labelFrame_4, textvariable=alpha_3)
                    value_c.grid(column=5, row=5)
                    label_d = ttk.Label(labelFrame_4, text='n: ')
                    label_d.grid(column=4, row=6)
                    en_3 = tk.StringVar(b1)
                    value_d = ttk.Entry(labelFrame_4, textvariable=en_3)
                    value_d.grid(column=5, row=6)
                    label_e = ttk.Label(labelFrame_4, text='Ks: ')
                    label_e.grid(column=4, row=7)
                    ks_3 = tk.StringVar(b1)
                    value_e = ttk.Entry(labelFrame_4, textvariable=ks_3)
                    value_e.grid(column=5, row=7)
                    label_f = ttk.Label(labelFrame_4, text='Density [kg/m3]: ')
                    label_f.grid(column=4, row=8)
                    v_dens_3 = tk.StringVar(b1)
                    value_f = ttk.Entry(labelFrame_4, textvariable=v_dens_3)
                    value_f.grid(column=5, row=8)
                    trash_3.append(label_a)
                    trash_3.append(label_b)
                    trash_3.append(label_c)
                    trash_3.append(label_d)
                    trash_3.append(label_e)
                    trash_3.append(label_f)
                    trash_3.append(value_a)
                    trash_3.append(value_b)
                    trash_3.append(value_c)
                    trash_3.append(value_d)
                    trash_3.append(value_e)
                    trash_3.append(value_f)

                for item in trash_3:
                    trash.append(item)

        global trash_1, trash_2, trash_3
        global modevars
        trash_1 = []
        trash_2 = []
        trash_3 = []
        modevars = []
        if group_count == 1:
            modevar_1 = tk.StringVar()
            modecombo_1 = ttk.Combobox(labelFrame_4, values=["Soil type", "Grain structure", "van Genuchten"],
                                     state='readonly', textvariable=modevar_1)
            modecombo_1.grid(column=1, row=2)
            modecombo_1.bind("<<ComboboxSelected>>", lambda eff: parameters(eff=None, ID=1, text=modevar_1.get()))
            trash.append(modecombo_1)
            modevars.append(modevar_1)

        if group_count == 2:
            modevar_1 = tk.StringVar()
            modecombo_1 = ttk.Combobox(labelFrame_4, values=["Soil type", "Grain structure", "van Genuchten"],
                                       state='readonly', textvariable=modevar_1)
            modecombo_1.grid(column=1, row=2)
            modecombo_1.bind("<<ComboboxSelected>>", lambda eff: parameters(eff=None, ID=1, text=modevar_1.get()))
            trash.append(modecombo_1)
            modevar_2 = tk.StringVar()
            modecombo_2 = ttk.Combobox(labelFrame_4, values=["Soil type", "Grain structure", "van Genuchten"],
                                       state='readonly', textvariable=modevar_2)
            modecombo_2.grid(column=3, row=2)
            modecombo_2.bind("<<ComboboxSelected>>", lambda eff: parameters(eff=None, ID=2, text=modevar_2.get()))
            trash.append(modecombo_2)
            modevars.append(modevar_1)
            modevars.append(modevar_2)

        if group_count == 3:
            modevar_1 = tk.StringVar()
            modecombo_1 = ttk.Combobox(labelFrame_4, values=["Soil type", "Grain structure", "van Genuchten"],
                                       state='readonly', textvariable=modevar_1)
            modecombo_1.grid(column=1, row=2)
            modecombo_1.bind("<<ComboboxSelected>>", lambda eff: parameters(eff=None, ID=1, text=modevar_1.get()))
            trash.append(modecombo_1)
            modevar_2 = tk.StringVar()
            modecombo_2 = ttk.Combobox(labelFrame_4, values=["Soil type", "Grain structure", "van Genuchten"],
                                       state='readonly', textvariable=modevar_2)
            modecombo_2.grid(column=3, row=2)
            modecombo_2.bind("<<ComboboxSelected>>", lambda eff: parameters(eff=None, ID=2, text=modevar_2.get()))
            trash.append(modecombo_2)
            modevar_3 = tk.StringVar()
            modecombo_3 = ttk.Combobox(labelFrame_4, values=["Soil type", "Grain structure", "van Genuchten"],
                                       state='readonly', textvariable=modevar_3)
            modecombo_3.grid(column=5, row=2)
            modecombo_3.bind("<<ComboboxSelected>>", lambda eff: parameters(eff=None, ID=3, text=modevar_3.get()))
            trash.append(modecombo_3)
            modevars.append(modevar_1)
            modevars.append(modevar_2)
            modevars.append(modevar_3)

    labelFrame_3 = ttk.LabelFrame(b1, text='Model geometry and mesh:')
    labelFrame_3.grid(column=0, row=2, padx=20, pady=20, sticky=tk.W)
    group_nr = tk.StringVar(b1)
    ttk.Label(labelFrame_3, text='Select number of geometry layers: ').grid(column=0, row=0)
    group_combo = ttk.Combobox(labelFrame_3, values=["1", "2", "3"], state='readonly', textvariable=group_nr)
    group_combo.grid(column=1, row=0)
    global trash
    trash = []
    group_combo.bind("<<ComboboxSelected>>", geometry)
    element_size = tk.StringVar(b1, value='0')
    ttk.Label(labelFrame_3, text='Element size: ').grid(column=0, row=1)
    el_s_entry = ttk.Entry(labelFrame_3, textvariable=element_size)
    el_s_entry.grid(column=1, row=1)

    def byebye():
        b1.destroy()

    def yaml_writer():
        file_path = (base_path / "../inputs/UZ_conf_file.yaml").resolve()
        global all_nucl_data
        global heights
        global modevars
        global sand_perc_1, sand_perc_2, sand_perc_3, silt_perc_1, silt_perc_2, silt_perc_3, clay_perc_1, \
            clay_perc_2, clay_perc_3, dens_1, dens_2, dens_3, soil_sel_1, soil_sel_2, soil_sel_3, theta_r_1, \
            theta_r_2, theta_r_3, theta_s_1, theta_s_2, theta_s_3, alpha_1, alpha_2, alpha_3, en_1, en_2, en_3, \
            ks_1, ks_2, ks_3, v_dens_1, v_dens_2, v_dens_3

        dictfile = {'simulation_parameters': {'simulation_time': float(period_e.get()), 'Dt': float(step_e.get()),
                                              'flow_iteration_count': float(iter_e.get()),
                                              'output_step_time': float(out_e.get())}}
        with open(file_path, 'w') as yaml_file:
            yaml.dump(dictfile, yaml_file, Dumper=MyDumper, default_flow_style=False, sort_keys=False)
        dictfile = {'units': {'length': length_c.get(), 'mass': mass_c.get(), 'time': time_c.get()}}
        with open(file_path, 'a') as yaml_file:
            yaml.dump(dictfile, yaml_file, Dumper=MyDumper, default_flow_style=False, sort_keys=False)

        dictfile = {'outputs': [{'entity': 'nodes', 'physical_quantity': 'c_water', 'file_format': 'gmesh_v2_ASCII',
                                 'file_name': "../outputs/UZ_results.msh"}]}
        with open(file_path, 'a') as yaml_file:
            yaml.dump(dictfile, yaml_file, Dumper=MyDumper, default_flow_style=False, sort_keys=False)

        h = 0
        for i in range(0, len(heights)):
            h = h + float(heights[i].get())
        dicts = []
        for i in range(0, len(heights)):
            if i == 0:
                if modevars[i].get() == 'Soil type':
                    dict = {'bottom': 0, 'parameters_mode': 'material', 'material': soil_sel_1.get()}
                if modevars[i].get() == 'Grain structure':
                    dict = {'bottom': 0, 'parameters_mode': 'granular_structure', 'sand': float(sand_perc_1.get()),
                            'silt': float(silt_perc_1.get()), 'clay': float(clay_perc_1.get()),
                            'density_kg_m3': float(dens_1.get())}
                if modevars[i].get() == 'van Genuchten':
                    dict = {'bottom': 0, 'parameters_mode': 'van_genuchten', 'theta_r': float(theta_r_1.get()),
                            'theta_s': float(theta_s_1.get()), 'alpha': float(alpha_1.get()), 'n': float(en_1.get()),
                            'Ks': float(ks_1.get()), 'density_kg_m3': float(v_dens_1.get())}
                dicts.append(dict)
            if i == 1:
                if modevars[i].get() == 'Soil type':
                    dict = {'bottom': float(heights[0].get()), 'parameters_mode': 'material',
                            'material': soil_sel_2.get()}
                if modevars[i].get() == 'Grain structure':
                    dict = {'bottom': float(heights[0].get()), 'parameters_mode': 'granular_structure',
                            'sand': float(sand_perc_2.get()), 'silt': float(silt_perc_2.get()),
                            'clay': float(clay_perc_2.get()), 'density_kg_m3': float(dens_2.get())}
                if modevars[i].get() == 'van Genuchten':
                    dict = {'bottom': float(heights[0].get()), 'parameters_mode': 'van_genuchten',
                            'theta_r': float(theta_r_2.get()), 'theta_s': float(theta_s_2.get()),
                            'alpha': float(alpha_2.get()), 'n': float(en_2.get()), 'Ks': float(ks_2.get()),
                            'density_kg_m3': float(v_dens_2.get())}
                dicts.append(dict)
            if i == 2:
                if modevars[i].get() == 'Soil type':
                    dict = {'bottom': float(heights[0].get()) + float(heights[1].get()),
                            'parameters_mode': 'material', 'material': soil_sel_3.get()}
                if modevars[i].get() == 'Grain structure':
                    dict = {'bottom': float(heights[0].get()) + float(heights[1].get()),
                            'parameters_mode': 'granular_structure', 'sand': float(sand_perc_3.get()),
                            'silt': float(silt_perc_3.get()), 'clay': float(clay_perc_3.get()),
                            'density_kg_m3': float(dens_3.get())}
                if modevars[i].get() == 'van Genuchten':
                    dict = {'bottom': float(heights[0].get()) + float(heights[1].get()),
                            'parameters_mode': 'van_genuchten', 'theta_r': float(theta_r_3.get()),
                            'theta_s': float(theta_s_3.get()), 'alpha': float(alpha_3.get()),
                            'n': float(en_3.get()), 'Ks': float(ks_3.get()), 'density_kg_m3': float(v_dens_3.get())}
                dicts.append(dict)

        dictfile = {'mesh': {'element_height': float(element_size.get()), 'height': h, 'horizons': dicts}}
        with open(file_path, 'a') as yaml_file:
            yaml.dump(dictfile, yaml_file, Dumper=MyDumper, default_flow_style=False, sort_keys=False)

        with open(top_bc_filename) as tbc_file:
            tbc_file.readline()
            tbcs = []
            for line in tbc_file:
                line = line.rstrip().split(';')
                if int(line[1]) == 0:
                    tbc = {'time': float(line[0]), 'type': 'dirichlet', 'head': float(line[2])}
                    tbcs.append(tbc)
                if int(line[1]) == 1:
                    tbc = {'time': float(line[0]), 'type': 'neumann', 'flux': float(line[2])}
                    tbcs.append(tbc)

        with open(bottom_bc_filename) as bbc_file:
            bbc_file.readline()
            bbcs = []
            for line in bbc_file:
                line = line.rstrip().split(';')
                if int(line[1]) == 0:
                    bbc = {'time': float(line[0]), 'type': 'dirichlet', 'head': float(line[2])}
                    bbcs.append(bbc)
                if int(line[1]) == 1:
                    bbc = {'time': float(line[0]), 'type': 'neumann', 'flux': float(line[2])}
                    bbcs.append(bbc)
        sources = [{'bottom': float('0'), 'flux_of_height_unit': float('0')}]
        dictfile = {'flow': {'top_boundary_conditions': tbcs, 'bottom_boundary_conditions': bbcs, 'initial_conditions':
            [{'top_head': float(ic_top.get())}, {'bottom': float('0'), 'head': float(ic_bottom.get())}],
                             'sources': sources}}
        with open(file_path, 'a') as yaml_file:
            yaml.dump(dictfile, yaml_file, Dumper=MyDumper, default_flow_style=False, sort_keys=False)

        sz_h = ((float(ic_bottom.get()))/abs(float(ic_top.get())-float(ic_bottom.get())))*h

        isotopes = []
        for i in range(0, len(all_nucl_data)):
            isotope = {'name': all_nucl_data[i][0], 'diff_coef_m2_s': float(all_nucl_data[i][1].get()),
                       'dist_coef_m3_kg': float(all_nucl_data[i][2].get())}
            isotopes.append(isotope)

        decays = []
        for i in range(0, len(all_nucl_data)-1):
            decay = {'isotope': all_nucl_data[i][0], 'new_isotope': all_nucl_data[i][4].get(),
                     'half_life': float(all_nucl_data[i][3].get())}
            decays.append(decay)

        trans_top_BCs = []
        for i in range(0, len(all_nucl_data)):
            trans_top_BC = {'isotope': all_nucl_data[i][0], 'time_function': [{'time': float('0'),
                                                                               'c_flux': float('0')}]}
            trans_top_BCs.append(trans_top_BC)

        trans_bottom_BCs = []
        sz_concs = []
        with open((base_path / "../inputs/saturated_concentrations.csv").resolve(), 'r') as conc_file:
            conc_file.readline()
            conc_file.readline()
            for line in conc_file:
                line = line.rstrip().split(';')
                line = [float(i) for i in line]
                sz_concs.append(line)
        for i in range(0, len(all_nucl_data)-1):
            conc_list = []
            for j in range(0, len(sz_concs)):
                conc = {'time': sz_concs[j][0], 'c_flux': sz_concs[j][i+1]}
                conc_list.append(conc)
            trans_bottom_BC = {'isotope': all_nucl_data[i][0], 'time_function': conc_list}
            trans_bottom_BCs.append(trans_bottom_BC)

        trans_bottom_BC = {'isotope': 'residual', 'time_function': [{'time': float('0'), 'c_flux': float('0')}]}
        trans_bottom_BCs.append(trans_bottom_BC)

        trans_ics = []
        for i in range(0, len(all_nucl_data) - 1):
            trans_ic = {'isotope': all_nucl_data[i][0], 'concentration_in_water': [{'bottom': float('0'),
                                                                                    'c': sz_concs[0][i+1]},
                                                                                   {'bottom': sz_h, 'c': float('0')}]}
            trans_ics.append(trans_ic)
        trans_ic = {'isotope': 'residual', 'concentration_in_water': [{'bottom': float('0'), 'c': float('0')}]}
        trans_ics.append(trans_ic)

        options = {'tortuosity': is_tort.get(), 'dispersivity': float(disp_val.get()),
                   'numerical_scheme': num_scheme.get(), 'isotopes': isotopes, 'isotopes_half_life': decays,
                   'top_boundary_conditions': trans_top_BCs, 'bottom_boundary_conditions': trans_bottom_BCs,
                   'saturated_zone_concentration': {'apply': 'yes', 'height': sz_h}, 'initial_conditions': trans_ics}
        dictfile = {'transport': options}
        with open(file_path, 'a') as yaml_file:
            yaml.dump(dictfile, yaml_file, Dumper=MyDumper, default_flow_style=False, sort_keys=False)
        label_text.set('Computing flow and transport in the unsaturated zone.\nThis may take a while, please wait...')
        main.focus_force()
        execute_Biorad1()
        proc_UZ_res(length_c.get(), mass_c.get(), h)
        b1.destroy()
        label_text.set('Unsaturated zone simulation complete.')

    labelFrame_5 = ttk.LabelFrame(b1, text='Actions:')
    labelFrame_5.grid(column=0, row=4, padx=20, pady=20, sticky=tk.W)
    ttk.Button(labelFrame_5, text='Quit', command=byebye).pack(side=tk.LEFT)
    ttk.Button(labelFrame_5, text='Write input YAML and run', command=yaml_writer).pack(side=tk.LEFT)

    def load_top_bc():
        global top_bc_filename
        top_bc_filename = askopenfilename(title="Select CSV file", filetypes=[("CSV", "*.csv")])
        top_BC_label.configure(text=top_bc_filename)
        b1.focus_force()

    def load_bottom_bc():
        global bottom_bc_filename
        bottom_bc_filename = askopenfilename(title="Select CSV file", filetypes=[("CSV", "*.csv")])
        bottom_BC_label.configure(text=bottom_bc_filename)
        b1.focus_force()

    labelFrame_6 = ttk.LabelFrame(b1, text='Flow boundary conditions (TOP):')
    labelFrame_6.grid(column=1, row=0, padx=20, pady=20, sticky=tk.W)
    load_button_1 = ttk.Button(labelFrame_6, text="Browse a File", command=load_top_bc)
    load_button_1.grid(column=0, row=0)
    top_BC_label = ttk.Label(labelFrame_6, text='')
    top_BC_label.grid(column=0, row=1)

    labelFrame_7 = ttk.LabelFrame(b1, text='Flow boundary conditions (BOTTOM):')
    labelFrame_7.grid(column=1, row=1, padx=20, pady=20, sticky=tk.W)
    load_button_2 = ttk.Button(labelFrame_7, text="Browse a File", command=load_bottom_bc)
    load_button_2.grid(column=0, row=0)
    bottom_BC_label = ttk.Label(labelFrame_7, text='')
    bottom_BC_label.grid(column=0, row=1)

    labelFrame_8 = ttk.LabelFrame(b1, text='Flow initial conditions:')
    labelFrame_8.grid(column=1, row=2, padx=20, pady=20, sticky=tk.W)
    ttk.Label(labelFrame_8, text='Pressure head at the TOP: ').grid(column=0, row=0)
    ttk.Label(labelFrame_8, text='(must be negative)').grid(column=0, row=1)
    ttk.Label(labelFrame_8, text='Pressure head at the BOTTOM: ').grid(column=0, row=2)
    ttk.Label(labelFrame_8, text='(must be positive)').grid(column=0, row=3)
    ic_top = tk.StringVar(b1, value='0')
    ic_bottom = tk.StringVar(b1, value='0')
    ict_entry = ttk.Entry(labelFrame_8, textvariable=ic_top)
    ict_entry.grid(column=1, row=0)
    icb_entry = ttk.Entry(labelFrame_8, textvariable=ic_bottom)
    icb_entry.grid(column=1, row=2)

    labelFrame_9 = ttk.LabelFrame(b1, text='Transport simulation options:')
    labelFrame_9.grid(column=2, row=0, padx=20, pady=20, sticky=tk.W)
    ttk.Label(labelFrame_9, text='Include tortuosity: ').grid(column=0, row=0)
    ttk.Label(labelFrame_9, text='Enter dispersivity: ').grid(column=0, row=1)
    ttk.Label(labelFrame_9, text='Select numerical scheme: ').grid(column=0, row=2)
    is_tort = tk.StringVar(b1, value='yes')
    chck_tort = ttk.Checkbutton(labelFrame_9, variable=is_tort, onvalue='yes', offvalue='no')
    chck_tort.grid(column=1, row=0)
    disp_val = tk.StringVar(b1, value='0')
    disp_entry = ttk.Entry(labelFrame_9, textvariable=disp_val)
    disp_entry.grid(column=1, row=1)
    num_scheme = tk.StringVar(b1, value='implicit')
    num_combo = ttk.Combobox(labelFrame_9, values=["implicit", "explicit", "crank_nicolson"], state='readonly',
                             textvariable=num_scheme)
    num_combo.grid(column=1, row=2)

    def param_win():
        global all_nucl_data

        def save_close():
            with open((base_path / "../inputs/nuclide_parameters.csv").resolve(), 'w') as nucl_file:
                nucl_file.write('Nuclide; De [m2/s]; Kd [m3/kg]; Half-life [units of time]; Product\n')
                for i in range(0, len(all_nucl_data)):
                    if i != len(all_nucl_data)-1:
                        nucl_file.write(all_nucl_data[i][0] + ';' + (all_nucl_data[i][1].get()) + ';' +
                                        (all_nucl_data[i][2].get()) + ';' + (all_nucl_data[i][3].get()) +
                                        ';' + all_nucl_data[i][4].get() + '\n')
                    else:
                        nucl_file.write(all_nucl_data[i][0] + ';' + (all_nucl_data[i][1].get()) + ';' +
                                        (all_nucl_data[i][2].get()) + '\n')
            nuclide_label.configure(text='Done.')
            b1.focus_force()
            b1_params.destroy()

        def load_ff():
            filename = askopenfilename(filetypes=[("CSV", "*.csv")])
            with open(filename, 'r') as loaded_f:
                loaded_data = []
                loaded_products = []
                loaded_f.readline()
                for line in loaded_f:
                    line = line.rstrip().split(';')
                    if len(line) == 5:
                        if line[4] not in products:
                            param_win_button_1_label.configure(text='File could not be loaded.\n(data not compatible)')
                            b1_params.focus_force()
                            return
                    loaded_products.append(line[0])
                    loaded_data.append(line)
            if loaded_products != products:
                param_win_button_1_label.configure(text='File could not be loaded.\n(data not compatible)')
                b1_params.focus_force()
                return
            for i in range(0, len(loaded_products)):
                item = loaded_data[i]
                for j in range(1, len(item)):
                    all_nucl_data[i][j].set(loaded_data[i][j])
            b1_params.focus_force()

        nuclides = []
        products = []
        with open((base_path / "../outputs/temp.txt").resolve(), 'r') as nucl_file:
            for line in nucl_file:
                line=line.rstrip()
                nuclides.append(line)
                products.append(line)
        b1_params = tk.Toplevel()
        b1_params.title('Unsaturated zone nuclide tranport parameters')
        b1_params.geometry("800x600")

        lF_1 = ttk.LabelFrame(b1_params)
        lF_1.grid(column=0, row=0, padx=20, pady=20, sticky=tk.W)
        param_win_button_1 = ttk.Button(lF_1, text='Load from file', command=load_ff)
        param_win_button_1.grid(row=0, column=0)
        param_win_button_1_label = ttk.Label(lF_1, text='')
        param_win_button_1_label.grid(row=1, column=0)
        param_win_button_2 = ttk.Button(lF_1, text='Save and close', command=save_close)
        param_win_button_2.grid(row=0, column=1)
        lF_2 = ttk.LabelFrame(b1_params)
        lF_2.grid(column=0, row=1, padx=20, pady=20, sticky=tk.W, columnspan=3)
        frame_1 = ScrollableFrame(lF_2)
        ttk.Label(frame_1.scrollable_frame, text='Nuclide').grid(column=0, row=0, padx=3, pady=3)
        ttk.Label(frame_1.scrollable_frame, text='De [m2/s]').grid(column=1, row=0, padx=3, pady=3)
        ttk.Label(frame_1.scrollable_frame, text='Kd [m3/kg]').grid(column=2, row=0, padx=3, pady=3)
        ttk.Label(frame_1.scrollable_frame, text='Half-life [units of time]').grid(column=3, row=0, padx=3, pady=3)
        ttk.Label(frame_1.scrollable_frame, text='Product').grid(column=4, row=0, padx=3, pady=3)
        all_nucl_data = []
        counter = 0
        products.append('residual')
        for item in nuclides:
            counter += 1
            nucl_data = []
            ttk.Label(frame_1.scrollable_frame, text=item).grid(column=0, row=counter, padx=3, pady=3)
            nucl_data.append(item)
            De = tk.StringVar(b1_params, value = '0')
            Kd = tk.StringVar(b1_params, value = '0')
            h_l = tk.StringVar(b1_params, value = '0')
            prod = tk.StringVar(b1_params, value = 'residual')
            e = ttk.Entry(frame_1.scrollable_frame, textvariable=De)
            e.grid(column=1, row=counter, padx=3, pady=3)
            e = ttk.Entry(frame_1.scrollable_frame, textvariable=Kd)
            e.grid(column=2, row=counter, padx=3, pady=3)
            e = ttk.Entry(frame_1.scrollable_frame, textvariable=h_l)
            e.grid(column=3, row=counter, padx=3, pady=3)
            c = ttk.Combobox(frame_1.scrollable_frame, values=[i for i in products if i != item],
                             state='readonly', textvariable=prod)
            c.grid(column=4, row=counter, padx=3, pady=3)
            nucl_data.append(De)
            nucl_data.append(Kd)
            nucl_data.append(h_l)
            nucl_data.append(prod)
            all_nucl_data.append(nucl_data)
        nucl_data = []
        ttk.Label(frame_1.scrollable_frame, text='residual').grid(column=0, row=counter+1, padx=3, pady=3)
        De = tk.StringVar(b1_params, value='0')
        Kd = tk.StringVar(b1_params, value='0')
        e = ttk.Entry(frame_1.scrollable_frame, textvariable=De)
        e.grid(column=1, row=counter+1, padx=3, pady=3)
        e = ttk.Entry(frame_1.scrollable_frame, textvariable=Kd)
        e.grid(column=2, row=counter+1, padx=3, pady=3)
        nucl_data.append('residual')
        nucl_data.append(De)
        nucl_data.append(Kd)
        all_nucl_data.append(nucl_data)
        frame_1.grid(column=0, row=0)

    labelFrame_10 = ttk.LabelFrame(b1, text='Radionuclide data:')
    labelFrame_10.grid(column=2, row=1, padx=20, pady=20, sticky=tk.W)
    nuclide_button = ttk.Button(labelFrame_10, text='Specify', command=param_win)
    nuclide_button.grid(column=0, row=0)
    nuclide_label = ttk.Label(labelFrame_10, text='')
    nuclide_label.grid(column=0, row=1)


def call_biorad_2():
    class B2(tk.Tk):
        def __init__(self):
            super(B2, self).__init__()
            global scenario_var
            self.title("Biosphere module user options.")
            self.minsize(800, 600)
            self.labelFrame_1 = ttk.LabelFrame(self, text="Open database file:")
            self.labelFrame_1.grid(column=0, row=1, padx=20, pady=20, sticky=tk.W)
            but_1 = ttk.Button(self.labelFrame_1, text="Browse A File", command=self.DB_fileDialog)
            but_1.grid(column=1, row=1)

            self.labelFrame_2 = ttk.LabelFrame(self, text="Select input scenario:")
            self.labelFrame_2.grid(column=0, row=2, padx=20, pady=20, sticky=tk.W)
            scenario_var = tk.IntVar(self, 0)
            tk.Radiobutton(self.labelFrame_2, text="Saturated zone only.", variable=scenario_var, value=0)\
                .grid(row=2, sticky=tk.W)
            tk.Radiobutton(self.labelFrame_2, text="Saturated and unsaturated zone.", variable=scenario_var, value=1)\
                .grid(row=3, sticky=tk.W)

            self.labelFrame_3 = ttk.LabelFrame(self, text="Select a basket type:")
            self.labelFrame_3.grid(column=0, row=3, padx=20, pady=20, sticky=tk.W)

            ttk.Button(self, text="Quit biosphere module", command=self.cl).grid(column=3, padx=20, pady=20, row=3)


        def DB_fileDialog(self):
            self.filename = askopenfilename(title="Select database file", filetypes=[("Database file", "*.db")])
            self.label = ttk.Label(self.labelFrame_1, text="")
            self.label.grid(column=1, row=2)
            self.label.configure(text=self.filename)
            DB_filename = self.filename
            print(DB_filename)
            conn = sq.connect(DB_filename)
            baskets = []
            for row in conn.execute("SELECT identificator FROM Basket_type"):
                baskets.append(row[0])
            basket_var = tk.StringVar(b2, baskets[0])
            for i in range(0, len(baskets)):
                tk.Radiobutton(b2.labelFrame_3, text=baskets[i], variable=basket_var, value=baskets[i], indicatoron=0)\
                    .pack(anchor=tk.CENTER)
            nucl_file = open((base_path / "../outputs/temp.txt").resolve(), 'r')
            fields = []
            for line in nucl_file:
                line = line.rstrip()
                fields.append(line)
            nucl_file.close()
            self.write_button = ttk.Button(b2, text="Write configuration file and run Biosphere module.",
                                           command=lambda: self.wr(fields, basket_var.get(), DB_filename))\
                .grid(column=3, padx=20, pady=20, row=2)
            b2.focus_force()

        def cl(self):
            b2.quit()
            b2.destroy()

        def wr(self, fields, bb, DB_file):
            dictfile = {'nuclides': fields}
            with open((base_path / "../inputs/Biosphere_configuration_file.yaml").resolve(), 'w') as yaml_file:
                yaml.dump(dictfile, yaml_file, Dumper=MyDumper, default_flow_style=False, sort_keys=False)
            with open((base_path / "../inputs/saturated_concentrations.csv").resolve(), 'r') as saturated_results:
                lines = []
                times = []
                concs = []
                msh_el = [int(saturated_results.readline().rstrip()[1:-1])]
                header = saturated_results.readline().rstrip().split(';')
                if len(header) != len(fields)+1:
                    print('WARNING! Nuclide counts do not match"')
                for line in saturated_results:
                    line = line.rstrip()
                    lines.append(line)
                for item in lines:
                    item = item.split(';')
                    item = [float(i) for i in item]
                    times.append(item[0])
                    conc_list = item[1:]
                    conc_list = [x * multi for x in conc_list]
                    concs.append(conc_list)
                dictfile = {'msh_el': msh_el}
                with open((base_path / "../inputs/Biosphere_configuration_file.yaml").resolve(), 'a') as yaml_file:
                    yaml.dump(dictfile, yaml_file, Dumper=MyDumper, default_flow_style=False, sort_keys=False)
            if scenario_var.get() == 1:
                with open((base_path / "../inputs/unsaturated_concentrations.csv").resolve(), 'r') as unsaturated_results:
                    lines = []
                    times_UZ = []
                    concs_UZ = []
                    header = unsaturated_results.readline().rstrip().split(';')
                    if len(header) != len(fields) + 1:
                        print('WARNING! Nuclide counts do not match"')
                    for line in unsaturated_results:
                        line = line.rstrip()
                        lines.append(line)
                    for item in lines:
                        item = item.split(';')
                        item = [float(i) for i in item]
                        times_UZ.append(item[0])
                        conc_list = item[1:]
                        conc_list = [x * multi for x in conc_list]
                        concs_UZ.append(conc_list)

            dictfile = {'baskets': [bb]}
            with open((base_path / "../inputs/Biosphere_configuration_file.yaml").resolve(), 'a') as yaml_file:
                yaml.dump(dictfile, yaml_file, Dumper=MyDumper, default_flow_style=False, sort_keys=False)
            if is_activities:
                dictfile = {'contamination_data_type': 'activity'}
            else:
                dictfile = {'contamination_data_type': 'concentration'}
            with open((base_path / "../inputs/Biosphere_configuration_file.yaml").resolve(), 'a') as yaml_file:
                yaml.dump(dictfile, yaml_file, Dumper=MyDumper, default_flow_style=False, sort_keys=False)
            if scenario_var.get()==0:
                dictfile = {'activity_type': 'gw'}
            else:
                dictfile = {'activity_type': 'water'}
            with open((base_path / "../inputs/Biosphere_configuration_file.yaml").resolve(), 'a') as yaml_file:
                yaml.dump(dictfile, yaml_file)
            dictfile = {'activity_file_name': ["../inputs/concentration_file.yaml"]}
            with open((base_path / "../inputs/Biosphere_configuration_file.yaml").resolve(), 'a') as yaml_file:
                yaml.dump(dictfile, yaml_file, Dumper=MyDumper, default_flow_style=False, sort_keys=False)
            dictfile = {'database_file_name': [DB_file]}
            with open((base_path / "../inputs/Biosphere_configuration_file.yaml").resolve(), 'a') as yaml_file:
                yaml.dump(dictfile, yaml_file, Dumper=MyDumper, default_flow_style=False, sort_keys=False)
            dictfile = {'results_file_name': ["../outputs/Biosphere_module_results.npz"]}
            with open((base_path / "../inputs/Biosphere_configuration_file.yaml").resolve(), 'a') as yaml_file:
                yaml.dump(dictfile, yaml_file, Dumper=MyDumper, default_flow_style=False, sort_keys=False)

            dictfile = {'times': times}
            with open((base_path / "../inputs/concentration_file.yaml").resolve(), 'w') as yaml_conc:
                yaml.dump(dictfile, yaml_conc, Dumper=MyDumper, default_flow_style=False, sort_keys=False)

            blocks = {}
            for i in range(0, len(fields)):
                conc = []
                for j in range(0, len(concs)):
                    conc.append(concs[j][i])
                blocks[fields[i]] = conc
            if scenario_var.get() == 0:
                dictfile = {'elements': {str(msh_el[0]): {'water_ground': blocks}}}
                with open((base_path / "../inputs/concentration_file.yaml").resolve(), 'a') as yaml_conc:
                    yaml.dump(dictfile, yaml_conc, Dumper=MyDumper, default_flow_style=False, sort_keys=False)
            else:
                blocks_UZ = {}
                concs_UZ_time_ax = []
                for t in times:
                    if t in times_UZ:
                        index = times_UZ.index(t)
                        concs_UZ_time_ax.append(concs_UZ[index])
                    elif t < times_UZ[0]:
                        zero_c = [0 for i in range(0, len(fields))]
                        concs_UZ_time_ax.append(zero_c)
                    elif t > times_UZ[-1]:
                        concs_UZ_time_ax.append(concs_UZ[-1])
                    else:
                        for i in range(0, len(times_UZ)):
                            if t < times_UZ[i]:
                                diff = t - times_UZ[i-1]
                                interpol_c = []
                                for j in range(0, len(fields)):
                                    interpol_c.append(concs_UZ[i-1][j]+(diff/(times_UZ[i]-times_UZ[i-1]))*
                                                      (concs_UZ[i][j]-concs_UZ[i-1][j]))
                                concs_UZ_time_ax.append(interpol_c)
                                break
                for i in range(0, len(fields)):
                    conc = []
                    for j in range(0, len(concs_UZ_time_ax)):
                        conc.append(concs_UZ_time_ax[j][i])
                    blocks_UZ[fields[i]] = conc
                dictfile = {'elements': {str(msh_el[0]): {'water_ground': blocks, 'water_surface': blocks_UZ}}}
                with open((base_path / "../inputs/concentration_file.yaml").resolve(), 'a') as yaml_conc:
                    yaml.dump(dictfile, yaml_conc, Dumper=MyDumper, default_flow_style=False, sort_keys=False)
            log_file = open((base_path / "../outputs/biosphere_module_log.txt").resolve(), 'w')
            brl.biorad2_v05((base_path / "../inputs/Biosphere_configuration_file.yaml").resolve(), log_file, 0)
            log_file.close()
            label_text.set('Biosphere module computations completed.\n')
            main.focus_force()
            visualize_button = ttk.Button(main, text="Process and visualize Biosphere module results.",
                                          command=show_me).pack()
            b2.quit()
            b2.destroy()

    def run_b2():
        global is_activities
        is_activities = activity_switch.get()
        mass_unit = mass_d.get()
        length_unit = length_d.get()
        global multi
        multi = 1
        if length_unit == 'mm':
            multi = multi * 1e9
        elif length_unit == 'cm':
            multi = multi * 1e6
        elif length_unit == 'dm':
            multi = multi * 1e3
        if mass_unit == 'ng':
            multi = multi / 1e12
        elif mass_unit == 'ug':
            multi = multi / 1e9
        elif mass_unit == 'mg':
            multi = multi / 1e6
        elif mass_unit == 'g':
            multi = multi / 1e3
        if is_activities:
            multi = 1
        global time_unit
        time_unit = time_d.get()
        units.quit()
        units.destroy()
        global b2
        b2 = B2()
        b2.mainloop()

    units = tk.Toplevel()
    units.title('Select units of input files:')
    units_label_text = tk.StringVar()
    units_label_text.set('This will take care of unit conversion. Select units of input files (concentration evolutions'
                         ' in saturated and unsaturated zones). Biosphere module computes in units of kg/m3 or '
                         'Bq/kg.\n\n'
                         'If you check the Inputs are in Bq/kg field you may ignore the selection of mass '
                         'and length units.')
    units_label = tk.Label(units, textvariable=units_label_text)
    units_label.grid(column=0, padx=20, pady=20, row=0, columnspan=2)
    labelFrame_1 = ttk.LabelFrame(units, text='Select units:')
    labelFrame_1.grid(column=0, row=1, padx=20, pady=20, sticky=tk.W)
    ttk.Label(labelFrame_1, text='Inputs are in Bq/kg: ').grid(column=0, row=0)
    activity_switch = tk.BooleanVar(units, value=False)
    x2 = ttk.Checkbutton(labelFrame_1, variable=activity_switch)
    x2.grid(column=1, row=0)
    ttk.Label(labelFrame_1, text='Length: ').grid(column=0, row=1)
    length_d = tk.StringVar(units, value="m")
    x3 = ttk.Combobox(labelFrame_1, values=["mm", "cm", "dm", "m"], state='readonly', textvariable=length_d)
    x3.grid(column=1, row=1)
    mass_d = tk.StringVar(units, value="kg")
    ttk.Label(labelFrame_1, text='Mass: ').grid(column=0, row=2)
    x4 = ttk.Combobox(labelFrame_1, values=["ng", "ug", "mg", "g", "kg"], state='readonly', textvariable=mass_d)
    x4.grid(column=1, row=2)
    time_d = tk.StringVar(units, value="s")
    ttk.Label(labelFrame_1, text='Time: ').grid(column=0, row=3)
    x5 = ttk.Combobox(labelFrame_1, values=["s", "h", "day", "year"], state='readonly', textvariable=time_d)
    x5.grid(column=1, row=3)
    x0 = ttk.Button(units, text='OK', command=run_b2)
    x0.grid(column=1, padx=20, pady=20, row=1)


main = tk.Tk()
main.geometry("500x500")
main.title('BIORAD_GUI')
label_text = tk.StringVar()
label_text.set("Messages will be shown here.")
label = tk.Label(main, textvariable=label_text).pack()

# this part creates SW menu
hlavniMenu = tk.Menu(main)

menuSoubor = tk.Menu(hlavniMenu, tearoff=0)
menuSoubor.add_command(label="Load Flow123d result", command=open_transport_file)
menuSoubor.add_command(label="Unsaturated zone module", command=call_biorad_1)
menuSoubor.add_command(label='Biosphere module', command=call_biorad_2)
menuSoubor.add_separator()
menuSoubor.add_command(label="Quit", command=main.quit)
hlavniMenu.add_cascade(label="File", menu=menuSoubor)

menuNapoveda = tk.Menu(hlavniMenu, tearoff=0)
menuNapoveda.add_command(label="About", command=what_up)
hlavniMenu.add_cascade(label="Help", menu=menuNapoveda)

main.config(menu=hlavniMenu)

main.mainloop()
