def proc_UZ_res(length, mass, height):
    mult = 1
    depth = 0.3
    if length == 'mm':
        mult = mult*1e9
        depth = 300
    elif length == 'cm':
        mult = mult*1e6
        depth = 30
    elif length == 'dm':
        mult = mult*1e3
        depth = 3
    if mass == 'ng':
        mult = mult/1e12
    elif mass == 'ug':
        mult = mult/1e9
    elif mass == 'mg':
        mult = mult/1e6
    elif mass == 'g':
        mult = mult/1e3
    with open('../outputs/UZ_results.msh', 'r') as input_file:
        for line in input_file:
            if line == '$Nodes\n':
                break
        node_count = int(input_file.readline().rstrip())
        nodes = []
        weights = []
        first = False
        for i in range(0, node_count):
            line = input_file.readline().rstrip().split()
            z = float(line[3])
            if z >= (height-depth):
                if not first:
                    first = True
                    nodes.append(id_prev)
                    weights.append((z-(height-depth))/(z-z_prev))
                nodes.append(int(line[0]))
                weights.append(1)
            z_prev = z
            id_prev = int(line[0])
        print(weights)
        print(nodes)
    output_file = open('../inputs/unsaturated_concentrations.csv', 'w')
    output_file.write('time')
    with open('../outputs/temp.txt', 'r') as tracer_file:
        for line in tracer_file:
            output_file.write(';')
            output_file.write(line.rstrip())
        output_file.write('\n')
    with open('../outputs/UZ_results.msh', 'r') as input_file:
        block = []
        save = False
        time = 0.0
        output_file.write(str(time))
        for line in input_file:
            if line.rstrip() == '$NodeData':
                save = True
                block = []
            if line.rstrip() == '$EndNodeData':
                save = False
                if block[2] != 'residual':
                    if float(block[4]) > time:
                        time = float(block[4])
                        output_file.write('\n')
                        output_file.write(str(time))
                    conc = 0
                    for i in range(0, len(nodes)-1):
                        conc += float(block[-(i+1)].split()[1])
                    if weights[0] > 0:
                        conc += float(block[-len(nodes)+1].split()[1]) + \
                                (float(block[-len(nodes)].split()[1])
                                 -float(block[-len(nodes)+1].split()[1]))*weights[0]
                        conc = conc/len(nodes)
                    else:
                        conc = conc/(len(nodes)-1)
                    if conc < 0:
                        conc = 0
                    output_file.write(';' + str(conc))

            if save:
                block.append(line.rstrip())


#proc_UZ_res('m', 'g', 10)
