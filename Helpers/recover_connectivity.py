class LigandPDB:
    def __init__(self, pdb_file, ligand_chain="L"):
        self.pdb_file = pdb_file
        self.ligand_chain = ligand_chain
        self.content = self.read_pdb()

    def read_pdb(self):
        with open(self.pdb_file) as pdbf:
            content = pdbf.read()
        return content

    def get_models_dict(self):
        model_list = self.content.split("MODEL")
        model_dict = {}
        for n, model in enumerate(model_list[1:]):
            model_dict[n] = model
        return model_dict

    def get_model_ligand_content(self, model):
        lines = self.get_models_dict()[model].split("\n")
        ligand_lines = []
        ligand_blank = []
        for line in lines:
            if line.startswith("HETATM"):
                if line[21] == " ":
                    ligand_blank.append(line)
                if line[21] == self.ligand_chain or line[21] == " ":
                    ligand_lines.append(line)
        if ligand_blank:
            print("{}\nWARNING: We found that the ligand chain is a blank space. "
                  "Previuos lines will be included in the ligand!".format("\n".join(ligand_blank)))

        return "\n".join(ligand_lines)

    def get_connects(self):
        lines = self.content.split("\n")
        connect_lines = []
        for line in lines:
            if line.startswith("CONECT"):
                connect_lines.append(line)
        return "\n".join(connect_lines)

    def get_model_ligand_name_index_dictionary(self, model):
        ligand_lines = self.get_model_ligand_content(model=model).split("\n")
        lig_dict = {}
        for line in ligand_lines:
            lig_dict[line[12:16].strip()] = line[6:11].strip()
        return lig_dict


def create_index_relation_between_connect_and_to_connect(pdb_connected, pdb2connect):
    pdb_connected_names_dict = pdb_connected.get_model_ligand_name_index_dictionary(0)
    pd2b_connected_names_dict = pdb2connect.get_model_ligand_name_index_dictionary(0)
    index_relations = {}
    for key, value in pd2b_connected_names_dict.items():
        index_relations[pdb_connected_names_dict[key]] = value
    return index_relations


def rebuild_connect_line(connections_list):
    new_connect_pattern = "{:5s}"*len(connections_list)
    new_connect_line = "CONECT " + new_connect_pattern.format(*connections_list)
    return new_connect_line


def recover_connectivity(pdb_connected, pdb_to_connect, ligand_chain_connected="L", ligand_chain_to_connect="L"):
    pdbc = LigandPDB(pdb_file=pdb_connected, ligand_chain=ligand_chain_connected)
    pdb2c = LigandPDB(pdb_file=pdb_to_connect, ligand_chain=ligand_chain_to_connect)
    index_relations = create_index_relation_between_connect_and_to_connect(pdbc, pdb2c)
    connects_pdb_connected = pdbc.get_connects().split("\n")
    connects_lines = []
    for line in connects_pdb_connected:
        connections_list = line.split()[1:]
        for n, index in enumerate(connections_list):
            new_index = index_relations[index]
            connections_list[n] = new_index
        new_connect = rebuild_connect_line(connections_list)
        connects_lines.append(new_connect)
    connects = "\n".join(connects_lines)






print(recover_connectivity(pdb_connected="/home/carles/Almirall/Data2Victor/ligands/CHEMBL1164264.pdb",
                           pdb_to_connect="/home/carles/Almirall/Data2Victor/trajectories/trajectory_1.pdb"))


