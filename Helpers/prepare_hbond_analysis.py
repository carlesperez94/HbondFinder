import glob
import templatizer
import string

ATOMS_DICT_LIST = ['{438: "H"}', '{436: "H"}',
                   '{438: "O"}']
ATOM_KEY = ['{438: "H"}']


def main(path_names_patt, template, resname="LIG", cpus=32, angle_original=2.0944, distance_original=0.25):
    names = glob.glob(path_names_patt)
    template = templatizer.read_file(template)
    counter_a = 1
    for atom in ATOMS_DICT_LIST:
        angle = angle_original
        distance = distance_original
        if atom in ATOM_KEY:
            angle = 2.44
            distance= 0.21
        counter_b = 1
        for i in range(2):
            if i == 1:
                pseudo = "False"
            else:
                pseudo = "True"
            for name in names:
                keys = {"NAME":name.split("/")[-1],
                        "ATOM": atom,
                        "PSEUDO": pseudo,
                        "RES": resname,
                        "ANGLE": angle,
                        "DIST": distance,
                        "CPUS": cpus
                       }
                t = string.Template(template)
                replaced_content = t.substitute(keys)
                templatizer.write_file(replaced_content, "hbonds_{}_{}{}.yaml".format(name.split("/")[-1], counter_a, counter_b))
            counter_b += 1
        counter_a += 1


main("/gpfs/projects/bsc72/vs/alm/ITK/validation/running/run/4PP9CHEMBL*", 
     "/gpfs/projects/bsc72/vs/alm/scripts/configurations/templates/template.yaml")
    
    
