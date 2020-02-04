import time
import os
import sys
import glob
import argparse
import multiprocessing as mp

import yaml
import numpy as np
import pandas as pd
import mdtraj as md
import hbond_mod as hm
import atoms_to_explore as a2e

__author__ = "Carles Perez Lopez"
__version__ = "1.0.0"
__maintainer__ = "Carles Perez Lopez"
__email__ = "carlesperez94@gmail.com"


def parseargs_yaml(args=[]):
    parser = argparse.ArgumentParser(description='Program to count hbonds from a list of PDB files.')
    parser.add_argument('input_file', type=str, help='Yaml input file')
    args = parser.parse_args(args) if args else parser.parse_args()
    return args


class YamlParser(object):

    def __init__(self, yamlfile):
        self.yamlfile = yamlfile
        self.parse()

    def parse_yaml(self):
        with open(self.yamlfile, 'r') as stream:
            try:
                data = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                raise(exc)
        return data

    def parse(self):
        data = self.parse_yaml()
        self.path = data.get("path", None)
        self.ligname = data.get("ligname", "LIG")
        self.specific = data.get("specific_hbonds", None)
        self.distance = data.get("distance", 0.25)
        self.angle = data.get("angle", 2.094)
        self.pseudo = data.get("pseudo", False)
        self.outpath = data.get("outpath", "hbond_analysis/")
        self.cpus = data.get("cpus", 8)
        self.top = data.get("topology", None)
        self.report_in = data.get("report_in", "report_")
        self.report_out = data.get("report_out", "report_hb_")


def parse_arguments():
    """
        Parse user arguments
        Output: list with all the user arguments
    """

    parser = argparse.ArgumentParser(description="""Description: Program to count hbonds from a list of PDB files.""")
    required_named = parser.add_argument_group('required named arguments')
    required_named.add_argument("path", type=str, help="Patter to PDB files.")
    parser.add_argument("-l", "--ligname", type=str, default="LIG",
                        help="Resname of the ligand.")
    parser.add_argument("-s", "--specific", type=dict, default=a2e.HBOND_ATOMS,
                        help="Dictionary with the specific atoms to compute H bonding events.")
    parser.add_argument("-d", "--distance", type=float, default=0.25,
                        help="H-X distance threshold (in nm). Default: 0.25nm.")
    parser.add_argument("-a", "--angle", type=float, default=2.0 * np.pi / 3.0,
                        help="Angle threshold (in rad). Default: 2pi/3 (120º).")
    parser.add_argument("-p", "--pseudo", default=False, action="store_true",
                        help="If set it also counts pseudo-hydrogen bonds (C-H)")
    parser.add_argument("-o", "--outpath", default="hbond_analysis",
                        help="Output path. By default: hbond_analysis/")
    parser.add_argument("-cpu", "--cpus", default=4, type=int,
                        help="Number of processors to paralelize the process")
    parser.add_argument("-t", "--top", default=None, type=str,
                        help="Topology path. By default it uses the standard Adaptive location.")
    parser.add_argument("-ro", "--report_out", default="report_hb_", type=str,
                        help="Report's output prefix. Default: 'report_hb_'")
    parser.add_argument("-ri", "--report_in", default="report_", type=str,
                        help="Report's output prefix. Default: 'report_'")

    args = parser.parse_args()

    return args.path, args.ligname, args.specific, args.distance, args.angle, args.pseudo, args.outpath, args.cpus, \
           args.top, args.report_in, args.report_out


def select_ligand(traj, resname="LIG"):
    """
    It selects the atoms with an specific residue name (Usually the ligand)"
    :param traj: mdtraj.Trajectory object.
    :param resname: residue name.
    :type resname: str
    :return: list with selection indexes.
    """
    lig = traj.topology.select("resname '{}'".format(resname))
    return lig


def select_specific_atoms(traj, specific_dict):
    """
    It selects atoms of certain residue number (key of dict) and atomnames (value: list of names).
    :param traj: mdtraj.Trajectory object.
    :param specific_dict: dictionary with residue numbers as keys and values are list of atom names for each residue.
    :return: list of indexes of the selection
    """
    indexes = set()
    for resseq, atoms in specific_dict.items():
        specifics = traj.topology.select("resSeq {} and name {}".format(resseq, " ".join(atoms)))
        for index in specifics:
            indexes.add(index)
    return list(indexes)


def find_hbond_in_snapshot(trajectory, n, hbonds_lig, ligand_indexes, specifics=None, distance=0.25, angle=2.0 * np.pi / 3.0,
                          pseudo=False):
    snapshot = trajectory[n]
    hbonds = hm.baker_hubbard(traj=snapshot, distance=distance, angle=angle, pseudo=pseudo)
    for hbond in hbonds:
        print(snapshot.topology.atom(hbond[0]), snapshot.topology.atom(hbond[2]))
        if not specifics:
            if any(atom in ligand_indexes for atom in hbond) and not all(atom in hbond for atom in ligand_indexes):

                hbonds_lig.setdefault(n, []).append(hbond)
        else:
            if any(atom in ligand_indexes for atom in hbond) and any(atom in specifics for atom in hbond) and \
                   not all(atom in hbond for atom in ligand_indexes):
                hbonds_lig.setdefault(n, []).append(hbond)
                print("SPECIFIC BOND: {} {}".format(snapshot.topology.atom(hbond[0]),snapshot.topology.atom(hbond[2])))
    try:
        print("{} hbonds found in snapshot {}".format(len(hbonds_lig[n]), n))
    except KeyError:
        print("0 hbonds found in snapshot {}".format(n))


def find_ligand_hbonds(trajectory, ligand_indexes, specifics=None, distance=0.25, angle=2.0 * np.pi / 3.0,
                       pseudo=False):
    hbonds_lig = {}
    for n in range(0, trajectory.n_frames):
        find_hbond_in_snapshot(trajectory, n, hbonds_lig, ligand_indexes, specifics, distance, angle, pseudo)
    return hbonds_lig


def get_column_name(atoms_to_explore_dict, pseudo=False):
    """
    It joins all residue numbers and atom names to create a unique column name for each search.
    :param atoms_to_explore_dict: dictionary with keys as residue numbers and values as list of atom names
    :return: column name
    """
    column_name = []
    for resid, atomnames_list in atoms_to_explore_dict.items():
        temp_list = []
        for atom in atomnames_list:
            temp_list.append(atom)
        column_name.append("{}{}".format(resid, "_".join(temp_list)))
    if pseudo:
        column_name.append("ps")
    return "-".join(column_name)


def read_and_rewrite_reports(report, trajectory, hbonds_in_traj, own_path, specific_hbonds,
                             report_in_pref="report_", report_out_pref="report_hb_", sep="   ", pseudo=False):
    report_number = report.split(report_in_pref)[-1]
    # Create a folder for each epoch
    report_epoch = report.split("/")[-2]
    if not os.path.exists(os.path.join(own_path, report_epoch)):
        os.mkdir(os.path.join(own_path, report_epoch))
    # Define output path of the report
    report_out = os.path.join(own_path, report_epoch, "{}{}".format(report_out_pref, report_number))
    if os.path.exists(report_out):
        print("REPORT {} Already exists. Extra data will be appended to the current file!".format(report_out))
        new_report = pd.read_csv(report_out, sep="\t", header=0, engine="python")
    else:
        new_report = pd.read_csv(report, sep=sep, header=0, engine="python")
    column_name = get_column_name(specific_hbonds, pseudo)  # Renaming automatically report columns
    for n in range(0, len(trajectory)):
        try:
            new_report.loc[n, column_name] = len(hbonds_in_traj[n])

        except KeyError:  # To avoid having empty gaps
            new_report.loc[n, column_name] = 0
    #  Writing output file
    new_report.to_csv(report_out, sep="\t", index=False)


def print_info(path_to_pdbs, resname, specifics, distance, angle, pseudo, outpath, proc, top, rep_in, rep_out):
    print("""Path to PDBs: {}\nResidue name: {}\nSpecific Hbonds: {}\nCutoff distance: {}\nCutoff angle: {}
Pseudo mode: {}\nOutput directory: {}\nProcessors: {}\nTopology: {}\nReport input prefix: {}
Report output prefix: {}\n================="""
          .format(str(path_to_pdbs), str(resname), str(specifics), str(distance), str(angle),
          str(pseudo), str(outpath), str(proc), str(top), str(rep_in), str(rep_out)))


def find_hbonds_in_pdb(pdb, resname, report, specifics=None, distance=0.25, angle=2.0 * np.pi / 3.0,
                       pseudo=False, own_path=None, top=None, report_in_pref="report_", report_out_pref="report_hb_"):

    path_pdb = os.path.abspath(pdb)
    if not top:
        top_standard = "../../topologies/topology_0.pdb"
        top = os.path.normpath(os.path.join(path_pdb, top_standard)) #  'top' variable is not currently used, but could
                                                                     #  be useful when reading xtc format
    traj = md.load(pdb)
    lig = select_ligand(traj, resname)
    if specifics:
        specific_index = select_specific_atoms(traj, specifics)
        hbonds_in_traj = find_ligand_hbonds(traj, lig, specific_index, distance=distance, angle=angle, pseudo=pseudo)
    else:
        hbonds_in_traj = find_ligand_hbonds(traj, lig, specifics=None, distance=distance, angle=angle,
                                            pseudo=pseudo)
    read_and_rewrite_reports(report, trajectory=traj, hbonds_in_traj=hbonds_in_traj,
                             specific_hbonds=specific, own_path=own_path, report_in_pref=report_in_pref,
                             report_out_pref=report_out_pref, pseudo=pseudo)


def main(path_to_pdbs, resname="LIG", specifics=None, distance=0.25, angle=2.0 * np.pi / 3.0, pseudo=False,
         outpath=".", proc=4, top="../topologies/topology_0.pdb", rep_in="report_", rep_out="report_hb_"):
    print_info(path_to_pdbs, resname, specifics, distance, angle, pseudo, outpath, proc, top, rep_in, rep_out)
    pdbs = sorted(glob.glob(path_to_pdbs))
    reports_paths_list = ["/".join(os.path.abspath(pdb).split("/")[:-1]) for pdb in pdbs]
    reports_list = []
    for pdb in pdbs:
        pdb_id = pdb.split(".pdb")[0].split("trajectory_")[-1]
        for report_path in reports_paths_list:
             reports_list.append(os.path.join(report_path, "{}{}".format(rep_in, pdb_id)))
    reports_list = sorted(list(set(reports_list)))
    own_path = os.path.join(outpath)
    if not os.path.exists(own_path):
        os.mkdir(own_path)
    p = mp.Pool(int(proc))
    multi = []
    for pdb, report in zip(pdbs, reports_list):
        multi.append(p.apply_async(find_hbonds_in_pdb, [pdb, resname, report, specifics, distance, angle,
                                                        pseudo, own_path, top, rep_in, rep_out]))
    for process in multi:
         process.get()
    p.close()
    p.join()


if __name__ == '__main__':
    if os.path.splitext(sys.argv[1])[-1] == ".yaml":
        print("YAML input\n=================")
        arguments = parseargs_yaml()
        arguments = YamlParser(arguments.input_file)
        args = arguments.path, arguments.ligname, arguments.specific, arguments.distance, arguments.angle, \
               arguments.pseudo, arguments.outpath, arguments.cpus, arguments.top, arguments.report_in, \
               arguments.report_out
        path, ligname, specific, distance, angle, pseudo, out, cpus, top, rep_in, rep_out = args
    else:
        print("Command-line input")
        path, ligname, specific, distance, angle, pseudo, out, cpus, top, rep_in, rep_out = parse_arguments()
    main(path, ligname, specific, distance, angle, pseudo, out, cpus, top, rep_in, rep_out)



