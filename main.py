import time
import os
import glob
import argparse
import multiprocessing as mp
from itertools import product

import numpy as np
import pandas as pd
import mdtraj as md
import hbond_mod as hm
import atoms_to_explore as a2e

__author__ = "Carles Perez Lopez"
__version__ = "1.0.0"
__maintainer__ = "Carles Perez Lopez"
__email__ = "carlesperez94@gmail.com"


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
    parser.add_argument("-a", "--angle", type=float, default=2.0 * np.pi / 3.0,
                        help="Angle threshold (in rad). Default: 2pi/3 (120º).")
    parser.add_argument("-d", "--distance", type=float, default=0.25,
                        help="H-X distance threshold (in nm). Default: 0.25nm.")
    parser.add_argument("-p", "--pseudo", default=False, action="store_true",
                        help="If set it also counts pseudo-hydrogen bonds (C-H)")
    parser.add_argument("-o", "--outpath", default="hbond_analysis",
                        help="Output path. By default: hbond_analysis/")
    parser.add_argument("-cpu", "--cpus", default=4, type=int,
                        help="Number of processors to paralelize the process")
    parser.add_argument("-t", "--top", default=None, type=str,
                        help="Topology path. By default it uses the standard Adaptive location.")
    parser.add_argument("-r", "--report", default="report_hb_", type=str,
                        help="Report's output prefix. Default: 'report_hb_'")

    args = parser.parse_args()

    return args.path, args.ligname, args.specific, args.angle, args.distance, args.pseudo, args.outpath, args.cpus, \
           args.top, args.report


def select_ligand(traj, resname="LIG"):
    lig = traj.topology.select("resname '{}'".format(resname))
    return lig


def select_specific_atoms(traj, specific_dict):
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


def get_column_name(atoms_to_explore_dict):
    """
    It joins all residue numbers and atom names to create a unique column name for each search.
    :param atoms_to_explore_dict: dictionary with keys as residue numbers and values as list of atom names
    :return: column name
    """
    column_name = []
    for resid, atomnames_list in atoms_to_explore_dict.items():
        column_name.append(resid)
        for atom in atomnames_list:
            column_name.append(atom)
    return "".join(column_name)


def find_hbonds_in_pdb(pdb, resname, path_to_in_report, specifics=None, distance=0.25, angle=2.0 * np.pi / 3.0, pseudo=False,
                       own_path=None, top=None, report_out_pref="report_hb_"):
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
    pdb_id = pdb.split(".pdb")[0].split("trajectory_")[-1]
    sep = "   "
    report = path_to_in_report + "/report_{}".format(pdb_id)
    new_report = pd.read_csv(report, sep=sep, header=0, engine="python")
    column_name = get_column_name(a2e.HBOND_ATOMS)
    for n in range(0, len(traj)):
        try:
            new_report.loc[n, column_name] = len(hbonds_in_traj[n])

        except KeyError:
            new_report.loc[n, column_name] = 0
    report_out = own_path + "/{}{}".format(report_out_pref, pdb_id)
    new_report.to_csv(report_out, sep="\t", index=False)


def main(path_to_pdbs, resname="LIG", specifics=None, distance=0.25, angle=2.0 * np.pi / 3.0, pseudo=False,
         outpath=".", proc=4, top="../topologies/topology_0.pdb", rep_out="report_hb_"):
    pdbs = sorted(glob.glob(path_to_pdbs))
    path_to_report = "/".join(os.path.abspath(path_to_pdbs).split("/")[:-1])
    own_path = os.path.join(outpath)
    if not os.path.exists(own_path):
        os.mkdir(own_path)
    p = mp.Pool(int(proc))
    multi = []
    for pdb in pdbs:
        multi.append(p.apply_async(find_hbonds_in_pdb, [pdb, resname, path_to_report, specifics, distance, angle,
                                                        pseudo, own_path, top, rep_out]))
    for process in multi:
         process.get()
    p.close()
    p.join()


if __name__ == '__main__':
    path, ligname, specific, angle, distance, pseudo, out, cpus, top, rep_out = parse_arguments()
    main(path, ligname, specific, distance, angle, pseudo, out, cpus, top, rep_out)



