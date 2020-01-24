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
                        help="Angle threshold (in rad).")
    parser.add_argument("-d", "--distance", type=float, default=0.25,
                        help="Distance threshold (in nm).")
    parser.add_argument("-p", "--pseudo", default=False, action="store_true",
                        help="If set it also counts pseudo-hydrogen bonds (C-H)")
    parser.add_argument("-o", "--outpath", default=None,
                        help="Output path")
    parser.add_argument("-cpu", "--cpus", default=4, type=int,
                        help="Number of processors to paralelize the process")
    parser.add_argument("-t", "--top", default=None, type=str,
                        help="Topology path. By default it uses the standard Adaptive location.")

    args = parser.parse_args()

    return args.path, args.ligname, args.specific, args.angle, args.distance, args.pseudo, args.outpath, args.cpus, args.top


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
#            common = [atom for atom in hbond if atom in ligand_indexes]
#            if len(common) > 0 and len(common) < 2:
                hbonds_lig.setdefault(n, []).append(hbond)
        else:
            if any(atom in ligand_indexes for atom in hbond) and any(atom in specifics for atom in hbond) and \
                   not all(atom in hbond for atom in ligand_indexes):
#            common = [atom for atom in hbond if atom in ligand_indexes]
#            specials = [atom for atom in hbond if atom in specifics]
#            if len(common) > 0:
#                for com in common:
#                    print("COMMON FOUND")
#                    print(snapshot.topology.atom(com))
#            if len(specials) > 0:
#                for spe in specials:
#                    print("SPECIAL ATOM FOUND")
#                    print(snapshot.topology.atom(spe))
#            if len(common) > 0 and len(specials) > 0 and len(common) < 3:
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


def find_hbonds_in_pdb(pdb, resname, path_to_report, specifics=None, distance=0.25, angle=2.0 * np.pi / 3.0, pseudo=False, 
                       own_path=None, top=None):
    path_pdb = os.path.abspath(pdb)
    if not top:
        top_standard = "../../topologies/topology_0.pdb"
        top = os.path.normpath(os.path.join(path_pdb, top_standard))
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
    report = path_to_report + "/report_{}".format(pdb_id)
    new_report = pd.read_csv(report, sep=sep, header=0, engine="python")
    for n in range(0, len(traj)):
        try:
            new_report.loc[n, a2e.COLUMN_NAME] = len(hbonds_in_traj[n])
            #print(n)
            #for hbond in hbonds_in_traj[n]:
                #print(traj.topology.atom(hbond[0]), traj.topology.atom(hbond[2]))
        except KeyError:
            new_report.loc[n, a2e.COLUMN_NAME] = 0
    if not own_path:
        report_out = outpath + "/report_{}".format(pdb_id)
    else:
        report_out = own_path + "/report_{}".format(pdb_id)
    new_report.to_csv(report_out, sep="\t", index=False)


def main(path_to_pdbs, resname="LIG", specifics=None, distance=0.25, angle=2.0 * np.pi / 3.0, pseudo=False,
         outpath=None, proc=4, top="../topologies/topology_0.pdb"):
    pdbs = sorted(glob.glob(path_to_pdbs))
    print(pdbs)
    path_to_report = "/".join(os.path.abspath(path_to_pdbs).split("/")[:-1])
    identifier = os.path.abspath(pdbs[0]).split("/")[-5]
    if not outpath:
        outpath = "/".join(os.path.abspath(path_to_pdbs).split("/")[:-1])
    else:
        own_path = os.path.join(outpath, identifier)
        if not os.path.exists(own_path):
            os.mkdir(own_path)
    p = mp.Pool(int(proc))
    multi = []
    for pdb in pdbs:
        multi.append(p.apply_async(find_hbonds_in_pdb, [pdb, resname, path_to_report, specifics, distance, angle, pseudo, own_path, top]))
    for process in multi:
         process.get()
    p.close()
    p.join()


if __name__ == '__main__':
    path, ligname, specific, angle, distance, pseudo, out, cpus, top = parse_arguments()
    main(path, ligname, specific, distance, angle, pseudo, out, cpus, top)



