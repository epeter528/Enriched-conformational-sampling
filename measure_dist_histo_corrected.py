#! /usr/bin/python

#################################################################################
# Calculates data matrix from the given pdb file
#################################################################################

import string, sys, math, optparse, mmLib.PDB, glob
from mmLib.FileIO import *
from mmLib.AtomMath import *

data_matrix_file_name = "data_matrix.csv"

# which atoms take place in the torsion angle
distance_def = {
    "dist1":["C","N"],
    "dist2":["C","CA"],
    "dist3":["C","CB"],
    "dist4":["C","C"],
    "dist5":["C","N"],
    "dist6":["N","CA"],
    "dist7":["N","CB"],
    "dist8":["N","C"],
    "dist9":["N","N"],
    "dist10":["C","CA"],
    "dist11":["C","CB"],
     }

nrdistances = len(distance_def)

# character string characterizing not available data
not_available = "-"

# list of residues I consider to be standard
std_resid_list = ["ALA"]

# return torsion angle of distance_name using residue1 and residue2 (b,g,d are calculated from atoms
# residing only in the actual residue1, no residue2 is needed, therefore its implicit value is None)
def calc_backbone_torsion(distance_name, residue1, residue2):

    for distnr in range(1, nrdistances+1):
	if distance_name == "dist" + str(distnr):
	    atom1 = residue1.get_atom(distance_def[distance_name][0])
	    atom2 = residue2.get_atom(distance_def[distance_name][1])

    if atom1 == None: return not_available
    if atom2 == None: return not_available

    distance =  calc_distance(atom1, atom2) #*180/math.pi

    if distance > 10.0:
	return not_available
    else:
	return "%.2f" % distance

def get_torsions(residue_list):
    row = {}

    residue1 = residue_list[0]
    residue2 = residue_list[1]

    for distnr in range(1, nrdistances+1):
	row["dist" + str(distnr)] = calc_backbone_torsion("dist" + str(distnr), residue1, residue2)

    return row

def get_data_matrix(pdb_files):
    data_matrix = []
    for pdb_file in pdb_files:
        #print pdb_file
	struct = LoadStructure(fil=pdb_file, format="PDB")

	pdb_code = "XXXX"

        # go through all chains in the structure
        for chain in struct:
	    nres = chain.count_fragments()
	    if nres >= 4:
		for fragment in chain:
                    for fragment2 in chain:             
                      
                      if fragment == fragment2:
                          continue

		      data_matrix_row = {}
		      residue_list = []

		      residue1 = fragment
		      residue2 = fragment2

		      try:
#			if residue1.res_name != "GLY":
#			    continue
		        	if residue1.res_name != "ALA":
			            continue
              			if residue2.res_name != "ALA":
	               		    continue
#			if residue4.res_name != "GLY":
#			    continue
		      except Exception, e1:
		        continue

		      residue_list.append(residue1)
		      residue_list.append(residue2)

		      my_id_q = pdb_code+"_"+chain.chain_id+"_"+residue1.res_name+residue1.fragment_id+"_"+residue2.res_name+residue2.fragment_id

		      data_matrix_row["my_id_q"] = my_id_q

		      t = get_torsions(residue_list)
		      data_matrix_row.update(t)

		      data_matrix.append(data_matrix_row)

    return data_matrix

def write_data_matrix(data_matrix):
    for record in data_matrix:
        print record["my_id_q"] + "," + ",".join(record["dist"+str(distnr)] for distnr in range(1, nrdistances+1))

#!!!!!!!!!!!!!!!!!!!!! MAIN BODY !!!!!!!!!!!!!!!!!!!!!

pdb_files = glob.glob(sys.argv[1])

data_matrix = get_data_matrix(pdb_files)
write_data_matrix(data_matrix)
