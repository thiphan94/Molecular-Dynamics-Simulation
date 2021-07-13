"""Calculation of atomic distances."""

import math
import logging
import linecache
import re
import numpy as np
import sys

# """FUnction to open the file."""
#
#
# def open_file(name_file):
#     try:
#         f = open(name_file, "r")
#         # lines = f.readlines()
#         # print(len(lines))
#         return f
#     except Exception as e:
#         logger.info(f"Error while opening file: {e}")
#         return 0


"""Function to calculate frames in file."""


def count_frame(file):
    with open(file) as f:
        count = 0
        for line in f:
            # line starts with "G"
            if line.startswith("G"):
                count += 1
        return count


"""Function to calculate the end-to-end distance."""


def count_distance(file, number, from_index, to_index):
    with open(file) as f:
        lines = f.readlines()
        for index, line in enumerate(lines[from_index:to_index]):
            if line == number:
                atom1 = lines[index + 1]
                # print(atom1)
                # print(index)
                atom2 = lines[index + int(number)]
                # print(atom2)
                list_atom1 = atom1.split()
                list_atom2 = atom2.split()
                list_atom1 = list(map(float, list_atom1[3:6]))
                list_atom2 = list(map(float, list_atom2[3:6]))
                distance = (
                    ((list_atom2[0] - list_atom1[0]) ** 2)
                    + ((list_atom2[1] - list_atom1[1]) ** 2)
                    + ((list_atom2[2] - list_atom1[2]) ** 2)
                ) ** 0.5

                print(distance)


"""Function to calculate dihedral angles."""


def dihedral_angle(file, number):
    with open(file) as f:
        last = int(number) + 2
        lines = f.readlines()
        # list of atoms to calculate angle phi
        list_angle_phi = []
        # list of atoms to calculate angle psi
        list_angle_psi = []
        # list of coordinates_atoms to calculate angle phi
        list_coordinates_phi = []
        # list of coordinates_atoms to calculate angle psi
        list_coordinates_psi = []
        # index of four atoms in chain phi
        phi = []
        # index of four atoms in chain psi
        psi = []
        phi_coordinates = []
        psi_coordinates = []
        vector_psi = []
        vector_phi = []
        for index, line in enumerate(lines[0:last]):
            if len(lines[index].split()) > 4 and lines[index].split()[1] != "by":
                if lines[index].split()[1] == "N":
                    psi.append(lines[index].split()[2])
                    psi_coordinates.append(lines[index].split()[3:6])
                    if len(phi) == 1:
                        phi.append(lines[index].split()[2])
                        phi_coordinates.append(lines[index].split()[3:6])
                if lines[index].split()[1] == "CA":
                    psi.append(lines[index].split()[2])
                    psi_coordinates.append(lines[index].split()[3:6])
                    if len(phi) == 2:
                        phi.append(lines[index].split()[2])
                        phi_coordinates.append(lines[index].split()[3:6])
                if lines[index].split()[1] == "C":
                    psi.append(lines[index].split()[2])
                    psi_coordinates.append(lines[index].split()[3:6])
                    if len(phi) == 0 or len(psi) == 3:
                        phi.append(lines[index].split()[2])
                        phi_coordinates.append(lines[index].split()[3:6])
                # print("psi:", psi)
                if len(phi) == 4:
                    list_angle_phi.append(phi[:])
                    del phi[:3]
                if len(psi) == 4:
                    list_angle_psi.append(psi[:])
                    del psi[:3]
                if len(phi_coordinates) == 4:
                    list_coordinates_phi.append(phi_coordinates[:])
                    del phi_coordinates[:3]
                if len(psi_coordinates) == 4:
                    list_coordinates_psi.append(psi_coordinates[:])
                    del psi_coordinates[:3]

        print(list_coordinates_psi)
        print(list_coordinates_phi)
        vector_psi, vector_phi = vector(list_coordinates_psi, list_coordinates_phi)

        arccos_angle(vector_psi, vector_phi)


"""Function to calculate vector between 4 atoms."""


def vector(list_coordinates_psi, list_coordinates_phi):
    ij_psi = []
    kj_psi = []
    kl_psi = []

    ij_phi = []
    kj_phi = []
    kl_phi = []
    vector_psi = []
    vector_phi = []

    for chain_psi in list_coordinates_psi:

        x_ij_psi = float(chain_psi[1][0]) - float(chain_psi[0][0])
        y_ij_psi = float(chain_psi[1][1]) - float(chain_psi[0][1])
        z_ij_psi = float(chain_psi[1][2]) - float(chain_psi[0][2])

        ij_psi.extend((x_ij_psi, y_ij_psi, z_ij_psi))

        x_kj_psi = float(chain_psi[2][0]) - float(chain_psi[1][0])
        y_kj_psi = float(chain_psi[2][1]) - float(chain_psi[1][1])
        z_kj_psi = float(chain_psi[2][2]) - float(chain_psi[1][2])

        kj_psi.extend((x_kj_psi, y_kj_psi, z_kj_psi))

        x_kl_psi = float(chain_psi[3][0]) - float(chain_psi[2][0])
        y_kl_psi = float(chain_psi[3][1]) - float(chain_psi[2][1])
        z_kl_psi = float(chain_psi[3][2]) - float(chain_psi[2][2])

        kl_psi.extend((x_kl_psi, y_kl_psi, z_kl_psi))

        vector_psi.append([ij_psi[:], kj_psi[:], kl_psi[:]])

        del ij_psi[:]
        del kj_psi[:]
        del kl_psi[:]

    for chain_phi in list_coordinates_phi:

        x_ij_phi = float(chain_phi[1][0]) - float(chain_phi[0][0])
        y_ij_phi = float(chain_phi[1][1]) - float(chain_phi[0][1])
        z_ij_phi = float(chain_phi[1][2]) - float(chain_phi[0][2])

        ij_phi.extend((x_ij_phi, y_ij_phi, z_ij_phi))

        x_kj_phi = float(chain_phi[2][0]) - float(chain_phi[1][0])
        y_kj_phi = float(chain_phi[2][1]) - float(chain_phi[1][1])
        z_kj_phi = float(chain_phi[2][2]) - float(chain_phi[1][2])

        kj_phi.extend((x_kj_phi, y_kj_phi, z_kj_phi))

        x_kl_phi = float(chain_phi[3][0]) - float(chain_phi[2][0])
        y_kl_phi = float(chain_phi[3][1]) - float(chain_phi[2][1])
        z_kl_phi = float(chain_phi[3][2]) - float(chain_phi[2][2])

        kl_phi.extend((x_kl_phi, y_kl_phi, z_kl_phi))

        vector_phi.append([ij_phi[:], kj_phi[:], kl_phi[:]])

        del ij_phi[:]
        del kj_phi[:]
        del kl_phi[:]

    # print(vector_psi)
    # print(vector_phi)
    return vector_psi, vector_phi


"""Function to calculate vector im and vector ln."""


def arccos_angle(vector_psi, vector_phi):

    for vector in vector_psi:
        vector_ij = np.array(vector[0])
        vector_kj = np.array(vector[1])
        vector_kl = np.array(vector[2])
        print(vector_ij)
        print(vector_kj)
        print(vector_kl)

        im = vector_ij - np.dot(
            (np.dot(vector_ij, vector_kj) / np.linalg.norm(vector_kj)), vector_kj
        )
        print("im =", im)

        ln = -vector_kl + np.dot(
            (np.dot(vector_kl, vector_kj) / np.linalg.norm(vector_kj)), vector_kj
        )

        print("ln =", ln)

        value_arccos = np.dot(im, ln) / np.dot(np.linalg.norm(im), np.linalg.norm(ln))
        print(value_arccos)
        # pass
        sign_angle = np.sign(np.dot(vector_ij, np.cross(vector_kj, vector_kl)))
        print(sign_angle)
        angle = math.degrees(np.arccos(value_arccos) * sign_angle)
        print(angle)


"""Main function."""


def Molecular_Dynamics_Simulation(file_input):
    # file_opened = open_file(file_input)
    try:
        with open(file_input) as f:
            number_frames = count_frame(file_input)
            print("Number of frames:", number_frames)
            number_atoms = linecache.getline(file_input, 2)
            print("Number of atoms:", number_atoms)
            from_index = 0
            to_index = int(number_atoms) + 2
            for i in range(0, int(number_frames)):
                print(from_index)
                print("to:  ", int(to_index))
                # to_index = int(number_atoms) + 3
                count_distance(file_input, number_atoms, from_index, to_index)

                from_index = from_index + 3 + int(number_atoms)
                to_index += int(number_atoms) + 3
            # dihedral_angle(file_input, number_atoms)

    except Exception as e:
        logging.info(f"Error while opening file: {e}")


# file_input = "data.gro"
file_input = str(sys.argv[1])
Molecular_Dynamics_Simulation(file_input)
