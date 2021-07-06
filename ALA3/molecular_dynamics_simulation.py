"""Calculation of atomic distances."""

import math
import logging
import linecache
import re

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


def count_distance(file, number):
    with open(file) as f:
        lines = f.readlines()
        # for index, line in enumerate(lines):
        for index, line in enumerate(lines[0:50]):
            if line == number:
                atom1 = lines[index + 1]
                atom2 = lines[index + int(number)]
                list_atom1 = atom1.split()
                list_atom2 = atom2.split()
                list_atom1 = list(map(float, list_atom1[3:6]))
                list_atom2 = list(map(float, list_atom2[3:6]))
                distance = (
                    ((list_atom2[0] - list_atom1[0]) ** 2)
                    + ((list_atom2[1] - list_atom1[1]) ** 2)
                    + ((list_atom2[2] - list_atom1[2]) ** 2)
                ) ** 0.5

                # print(distance)


"""Function to calculate dihedral angles."""


def dihedral_angle(file, number):
    with open(file) as f:
        last = int(number) + 2
        lines = f.readlines()
        # list of atoms to calculate angle phi
        list_angle_phi = []
        # list of atoms to calculate angle psi
        list_angle_psi = []
        # index of four atoms in chain phi
        phi = []
        # index of four atoms in chain psi
        psi = []
        for index, line in enumerate(lines[0:last]):
            if len(lines[index].split()) > 4 and lines[index].split()[1] != "by":
                if lines[index].split()[1] == "N":
                    phi.append(lines[index].split()[2])
                    if len(psi) == 1:
                        psi.append(lines[index].split()[2])
                if lines[index].split()[1] == "CA":
                    phi.append(lines[index].split()[2])
                    if len(psi) == 2:
                        psi.append(lines[index].split()[2])
                if lines[index].split()[1] == "C":
                    phi.append(lines[index].split()[2])
                    if len(psi) == 0 or len(psi) == 3:
                        psi.append(lines[index].split()[2])
                # print("psi:", psi)
                if len(phi) == 4:
                    list_angle_phi.append(phi[:])
                    del phi[:3]
                if len(psi) == 4:
                    list_angle_psi.append(psi[:])
                    del psi[:3]

        print(list_angle_phi)
        print(list_angle_psi)

        coordinates_atom(file, number, list_angle_phi, list_angle_psi)


def coordinates_atom(file, number, list_phi, list_psi):
    with open(file) as f:
        last = int(number) + 2
        lines = f.readlines()
        # list of coordinates_atoms to calculate angle phi
        list_coordinates_phi = []
        # list of coordinates_atoms to calculate angle psi
        list_coordinates_psi = []
        # coordinates of atom in chain phi
        phi = []
        # coordinates of atom in chain psi
        psi = []
        i = 0
        #     print("ok")
        for i in range(len(list_phi)):
            for index, line in enumerate(lines[0:last]):
                # print(line)
                if len(lines[index].split()) > 4 and lines[index].split()[1] != "by":

                    if lines[index].split()[2] == list_phi[i][0]:
                        phi.append(lines[index].split()[3:6])
                        print(phi)
                    # psi.append(lines[index].split()[3:6])

                    if (
                        lines[index].split()[2]
                        == list_phi[i][1]
                        # or lines[index].split()[2] == list[1][1]
                    ):
                        phi.append(lines[index].split()[3:6])
                        # psi.append(lines[index].split()[3:6])

                    if (
                        lines[index].split()[2]
                        == list_phi[i][2]
                        # or lines[index].split()[2] == list[1][2]
                    ):
                        phi.append(lines[index].split()[3:6])
                        # psi.append(lines[index].split()[3:6])

                    if (
                        lines[index].split()[2]
                        == list_phi[i][3]
                        # or lines[index].split()[2] == list[1][3]
                    ):
                        phi.append(lines[index].split()[3:6])
                        # psi.append(lines[index].split()[3:6])

                    if len(phi) == 4:
                        list_coordinates_phi.append(phi[:])
                        del phi[:]
                    if len(psi) == 4:
                        list_coordinates_psi.append(psi[:])
                        del psi[:3]
        print(list_coordinates_phi)
        print(list_coordinates_psi)


"""Main function."""


def Molecular_Dynamics_Simulation(file_input):
    # file_opened = open_file(file_input)
    try:
        with open(file_input) as f:
            number_frames = count_frame(file_input)
            print("Number of frames:", number_frames)
            number_atoms = linecache.getline(file_input, 2)
            print("Number of atoms:", number_atoms)
            count_distance(file_input, number_atoms)

            dihedral_angle(file_input, number_atoms)

    except Exception as e:
        logging.info(f"Error while opening file: {e}")


file_input = "data.gro"
Molecular_Dynamics_Simulation(file_input)
