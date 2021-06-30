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

                print(distance)


"""Function to calculate dihedral angles."""


def dihedral_angle(file):
    with open(file) as f:
        lines = f.readlines()
        # for index, line in enumerate(lines):
        for index, line in enumerate(lines[0:20]):
            if len(line.split()) > 2 and line.split()[1] != "by":
                print(line.split()[1])
                if line.split()[1] in ["C", "N"]:
                    list_alpha.append(line.split()[1])
        # print(lines[2].split()[1])
        # #     print("ok")


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
            dihedral_angle(file_input)

    except Exception as e:
        logging.info(f"Error while opening file: {e}")


file_input = "data.gro"
Molecular_Dynamics_Simulation(file_input)
