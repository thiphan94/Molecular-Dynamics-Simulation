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
        for index, line in enumerate(lines):
            if line == number:
                atom1 = lines[index + 1]
                atom2 = lines[index + int(number)]
                list_atom1 = atom1.split()
                print(type(list_atom1[3]))
                list_atom2 = atom2.split()
                x = int(list_atom1[4])
                print(x)
                # print(list_atom1[3])
                # distance = sqrt(
                #     pow((list_atom1[3] - list_atom2[3]), 2)
                #     + pow((list_atom1[4] - list_atom2[4]), 2)
                #     + pow((list_atom1[5] - list_atom2[5]), 2)
                # )
                #
                # print(distance)


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
            # line_3 = linecache.getline(file_input, 3)
            # print(line_3.split()[3])

    except Exception as e:
        logging.info(f"Error while opening file: {e}")
        # return 0


file_input = "data.gro"
Molecular_Dynamics_Simulation(file_input)
