"""Calculation of atomic distances."""

import math


"""FUnction to open the file."""


def open_file(name_file):
    try:
        f = open(name_file, "r")
        return f
    except IOError:
        print("Error: File does not appear to exist !")
        return 0


"""Function to count frames in file."""


def count_frame(file):
    count = 0
    for line in file:
        # line starts with "G"
        if line.startswith("G"):
            count += 1
    return count


"""Main function."""


def Molecular_Dynamics_Simulation(file_input):
    file_opened = open_file(file_input)
    number_frames = count_frame(file_opened)
    print("Number of frames in file:", number_frames)


file_input = "data.gro"
Molecular_Dynamics_Simulation(file_input)
