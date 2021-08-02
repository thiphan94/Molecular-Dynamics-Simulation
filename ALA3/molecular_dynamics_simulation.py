"""Calculation of atomic distances."""

import math
import logging
import linecache
import re
import numpy as np
import sys
import pandas as pd


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
        lines = lines[from_index:to_index]
        for index, line in enumerate(lines):

            pattern = "=(.+?)step"
            if line == number:
                t = re.search("=(.+?)step", lines[index - 1]).group(1)

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

        return float(t), distance


def count_angle(file, number):
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
                    psi.append(lines[index].split()[2])

                    if len(phi) == 1:
                        phi.append(lines[index].split()[2])

                if lines[index].split()[1] == "CA":
                    psi.append(lines[index].split()[2])

                    if len(phi) == 2:
                        phi.append(lines[index].split()[2])

                if lines[index].split()[1] == "C":
                    psi.append(lines[index].split()[2])

                    if len(phi) == 0 or len(psi) == 3:
                        phi.append(lines[index].split()[2])

                if len(phi) == 4:
                    list_angle_phi.append(phi[:])
                    del phi[:3]
                if len(psi) == 4:
                    list_angle_psi.append(psi[:])
                    del psi[:3]

        return (list_angle_psi, list_angle_phi)


"""Function to calculate dihedral angles."""


def dihedral_angle(file, number, list_angle_psi, list_angle_phi, from_index, to_index):
    with open(file) as f:
        lines = f.readlines()
        lines = lines[from_index:to_index]

        list_coordinates_psi = []
        list_coordinates_phi = []
        phi_coordinates = []
        psi_coordinates = []
        vector_psi = []
        vector_phi = []

        for index, line in enumerate(lines):
            if len(lines[index].split()) > 4 and lines[index].split()[1] != "by":
                if any(
                    lines[index].split()[2] in sublist for sublist in list_angle_psi
                ):

                    psi_coordinates.append(lines[index].split()[3:6])

                if any(
                    lines[index].split()[2] in sublist for sublist in list_angle_phi
                ):

                    phi_coordinates.append(lines[index].split()[3:6])

                if len(phi_coordinates) == 4:
                    list_coordinates_phi.append(phi_coordinates[:])
                    del phi_coordinates[:3]

                if len(psi_coordinates) == 4:
                    list_coordinates_psi.append(psi_coordinates[:])
                    del psi_coordinates[:3]

        vector_psi, vector_phi = vector(list_coordinates_psi, list_coordinates_phi)

        angle_psi, angle_phi = value_angle(vector_psi, vector_phi)

        return angle_psi, angle_phi


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

    return vector_psi, vector_phi


"""Function to calculate vector im and vector ln."""


def value_angle(vector_psi, vector_phi):
    list_angle_psi = []
    list_angle_phi = []
    for vector in vector_psi:
        vector_ij_psi = np.array(vector[0])
        vector_kj_psi = np.array(vector[1])
        vector_kl_psi = np.array(vector[2])

        im_psi = vector_ij_psi - np.dot(
            (np.dot(vector_ij_psi, vector_kj_psi) / np.linalg.norm(vector_kj_psi)),
            vector_kj_psi,
        )

        ln_psi = -vector_kl_psi + np.dot(
            (np.dot(vector_kl_psi, vector_kj_psi) / np.linalg.norm(vector_kj_psi)),
            vector_kj_psi,
        )

        value_arccos_psi = np.dot(im_psi, ln_psi) / np.dot(
            np.linalg.norm(im_psi), np.linalg.norm(ln_psi)
        )

        sign_angle_psi = np.sign(
            np.dot(vector_ij_psi, np.cross(vector_kj_psi, vector_kl_psi))
        )

        angle_psi = math.degrees(np.arccos(value_arccos_psi) * sign_angle_psi)
        list_angle_psi.append(angle_psi)

    for vector in vector_phi:
        vector_ij_phi = np.array(vector[0])
        vector_kj_phi = np.array(vector[1])
        vector_kl_phi = np.array(vector[2])

        im_phi = vector_ij_phi - np.dot(
            (np.dot(vector_ij_phi, vector_kj_phi) / np.linalg.norm(vector_kj_phi)),
            vector_kj_phi,
        )

        ln_phi = -vector_kl_phi + np.dot(
            (np.dot(vector_kl_phi, vector_kj_phi) / np.linalg.norm(vector_kj_phi)),
            vector_kj_phi,
        )

        value_arccos_phi = np.dot(im_phi, ln_phi) / np.dot(
            np.linalg.norm(im_phi), np.linalg.norm(ln_phi)
        )

        sign_angle_phi = np.sign(
            np.dot(vector_ij_phi, np.cross(vector_kj_phi, vector_kl_phi))
        )

        angle_phi = math.degrees(np.arccos(value_arccos_phi) * sign_angle_phi)
        list_angle_phi.append(angle_phi)
    return list_angle_psi, list_angle_phi


"""Main function."""


def Molecular_Dynamics_Simulation(file_input, file_data):
    try:
        with open(file_input) as f:
            number_frames = count_frame(file_input)
            print("Number of frames:", number_frames)
            number_atoms = linecache.getline(file_input, 2)
            print("Number of atoms:", number_atoms)

            from_index = 0
            to_index = int(number_atoms) + 2
            # calculate number of psi angle and phi angle
            list_angle_psi, list_angle_phi = count_angle(file_input, number_atoms)
            number_angle_psi = len(list_angle_psi)
            number_angle_phi = len(list_angle_phi)

            # initialize size of matrix
            cols = 2 + number_angle_phi + number_angle_psi

            rows = number_frames

            values = [[0 for i in range(cols)] for j in range(rows)]
            index_psi = 1 + number_angle_psi

            columns = []

            for i in range(0, cols):
                if i == 0:
                    columns.append("t")
                elif i == 1:
                    columns.append("distance")
                elif 1 < i <= index_psi:
                    columns.append("psi" + str(i - 1))
                elif index_psi < i:
                    columns.append("phi" + str(i - number_angle_psi - 1))

            for i in range(0, number_frames):
                t, distance = count_distance(
                    file_input, number_atoms, from_index, to_index
                )
                # write t to matrix
                values[i][0] = t
                # write distance to matrix
                values[i][1] = distance

                angle_psi, angle_phi = dihedral_angle(
                    file_input,
                    number_atoms,
                    list_angle_psi,
                    list_angle_phi,
                    from_index,
                    to_index,
                )

                for j in range(len(values[0]) - 2):
                    if j < len(angle_psi):
                        values[i][j + 2] = angle_psi[j]
                    else:
                        values[i][j + 2] = angle_phi[j - 4]

                from_index = from_index + 3 + int(number_atoms)
                to_index += int(number_atoms) + 3

            df = pd.DataFrame(values, columns=columns)

            df.to_csv(file_data, index=False)

    except Exception as e:
        logging.info(f"Error while opening file gro: {e}")


file_input = str(sys.argv[1])
file_data = str(sys.argv[2])
Molecular_Dynamics_Simulation(file_input, file_data)
