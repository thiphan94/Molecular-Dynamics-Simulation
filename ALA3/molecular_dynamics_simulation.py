"""Calculation of atomic distances."""
import math
import logging
import re
import numpy as np
import sys
import pandas as pd
import time


def count_frame(lines):
    """Function to calculate frames in file."""
    return sum(1 for line in lines if line.startswith("G"))


def count_distance(lines, number, from_index, to_index):
    """Function to calculate the end-to-end distance."""
    sublines = lines[from_index:to_index]
    for index, line in enumerate(sublines):
        if line == number:
            t = re.search("=(.+?)step", sublines[index - 1]).group(1)

            atom1 = sublines[index + 1]
            atom2 = sublines[index + int(number)]
            list_atom1 = atom1.split()
            list_atom2 = atom2.split()
            list_atom1 = list(map(float, list_atom1[3:6]))
            list_atom2 = list(map(float, list_atom2[3:6]))
            distance = np.linalg.norm(np.array(list_atom2) - np.array(list_atom1))
    return float(t), distance


def count_angle(lines, number):
    last = int(number) + 2
    # list of atoms to calculate angle phi
    list_angle_phi = []
    # list of atoms to calculate angle psi
    list_angle_psi = []

    # index of four atoms in chain phi
    phi = []
    # index of four atoms in chain psi
    psi = []

    for line in lines[0:last]:
        if len(line.split()) > 4 and line.split()[1] != "by":
            if line.split()[1] == "N":
                psi.append(line.split()[2])

                if len(phi) == 1:
                    phi.append(line.split()[2])

            if line.split()[1] == "CA":
                psi.append(line.split()[2])

                if len(phi) == 2:
                    phi.append(line.split()[2])

            if line.split()[1] == "C":
                psi.append(line.split()[2])

                if len(phi) == 0 or len(psi) == 3:
                    phi.append(line.split()[2])

            if len(phi) == 4:
                list_angle_phi.append(phi[:])
                del phi[:3]
            if len(psi) == 4:
                list_angle_psi.append(psi[:])
                del psi[:3]

    return (list_angle_psi, list_angle_phi)


def dihedral_angle(lines, list_angle_psi, list_angle_phi, from_index, to_index):
    """Function to calculate dihedral angles."""
    sublines = lines[from_index:to_index]

    list_coordinates_psi = []
    list_coordinates_phi = []
    phi_coordinates = []
    psi_coordinates = []
    vector_psi = []
    vector_phi = []

    for line in sublines:
        if len(line.split()) > 4 and line.split()[1] != "by":
            if any(line.split()[2] in sublist for sublist in list_angle_psi):
                psi_coordinates.append(line.split()[3:6])

            if any(line.split()[2] in sublist for sublist in list_angle_phi):
                phi_coordinates.append(line.split()[3:6])

            if len(phi_coordinates) == 4:
                list_coordinates_phi.append(phi_coordinates[:])
                del phi_coordinates[:3]

            if len(psi_coordinates) == 4:
                list_coordinates_psi.append(psi_coordinates[:])
                del psi_coordinates[:3]

    vector_psi, vector_phi = vector(list_coordinates_psi, list_coordinates_phi)

    angle_psi, angle_phi = value_angle(vector_psi, vector_phi)

    return angle_psi, angle_phi


def vector(coordinates):
    coordinates_vector = np.array(coordinates, dtype=float)
    vectors = np.array(
        [
            coordinates_vector[1, :] - coordinates_vector[0, :],
            coordinates_vector[1, :] - coordinates_vector[2, :],
            coordinates_vector[3, :] - coordinates_vector[2, :],
        ]
    )
    return vectors


# def vector(list_coordinates_psi, list_coordinates_phi):
#     """Function to calculate vector between 4 atoms."""
#     print("Vector")
#     print(list_coordinates_psi)
#     print(list_coordinates_phi)
#     ij_psi = []
#     kj_psi = []
#     kl_psi = []

#     ij_phi = []
#     kj_phi = []
#     kl_phi = []
#     vector_psi = []
#     vector_phi = []

#     for chain_psi in list_coordinates_psi:
#         for i in range(3):
#             ij = float(chain_psi[1][i]) - float(chain_psi[0][i])
#             kj = float(chain_psi[2][i]) - float(chain_psi[1][i])
#             kl = float(chain_psi[3][i]) - float(chain_psi[2][i])
#             ij_psi.append(ij)
#             kj_psi.append(kj)
#             kl_psi.append(kl)

#         vector_psi.append([ij_psi[:], kj_psi[:], kl_psi[:]])

#         del ij_psi[:]
#         del kj_psi[:]
#         del kl_psi[:]

#     for chain_phi in list_coordinates_phi:
#         for i in range(3):
#             ij = float(chain_phi[1][i]) - float(chain_phi[0][i])
#             kj = float(chain_phi[2][i]) - float(chain_phi[1][i])
#             kl = float(chain_phi[3][i]) - float(chain_phi[2][i])
#             ij_phi.append(ij)
#             kj_phi.append(kj)
#             kl_phi.append(kl)

#         vector_phi.append([ij_phi[:], kj_phi[:], kl_phi[:]])

#         del ij_phi[:]
#         del kj_phi[:]
#         del kl_phi[:]
#     print(vector_psi, vector_phi)
#     print("End vector")
#     return vector_psi, vector_phi


def calculate_angle(vector_ij, vector_kj, vector_kl):
    im = vector_ij - np.dot(
        (np.dot(vector_ij, vector_kj) / np.linalg.norm(vector_kj) ** 2),
        vector_kj,
    )

    ln = -vector_kl + np.dot(
        (np.dot(vector_kl, vector_kj) / np.linalg.norm(vector_kj) ** 2),
        vector_kj,
    )

    value_arccos = np.dot(im, ln) / (np.linalg.norm(im) * np.linalg.norm(ln))

    sign_angle = np.sign(np.dot(vector_ij, np.cross(vector_kj, vector_kl)))

    angle = math.degrees(np.arccos(value_arccos) * sign_angle)
    return angle


def value_angle(vector_psi, vector_phi):
    """Function to calculate vector im and vector ln."""
    list_angle_psi = []
    list_angle_phi = []
    for vector in vector_psi:
        vector_ij_psi = np.array(vector[0])
        vector_kj_psi = np.array(vector[1])
        vector_kl_psi = np.array(vector[2])

        angle_psi = calculate_angle(vector_ij_psi, vector_kj_psi, vector_kl_psi)
        list_angle_psi.append(angle_psi)

    for vector in vector_phi:
        vector_ij_phi = np.array(vector[0])
        vector_kj_phi = np.array(vector[1])
        vector_kl_phi = np.array(vector[2])

        angle_phi = calculate_angle(vector_ij_phi, vector_kj_phi, vector_kl_phi)
        list_angle_phi.append(angle_phi)
    return list_angle_psi, list_angle_phi


def read_frames(
    lines, number_atoms, from_index, to_index, list_angle_psi, list_angle_phi
):
    t, distance = count_distance(lines, number_atoms, from_index, to_index)
    values = [t, distance]

    angle_psi, angle_phi = dihedral_angle(
        lines,
        list_angle_psi,
        list_angle_phi,
        from_index,
        to_index,
    )
    for psi, phi in zip(angle_psi, angle_phi):
        values = values + [psi, phi]

    return values


def Molecular_Dynamics_Simulation(file_input, file_data):
    """Main function."""
    try:
        with open(file_input) as f:
            lines = f.readlines()
    except Exception as e:
        logging.info(f"Error while opening file gro: {e}")

    number_frames = count_frame(lines)
    print("Number of frames:", number_frames)
    number_atoms = lines[1]
    print("Number of atoms:", number_atoms)
    from_indices = [(3 + int(number_atoms)) * i for i in range(number_frames)]
    to_indices = [
        int(number_atoms) + 2 + (3 + int(number_atoms)) * i
        for i in range(number_frames)
    ]
    from_indices = [from_indices[0]]
    to_indices = [to_indices[0]]
    # calculate number of psi angle and phi angle
    list_angle_psi, list_angle_phi = count_angle(lines, number_atoms)
    number_angle_psi = len(list_angle_psi)
    number_angle_phi = len(list_angle_phi)

    # initialize size of matrix
    cols = 2 + number_angle_phi + number_angle_psi

    index_psi = 1 + number_angle_psi

    columns = []
    index_psi = 0
    index_phi = 1

    for i in range(0, cols):
        if i == 0:
            columns.append("t")
        elif i == 1:
            columns.append("distance")
        elif i % 2 == 0 and i != 0:
            columns.append("psi" + str(index_psi))
            index_psi += 1
        else:
            columns.append("phi" + str(index_phi))
            index_phi += 1

    frames = [
        read_frames(
            lines, number_atoms, from_index, to_index, list_angle_psi, list_angle_phi
        )
        for from_index, to_index in zip(from_indices, to_indices)
    ]
    df = pd.DataFrame(np.array(frames), columns=columns)
    df.to_csv(file_data, index=False)


if __name__ == "__main__":
    start_time = time.time()
    file_input = str(sys.argv[1])
    file_data = str(sys.argv[2])
    Molecular_Dynamics_Simulation(file_input, file_data)
    print("--- %s seconds ---" % (time.time() - start_time))
