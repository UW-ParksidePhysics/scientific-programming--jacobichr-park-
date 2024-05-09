"""
This module handles reading, analyzing, and plotting data related to energy and volume calculations
from density functional theory simulations.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
from convert_units import convert_units
from read_two_columns_text import read_two_columns_text
from calculate_bivariate_statistics import calculate_bivariate_statistics
from calculate_quadratic_fit import calculate_quadratic_fit
from equations_of_state import fit_eos
from annotate_plot import annotate_plot
from generate_matrix import generate_matrix
from calculate_lowest_eigenvectors import calculate_lowest_eigenvectors


def parse_file_name(filepath):
    """
    Extracts metadata from a formatted filename.

    Args:
        filepath (str): The file path to parse.

    Returns:
        tuple: Containing three strings (chemical_symbol, crystal_symmetry_symbol,
               density_functional_acronym).
    """
    filename = os.path.basename(filepath)
    parts = filename.split('.')
    if len(parts) < 4:
        raise ValueError("Filename does not contain enough parts to extract the necessary information.")

    return parts[0], parts[1], parts[2]


def plot_data_and_fit(data, eos_fit_curve, eos_parameters, chemical_symbol, crystal_symmetry_symbol, display_graph):
    """
    Plots data points and fit curve, annotates the plot, and displays or saves the figure.

    Args:
        data (tuple): Tuple of numpy arrays (volumes, energies).
        eos_fit_curve (numpy array): Energy values from EOS fit.
        eos_parameters (numpy array): Parameters from EOS fit, contains equilibrium_volume etc.
        chemical_symbol (str): Chemical symbol.
        crystal_symmetry_symbol (str): Crystal symmetry symbol.
        display_graph (bool): If True, display the graph; otherwise, save the graph to a file.
    """
    volumes, energies = data

    fit_volumes = np.linspace(min(volumes), max(volumes), num=len(eos_fit_curve))
    volumes_angstrom = np.array([convert_units(volume, 'cubic_bohr_per_atom', 'cubic_angstroms_per_atom') for volume in volumes])
    fit_volumes_angstrom = np.array([convert_units(volume, 'cubic_bohr_per_atom', 'cubic_angstroms_per_atom') for volume in fit_volumes])
    energies_eV = np.array([convert_units(energy, 'rydberg_per_atom', 'electron_volts_per_atom') for energy in energies])
    eos_fit_curve_eV = np.array([convert_units(energy, 'rydberg_per_atom', 'electron_volts_per_atom') for energy in eos_fit_curve])

    fig, ax = plt.subplots()
    ax.scatter(volumes_angstrom, energies_eV, color='blue', marker='o', label='Data Points')
    ax.plot(fit_volumes_angstrom, eos_fit_curve_eV, 'k-', label='Fit Curve')

    x_buffer = 0.1 * (np.max(fit_volumes_angstrom) - np.min(fit_volumes_angstrom))
    y_buffer = 0.1 * (np.max(energies_eV) - np.min(energies_eV))
    ax.set_xlim([np.min(fit_volumes_angstrom) - x_buffer, np.max(fit_volumes_angstrom) + x_buffer])
    ax.set_ylim([np.min(energies_eV) - y_buffer, np.max(energies_eV) + y_buffer])

    ax.set_xlabel(r'$V$ ($\mathit{\AA}^3$/atom)', fontsize=12)
    ax.set_ylabel(r'$E$ (eV/atom)', fontsize=12)
    ax.set_title(f"Murnaghan Equation of State for {chemical_symbol} in DFT {crystal_symmetry_symbol}")

    annotations = {
        chemical_symbol: {'position': np.array([0.05, 0.95]), 'alignment': ['left', 'top'], 'fontsize': 12},
        crystal_symmetry_symbol: {'position': np.array([0.5, 0.9]), 'alignment': ['center', 'top'], 'fontsize': 12, 'style': 'italic'},
        f'K0 = {eos_parameters[1]:.1f} GPa': {'position': np.array([0.5, 0.85]), 'alignment': ['center', 'top'], 'fontsize': 12},
        f'Created by Jacob Christensen, {datetime.now().date().isoformat()}': {'position': np.array([0.05, 0.01]), 'alignment': ['left', 'bottom'], 'fontsize': 10}
    }
    annotate_plot(ax, annotations)

    ax.legend()
    ax.grid(True)

    if display_graph:
        plt.show()
    else:
        plt.savefig('Christensen.Al.Fm-3m.GGA-PBE.MurnaghanEquationOfState.png')


def fit_an_equation_of_state(filename, display_graph=True):
    """
    Processes the file and plots results.

    Args:
        filename (str): Path to the data file.
        display_graph (bool): Controls whether to display or save the graph.
    """
    chemical_symbol, crystal_symmetry_symbol, density_functional_acronym = parse_file_name(filename)
    print(f"Processed file for: {chemical_symbol} in {crystal_symmetry_symbol} structure using {density_functional_acronym} approximation.")

    try:
        data = read_two_columns_text(filename)
    except OSError as e:
        print(e)
        return

    volumes = data[0]
    energies = data[1]
    statistics = calculate_bivariate_statistics(data)
    print(f"Statistics:\nMean of Y: {statistics[0]}\nStandard Deviation of Y: {statistics[1]}\nMin X: {statistics[2]}, Max X: {statistics[3]}\nMin Y: {statistics[4]}, Max Y: {statistics[5]}")

    quadratic_coefficients = calculate_quadratic_fit(data)
    print(f"Quadratic Coefficients: {quadratic_coefficients}")

    eos_fit_curve, eos_parameters = fit_eos(volumes, energies, quadratic_coefficients, eos='murnaghan', number_of_points=120)
    print(f"EOS Fit Curve: {eos_fit_curve}")
    print(f"EOS Parameters: {eos_parameters}")

    plot_data_and_fit(data, eos_fit_curve, eos_parameters, chemical_symbol, crystal_symmetry_symbol, display_graph)


def visualize_vectors_in_space(minimum_x, maximum_x, number_of_dimensions, potential_name, potential_parameter, display_graph=True):
    """
    Generates spatial grid and matrix, calculates eigenvectors and eigenvalues, and plots the results.

    Args:
        minimum_x (float): Minimum x value for spatial grid.
        maximum_x (float): Maximum x value for spatial grid.
        number_of_dimensions (int): Number of spatial grid points.
        potential_name (str): Name of the potential used.
        potential_parameter (float): Parameter for the potential.
        display_graph (bool): If True, displays the graph; otherwise, saves it.
    """
    spatial_grid = np.linspace(minimum_x, maximum_x, num=number_of_dimensions)
    matrix = generate_matrix(minimum_x, maximum_x, number_of_dimensions, potential_name, potential_parameter)

    eigenvalues, eigenvectors = calculate_lowest_eigenvectors(matrix, number_of_eigenvectors=3)

    if np.any(eigenvectors[0] < 0):
        eigenvectors[0] *= -1

    fig, ax = plt.subplots()
    for i in range(len(eigenvectors)):
        ax.plot(spatial_grid, eigenvectors[i], label=f"$\psi_{i+1}, E_{i+1} = {eigenvalues[i]:.4f}$ a.u.")

    ax.set_xlabel("x [a.u.]")
    ax.set_ylabel("$\psi_n (x) [a.u.]")
    title_text = f"Select Wavefunctions for a {potential_name} Potential on a Spatial Grid of {number_of_dimensions} Points"
    ax.set_title(title_text)

    max_eigenvalue_component = np.max(np.abs(eigenvectors))
    ax.set_ylim(-2 * max_eigenvalue_component, 2 * max_eigenvalue_component)

    ax.axhline(0, color='black', linewidth=1)

    signature_text = f"Created by Jacob Christensen, {datetime.now().date().isoformat()}"
    annotations = {
        signature_text: {'position': [0.01, 0.01], 'alignment': ['left', 'bottom'], 'fontsize': 10}
    }
    annotate_plot(ax, annotations)

    ax.legend()
    ax.grid(True)

    if display_graph:
        plt.show()
    else:
        plt.savefig('Christensen.harmonic.Eigenvector(1_2_3).png')

if __name__ == "__main__":
    # Visualize Vectors in Space
    visualize_vectors_in_space(minimum_x=-10, maximum_x=10, number_of_dimensions=120, potential_name='harmonic', potential_parameter=1.0, display_graph=True)

    # Fit an Equation of State
    current_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(current_dir, 'Al.Fm-3m.GGA-PBE.volumes_energies.dat')
    fit_an_equation_of_state(file_path, display_graph=True)