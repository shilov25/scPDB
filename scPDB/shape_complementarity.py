from scPDB.pdb_parser import split_complex_into_binders
import sys
import os
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import seaborn as sns
import matplotlib.pyplot as plt
import pymol
from pymol import cmd
import trimesh

def plot_sc_histogram(run_name, binder1_chains, binder2_chains, sc_b1_to_b2, sc_b2_to_b1, sc, title, save_path):
    """
    Plot histogram of Sc values for each binder.
    """
    b1_scores = pd.DataFrame({"Binder": f"{run_name}_{binder1_chains}", "sc_score": sc_b1_to_b2})
    b2_scores = pd.DataFrame({"Binder": f"{run_name}_{binder2_chains}", "sc_score": sc_b2_to_b1})
    sc_scores = pd.concat([b1_scores, b2_scores], ignore_index=True)

    ax = sns.histplot(data=sc_scores,
                binrange=[-1,1],
                x='sc_score',
                stat='count', 
                bins=100, 
                element='step',
                kde=True,
                hue='Binder',
                palette=['goldenrod','darkslategray'],
                alpha=0.3,
                )
    
    #show mean and median lines
    ax.axvline(np.median(sc_b1_to_b2), color='goldenrod', linestyle='--', linewidth=1.5, 
               label=f'{run_name}_{binder1_chains} Median: {np.median(sc_b1_to_b2):.2f}')

    ax.axvline(np.median(sc_b2_to_b1), color='darkslategray', linestyle='--', linewidth=1.5,
               label=f'{run_name}_{binder2_chains} Median: {np.median(sc_b2_to_b1):.2f}')
    
    ax.axvline(np.median(sc), color='red', linestyle='solid', linewidth=2,
               label=f'Sc Score: {sc:.2f}')

    plt.tight_layout()
    plt.title(title)
    plt.legend()
    plt.xlabel('sc scores')
    plt.ylabel('count')
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Saved sc histogram for {run_name} to {save_path}.")
    plt.close()

def create_binder_surface(binder1_chains, binder2_chains, binder1_pdb, binder2_pdb, run_name, output_dir):
    """
    Create full binder mesh using pyMOL. Mesh provides a
    higher density of points to sample from for Sc calculations.
    Meshes are saved as .obj files in output directory.
    """
    binder1_filepath = os.path.join(output_dir, run_name, f'{run_name}_{binder1_chains}.obj')
    binder2_filepath = os.path.join(output_dir, run_name, f'{run_name}_{binder2_chains}.obj')

    pymol.finish_launching(['pymol', '-c', '-q'])

    cmd.load(binder1_pdb, 'binder1')
    cmd.load(binder2_pdb, 'binder2')
    cmd.remove('solvent')

    cmd.set('mesh_quality', 2)
    cmd.set('mesh_solvent', 1)
    cmd.set('solvent_radius', 1.4)
    cmd.set('mesh_skip', 0)

    cmd.hide('everything')
    cmd.delete('binder2')
    cmd.show('surface', 'binder1')
    cmd.save(binder1_filepath, 'binder1')

    cmd.load(binder2_pdb, 'binder2')
    cmd.remove('solvent')

    cmd.set('mesh_quality', 2)
    cmd.set('mesh_solvent', 1)
    cmd.set('solvent_radius', 1.4)
    cmd.set('mesh_skip', 0)

    cmd.hide('everything')
    cmd.delete('binder1')
    cmd.show('surface', 'binder2')
    cmd.save(binder2_filepath, 'binder2')

    cmd.reinitialize('everything')

    binder1_mesh = trimesh.load(binder1_filepath)
    binder2_mesh = trimesh.load(binder2_filepath)

    binder1_surface_area = binder1_mesh.area
    binder2_surface_area = binder2_mesh.area 

    return binder1_mesh, binder2_mesh, binder1_surface_area, binder2_surface_area, binder1_filepath, binder2_filepath

def create_binder_interfaces(binder1_mesh, binder2_mesh, interface_thresh):
    """
    Creates interface surfaces P_a and P_b for Binder 1 and 
    Binder 2, respectively, using cKDTree. Default threshold 
    for defining interface points is 2 Angstroms.
    """
    b1_tree = cKDTree(binder1_mesh.vertices)
    b2_tree = cKDTree(binder2_mesh.vertices)

    #find interface vertices for binder 1
    b1_neighbors = b2_tree.query_ball_point(binder1_mesh.vertices, r=interface_thresh)
    b1_vertices_mask = np.array([len(neighbors) > 0 for neighbors in b1_neighbors])

    #find interface faces for binder 1
    b1_faces = []
    for i, face in enumerate(binder1_mesh.faces):
        if np.any(b1_vertices_mask[face]):
            b1_faces.append(i)

    #compile interface mesh for binder 1
    b1_interface_mesh = binder1_mesh.submesh([b1_faces])[0]
    
    #find interface vertices for binder 2
    b2_neighbors = b1_tree.query_ball_point(binder2_mesh.vertices, r=interface_thresh)
    b2_vertices_mask = np.array([len(neighbors) > 0 for neighbors in b2_neighbors])

    #find interface faces for binder 2
    b2_faces = []
    for i, face in enumerate(binder2_mesh.faces):
        if np.any(b2_vertices_mask[face]):
            b2_faces.append(i)

    #compile interface mesh for binder 2
    b2_interface_mesh = binder2_mesh.submesh([b2_faces])[0]
    
    b2_interface_area = b2_interface_mesh.area
    b1_interface_area = b1_interface_mesh.area

    return b1_interface_mesh, b2_interface_mesh, b1_interface_area, b2_interface_area

def sample_surface_points(interface_mesh, interface_area, surface_point_density):
    """
    Samples # of surface points/Angstrom^2 over the interface area.
    """
    num_surface_points = round(interface_area * surface_point_density)
    sampled_coords, indices = trimesh.sample.sample_surface_even(interface_mesh, num_surface_points)

    return sampled_coords, indices

def find_nearest_neighbor(b1_sampled_coords, b2_sampled_coords):
    """
    For each sampled surface point on Binder 1, find it's nearest neighbor
    on Binder 2, store its index and distance. Repeat vice versa for Binder 2.
    We will use this for calculating the distance component of the Sc.
    """
    nearest_neighbors = {}

    #collect nearest neighbors for points on binder 2
    b1_point_tree = cKDTree(b1_sampled_coords)
    b2_neighbor_distances, b2_neighbor_indices = b1_point_tree.query(b2_sampled_coords, k=1)
    b2_neighbor_coords = b1_sampled_coords[b2_neighbor_indices]

    nearest_neighbors['b2_neighbor_indices'] = b2_neighbor_indices
    nearest_neighbors['b2_neighbor_distances'] = b2_neighbor_distances
    nearest_neighbors['b2_neighbor_coords'] = b2_neighbor_coords

    b2_point_tree = cKDTree(b2_sampled_coords)
    b1_neighbor_distances, b1_neighbor_indices = b2_point_tree.query(b1_sampled_coords, k=1)
    b1_neighbor_coords = b2_sampled_coords[b1_neighbor_indices]

    nearest_neighbors['b1_neighbor_indices'] = b1_neighbor_indices
    nearest_neighbors['b1_neighbor_distances'] = b1_neighbor_distances
    nearest_neighbors['b1_neighbor_coords'] = b1_neighbor_coords

    return nearest_neighbors

def calculate_surface_normal_vectors(interface_mesh, indices):
    """
    For each point, a normal vector is calculated. The normal 
    vector points OUT from the face at the sampled point.
    """
    binder_normals = []
    for i in indices:
        face_normal = interface_mesh.face_normals[i]
        binder_normals.append(face_normal)

    return np.array(binder_normals)

def calculate_neighbor_normal_vectors(nearest_neighbors, binder1_normals, binder2_normals):
    """
    Rather than recalculate normals, we can search
    Binder 1 and Binder 2 normals by index and flip their signs
    so that neighbor normals point INTO the surface.
    """
    #get normals for binder 1 neighbors
    b1_neighbor_indices = nearest_neighbors['b1_neighbor_indices']
    b1_neighbor_normals = (-1) * binder2_normals[b1_neighbor_indices]
        
    #get normals for binder 2 neighbors
    b2_neighbor_indices = nearest_neighbors['b2_neighbor_indices']
    b2_neighbor_normals = (-1) * binder1_normals[b2_neighbor_indices]

    return np.array(b1_neighbor_normals), np.array(b2_neighbor_normals)

def remove_obj_files(binder_filepath, save_obj_files):
    if not save_obj_files:
        os.remove(binder_filepath)

def remove_binder_files(binder_pdb, save_binder_files):
    if not save_binder_files:
        os.remove(binder_pdb)

def calculate_sc(binder1_chains, binder2_chains, binder1_pdb, binder2_pdb, run_name, output_dir, interface_thresh, surface_point_density, w_constant, save_obj_files, save_binder_files):
    """
    Calculates Lawrence & Coleman shape complementarity score between
    two points nearest each other: x_a (on surface P_a) and x'_a 
    (on surface P_b). w is a scalar of 4 Angstroms, controls rate 
    of decay of Sc with distance. The set of formulas are as follows:

    S(A->B)(x_a) = (n_a . n'_a)exp[-w(|x_a - x'_a|)^2]
    S(B->A)(x_b) = (n_b . n'_b)exp[-w(|x_b - x'_b|)^2]
    Sc = 0.5 * {S(A->B)(x_a)} + {S(B->A)(x_b)}

    where the curly braces {} denote the MEDIAN of the distribution.
    In terms of this program, binder1 is A and binder2 is B.
    """
    #create surface mesh for binder 1 and binder 2
    binder1_mesh, binder2_mesh, _, _, binder1_filepath, binder2_filepath = create_binder_surface(
        binder1_chains,
        binder2_chains,
        binder1_pdb, 
        binder2_pdb, 
        run_name,
        output_dir,
        )
    #get interfaces and their areas for binder 1 and binder 2
    b1_interface_mesh, b2_interface_mesh, b1_interface_area, b2_interface_area = create_binder_interfaces(
        binder1_mesh,
        binder2_mesh,
        interface_thresh,
        )
    #get point samples across interface surface for binder 1 and binder 2
    b1_sampled_coords, b1_indices = sample_surface_points(
        b1_interface_mesh, 
        b1_interface_area, 
        surface_point_density,
        )
    b2_sampled_coords, b2_indices = sample_surface_points(
        b2_interface_mesh,
        b2_interface_area,
        surface_point_density,
        )
    #find nearest neighbors on binder 2 for sampled points on binder 1 and vice versa
    nearest_neighbors = find_nearest_neighbor(
        b1_sampled_coords, 
        b2_sampled_coords,
        )
    #get normals for sampled points on binders 1 and 2
    binder1_normals = calculate_surface_normal_vectors(
        b1_interface_mesh,
        b1_indices,
        )
    binder2_normals = calculate_surface_normal_vectors(
        b2_interface_mesh,
        b2_indices,
        )
    #get normals for neighbors of sampled points on binders 1 and 2
    b1_neighbor_normals, b2_neighbor_normals = calculate_neighbor_normal_vectors(
        nearest_neighbors, 
        binder1_normals, 
        binder2_normals,
        )
    #remove files, unless flagged
    remove_obj_files(binder1_filepath, save_obj_files)
    remove_obj_files(binder2_filepath, save_obj_files)

    remove_binder_files(binder1_pdb, save_binder_files)
    remove_binder_files(binder2_pdb, save_binder_files)

    #Calculate S(b1->b2)(x_b1)
    b1_orientation_component = np.sum(binder1_normals * b1_neighbor_normals, axis=1) #(n_a . n'_a)
    b1_distance_component = np.exp((-w_constant)*(nearest_neighbors['b1_neighbor_distances']**2))
    sc_b1_to_b2 = b1_orientation_component * b1_distance_component

    #Calculate S(b2->b1)(x_b2)
    b2_orientation_component = np.sum(binder2_normals * b2_neighbor_normals, axis=1) #(n_b . n'_b)
    b2_distance_component = np.exp((-w_constant)*(nearest_neighbors['b2_neighbor_distances']**2))
    sc_b2_to_b1 = b2_orientation_component * b2_distance_component

    sc = 0.5 * (np.median(sc_b1_to_b2) + np.median(sc_b2_to_b1))

    plot_sc_histogram(run_name, binder1_chains, binder2_chains, sc_b1_to_b2, sc_b2_to_b1, sc, 
                      title=f"Sc Score Distribution for {run_name}", 
                      save_path=f"{output_dir}{run_name}/{run_name}_sc_histogram.svg")

    sc = 0.5 * (np.median(sc_b1_to_b2) + np.median(sc_b2_to_b1))

    return np.median(sc_b1_to_b2), np.median(sc_b2_to_b1), sc, b1_interface_area, b2_interface_area