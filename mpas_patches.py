import os
import sys
import time
import pickle as pkle

import numpy as np
import matplotlib.collections as mplcollections
import matplotlib.patches as patches
import matplotlib.path as path

from netCDF4 import Dataset

''' This module creates or retrives a collection of MPL Path Patches for an
MPAS unstructured mesh.

Given an MPAS mesh file, `get_mpas_patches` will create a Path Patch for each
MPAS grid, by looping over a Cell's vertices. Because this operation is a nCell
* nEdge operation, it will take some quite some time.

However, once a patch collection is created it is saved (using Python's Pickle
module) as a 'patch' file. This patch file can be loaded for furture plots on
that mesh, which will speed up future plots creation.

This module was created with much help and guidence from the following
repository:

* https://github.com/lmadaus/mpas_python

'''
def update_progress(job_title, progress):
    length = 40
    block = int(round(length*progress))
    msg = "\r{0}: [{1}] {2}%".format(job_title, "#"*block + "-"*(length-block),
    round(progress*100, 2))
    if progress >= 1: msg += " DONE\r\n"
    sys.stdout.write(msg)
    sys.stdout.flush()

def generate_mesh_patch_fname(mesh):
    ''' Generate the default mesh patch filename. eg: x1.10242.patches. '''
    nCells = len(mesh.dimensions['nCells'])
    mesh_patches_fname = str(nCells)+'.'+'patches'
    return mesh_patches_fname

def get_mpas_patches(mesh, pickle=True, pickleFile=None, **kwargs):

    force = kwargs.get('force', False)

    nCells = mesh.dimensions['nCells'].size
    maxEdges = mesh.dimensions['maxEdges'].size
    verticesOnCell = mesh.variables['verticesOnCell'][:]
    nEdgesOnCell = mesh.variables['nEdgesOnCell'][:]
    latVertex = mesh.variables['latVertex'][:]
    lonVertex = mesh.variables['lonVertex'][:]

    if pickleFile:
        pickle_fname = pickleFile
    else:
        pickle_fname = generate_mesh_patch_fname(mesh)

    print("Creating a patch file: ", pickle_fname)

    # Try to load the pickle file
    if(os.path.isfile(pickle_fname) and not force):
        pickled_patches = open(pickle_fname,'rb')
        try:
            patch_collection = pkle.load(pickled_patches)
            pickled_patches.close()
            print("Pickle file (", pickle_fname, ") loaded succsfully")
            return patch_collection
        except:
            print("ERROR: Error while trying to read the pickled patches")
            print("ERROR: The pickle file may be corrupted and was probably not created")
            print("ERROR: succesfully")
            print("ERROR: Please delete it and try again")
            sys.exit(-1)

    if force:
        print("\nForce flag set to True, overwriting pickle file if it exists...")
    else:
        print("\nNo pickle file found, creating patches...")

    print("If this is a large mesh, then this proccess will take a while...")

    mesh_patches = [None] * nCells

    # make a patch for each vertex
    vertices = verticesOnCell[:,:]
    latVertex = latVertex[vertices - 1]
    lonVertex = lonVertex[vertices - 1]
    del(vertices)
    del(verticesOnCell)

    # Normalize latitude and longitude
    lonVertex *= (180.0/np.pi)
    latVertex *= (180.0/np.pi)


    llVertexOnCell = np.stack((lonVertex, latVertex))
    del(lonVertex)
    del(latVertex)

    llVertexOnCell = np.moveaxis(llVertexOnCell, 0, 1)
    llVertexOnCell = np.moveaxis(llVertexOnCell, 1, 2)


    for i in range(llVertexOnCell.shape[0]):
        diff = np.subtract(llVertexOnCell[i,:,0], llVertexOnCell[i,0,0])
        llVertexOnCell[i,diff > 180.0,0] -= 360
        llVertexOnCell[i,diff < -180.0,0] += 360

        cell_patch = path.Path(llVertexOnCell[i,0:nEdgesOnCell[i]+1,:],
                                 closed=True,
                                 readonly=True)
        mesh_patches[i] = patches.PathPatch(cell_patch)
        if i % 250 == 0:
            update_progress("Creating Patch file: "+pickle_fname, i/nCells)

    print("\n")

    # Create patch collection
    patch_collection = mplcollections.PatchCollection(mesh_patches)

    # Pickle the patch collection
    if pickle:
        pickle_file = open(pickle_fname, 'wb')
        pkle.dump(patch_collection, pickle_file)
        pickle_file.close()
        print("\nCreated a patch file for mesh: ", pickle_fname)
    else:
        print("\nPatch Collection created, but was not saved")
        print("Enable Pickling by calling this function with pickle=True")

    return patch_collection


if __name__ == "__main__":
    import argparse

    description = """ Create a MatPltLib Patch Collection of an MPAS-A Mesh and
    save it as a Pickle File.

    Multi-threaded currently not available.
    """

    parser = argparse.ArgumentParser(description=description)
    
    parser.add_argument('mesh',
                        help='Path to NetCDF MPAS Mesh',
                        type=str) 
    parser.add_argument('-o', '--output',
                        help='Location and name of the output file',
                        type=str,
                        default=None) 
    parser.add_argument('-n', '--nThreads',
                        help='Number of Threads',
                        type=int,
                        default=1) 
    parser.add_argument('-f', '--force',
                        help='Overwrite existing patch file',
                        action='store_true',
                        default=False) 

    args = parser.parse_args()

    force = args.force
    mesh = args.mesh 
    output = args.output
    nThreads = args.nThreads

    if not os.path.isfile(mesh):
        print("ERROR: ", mesh)
        print("ERROR: Mesh file was not found")
        sys.exit(-1)

    # Check to see that the file is a dataset
    try:
        mesh = Dataset(mesh)
    except:
        print("ERROR: ", mesh)
        print("ERROR: Was not a valid NetCDF file")

    # See if the specified output file already exists...
    if not output:
        output = generate_mesh_patch_fname(mesh)

    # If the patch file exists and force is not specified, exit
    if os.path.isfile(output) and not force:
        print("ERROR: '", output, "' patch file already exists")
        sys.exit(-1)

    get_mpas_patches(mesh, pickleFile=output, force=force)
