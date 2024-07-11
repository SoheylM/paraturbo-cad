import cadquery as cq
import numpy as np
import os
import pickle
import trimesh
import time

class HELPER():
    def __init__(self, base_path=None):
        if base_path is None:
            self.cwf = os.getcwd().replace("\\", "/") + '/../paraturbo-cad/'
        else:
            self.cwf = base_path.replace("\\", "/")  # Ensure it's always in the correct format
        #self.cwf = os.getcwd().replace("\\", "/") + '/../paraturbo-cad/'

    def importpickle(self,filename):
        # Importing pickle file and closing the file
        file = open(self.cwf  + filename + '.pickle', 'rb')
        Element = pickle.load(file)
        file.close

        return Element
    
    def assemble(self,files,filename, *settings):
        # Assembling the parts that are given and saving as step
        assembly = cq.Assembly(name=filename)
        for i in range(0,len(files)):
            assembly.add(files[i],name='Subassembly '+str(i+1))

        assembly.save(self.cwf  + '/STEP/' + filename + '.step')
        print('Pausing 18 seconds for writing '+ filename +' STEP.')
        time.sleep(18)
        #cq.exporters.export(assembly.toCompound(), self.cwf + '/STL/Turbocompressor.stl')
        if 'stl' or 'STL' in settings:
            self.convert_step_to_stl(self.cwf + '/STEP/Turbocompressor', self.cwf + '/STL/Turbocompressor')

    def convert_step_to_stl(self, input_filename, output_filename, gmsh_args=None):
        # Example usage:
        #convert_step_to_stl('Turbocompressor.STEP', 'Turbocompressor.STL')

        if gmsh_args is None:
            num_threads = os.cpu_count() or 1  # Get the number of CPUs, default to 1 if unable to detect
            gmsh_args = [
                ("Mesh.Algorithm", 1),  # Different algorithm types
                ("Mesh.CharacteristicLengthFromCurvature", 50),  # Tuning the smoothness
                ("General.NumThreads", num_threads),  # Multithreading capability
                ("Mesh.MinimumCirclePoints", 50)
            ]

        mesh = trimesh.Trimesh(**trimesh.interfaces.gmsh.load_gmsh(
            file_name = input_filename + '.step', 
            gmsh_args=gmsh_args
        )) 

        # Get and print some properties of the formed mesh (default output in mm)
        print("Mesh Volume: ", mesh.volume)
        print("Mesh Bounding Box Volume: ", mesh.bounding_box_oriented.volume)
        print("Mesh Area: ", mesh.area)

        # Export the new mesh in the STL format
        mesh.export(output_filename + '.stl')

