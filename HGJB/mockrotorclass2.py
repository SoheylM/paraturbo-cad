import cadquery as cq
import numpy as np
import os
import pickle
from math import sin, cos, pi, tan
import timeit
from enum import Enum, auto
import cq_warehouse.extensions




class Rotor():
    def __init__(self):
        #open file
        self.file = open('Element_23_08_19.pickle', 'rb')

        # dump info to that file
        self.Element = pickle.load(self.file)

        #close file
        self.file.close()

        self.Laenge = self.Element['Laenge']
        #change list into numpy array
        self.Laenge = np.array(self.Laenge)
        #change units from m to mm to avoid hollow visual effect

        self.Laenge = 1000*self.Laenge 
        #multiplying Laenge streeeetches the part out, so don't
        self.DI1 = 1000*np.array(self.Element['DI1'])
        self.DI2 = 1000*np.array(self.Element['DI2'])
        self.DI3 = 1000*np.array(self.Element['DI3'])
        self.DA1 = 1000*np.array(self.Element['DA1'])
        self.DA2 = 1000*np.array(self.Element['DA2'])
        self.DA3 = 1000*np.array(self.Element['DA3'])
        #print('type DI1', type(DI1))

        self.sys_pos = self.Element['sys_pos']
        self.pos_hgjb1 = self.sys_pos['pos_hgjb1']
        self.pos_hgjb2 =self.sys_pos['pos_hgjb2']
        self.cwf = os.getcwd().replace("\\", "/")
        #Begin code adapted from Christophe's Matlab
        #Parameters 
        self.N_HG = 28 #number of grooves generally between 26 - 30
        self.alpha_HG = Element['parameters']['hgjb1']['alpha']#0.68 #given
        self.beta_HG = Element['parameters']['hgjb1']['beta'] #-135 #given 
        self.beta_HG = self.beta_HG*pi/180
        self.gamma_HG = Element['parameters']['hgjb1']['gamma'] #0.89 #given
        self.h_gr = Element['parameters']['hgjb1']['hg'] #16 #groove depth given in micrometers
        self.h_rr = 9 #clearance on radiu given in micrometers
        self.D = 16 #on drawing [mm]
        self.L = self.Laenge[self.pos_hgjb1]#28 #length of HGJB on drawing [mm]
        self.L_land=self.L-(self.gamma_HG*self.L) #Value for CAD
        self.L=self.L+0.8 #oversized length for safety generally between 0.6 - 1
        self.Spiral_step = pi*self.D*tan(self.beta_HG)
        self.Spiral_height = self.L/2
        self.a_HG = (pi*self.D*self.alpha_HG)/self.N_HG #mm
        self.a_HG_plus_b_HG = self.a_HG/self.alpha_HG #mm
        self.h_rr_tot = self.h_rr*2 #diametral clearance given in micrometers
        #end code adapted from Christophe's Matlab
        
        self.dist1 = 0
        for d in range(self.pos_hgjb1):
            self.dist1=self.dist1+self.Laenge[d]
            
        #find distance to first HGJB
        self.dist2 = 0
        for d in range(self.pos_hgjb2):
            self.dist2=self.dist2+self.Laenge[d]
            
        #find distance to center of first HGJB
        self.DistCenter1=self.dist1+self.L/2

        #find distance to center of second HGJB
        self.DistCenter2=self.dist2+self.L/2

        #define separation angle between grooves
        self.sepang=360/self.N_HG

        #length between parallelogram verticals
        self.LenBetwVert = self.L/2 - self.L_land/2

        #Angle between parallelogram diagonal and reference horizontal
        self.Betaprime = abs(self.beta_HG) - 90

        #gap between horizontal leaving from parallelogram lowest 
        #corner and parallelogram vertical

        self.gap = self.LenBetwVert*tan(self.Betaprime)

        #create removal cylinder1
        self.CylLen1 = self.Laenge[self.pos_hgjb1]
        self.CylRadOut1= self.DA3[self.pos_hgjb1]/2

        #create removal cylinder2
        self.CylLen2 = self.Laenge[self.pos_hgjb2]
        self.CylRadOut2= self.DA3[self.pos_hgjb2]/2
        self.rot = cq.importers.importStep(self.cwf + '/Rotor.stp')
        
    def create_HGJB(self): 
        
        self.rot = self.rot.rotate((0,0,0),(1,0,0),270)
        #find distance to first HGJB
        

        #first cylinder to be projected onto
        cylinder1 = cq.Solid.makeCylinder(
             self.CylRadOut1, self.CylLen1+2, pnt=cq.Vector(0, self.DistCenter1+self.L/2+1, 0), dir=cq.Vector(0, -1, 0)
        )
        #first removal cylinder
        removalcylinder1 = cq.Solid.makeCylinder(
              self.CylRadOut1*2, self.CylLen1, pnt=cq.Vector(0, self.DistCenter1+self.L/2, -self.DA3[self.pos_hgjb1]*0.8), dir=cq.Vector(0, -1, 0)
        )

        #second cylinder to be projected onto
        cylinder2 = cq.Solid.makeCylinder(
             self.CylRadOut2, self.CylLen2+2, pnt=cq.Vector(0, self.DistCenter2+self.L/2+1, 0), dir=cq.Vector(0, -1, 0)
        )
        #second removal cylinder
        removalcylinder2 = cq.Solid.makeCylinder(
              self.CylRadOut2*2, self.CylLen2, pnt=cq.Vector(0, self.DistCenter2+self.L/2, -self.DA3[self.pos_hgjb1]*0.8), dir=cq.Vector(0, -1, 0)
        )

        #direction of projection 
        projection_direction = cq.Vector(0, 0, 1)

        #global coordinates of origin of parallelogram1
        # yp1 = DistCenter1+L/2
        xp1 = self.LenBetwVert*tan(self.Betaprime)/2 
        yp2 = -self.LenBetwVert*tan(self.Betaprime)
        xp2 = -self.LenBetwVert


        parallelogram1a = (
              cq.Workplane("YX", origin=((self.gap+self.a_HG)/2, self.DistCenter1+self.L/2, -2*self.CylRadOut1))
              #when viewed / \, y to the right and x up, z into screen
              #origin at upper outside corner
              #points below in standard x and y coordinates
              .lineTo(-self.LenBetwVert,-self.gap) #upper inside corner 
              .lineTo(-self.LenBetwVert,-self.gap-self.a_HG) #lower inside corner
              .lineTo(0,-self.a_HG) #lower outside corner
              .close()
              .extrude(1)
              .faces("<Z")
              .val()
          )

        parallelogram2a = (
              cq.Workplane("YX", origin=((self.gap+self.a_HG)/2, self.DistCenter1-self.L/2, -2*self.CylRadOut1))
              #when viewed / \, y to the right and x up, z into screen
              #origin at upper outside corner
              #points below in standard x and y coordinates
              .lineTo(self.LenBetwVert,-self.gap) #upper inside corner 
              .lineTo(self.LenBetwVert,-self.gap-self.a_HG) #lower inside corner
              .lineTo(0,-self.a_HG) #lower outside corner
              .close()
              .extrude(1)
              .faces("<Z")
              .val()
          )

        parallelogram1b = (
              cq.Workplane("YX", origin=((self.gap+self.a_HG)/2, self.DistCenter2+self.L/2, -2*self.CylRadOut2))
              #when viewed / \, y to the right and x up, z into screen
              #origin at upper outside corner
              #points below in standard x and y coordinates
              .lineTo(-self.LenBetwVert,-self.gap) #upper inside corner 
              .lineTo(-self.LenBetwVert,-self.gap-self.a_HG) #lower inside corner
              .lineTo(0,-self.a_HG) #lower outside corner
              .close()
              .extrude(1)
              .faces("<Z")
              .val()
          )

        parallelogram2b = (
              cq.Workplane("YX", origin=((self.gap+self.a_HG)/2, self.DistCenter2-self.L/2, -2*self.CylRadOut2))
              #when viewed / \, y to the right and x up, z into screen
              #origin at upper outside corner
              #points below in standard x and y coordinates
              .lineTo(self.LenBetwVert,-self.gap) #upper inside corner 
              .lineTo(self.LenBetwVert,-self.gap-self.a_HG) #lower inside corner
              .lineTo(0,-self.a_HG) #lower outside corner
              .close()
              .extrude(1)
              .faces("<Z")
              .val()
          )

        #fist HGJB
        #project first parallelogram onto cylinder
        parallelogram1a_projected = parallelogram1a.projectToShape(cylinder1, projection_direction)

        #turn first parallelogram into 3D shape on cylinder surface
        parallelogram1a_solids = cq.Compound.makeCompound(
              [f.thicken(self.h_gr/1000) for f in parallelogram1a_projected]
          )
        parallelogram1a_solids = parallelogram1a_solids.cut(removalcylinder1)

        #project second parallelogram onto cylinder
        parallelogram2a_projected = parallelogram2a.projectToShape(cylinder1, projection_direction)

        #turn first parallelogram into 3D shape on cylinder surface
        parallelogram2a_solids = cq.Compound.makeCompound(
              [f.thicken(self.h_gr/1000) for f in parallelogram2a_projected]
          )
        parallelogram2a_solids = parallelogram2a_solids.cut(removalcylinder1)


        #second HGJB
        #project first parallelogram onto cylinder
        parallelogram1b_projected = parallelogram1b.projectToShape(cylinder2, projection_direction)

        #turn first parallelogram into 3D shape on cylinder surface
        parallelogram1b_solids = cq.Compound.makeCompound(
              [f.thicken(self.h_gr/1000) for f in parallelogram1b_projected]
          )
        parallelogram1b_solids = parallelogram1b_solids.cut(removalcylinder2)

        #project second parallelogram onto cylinder
        parallelogram2b_projected = parallelogram2b.projectToShape(cylinder2, projection_direction)

        #turn first parallelogram into 3D shape on cylinder surface
        parallelogram2b_solids = cq.Compound.makeCompound(
              [f.thicken(self.h_gr/1000) for f in parallelogram2b_projected]
          )
        parallelogram2b_solids = parallelogram2b_solids.cut(removalcylinder2)

        # for i in range(N_HG):
        #     cylinder = cylinder.cut(parallelogram1_solids)
        #     #cylinder = cylinder.cut(parallelogram2_solids)
        #     cylinder = cylinder.transformed(rotate=(0,sepang,0))

        for i in range(self.N_HG):
            self.rot = self.rot.cut(parallelogram1a_solids)
            self.rot = self.rot.cut(parallelogram2a_solids)
            self.rot = self.rot.cut(parallelogram1b_solids)
            self.rot = self.rot.cut(parallelogram2b_solids)
            self.rot = self.rot.rotate((0,0,0),(0,1,0),self.sepang)

        self.rot = self.rot.rotate((0,0,0),(1,0,0),-270)
        
        
        
        return self.rot
        
        


t = Rotor().create_HGJB()


show_object(t)
