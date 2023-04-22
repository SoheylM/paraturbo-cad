# SGTB Code

import cadquery as cq
import numpy as np
import os
import pickle

class SGTB():
    def __init__(self):
        self.n_grooves_SGTB = 28
        self.color = True
        self.sectionview = False
        self.cwf = os.path.dirname(os.path.abspath(__file__))

    def importpickle(self,filename):
        file = open(self.cwf  + '/' + filename, 'rb')
        Element = pickle.load(file)
        file.close

        return Element

    def parameters(self,Element):
        if type(Element) == dict:
        # Check for information exist in dictionary
            if all(x in Element.keys() for x in ['Laenge','DA3','sys_pos','alpha_SGTB','beta_SGTB','gamma_SGTB','hg_SGTB','hr_SGTB']):
                    if len(Element['Laenge']) == len(Element['DA3']):
                        # Taking variables from Element
                        self.length = 1000*np.array(Element['Laenge'])
                        self.DO3 = 1000*np.array(Element['DA3'])
                        self.pos_SGTB = np.int64(Element['sys_pos']['pos_sgtb'])
                        self.alpha_SGTB = np.float64(Element['alpha_SGTB'][0])
                        self.beta_SGTB = np.float64(Element['beta_SGTB'][0])
                        self.beta_SGTB = (np.pi*self.beta_SGTB)/180 
                        self.gamma_SGTB = np.float64(Element['gamma_SGTB'][0])
                        self.hg_SGTB = 1000*np.float64(Element['hg_SGTB'][0])
                        self.hr_SGTB = 1000*np.float64(Element['hr_SGTB'][0])        
                    else:
                        print('SGTB.parameters: Size of the needed dictionary values are not equal.')
                        return
            else:
                print('SGTB.parameters: Element dictionary does not include all the needed keys.')
                return
        else:
            print('SGTB.parameters: Element type is not dictionary.')
            return
        
    def parameters_manual(self,Length,DO3,pos_SGTB,alpha_SGTB,beta_SGTB,gamma_SGTB,hg_SGTB,hr_SGTB):
        if type(Length) == list and type(DO3) == list:
                if len(Length) == len(DO3):
                    # Taking variables from user
                    self.length = np.array(Length)
                    self.DO3 = np.array(DO3)
                    self.pos_SGTB = np.int64(pos_SGTB)
                    self.alpha_SGTB = np.float64(alpha_SGTB)
                    self.beta_SGTB = np.float64(beta_SGTB)
                    self.gamma_SGTB = np.float64(gamma_SGTB)
                    self.hg_SGTB = np.float64(hg_SGTB)
                    self.hr_SGTB = np.float64(hr_SGTB)
                else:
                    print('SGTB.parameters_manual: Size of the given list values are not equal.')
                    return
        else:
            print('SGTB.parameters_manual: The type of the given variables are not suitable.')
            return
        
    def settings(self,color,sectionview):
        self.color = color
        self.sectionview = sectionview

    def grooves(self,n_grooves):
        self.n_grooves_SGTB = n_grooves

    def CAD(self):
        # Parameters
        Ri  = np.round((self.DO3[self.pos_SGTB+1]/2),1)+0.5       # on drawing - SGTB inner diameter
        R0  = np.round((self.DO3[self.pos_SGTB]/2),1)             # on drawing - Rotor axial stop
        n_SG = self.n_grooves_SGTB                                # number of grooves generally between 28 - 30
        Rg = self.gamma_SGTB*(R0-Ri)+Ri  
        a_SG = (2*np.pi*Rg*self.alpha_SGTB)/n_SG 
        phi_lag_SG = a_SG/Rg

        # Building the grooves
        R0 = R0 + 0.3                                 # oversized radius for safety between 0.3 - 0.5
        n = 10                                        # nb of points for each log spiral splines
        alpha = np.linspace(-0.2, 1, num=n)
        alpha = np.linspace(0, np.log(R0/Rg)/np.tan(self.beta_SGTB), num=n)
        
        # Generating the spiral
        r_SGTB = Rg * np.exp(np.tan(self.beta_SGTB)*alpha)
        alpha_lag = alpha + phi_lag_SG

        # Building the data vector for CAD
        Groovy_up_x = r_SGTB*np.cos(alpha)
        Groovy_up_y = r_SGTB*np.sin(alpha)
        Groovy_down_x = r_SGTB*np.cos(alpha_lag)
        Groovy_down_y = r_SGTB*np.sin(alpha_lag)

        # Reversing the data for the lower spline
        Groovy_down_x = np.flip(Groovy_down_x)
        Groovy_down_y = np.flip(Groovy_down_y)

        # Creating arcs for the sketch
        n_circle = 5                                 # nb of points for each arc
        angle1 = np.arctan(Groovy_up_y[-1]/Groovy_up_x[-1])
        angle2 = np.arctan(Groovy_down_y[0]/Groovy_down_x[0])
        circle1angles = np.linspace(angle1,angle2, num=n_circle)
        Upper_arc_x = R0*np.cos(circle1angles)
        Upper_arc_y = R0*np.sin(circle1angles)

        angle3 = np.arctan(Groovy_down_y[-1]/Groovy_down_x[-1])
        angle4 = np.arctan(Groovy_up_y[0]/Groovy_up_x[0])
        circle2angles = np.linspace(angle3,angle4, num=n_circle)
        Lower_arc_x = Rg*np.cos(circle2angles)
        Lower_arc_y = Rg*np.sin(circle2angles)

        # Making list of coordinates for Cadquery spline function
        firstcurve = []; secondcurve = []; firstcircle = []; secondcircle = []

        for i in range(0,n):
            firstcurve.append((Groovy_up_x[i],Groovy_up_y[i]))
            secondcurve.append((Groovy_down_x[i],Groovy_down_y[i]))

        for i in range(0,n_circle):
            firstcircle.append((Upper_arc_x[i], Upper_arc_y[i]))
            secondcircle.append((Lower_arc_x[i], Lower_arc_y[i]))

        # CAD of the bearing
        wp = cq.Workplane('XY')
        b_di = Ri*2
        b_do = np.round(R0*2*1.4,1)
        b_l = np.round(b_do*0.09,1)

        thrust_right = wp.cylinder(b_l,b_do/2,
        direct=(0,0,1),angle=360,centered=(True,True,False),
        combine=False,clean=True).faces('>Z').hole(b_di,depth=b_l,clean=True)
        thrust_right = thrust_right.edges(cq.selectors.RadiusNthSelector(1)).chamfer(np.round(b_l/7,1))

        # Creating an empty curves dictionary and rotation angles
        curves = {}

        # Rotating, drawing and opening the grooves
        for i in range(0,n_SG):
            wp = wp.transformed(rotate=(0,0,360/n_SG))
            curves['curve'+str(i)] = wp.spline(firstcurve).spline(firstcircle)\
                .spline(secondcurve).spline(secondcircle).close()
            locals().update(curves)
            thrust_right = thrust_right.spline(firstcurve).spline(firstcircle)\
                .spline(secondcurve).spline(secondcircle).close().cutBlind(self.hg_SGTB,clean=True)

        if self.sectionview == True:
            thrust_right = thrust_right.rect(80,40,(-40,0)).cutThruAll()

        self.thrust_right = thrust_right

    def mirror(self):
        self.thrust_left = self.thrust_right.mirror(mirrorPlane = 'XY')

    def position(self):
        # Positioning the bearings
        pos_right = self.hr_SGTB
        for i in range(0,self.pos_SGTB+1):
            pos_right += self.length[i]

        pos_left = -self.hr_SGTB
        for i in range(0,self.pos_SGTB):
            pos_left += self.length[i]

        self.pos_right = pos_right
        self.pos_left = pos_left

    def combined(self):
        color = ('magenta4','gray50')
        assembly = cq.Assembly()
        if self.color == False:
            assembly.add(self.thrust_right,loc = cq.Location((0,0,self.pos_right),(1,0,0),0),
                name='rightSGTB',color=cq.Color(color[1]))
            assembly.add(self.thrust_left,loc = cq.Location((0,0,self.pos_left),(1,0,0),0),
                name='leftSGTB',color=cq.Color(color[1]))
        elif self.color == True:
            assembly.add(self.thrust_right,loc = cq.Location((0,0,self.pos_right),(1,0,0),0),
                name='rightSGTB',color=cq.Color(color[0]))
            assembly.add(self.thrust_left,loc = cq.Location((0,0,self.pos_left),(1,0,0),0),
                name='leftSGTB',color=cq.Color(color[0]))

        assembly.save(self.cwf  + '/SGTBs.step')

        return assembly

    def right(self):
        color = ('magenta4','gray50')
        assembly = cq.Assembly()
        if self.color == False:
            assembly.add(self.thrust_right,
                name='rightSGTB',color=cq.Color(color[1]))
        elif self.color == True:
            assembly.add(self.thrust_right,
                name='rightSGTB',color=cq.Color(color[0]))

        assembly.save(self.cwf  + '/SGTB Right.step')

        return assembly
    
    def left(self):
        color = ('magenta4','gray50')
        assembly = cq.Assembly()
        if self.color == False:
            assembly.add(self.thrust_left,
                name='leftSGTB',color=cq.Color(color[1]))
        elif self.color == True:
            assembly.add(self.thrust_left,
                name='leftSGTB',color=cq.Color(color[0]))

        assembly.save(self.cwf  + '/SGTB Left.step')

        return assembly