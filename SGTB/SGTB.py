import cadquery as cq
import numpy as np
import os

class SGTB():
    def __init__(self):
        self.n_grooves = 28
        self.cwf = os.getcwd().replace("\\", "/")

    def parameters(self,Element):
        # Checks for information exist in dictionary
        # Checks for are all needed keys are in Element, is Element a dictionary
        # If not raises an error
        if type(Element) == dict:
            if all(x in Element.keys() for x in ['Laenge','sys_pos','parameters']) and \
                all(x in Element['parameters'].keys() for x in ['sgtb']) and \
                    all(x in Element['parameters']['sgtb'].keys() for x in ['alpha','beta', 'hg', 'hr', 'Ri', 'Rg', 'Ro', 'gamma', 'L']):
                    # Taking variables from Element
                    self.length = 1000*np.array(Element['Laenge'])
                    self.pos = np.int64(Element['sys_pos']['pos_sgtb'])
                    self.alpha = np.float64(Element['parameters']['sgtb']['alpha'])
                    self.beta = np.float64(Element['parameters']['sgtb']['beta'])
                    self.beta = (np.pi*self.beta)/180 
                    self.gamma = np.float64(Element['parameters']['sgtb']['gamma'])
                    self.hg = 1000*np.float64(Element['parameters']['sgtb']['hg'])
                    self.hr = 1000*np.float64(Element['parameters']['sgtb']['hr'])
                    self.Ri= 1000*np.float64(Element['parameters']['sgtb']['Ri'])
                    self.Rg = 1000*np.float64(Element['parameters']['sgtb']['Rg'])
                    self.Ro = 1000*np.float64(Element['parameters']['sgtb']['Ro'])
                    self.L = 1000*np.float64(Element['parameters']['sgtb']['L'])                 
            else:
                raise KeyError('SGTB.parameters: Element dictionary does not include all the needed keys.')
        else:
            raise TypeError('SGTB.parameters: Element type is not dictionary.')
        
    def parameters_manual(self,Length,pos_SGTB,alpha_SGTB,beta_SGTB,gamma_SGTB,hg_SGTB,hr_SGTB,Ri_SGTB,Rg_SGTB,Ro_SGTB,L_SGTB):
        # Checks for information exist in given input
        # Checks for the type for some inputs
        # If not raises an error
        if type(Length) == list:
                # Taking variables from user
                self.length = np.array(Length)
                self.pos = np.int64(pos_SGTB)
                self.alpha = np.float64(alpha_SGTB)
                self.beta = np.float64(beta_SGTB)
                self.beta = (np.pi*self.beta)/180
                self.gamma = np.float64(gamma_SGTB)
                self.hg = np.float64(hg_SGTB)
                self.hr = np.float64(hr_SGTB)
                self.Ri = np.float64(Ri_SGTB)
                self.Rg = np.float64(Rg_SGTB)
                self.Ro = np.float64(Ro_SGTB)
                self.L = np.float64(L_SGTB)
        else:
            raise TypeError('SGTB.parameters_manual: The type of the given variables are not suitable.')
        
    def grooves(self,n_grooves):
        self.n_grooves = n_grooves

    def CAD(self,*settings,**times):
        # Parameters for constructing the SGTB
        Ri  = self.Ri                                        # On drawing - SGTB inner diameter
        Ro  = self.Ro                                        # On drawing - Rotor axial stop
        n_SG = self.n_grooves                                # Number of grooves generally between 28 - 30
        Rg = self.Rg  
        a_SG = (2*np.pi*Rg*self.alpha)/n_SG 
        phi_lag_SG = a_SG/Rg

        # Building the grooves
        Ro = Ro + 0.3                                        # Oversized radius for safety between 0.3 - 0.5
        n = 10                                               # Number of points for each log spiral splines
        alpha = np.linspace(0, np.log(Ro/Rg)/np.tan(self.beta), num=n)
        
        # Generating the spiral
        r_SGTB = Rg * np.exp(np.tan(self.beta)*alpha)
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
        n_circle = 5                                        # Number of points for each arc
        angle1 = np.arctan(Groovy_up_y[-1]/Groovy_up_x[-1])
        angle2 = np.arctan(Groovy_down_y[0]/Groovy_down_x[0])
        circle1angles = np.linspace(angle1,angle2, num=n_circle)
        Upper_arc_x = Ro*np.cos(circle1angles)
        Upper_arc_y = Ro*np.sin(circle1angles)

        angle3 = np.arctan(Groovy_down_y[-1]/Groovy_down_x[-1])
        angle4 = np.arctan(Groovy_up_y[0]/Groovy_up_x[0])
        circle2angles = np.linspace(angle3,angle4, num=n_circle)
        Lower_arc_x = Rg*np.cos(circle2angles)
        Lower_arc_y = Rg*np.sin(circle2angles)

        # Making list of coordinates for Cadquery spline function
        firstcurve = []; secondcurve = []; firstcircle = []; secondcircle = []

        # Adding coordinates to a list
        for i in range(0,n):
            firstcurve.append((Groovy_up_x[i],Groovy_up_y[i]))
            secondcurve.append((Groovy_down_x[i],Groovy_down_y[i]))

        for i in range(0,n_circle):
            firstcircle.append((Upper_arc_x[i], Upper_arc_y[i]))
            secondcircle.append((Lower_arc_x[i], Lower_arc_y[i]))

        # CAD of the bearing
        # Creating the workplane
        wp = cq.Workplane('XY')

        b_di = Ri*2
        b_do = np.round(Ro*2*1.4,1)
        b_l = self.L

        # Building the SGTB without grooves
        thrust_right = wp.cylinder(b_l,b_do/2,
        direct=(0,0,1),angle=360,centered=(True,True,False),
        combine=False,clean=True).faces('>Z').hole(b_di,depth=b_l,clean=True)
        thrust_right = thrust_right.edges(cq.selectors.RadiusNthSelector(1)).chamfer(np.round(b_l/7,1))

        # Creating an empty curves dictionary and rotation angles
        curves = {}

        if 'dramatize' in settings:
            if 'd' in times:
                multiplier = times['d']
                self.hg = self.hg * multiplier
            else:
                raise ValueError('SGTB.CAD: Multiplier is not specified correctly.')
            
        # Rotating, drawing and opening the grooves
        for i in range(0,n_SG):
            wp = wp.transformed(rotate=(0,0,360/n_SG))
            curves['curve'+str(i)] = wp.spline(firstcurve).spline(firstcircle)\
                .spline(secondcurve).spline(secondcircle).close()
            locals().update(curves)
            thrust_right = thrust_right.spline(firstcurve).spline(firstcircle)\
                .spline(secondcurve).spline(secondcircle).close().cutBlind(self.hg,clean=True)

        # Checks for section view
        if 'section view' in settings:
            thrust_right = thrust_right.rect(80,40,(-40,0)).cutThruAll()

        self.thrust_right = thrust_right

        # Specifying color
        if 'color' in settings:
            self.color = True
        else:
            self.color = False

        # Positioning the bearings
        pos_right = self.hr
        for i in range(0,self.pos+1):
            pos_right += self.length[i]

        pos_left = -self.hr
        for i in range(0,self.pos):
            pos_left += self.length[i]

        self.pos_right = pos_right
        self.pos_left = pos_left
    
    def mirror(self):
        # Mirros the right SGTB to get the left SGTB
        self.thrust_left = self.thrust_right.mirror(mirrorPlane = 'XY')

    def combined(self,*settings):
        # Creating an assembly file and adding SGTBs into it with or without color
        assembly = cq.Assembly(name='SGTB')
        if self.color == True:
            assembly.add(self.thrust_right,loc = cq.Location((0,0,self.pos_right),(1,0,0),0),
                name='Right SGTB',color=cq.Color('magenta4'))
            assembly.add(self.thrust_left,loc = cq.Location((0,0,self.pos_left),(1,0,0),0),
                name='Left SGTB',color=cq.Color('magenta4'))
        elif self.color == False:
            assembly.add(self.thrust_right,loc = cq.Location((0,0,self.pos_right),(1,0,0),0),
                name='Right SGTB',color=cq.Color('gray50'))
            assembly.add(self.thrust_left,loc = cq.Location((0,0,self.pos_left),(1,0,0),0),
                name='Left SGTB',color=cq.Color('gray50'))
            
        # Saves as stl if given in arguments
        if 'stl' or 'STL' in settings:
            cq.exporters.export(self.thrust_right, self.cwf + '/STL/SGTB Right.stl')
            cq.exporters.export(self.thrust_left, self.cwf + '/STL/SGTB Left.stl')
        
        # Saves as step
        assembly.save(self.cwf  + '/STEP/SGTBs.step')

        return assembly

    def right(self,*settings):
        # Creating an assembly file and adding right SGTB into it with or without color
        assembly = cq.Assembly(name='SGTB')
        if self.color == True:
            assembly.add(self.thrust_right,
                name='Right SGTB',color=cq.Color('magenta4'))
        elif self.color == False:
            assembly.add(self.thrust_right,
                name='Right SGTB',color=cq.Color('gray50'))
        
        # Saves as stl if given in arguments
        if 'stl' or 'STL' in settings:
            cq.exporters.export(self.thrust_right, self.cwf + '/STL/SGTB Right.stl')
        
        # Saves as step
        assembly.save(self.cwf  + '/STEP/SGTB Right.step')

        return assembly
    
    def left(self,*settings):
        # Creating an assembly file and adding right SGTB into it with or without color
        assembly = cq.Assembly(name='SGTB')
        if self.color == True:
            assembly.add(self.thrust_left,
                name='Left SGTB',color=cq.Color('magenta4'))
        elif self.color == False:
            assembly.add(self.thrust_left,
                name='Left SGTB',color=cq.Color('gray50'))
        
        # Saves as stl if given in arguments
        if 'stl' or 'STL' in settings:
            cq.exporters.export(self.thrust_left, self.cwf + '/STL/SGTB Left.stl')
       
        # Saves as step
        assembly.save(self.cwf  + '/STEP/SGTB Left.step')

        return assembly