import pickle
import numpy as np
import cadquery as cq

#open file
file = open('best_solution_Element.pickle', 'rb')

# dump info to that file
Element = pickle.load(file)

#close file
file.close()

#print('pos_hgjb1',Element.pos_hgjb1)

Laenge = Element['Laenge']
#change list into numpy array
Laenge = np.array(Laenge)
#change units from m to mm to avoid hollow visual effect

Laenge = 1000*Laenge 
#multiplying Laenge streeeetches the part out, so don't
DI1 = 1000*np.array(Element['DI1'])
DI2 = 1000*np.array(Element['DI2'])
DI3 = 1000*np.array(Element['DI3'])
DA1 = 1000*np.array(Element['DA1'])
DA2 = 1000*np.array(Element['DA2'])
DA3 = 1000*np.array(Element['DA3'])
#print('type DI1', type(DI1))

sys_pos = Element['sys_pos']
pos_hgjb1 = sys_pos['pos_hgjb1']
pos_hgjb2 =sys_pos['pos_hgjb2']
print(pos_hgjb1)
print(pos_hgjb2)

#NOTE: circle is defined by radius
#NOTE: hole is defined by diameter

WP = cq.Workplane('XY')

#part=WP.cylinder(10,2, centered =(True,True,False)) centers the cylinder on x and y axis but not on z axis

#initialize part
#start outside in: draw outermost cylinder, then make a hole
#this avoids drawing a cylinder, drawing another cylinder around it and putting a hole through both

#test presence of outermost cylinder, if found, draw and test presence of middle cylinder and inner cylinder, draw if present
if DA3[0] != DI3[0]:
    part3 = WP.cylinder(Laenge[0],DA3[0]/2, direct=(0,0,1), angle = 360, centered=(True,True,False),combine=False).hole(DI3[0],Laenge[0])
    show_object(part3)
    if DA2[0] != DI2[0]:
        part2=WP.cylinder(Laenge[0],DA2[0]/2, direct=(0,0,1), angle = 360,centered =(True,True,False),combine=False).hole(DI2[0],Laenge[0])    
        show_object(part2)
        if DA1[0] != DI1[0]:
                part1=WP.cylinder(Laenge[0],DA1[0]/2,direct=(0,0,1), angle = 360, centered =(True,True,False),combine=False).hole(DI1[0],Laenge[0])
                show_object(part1)
    elif DA1[0] != DI1[0]:
            part1=WP.cylinder(Laenge[0],DA1[0]/2, direct=(0,0,1), angle = 360,centered =(True,True,False),combine=False).hole(DI1[0],Laenge[0])
            show_object(part1)
#if no outermost cylinder, test presence of middle cyinder, draw if present, test for innermost cylinder and draw if present
elif DA2[0] != DI2[0]:
    part2=WP.cylinder(Laenge[0],DA2[0]/2, direct=(0,0,1), angle = 360,centered =(True,True,False),combine=False).hole(DI2[0],Laenge[0])
    show_object(part2)
    if DA1[0] != DI1[0]:
        part1=WP.cylinder(Laenge[0],DA1[0]/2,direct=(0,0,1), angle = 360, centered =(True,True,False),combine=False).hole(DI1[0],Laenge[0])
        show_object(part1)
#if no outermost and middle cylinder, test for inner cylinder and draw if present
elif DA1[0] != DI1[0]:
    part1 =WP.cylinder(Laenge[0],DA1[0]/2, direct=(0,0,1), angle = 360,centered =(True,True,False),combine=False).hole(DI1[0],Laenge[0])
    show_object(part1)
#if no cylinders present, print informational message
else:
    print("First section empty")
#part=WP.cylinder(Laenge[0],DA1[0]/2, centered=(True,True,False)).hole(DI1[0],Laenge[0])

#build successive cylinder stacks
#for outermost cylinder present, move to uppermost z face and add next cylinder
#check that the right face is being used at each stacking
for i in range(21):
    if DA3[i] != DI3[i]:
        WP = part3.faces('>Z')
    elif DA2[i] != DI2[i]:
        WP = part2.faces('>Z')
    elif DA1[i] != DI1[i]:
        WP = part1.faces('>Z')

    if DA3[i+1] != DI3[i+1]:
        part3=WP.cylinder(Laenge[i+1],DA3[i+1]/2, direct=(0,0,1), angle = 360,centered =(True,True,False),combine=False).faces('>Z').hole(DI3[i+1],Laenge[i+1])
        show_object(part3)
        if DA2[i+1] != DI2[i+1]:
            part2=WP.cylinder(Laenge[i+1],DA2[i+1]/2, direct=(0,0,1), angle = 360,centered =(True,True,False),combine=False).faces('>Z').hole(DI2[i+1],Laenge[i+1])
            show_object(part2)
        elif DA1[i+1] != DI1[i+1]:
            part1=WP.cylinder(Laenge[i+1],DA1[i+1]/2, direct=(0,0,1), angle = 360,centered =(True,True,False),combine=False).faces('>Z').hole(DI1[i+1],Laenge[i+1])
            show_object(part1)
    elif DA2[i+1] != DI2[i+1]:
        #print('DA2[i+1] != DI2[i+1]',i+1)
        part2=WP.cylinder(Laenge[i+1],DA2[i+1]/2, direct=(0,0,1), angle = 360,centered =(True,True,False),combine=False).faces('>Z').hole(DI2[i+1],Laenge[i+1])
        show_object(part2)
        if DA1[i+1] != DI1[i+1]:
            #print('DA1[i+1] != DI1[i+1]',i+1)
            part1=WP.cylinder(Laenge[i+1],DA1[i+1]/2, direct=(0,0,1), angle = 360,centered =(True,True,False),combine=False).faces('>Z').hole(DI1[i+1],Laenge[i+1])
            show_object(part1)
    elif DA1[i+1] != DI1[i+1]:
        part1=WP.cylinder(Laenge[i+1],DA1[i+1]/2, direct=(0,0,1), angle = 360,centered =(True,True,False),combine=False).faces('>Z').hole(DI1[i+1],Laenge[i+1])
        show_object(part1)

#try moving the workplane back to the front or making cylinder come back
#try cutting the parallelogram and then filling the hole left with another cylinder



