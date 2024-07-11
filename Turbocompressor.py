import cadquery as cq
import numpy as np
import pandas as pd
import time
import os
import sys

cwf = os.getcwd().replace("\\", "/")
sys.path.append(cwf  + '/../paraturbo-cad')
sys.path.append(cwf  + '/../paraturbo-cad/SGTB')
sys.path.append(cwf  + '/../paraturbo-cad/ROT')
sys.path.append(cwf  + '/../paraturbo-cad/COMP')
sys.path.append(cwf  + '/../paraturbo-cad/HGJB')

from helper_CAD import HELPER
from rotor_CAD import ROTOR
from sgtb_CAD import SGTB
#from impeller_old_but_correct_cut import IMPELLER
from Impeller_CAD import IMPELLER


# Addition 24/05/13
def main(case_folder=None):
    case_folder_default = False
    if case_folder is None:
        case_folder_default = True
        # Define a default folder or raise an error
        case_folder = os.getcwd().replace("\\", "/") + '/../paraturbo-cad/'
        print(f"No case folder specified. Using default: {case_folder}")

    # Ensure the default or specified directory exists
    if not os.path.exists(case_folder):
        os.makedirs(case_folder)
        print(f"Created directory: {case_folder}")
    # End addition 24/05/13
    
    t0 = time.time()
    # Example data for manual construction
    Length = [1.33373498, 1.33373498, 1.33373498, 1.33373498, 1.33373498, 1.33373498,  1.33373498, 1.33373498, 1.33373498, 2.26216371, 2.26216371, 12.22484027, 3.55328257, 4.05213368, 18.13412215, 39.15096023, 18.13412215, 14.33896901, 24.13178693]
    DI1 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0, 0.0, 0.0, 0.0, 7.37090243, 7.37090243, 7.37090243, 7.37090243, 7.37090243, 7.37090243, 0.0]
    DI2 = [4.57184857, 5.25196219, 5.60449431, 6.66769992, 8.61739416,11.98803844, 21.42243553, 16.09357927, 14.12109218, 13.57298226, 13.57298226, 12.22484027, 7.37090243, 7.37090243, 7.37090243, 7.37090243, 7.37090243, 7.37090243, 12.12462737]
    DI3 = [4.57184857, 10.38005986, 10.68632443, 11.62252172, 13.40640487, 16.95932086, 21.42243553, 16.09357927, 14.12109218, 13.57298226, 13.57298226, 13.57298226, 31.22649767, 13.57298226, 13.57298226, 13.57298226, 13.57298226, 13.57298226, 13.57298226]
    DO1 = [4.57184857, 5.25196219, 5.60449431, 6.66769992, 8.61739416, 11.98803844, 21.42243553, 16.09357927, 14.12109218, 13.57298226, 13.57298226, 12.22484027, 7.37090243, 7.37090243, 7.37090243, 7.37090243, 7.37090243, 7.37090243, 12.12462737]
    DO2 = [4.57184857, 10.38005986, 10.68632443, 11.62252172, 13.40640487, 16.95932086, 21.42243553, 16.09357927, 14.12109218, 13.57298226, 13.57298226, 13.57298226, 31.22649767, 13.57298226, 13.57298226, 13.57298226, 13.57298226, 13.57298226, 13.57298226]
    DO3 = [4.57184857, 10.38005986, 10.68632443, 11.62252172, 13.40640487, 16.95932086, 21.42243553, 16.09357927, 14.12109218, 13.57298226, 13.57298226, 13.57298226, 31.22649767, 13.57298226, 13.57298226, 13.57298226, 13.57298226, 13.57298226, 13.57298226]
    types1 = ['COMP1', 'COMP1', 'COMP1', 'COMP1', 'COMP1', 'COMP1', 'COMP1', 'COMP1', 'COMP1', 'COMP1', 'PLUG', 'PLUG', 'ROT', 'ROT', 'ROT', 'ROT', 'ROT', 'ROT', 'MAG']
    types2 = ['COMP1', 'COMP1', 'COMP1', 'COMP1', 'COMP1', 'COMP1', 'COMP1', 'COMP1', 'COMP1', 'COMP1', 'PLUG', 'ROT', 'ROT', 'ROT', 'ROT', 'ROT', 'ROT', 'ROT', 'ROT']
    types3 = ['COMP1', 'COMP1', 'COMP1', 'COMP1', 'COMP1', 'COMP1', 'COMP1', 'COMP1', 'COMP1', 'COMP1', 'PLUG', 'ROT', 'ROT', 'ROT', 'ROT', 'ROT', 'ROT', 'ROT', 'ROT']

    pos_SGTB = 12; alpha_SGTB = 0.52; beta_SGTB = -139.80; gamma_SGTB = 0.51; hg_SGTB = 0.0101; hr_SGTB = 0.0137; Ri_SGTB = 7.29; Rg_SGTB = 11.53; Ro_SGTB = 15.61; L_SGTB = 3.55
    pos_HGJB1 = 14; pos_HGJB2 = 16; alpha_HGJB = 0.68; beta_HGJB = -135.22; gamma_HGJB = 0.89; hg_HGJB = 0.016; hr_HGJB = 0.009

    Phi = 1.618; r_4 = 10; r_2h = r_4/(Phi**3); r_2s = r_4/Phi; R_rot = r_4/(Phi**2); b_4 = 0.5; b_6 = b_4/Phi; c = r_4/Phi; d = r_4-r_2h
    e = r_4/(Phi**2)-b_6; f = r_4 - R_rot; L_imp = r_2h+c+b_6+e; a = c-b_4; b = r_4-r_2s; delta_x = b_6/3; N_bld = 9; e_bld = 0.1
    R_rot = 5; L_ind = 40; beta_4 = -45; beta_2 = -56; beta_2s = -60; r_1 = 20; r_5 = 15; e_bld = 0.25; e_tip = 0.01; e_back = 0.01

    '''
    Use of Helper
    .importpickle:
        - input file name(string)
        - do not add the .pickle to the end of the file name

    .assemble:
        - input files(tuple), name of the output file(string)
        - combines all the files given in the tuple in a single assembly
    '''

    DesignTurbocompressor = HELPER(base_path=case_folder)

    if case_folder_default:
        Element = DesignTurbocompressor.importpickle('/../paraturbo-cad/ELEMENT/Element_to_CAD')
    else:
        # Path where the rotor's pickled element will be saved
        Element = DesignTurbocompressor.importpickle('/Element_to_CAD')

    '''
    SGTB Construction
    .parameters:
        - input Element(dictionary)

    .parameters.manual:
        - input Length(list), position(integer), alpha(float),
        beta(float), gamma(float), hg(float), hr(float), Ri(float),
        Rg(float), R0(float), L(float)
        - should not be used together with .parameters

    .grooves:
        - input number of grooves(integer)
        - by default, it is 28

    .CAD:
        - generates the CAD of SGTB for grooves towards the right
        - for color input 'color' in the end
        - for section view input 'section view' at the end
        - for dramatized groove depth for applications such as 3D printing input 'dramatize' in the end
        and d = (integer)

    .mirror:
        - mirrors the right SGTB to get the left SGTB

    .combined:
        - should be used together with .mirror
        - for saving as stl input 'stl' in the end

    .right:
        - for saving as stl input 'stl' in the end

    .left:
        - should be used together with .mirror
        - for saving as stl input 'stl' in the end
    '''
    start_time = time.time()
    DesignSGTB = SGTB(helper_instance=DesignTurbocompressor, base_path=case_folder)

    DesignSGTB.parameters(Element)
    # DesignSGTB.parameters_manual(Length,pos_SGTB,alpha_SGTB,beta_SGTB,gamma_SGTB,hg_SGTB,hr_SGTB,Ri_SGTB,Rg_SGTB,Ro_SGTB,L_SGTB)
    DesignSGTB.grooves(28)
    DesignSGTB.CAD('color')
    DesignSGTB.mirror()
    SGTBs = DesignSGTB.combined()
    SGTB_right = DesignSGTB.right()
    SGTB_left = DesignSGTB.left()
    end_time = time.time()
    time_SGTB = np.round(end_time - start_time,2)
    print("SGTB Execution Time: ", time_SGTB, "seconds")

    # show_object(SGTBs, name='SGTBs')
    # show_object(SGTB_right, name='SGTB Right')
    # show_object(SGTB_left, name='SGTB Left')

    '''
    Rotor Construction
    .parameters:
        - input Element(dictionary)

    .parameters.manual:
        - input Length(list), DI1(list), DI2(list), DI3(list),
        DO1(list), DO2(list), DO3(list)
        - should not be used together with .parameters
        - element types can be given in the end by specifying them as
        elem_type1 = (list), elem_type2 = (list), elem_type3=(list)

    .CAD:
        - generates the CAD of the rotor and returns it
        - for section view input 'section view' at the end

    .HGJB:
        - makes the parametrization for HGJB

    .HGJB_CAD:
        - input the CAD of the rotor from .CAD
        - generates grooves on the rotor

    .assemble:
        - assembles the rotor
        - for color input 'color' in the end
        - for saving as stl input 'stl' in the end
    '''

    start_time = time.time()
    DesignRotor = ROTOR('Joseph', helper_instance=DesignTurbocompressor,base_path=case_folder) # ROTOR('Joseph') ROTOR()

    DesignRotor.parameters(Element)
    # DesignRotor.parameters_manual(Length,DI1,DI2,DI3,DO1,DO2,DO3,pos_HGJB1,pos_HGJB2,alpha_HGJB,beta_HGJB,gamma_HGJB,hg_HGJB,hr_HGJB,elem_type1=types1,elem_type2=types2,elem_type3=types3)
    ROT = DesignRotor.CAD('color')
    DesignRotor.HGJB()
    DesignRotor.HGJB_CAD(ROT)
    print('grooves done.')
    print('ROT',ROT)
    Rot = DesignRotor.assemble()
    print('rotor assembled.')
    end_time = time.time()
    time_ROTOR_HGJB = np.round(end_time - start_time,2)
    print("ROTOR+HGJB Execution Time: ", time_ROTOR_HGJB, "seconds")

    # show_object(Rot, name='Rotor')

    '''
    Impeller Construction
    .parameters_impeller:
        - extracts the parameters from a pickle file
        - input Element(dictionary)

    .manualparams_impeller:
        - user defines parameters manually
        - input Element(dictionary), r_4(float), r_2s(float) ,beta_4(float) ,b_4(float) ,r_1(float), r_2h(float),
        r_5(float), e_bld(float), e_tip(float), e_back(float), L_ind(float), beta_2(float), beta_2s(float), N_bld(integer), R_rot(float)
        - for a custom rotor radius input 'manual_rotor' in the end
        - for using the rotor radius from the pickle file input 'auto_rotor' in the end
        - should not be used together with .parameters_impeller

    .hub:
        - models the hub
        - for section view input 'section view' in the end

    .blades_excel:
        - retrieves the coordinates of the blades from an excel file
        - input excel 'filename'(string)
        - should not be used together with .blades_coords

    .blades_coords:
        - extracts the coordinates of the blade curves from a pickle file
        - input Element(dictionary)
        - should not be used together with .blades_excel

    .model_blades:
        - models the blades
        - used after .blades_excel
        - input result of .blades_excel

    .rotate_blade:
        - patterns the blades
        - input result of .model_blades, 'bladename'(string)

    .assemble:
        - combines the impeller components in a common assembly and exports
        - input result of .hub, results of .rotate_blade
        - for exporting  as a stl file input 'stl' or 'STL' in the end
        - exports by default as step
    '''

    start_time = time.time()
    DesignImpeller = IMPELLER(helper_instance=DesignTurbocompressor,base_path=case_folder)
    DesignImpeller.parameters_impeller(Element)
    print('parameters_impeller.')

    # Optimizing the blades to match outlet angle under a reasonable sweep
    # This is the part where one ma want to introduce 2D code optimization,
    # blade generation rules, CFD based optimization etc.
    DesignImpeller.blade_opti()

    # To be used to only model the impeller independently
    # Imp.manualparams_impeller(Element,r_4,r_2s,beta_4,b_4,r_1,r_2h,r_5,e_bld,e_tip,e_back,L_ind,beta_2,beta_2s,N_bld,R_rot,'auto_rotor')

    Hub, Hub_solid = DesignImpeller.hub_v7()
    print('hub.')

    # Coords_mainblades = Imp.blades_excel('coordinates_blade_python.xlsx')
    # Coords_splitterblades = Imp.blades_excel('coordinates_splitter_python.xlsx')

    Coords_mainblades, Coords_splitterblades = DesignImpeller.blades_coords(Element)
    print('blades_coords.')

    Mainblade = DesignImpeller.model_blades_v8(Coords_mainblades)
    #Mainblade = DesignImpeller.trim_blades_old(Mainblade) # bugged
    Mainblades, Mainblades_solid_list = DesignImpeller.rotate_blade(Mainblade,'Main Blade')

    print('model_blades rotate_blade.')

    Splitterblade = DesignImpeller.model_blades_v8(Coords_splitterblades)
    #Splitterblade = DesignImpeller.trim_blades_old(Splitterblade) #bugged
    Splitterblades, Splitterblades_solid_list = DesignImpeller.rotate_blade(Splitterblade,'Splitter Blade')

    print('model_blades rotate_blade.')


    Compressor = DesignImpeller.assemble_v9((Hub_solid,Mainblades_solid_list,Splitterblades_solid_list))
    print('impeller assemble.')
    end_time = time.time()
    time_COMP = np.round(end_time - start_time,2)
    print("COMP Execution Time: ", time_COMP, "seconds")

    # show_object(Compressor, name = 'Compressor')

    start_time = time.time()
    DesignTurbocompressor.assemble((Rot,SGTBs,Compressor),'Turbocompressor', 'stl')
    print('turbocompressor assemble.')
    end_time = time.time()
    time_ASSEMBLY = np.round(end_time - start_time,2)
    print("Full ASSEMBLY Execution Time: ", time_ASSEMBLY, "seconds")

    time_TOTAL = np.round((time.time()-t0),2)
    print('Time: ' + str(time_TOTAL) + ' seconds')

     # Create a pandas Series for execution times
    times_data = pd.Series({
        "Time SGTB": time_SGTB,
        "Time ROTOR AND HGJB": time_ROTOR_HGJB,
        "Time COMP": time_COMP,
        "Time ASSEMBLY": time_ASSEMBLY,
        "Time TOTAL": time_TOTAL
    })
    # Path where to save the execution times CSV
    times_path = os.path.join(case_folder, "execution_times.csv")
    times_data.to_csv(times_path)

    print(f"Execution times saved to {times_path}")

# Addition 24/05/13
if __name__ == "__main__":
    # Check if the script was run with an argument, use it if present
    case_folder = sys.argv[1] if len(sys.argv) > 1 else None
    main(case_folder)
    # End addition 24/05/13