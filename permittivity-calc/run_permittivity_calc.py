# -*- coding: utf-8 -*-
"""Contains useful functions and examples that use AirlineData"""

#File input
from helper_functions import _get_file, get_METAS_data, perm_compare
# Make relative path
import os
# Dataclass
from sparam_data import AirlineData

# Get data folder path for example files
DATAPATH = os.path.abspath('..') + '/data/'     

        
def run_default(airline_name='VAL',**kwargs):
    """
    Run AirlineData on get_METAS_data with all the prompts and return the \
        instance.
    """
    return AirlineData(*get_METAS_data(airline=airline_name),**kwargs)

def run_example():
    rexolite_example = AirlineData(*get_METAS_data(airline='PAL',\
        file_path=DATAPATH + 'rexolite_PAL.txt'),name='Rexolite')
    serpentine_example = AirlineData(*get_METAS_data(airline='VAL',\
        file_path=DATAPATH + 'serpentine_dry.txt'),bulk_density=1.6,\
        name='Serpentine - 25$^\\circ$C',normalize_density=True,norm_eqn='LI')
    return rexolite_example, serpentine_example

def multiple_meas(file_path=None,airline=None):
    """
    Generate an instance of AirlineData for every file in a directory. Store \
        the intances in a list, and plot them all using perm_compare.
        
    Arguments
    ---------
    file_path (str): Full path of any file in the source directory. \
        (Optional - will produce file dialog box if not provided.)
    
    airlne (str): Name of airline used. Every measurement must have been made \
        in the same airline. (Optional - will prompt if not provided.)
        
    Return
    ------
    class_list (lst): List of generated class instances of AirlineData.
    """
    # Use _get_file to get the filepath and airline name if not provided
    if file_path:   # If file path provided as argument
        file = file_path
    elif not file_path:   # If file path not provided
        print('\n')
        print("Select any data file in the source folder. All .txt "+\
              "files in the source folder must be METAS data tables.")
        # Get the file path and the airline name
        airline, file, L_in = _get_file(airline,file_path)
    elif not airline:   # If file path is given but airline is not
        airline, file, L_in = _get_file(airline,file_path)
        
    # Get directory path    
    directory = os.path.dirname(file)
    # Use a list to maintain order for plotting
    class_list = []
    # Iterate through all .txt files in the directory and run AirlineData
    for file in os.listdir(directory):
        if file.endswith(".txt"):
            filename = os.path.splitext(file)[0]    # Use file name as plot label
            # Append each new instance to class list
            class_list.append(AirlineData(*get_METAS_data(airline,\
                                os.path.join(directory, file)),name=filename))
    
    # Plot all files        
    perm_compare(class_list)
    
    return class_list 
    
#%% MAIN
def main():
    ## Run Examples
    global rexolite_example
    global serpentine_example
    rexolite_example, serpentine_example = run_example()
    
if __name__ == '__main__':
    main()  # Comment to supress example
#    pass    # Uncomment to supress example