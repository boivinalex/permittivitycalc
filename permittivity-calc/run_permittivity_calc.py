# -*- coding: utf-8 -*-
"""Contains useful functions and examples that use AirlineData"""

#File input
from helper_functions import get_METAS_data, _get_file
# Plotting
import permittivity_plot as pplot
# Make relative path
import os
# Dataclass
from sparam_data import AirlineData

DATAPATH = os.path.dirname(__file__) + '/data/'


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

def perm_compare(classlist,allplots=False,**kwargs):
    """
    Given a list of AirlineData instances, plot their permittivity results \
        together using permittivity_plot_V1.py
        
    Arguments
    ---------
    classlist (list): List of instances of AirlineData
    
    allplots (bool): If True plot all of dielectric constant, loss factor and \
        loss tangent. If Flase plot only the dielectric constant and the loss \
        tangent. Default: False
    """
    freq = []
    dielec = []
    losstan = []
    labels = []
    for item in classlist:
        freq.append(item.freq)
        if item.normalize_density: # Check for normalize_density
            dielec.append(item.norm_dielec)
            losstan.append(item.norm_losstan)
        else:
            dielec.append(item.avg_dielec)
            losstan.append(item.avg_losstan)
        labels.append(item.name)
    kwargs["legend_label"] = labels
    if allplots:
        lossfac = []
        for item in classlist:
            lossfac.append(item.avg_lossfac)
        pplot.make_plot(freq,dielec,'d',**kwargs)
        pplot.make_plot(freq,lossfac,'lf',**kwargs)
        pplot.make_plot(freq,losstan,'lt',**kwargs)
    else:
        pplot.make_plot(freq,dielec,'d',**kwargs)
        pplot.make_plot(freq,losstan,'lt',**kwargs)
        
def run_default(airline_name='VAL',**kwargs):
    """
    Run AirlineData on get_METAS_data with all the prompts and return the \
        instance.
    """
    return AirlineData(*get_METAS_data(airline=airline_name),**kwargs)

def run_example(flag='single'):
    test = AirlineData(*get_METAS_data(airline='GAL',file_path=DATAPATH + \
                        '2.5hrs.txt'),bulk_density=2.0,temperature=None,\
                         name='Alumina Vac 2.5hrs',date='2017/04/07')
    if flag == 'single':
        atm = AirlineData(*get_METAS_data(airline='VAL',\
            file_path=DATAPATH + 'atm.txt'),bulk_density=None,\
            temperature=None,name='Alumina atm',date=None,corr=True,\
            solid_dielec=None,solid_losstan=None,particle_diameter=None,\
            particle_density=None,nrw=False)
        return test, atm
    elif flag == 'multiple':
        test2 = AirlineData(*get_METAS_data(),name='TRM')
        classlist = [test,test2]
        perm_compare(classlist)
        return test, test2, classlist
    
#%% MAIN
def main():
    ## Single file example:
    global test
    global atm
    test, atm = run_example()
    ## Multiple file example:
    #global test, test2, classlist
    #test, test2, classlist = run_example(flag='multiple')
    #pass    # Comment to run example
    
if __name__ == '__main__':
    main()
    #pass    # Comment to run example