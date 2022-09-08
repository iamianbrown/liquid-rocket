import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('TkAgg')

def read_data(datafile):
    data = pd.read_csv(datafile, delim_whitespace=True, skiprows=[0,1])
    data.columns = ['POSN (mm)','LOAD (kN)', 'STRAIN1 %', 'STRAIN2 %', 'STRESS (MPA)', 'TIME (sec)']

    global dis
    dis = find_discontinuities(data['STRAIN2 %'], 2)
    d = dis[0]
    # edit strain column to include strain vs. position data

    m , b = calc_strain_vs_pos(data, dis[0])
    data['STRAIN2 %'][d:] = data['POSN (mm)'][d:] * m + b
    # data = data.iloc[:d + 1] # fix concatenation
    return data

def calc_strain_vs_pos(data, discontinuity_index):
    '''
    Returns a linear interpolation of the strain vs. position graph to fill in strain data after strain meter falls apart
    Inputs:
    data - data containing position and strain info
    Outputs:
    m, b - slope and intercept of linear fit
    '''
    fig, ax = plt.subplots(1)
    ax.plot(data['POSN (mm)'], data['STRAIN2 %'])

    m = (data['STRAIN2 %'][discontinuity_index] - data['STRAIN2 %'][discontinuity_index - 5])/(data['POSN (mm)'][discontinuity_index] - data['POSN (mm)'][discontinuity_index - 5])
    b = data['STRAIN2 %'][discontinuity_index] - m * data['POSN (mm)'][discontinuity_index - 5]
    ax.set_xlabel('Position, mm')
    ax.set_ylabel('Strain, %')
    ax.plot(data['POSN (mm)'], m * data['POSN (mm)'] + b)
    ax.axhline(data['STRAIN2 %'][dis[0]], color='red', linestyle='dashed')
    return m, b

def calc_yield(SS_data):
    '''
    Returns the yield strength of the material represented by SS_data
    Inputs: 
    SS_data - pandas dataframe containing stress-strain data
    Outputs:
    yield_strength - yield strength of the material
    '''
    E = calculate_E(SS_data)
    b = SS_data['STRESS (MPA)'][50] - SS_data['STRAIN2 %'][50] * E

    # y = 0 => mx = -b, x=-b/m
    # calculate x-intercept
    x_int = -1 * b / E

    # shift curve 0.2% to the right
    x_int_new = x_int + 0.2
    b_new = -1 * x_int_new * E

    # find intersection of shifted curve and data
    for index, row in SS_data.iterrows():
        if row['STRESS (MPA)'] < E * row['STRAIN2 %'] + b_new:
            return row['STRESS (MPA)']

def plot_SS(SS_data, figsize=(3,3)):
    '''
    Returns a matplotlib figure containing the stress-strain curve based on sample data.
    Input:
    SS_data - sample stress-strain data
    Outputs:
    fig, ax - matplotlib figure and axes containing the stress-strain plot
    '''
    fig, ax = plt.subplots(1, figsize=figsize)

    ax.plot(SS_data['STRAIN2 %'], SS_data['STRESS (MPA)'])
    # E = calculate_E(SS_data)
    # # plot young's modulus
    # strain = np.linspace(0, 0.2, 100)
    # stress = strain * E
    # ax.plot(strain, stress, color='green', linestyle='dashed')
    ax.set_xlabel('Strain, %')
    ax.set_ylabel('Stress, MPa')

    # ax.axvline(SS_data['STRAIN2 %'][dis[0]], color='red') # point where strain gauge fails
    # fig.show()
    return fig, ax

def calculate_E(SS_data):
    '''
    Calculate the Young's Modulus given a set of SS data
    Inputs:
    SS_data - a pandas dataframe containing the stress-strain data
    Ouputs:
    E - Young's Modulus
    '''
    E = (SS_data['STRESS (MPA)'][50] - SS_data['STRESS (MPA)'][0])/(SS_data['STRAIN2 %'][50] - SS_data['STRAIN2 %'][0])
    return E

def calculate_ultimate(SS_data):
    '''
    Calculate the ultimate strength by maximizing the SS data
    Inputs:
    SS_data - stress-strain data
    Output:
    ultimate_strength - ultimate strength of the material
    '''
    return max(SS_data['STRESS (MPA)'])

def calculate_fracture(SS_data):
    '''
    Calculate the fracture strength of the material represented by SS_data
    Inputs:
    SS_data - pandas dataframe containing stress-strain data
    Output:
    fracture_stress - stress at which fracture occurs
    '''
    for index, row in SS_data.iterrows():
        if row['STRESS (MPA)'] <= 10 and row['STRAIN2 %'] > 2:
            return SS_data['STRESS (MPA)'][index - 10]

def find_discontinuities(array, threshold):
    '''
    Finds the indicies of the discontinuities in array, defined as any jump larger than threshold

    Inputs:
    array - array like variable containing numerical values
    threshold - quantity by which a discontinuity is defined

    Outputs:
    discontinuity_indices - indices of array which precede discontinuities, as defined by threshold
    '''
    differences = np.zeros(len(array) - 1)
    for i in range(len(array) - 1):
        differences[i] = array[i + 1] - array[i]
    
    discontinuity_indices = []

    for d in range(len(differences)):
        if differences[d] > threshold:
            discontinuity_indices.append(d)
    
    return discontinuity_indices

if __name__ == '__main__':
    file = 'SS_Data/316/316-1d1.txt' # get from dialog
    data = read_data(file)
    fig, ax = plot_SS(data)
    print(data['STRAIN2 %'])
    print(data.head)
    print(calc_strain_vs_pos(data))

    plt.show()
