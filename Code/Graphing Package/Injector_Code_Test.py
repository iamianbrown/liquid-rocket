import numpy as np
import matplotlib.pyplot as plt


class Injector:
    # -------------------------------------------Initial Data-------------------------------------------------#

    '''

    _r designates radial propellant
    _z designates axial propellant

    '''

    def __init__(self, chamber, d_c, rho_r, rho_z, d1, d2, C_d, delta_P, delta_P_o):  
        self.p_c = chamber.combustion_props['p_c']
        self.mdot_t = chamber.dict_values()['m_dot']
        self.OF_ratio = chamber.combustion_props['OF_ratio']
        
        self.injector_props = {
            'd_c': d_c,
            'rho_r': rho_r,
            'rho_z': rho_z,
            'd1': d1,
            'd2': d2,
            'C_d': C_d,
            'delta_P': delta_P,
            'delta_P_o': delta_P_o
        }
        
    def pintleParams(cea_data, gap_size, d_o, P_c):
        # cea_data is a list containing all relevant CEARUN parameters (T, OF, V_exit)
        # calculate the following parameters
        return theta, delta_p, BF, mdot_t

    P_c = 300 * 6894.76  # CHAMBER PRESSURE (Pascals)
    mdot_t = 1.989416305  # TOTAL MASS FLOW RATE (kg/s)
    OF = 3  # OXIDIZER TO FUEL RATIO
    d_c = .13  # CHAMBER DIAMETER (m)
    rho_r = 1141  # RADIAL PROPELLANT DENSITY (kg/m^3)
    rho_z = 810  # AXIAL PROPELLANT DENSITY (kg/m^3)
    d1 = 1.3  # ORIFICE ROW 1 DIAMETER (mm)
    d2 = 1.3  # ORIFICE ROW 2 DIAMETER (mm)
    C_d = 0.75  # DISCHARGE COEFFICIENT
    delta_P = 0.25 * P_c  # FUEL PRESSURE DROP
    delta_P_o = 0.15 * P_c  # OXIDIZER PRESSURE DROP

    # Calculate radial and axial propellant mass flow rates
    mdot_r = mdot_t / (1 + OF) * OF
    mdot_z = mdot_t / (1 + OF)

    # Calculate pintle diameter and radius using chamber-to-pintle ratio
    d_p = d_c / 4
    r_p = d_p / 2

    # Gap tolerance as provided by machine shop
    gap_tolerance = 0.05 / 1000

    # Constants needed for LMR and Spray Angle relation
    alpha = 0.7
    beta = 2.0

    # ---------------------------------------ORIFICE SIZING-------------------------------------------#

    # Orifice propellant discharge area

    '''
    _m designates value in meters
    _mm designates value in millimeters

    '''

    A_r_m = mdot_r / (C_d * np.sqrt(2 * rho_r * delta_P_o))  # m^2
    A_r_mm = A_r_m * 1000000  # mm^2

    a1 = np.pi * ((d1 * 0.5) ** 2)  # area of each orifice in row 1
    a2 = np.pi * ((d2 * 0.5) ** 2)  # area of "" row 2

    orifice_pairs = a1 + a2  # area per pair of orifices

    n_orifice_pairs = round(A_r_mm / orifice_pairs)

    percent_error = 100 * (n_orifice_pairs - (A_r_mm / orifice_pairs)) / (A_r_mm / orifice_pairs)

    BF = (n_orifice_pairs * (d1 + d2)) / (2 * np.pi * ((d_p * 1000) * 0.5))

    # -------------------------------------ANNULAR GAP SIZING-----------------------------------------#

    A_z_m = mdot_z / (C_d * np.sqrt(2 * rho_z * delta_P))
    A_z_mm = A_z_m * 1000000

    gap = np.sqrt((A_z_mm / np.pi) + ((r_p * 1000) ** 2)) - (r_p * 1000)
    gap_low = gap - (gap_tolerance * 1000)
    gap_high = gap + (gap_tolerance * 1000)

    # ---------------------------------------CALCULATE LMR AND THETA---------------------------------------#

    A_lr_c = np.pi * ((d1 / 2000) ** 2 + (d2 / 2000) ** 2)  # cross-sectional area for one orifice pair
    A_lz_c = (gap / 1000) * (
                d1 / 1000 + d2 / 1000)  # cross-sectional area of annular stream that impinges on the orifice pair

    U_r_c = mdot_r / (rho_r * A_r_m)
    U_z_c = mdot_z / (rho_z * A_z_m)

    LMR_c = (rho_r * (U_r_c ** 2) * A_lr_c) / (rho_z * (U_z_c ** 2) * A_lz_c)
    theta_c = alpha * np.arctan(beta * LMR_c) * (180 / np.pi) + 20 #theta_c = converging angle of nozzle

    # -------------------------------------------SETTING UP PLOTS-------------------------------------------------#

    # Array of gap sizes between 0.05 mm and 0.8 mm
    a = np.linspace(0.05 / 1000, 0.8 / 1000, 100000)

    A_lr = np.pi * ((d1 / 2 / 1000) ** 2)
    A_lz = (gap / 1000) * (d1 / 1000)
    U_r = mdot_r / (rho_r * A_r_m)
    U_z = mdot_z / (rho_z * np.pi * ((r_p + a) ** 2 - (r_p) ** 2))

    LMR = (rho_r * (U_r ** 2) * A_lr) / (rho_z * (U_z ** 2) * A_lz)
    theta = (alpha * np.arctan(beta * LMR) * (180 / np.pi)) + 20

    fig, axs = plt.subplots(2, 1, figsize=(12, 12))
    fig.suptitle('LMR and Theta vs annular gap')

    axs[0].plot(a, LMR)
    axs[1].plot(a, theta)

    axs[0].set_xlabel('annular gap size')
    axs[0].set_ylabel('LMR')
    axs[1].set_xlabel('annular gap size')
    axs[1].set_ylabel('spray angle')
    axs[0].locator_params(axis="x", nbins=15)
    axs[0].locator_params(axis="y", nbins=15)
    axs[1].locator_params(axis="x", nbins=15)
    axs[1].locator_params(axis="y", nbins=20)
    axs[0].set_xlim(a[0])
    axs[1].set_xlim(a[0])
    axs[0].set_ylim(0)
    axs[1].set_ylim(0)

    # --------------------------------LMR AND THETA TOLERANCED VALES---------------------------------------#

    LMR_low = np.interp(((gap / 1000) - gap_tolerance), a, LMR)
    LMR_high = np.interp(((gap / 1000) + gap_tolerance), a, LMR)

    theta_low = alpha * np.arctan(beta * LMR_low) * (180 / np.pi) + 20
    theta_high = alpha * np.arctan(beta * LMR_high) * (180 / np.pi) + 20

    axs[0].vlines(x=gap / 1000, label='With gap determined by DA eqn', linewidth=1, linestyle='dashed', color='red',
                  ymin=0, ymax=LMR_c)
    axs[0].vlines(x=(gap / 1000) - gap_tolerance, linewidth=1, linestyle='dashed', color='red', ymin=0, ymax=LMR_low)
    axs[0].vlines(x=(gap / 1000) + gap_tolerance, linewidth=1, linestyle='dashed', color='red', ymin=0, ymax=LMR_high)

    axs[1].vlines(x=gap / 1000, label='With gap determined by DA eqn', linewidth=1, linestyle='dashed', color='red',
                  ymin=0, ymax=theta_c)
    axs[1].vlines(x=(gap / 1000) - gap_tolerance, linewidth=1, linestyle='dashed', color='red', ymin=0, ymax=theta_low)
    axs[1].vlines(x=(gap / 1000) + gap_tolerance, linewidth=1, linestyle='dashed', color='red', ymin=0, ymax=theta_high)

    axs[0].hlines(y=LMR_low, xmin=0, xmax=(gap / 1000) - gap_tolerance, linewidth=1, linestyle='dashed', color='red')
    axs[0].hlines(y=LMR_high, xmin=0, xmax=(gap / 1000) + gap_tolerance, linewidth=1, linestyle='dashed', color='red')
    axs[1].hlines(y=theta_low, xmin=0, xmax=(gap / 1000) - gap_tolerance, linewidth=1, linestyle='dashed', color='red')
    axs[1].hlines(y=theta_high, xmin=0, xmax=(gap / 1000) + gap_tolerance, linewidth=1, linestyle='dashed', color='red')
    axs[0].hlines(y=LMR_c, xmin=0, xmax=(gap / 1000), linewidth=1, linestyle='dashed', color='red')
    axs[1].hlines(y=theta_c, xmin=0, xmax=(gap / 1000), linewidth=1, linestyle='dashed', color='red')



    # # ----------------------------WHAT WE WANT-----------------------------#
    #
    # gap_opt = np.interp(65, theta, a)
    # gap_opt_low = gap_opt - gap_tolerance
    # gap_opt_high = gap_opt + gap_tolerance
    #
    # LMR_opt = np.interp(gap_opt, a, LMR)
    # LMR_opt_low = np.interp(gap_opt_low, a, LMR)
    # LMR_opt_high = np.interp(gap_opt_high, a, LMR)
    #
    # theta_opt = np.interp(gap_opt, a, theta)
    # theta_opt_low = np.interp(gap_opt_low, a, theta)
    # theta_opt_high = np.interp(gap_opt_high, a, theta)
    #
    # axs[0].vlines(x=gap_opt, linewidth=1, linestyle='dashed', color='green', ymin=0, ymax=LMR_opt)
    # axs[0].vlines(x=gap_opt_low, linewidth=1, linestyle='dashed', color='green', ymin=0, ymax=LMR_opt_low)
    # axs[0].vlines(x=gap_opt_high, linewidth=1, linestyle='dashed', color='green', ymin=0, ymax=LMR_opt_high)
    # axs[0].hlines(y=LMR_opt, label='With gap adjusted for desired spray angle', xmin=0, xmax=gap_opt, linewidth=1,
    #               linestyle='dashed', color='green')
    # axs[0].hlines(y=LMR_opt_low, xmin=0, xmax=gap_opt_low, linewidth=1, linestyle='dashed', color='green')
    # axs[0].hlines(y=LMR_opt_high, xmin=0, xmax=gap_opt_high, linewidth=1, linestyle='dashed', color='green')
    #
    # axs[1].vlines(x=gap_opt, linewidth=1, linestyle='dashed', color='green', ymin=0, ymax=theta_opt)
    # axs[1].vlines(x=gap_opt_low, linewidth=1, linestyle='dashed', color='green', ymin=0, ymax=theta_opt_low)
    # axs[1].vlines(x=gap_opt_high, linewidth=1, linestyle='dashed', color='green', ymin=0, ymax=theta_opt_high)
    # axs[1].hlines(y=theta_opt, label='With gap adjusted for desired spray angle', xmin=0, xmax=gap_opt, linewidth=1,
    #               linestyle='dashed', color='green')
    # axs[1].hlines(y=theta_opt_low, xmin=0, xmax=gap_opt_low, linewidth=1, linestyle='dashed', color='green')
    # axs[1].hlines(y=theta_opt_high, xmin=0, xmax=gap_opt_high, linewidth=1, linestyle='dashed', color='green')

    axs[0].legend(loc='lower right')
    axs[1].legend(loc='lower right')

    plt.show()

    # ------------------------------------------OUTPUTS------------------------------------------------#

    print('Radial discharge area = {} mm^2'.format(round(A_r_mm, 2)))
    print('Axial discharge area = {} mm^2'.format(round(A_z_mm, 2)))
    print('Number of orifice pairs = {}'.format(n_orifice_pairs))
    print('Orifice area percent error = {}%'.format(round(percent_error, 2)))
    print('BF = {}'.format(round(BF, 2)))
    print('------------------------------------------')
    print('USING THE DISCHARGE AREA EQUATION:')
    print('Gap size = {}mm < {}mm < {}mm'.format(round(gap_low, 3), round(gap, 3), round(gap_high, 3)))
    print('LMR = {} < {} < {}'.format(round(LMR_low, 3), round(LMR_c, 3), round(LMR_high, 3)))
    print('Spray Angle = {}deg < {}deg < {}deg'.format(round(theta_low, 2), round(theta_c, 2), round(theta_high, 2)))
    # print('------------------------------------------')
    # print('FOR CHOICE OF SPRAY ANGLE')
    # print('Gap size = {}mm < {}mm <{}mm'.format(round(gap_opt_low * 1000, 3), round(gap_opt * 1000, 3),
    #                                             round(gap_opt_high * 1000, 3)))
    # print('LMR = {} < {} < {}'.format(round(LMR_opt_low, 3), round(LMR_opt, 3), round(LMR_opt_high, 3)))
    # print('Spray Angle = {}deg < {}deg < {}deg'.format(round(theta_opt_low, 2), round(theta_opt, 2),
    #                                                    round(theta_opt_high, 2)))