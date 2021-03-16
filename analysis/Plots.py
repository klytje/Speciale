import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.figure import figaspect
from pathlib import Path
from scipy.stats import chi2, sem
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes, mark_inset
from Functions import hist, exp_fit, exp_plot
from scipy.optimize import curve_fit

### GLOBALS ###
save_fig = True
step = 60
target = ['manganese', '1']
method = ['chi2', # standard minimal chi2 fit
          'chi2_poisson', # modified chi2 for a poisson distribution (which exponential decay follows)
          'time-interval'] # not implemented
method = method[1] # choose one of the above

model = [None, 'extending', 'non-extending'] # dead time model to use
model = model[0]
factor = 50 # dead time factor
folder = f'Data/{method}/{model}'

# select a plot from the list below
plots = ['snippet', # 0
         'dead time', # 1
         'dead time grouped', # 2
         'simulation', # 3
         'simulation distributions grouped', # 4
         'dead time rate dependence', # 5
         'pulser tests', # 6
         'spectra'] # 7
plots = plots[2]
target = ['simulated', int(9284.04/np.log(2))] # simulation

dt = 9.116e-06
interval = 0.5 # interval for dead time plots
cmax = 1000 # interval length for dead time plots
## SNIPPET OF DATA ###
if plots == 'snippet':
    data = np.array(pd.read_csv(Path(f'Data/{target[0]}_{target[1]}.txt'), delim_whitespace=True, skiprows=4))
    data = data[data[:, 1] > 0, :] # Ignore all counts with zero energy
    ch = np.array(data[:, 1])
    time = np.array(data[:, 0])*1e-8 # Converting to numpy arrays and seconds (MC2Analyzer outputs in units of 10ns)

    channels = []
    for i in range(len(time)):
        if time[i] < 600:
            channels.append(ch[i])
        else:
            break

    plt.figure()
    plt.hist(channels, bins=np.arange(0, 3000, 1))

### DEAD TIME PLOTS ###
elif plots == 'dead time':
    data = np.array(pd.read_csv(Path(f'Data/dead time distributions/correction_deadtime_test_{interval}_{cmax}_{target[0]}_{target[1]}_{method}_{model}.txt'), delim_whitespace=True))
    tau, chi = data[:, 0], data[:, 1]
    x = np.arange(interval, len(tau)*interval+interval, interval)
    tau = tau[x[:] > 0.5]
    chi = chi[x[:] > 0.5]
    x = x[x[:] > 0.5]

    plt.figure()
    plt.title(r'$\tau$')
    plt.plot(x, tau, '.r', label=r'$\tau$ measured')
    plt.xlabel('times the deadtime')
    plt.ylabel('value')
    plt.grid()
    if save_fig:
        plt.savefig(Path(f'Figures/dead time distribution/corrected_chi_distribution_{interval}_{cmax}_{model}_{target[0]}_{target[1]}.pdf'), dpi=300, format='pdf')

    plt.figure()
    plt.title(r'$\chi^2$')
    plt.plot(x, chi, '.r', label=r'$\chi^2$ measured')
    plt.xlabel('times the deadtime')
    plt.ylabel('value')
    plt.grid()
    if save_fig:
        plt.savefig(Path(f'Figures/dead time distribution/corrected_tau_distribution_{interval}_{cmax}_{model}_{target[0]}_{target[1]}.pdf'), dpi=300, format='pdf')

elif plots == 'dead time grouped':
    # data = np.loadtxt(Path(f'Data/{method}/{model}/info_{factor}_60_{target[0]}_{target[1]}.txt'), usecols=[1])  # Loading supplementary information
    # timestep, dof, N, duration, A, tau, c, gof = int(data[0]), int(data[1]), int(data[2]), data[3], data[4], data[5], data[6], data[7]  # Extracting information to base the fit on
    data = np.array(pd.read_csv(Path(f'Data/dead time distributions/{target[0]}_{target[1]}_chi2_extending.txt'), delim_whitespace=True))
    tau_e_u, tau_e, tau_e_l, chi_e, nfree_e, x_e = data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4], data[:, 5]*1e6 # tau is converted to seconds
    data = np.array(pd.read_csv(Path(f'Data/dead time distributions/{target[0]}_{target[1]}_chi2_non-extending.txt'), delim_whitespace=True))
    tau_ne_u, tau_ne, tau_ne_l, chi_ne, nfree_ne, x_ne = data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4], data[:, 5]*1e6 # tau is converted to seconds

    data = np.array(pd.read_csv(Path(f'Data/dead time distributions/{target[0]}_{target[1]}_chi2_poisson_extending.txt'), delim_whitespace=True))
    tau2_e_u, tau2_e, tau2_e_l, chi2_e, nfree2_e, x2_e = data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4], data[:, 5]*1e6 # tau is converted to seconds
    data = np.array(pd.read_csv(Path(f'Data/dead time distributions/{target[0]}_{target[1]}_chi2_poisson_non-extending.txt'), delim_whitespace=True))
    tau2_ne_u, tau2_ne, tau2_ne_l, chi2_ne, nfree2_ne, x2_ne = data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4], data[:, 5]*1e6 # tau is converted to seconds

    data = np.array(pd.read_csv(Path(f'Data/dead time distributions/{target[0]}_{target[1]}_time-interval_extending.txt'), delim_whitespace=True))
    tau_ti_u, tau_ti, tau_ti_l, chi_ti, nfree_ti, x_ti = data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4], data[:, 5]*1e6 # tau is converted to seconds
    # data = np.array(pd.read_csv(Path(f'Data/dead time distributions/{sim}correction_deadtime_test_5_100_{target[0]}_{target[1]}_time-interval_extending.txt'), delim_whitespace=True))
    # tau2_ti, unc2_ti, chi2_ti, x2_ti = data[:, 0]*60*np.log(2), data[:, 1]*60*np.log(2), data[:, 2], data[:, 3]*1e6 # tau is converted to seconds

    if target[0] == 'simulated':
        lim = 400
        lim2 = 80
    else:
        lim = 250
        lim2 = 40
    lim = 1000
    lim2 = 100
    tau_e_u, tau_e, tau_e_l, chi_e, nfree_e, x_e = tau_e_u[:lim], tau_e[:lim], tau_e_l[:lim], chi_e[:lim], nfree_e[:lim], x_e[:lim]
    tau_ne_u, tau_ne, tau_ne_l, chi_ne, nfree_ne, x_ne = tau_ne_u[:lim], tau_ne[:lim], tau_ne_l[:lim], chi_ne[:lim], nfree_ne[:lim], x_ne[:lim]
    tau2_e_u, tau2_e, tau2_e_l, chi2_e, nfree2_e, x2_e = tau2_e_u[:lim], tau2_e[:lim], tau2_e_l[:lim], chi2_e[:lim], nfree2_e[:lim], x2_e[:lim]
    tau2_ne_u, tau2_ne, tau2_ne_l, chi2_ne, nfree2_ne, x2_ne = tau2_ne_u[:lim], tau2_ne[:lim], tau2_ne_l[:lim], chi2_ne[:lim], nfree2_ne[:lim], x2_ne[:lim]
    tau_ti_u, tau_ti, tau_ti_l, chi_ti, nfree_ti, x_ti = tau_ti_u[:lim2], tau_ti[:lim2], tau_ti_l[:lim2], chi_ti[:lim2], nfree_ti[:lim2], x_ti[:lim2]

    dt = dt*1e6 # scales with the axes
    if target[0] == 'manganese':
        tau_ref = 9284.04
        tau_ref_std = 0.36
    elif target[0] == 'copper':
        tau_ref = 12.700*60*60
        tau_ref_std = 0.002*60*60
    else: # simulation
        tau_ref = 9284.04
        tau_ref_std = 0.36

    def plot_ref():
        plt.axhline(y=tau_ref, label='Reference', color='k', zorder=100)
        if target[0] != 'simulated':
            plt.axvline(dt, color='k', label=r'detector dead time')

    ### chi2 tau ###
    fig, ax = plt.subplots(3, 1, sharex=True)
    plt.sca(ax[0])
    plt.title(r'$\chi^2$')
    plt.errorbar(x_e, tau_e, fmt='.', yerr=[tau_e_u, -tau_e_l], alpha=0.5, label=r'$t_{1/2, o}$, $\chi^2_E$')
    plot_ref()
    plt.ylabel(r'$t_{1/2}$ [s]')
    plt.grid(axis='x')
    ax2 = ax[0].twinx()
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)  # hides all info of the twinned axis
    plt.ylabel('extending')

    plt.sca(ax[1])
    plt.errorbar(x_ne, tau_ne, fmt='.', yerr=[tau_ne_u, -tau_ne_l], alpha=0.5, label=r'$t_{1/2, o}$, $\chi^2_{NE}$')
    plot_ref()
    plt.ylabel(r'$t_{1/2}$ [s]')
    plt.grid(axis='x')
    ax2 = ax[1].twinx()
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)  # hides all info of the twinned axis
    plt.ylabel('non-extending')

    plt.sca(ax[2])
    plt.plot(x_e, chi_e / nfree_e, '.', label=r'$\chi^2_E$')
    plt.plot(x_ne, chi_ne / nfree_ne, '.', label=r'$\chi^2_{NE}$')
    plt.axhline(y=1, color='k')
    plt.xlabel(r'$\mu s$', horizontalalignment='right', x=1.0)
    plt.ylabel(r'$\chi^2$')
    plt.grid(axis='x')
    plt.legend()

    if save_fig:
        if target[0] == 'simulated':
            plt.savefig(Path(f'Figures/dead time distribution/chi2_sim.pdf'), dpi=300, format='pdf')
        else:
            plt.savefig(Path(f'Figures/dead time distribution/chi2.pdf'), dpi=300, format='pdf')

    ### poisson tau ###
    fig, ax = plt.subplots(3, 1, sharex=True)
    plt.sca(ax[0])
    plt.title(r'$\chi^2_\lambda$')
    plt.errorbar(x2_e, tau2_e, fmt='.', yerr=[tau2_e_u, -tau2_e_l], alpha=0.5, label=r'$t_{1/2, o}$, $\chi^2_{\lambda, E}$')
    plot_ref()
    plt.ylabel(r'$t_{1/2}$ [s]')
    plt.grid(axis='x')
    ax2 = ax[0].twinx()
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)  # hides all info of the twinned axis
    plt.ylabel('extending')

    plt.sca(ax[1])
    plt.errorbar(x2_ne, tau2_ne, fmt='.', yerr=[tau2_ne_u, -tau2_ne_l], alpha=0.5, label=r'$t_{1/2, o}$, $\chi^2_{\lambda, NE}$')
    plot_ref()
    plt.ylabel(r'$t_{1/2}$ [s]')
    plt.grid(axis='x')
    ax2 = ax[1].twinx()
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)  # hides all info of the twinned axis
    plt.ylabel('non-extending')
    plt.grid(axis='x')

    plt.sca(ax[2])
    plt.plot(x2_e, chi2_e/nfree2_e, '.', label=r'$\chi^2_{\lambda, E}$')
    plt.plot(x2_ne, chi2_ne/nfree2_ne, '.', label=r'$\chi^2_{\lambda, NE}$')
    plt.axhline(y=1, color='k')
    plt.xlabel(r'$\mu s$', horizontalalignment='right', x=1.0)
    plt.ylabel(r'$\chi^2$')
    plt.grid(axis='x')
    plt.legend()
    if save_fig:
        if target[0] == 'simulated':
            plt.savefig(Path(f'Figures/dead time distribution/chi2_poisson_sim.pdf'), dpi=300, format='pdf')
        else:
            plt.savefig(Path(f'Figures/dead time distribution/chi2_poisson.pdf'), dpi=300, format='pdf')

    ### time interval tau ###
    fig, ax = plt.subplots(2, 1, sharex=True)
    plt.sca(ax[0])
    plt.errorbar(x_ti, tau_ti, fmt='.', alpha=0.5, yerr=[tau_ti_u, -tau_ti_l], label=r'time interval')
    plot_ref()
    plt.title(r'time-interval')
    plt.ylabel(r'$t_{1/2}$')
    plt.grid(axis='x')

    plt.sca(ax[1])
    plt.plot(x_ti, chi_ti/nfree_ti, '.')
    plt.axhline(y=2*np.euler_gamma, color='k')
    plt.axhline(y=2*np.euler_gamma + 1.606/np.sqrt(33393147), color='k', linestyle='--')
    plt.axhline(y=2*np.euler_gamma - 1.606/np.sqrt(33393147), color='k', linestyle='--')
    plt.xlabel(r'$\mu s$', horizontalalignment='right', x=1.0)
    plt.ylabel(r'$\chi^2$')
    plt.grid(axis='x')
    if save_fig:
        if target[0] == 'simulated':
            plt.savefig(Path(f'Figures/dead time distribution/time-interval_sim.pdf'), dpi=300, format='pdf')
        else:
            plt.savefig(Path(f'Figures/dead time distribution/time-interval.pdf'), dpi=300, format='pdf')

    ### time interval degrees of freedom ###
    plt.figure()
    plt.plot(x_ti, nfree_ti, '.')
    plt.xlabel(r'$\mu s$', horizontalalignment='right', x=1.0)
    plt.ylabel(r'Degrees of freedom')

    plt.show()

### SIMULATION PLOTS ###
elif plots == 'simulation':
    # data = np.loadtxt(Path(f'{folder}/info_{factor}_{target[0]}_{target[1]}.txt'), usecols=[1])  # Loading supplementary information
    results = np.stack(np.loadtxt(Path(f'{folder}/simulated_fits_{step}_{target[0]}_{target[1]}.txt')), axis=1) # {factor}_{target[0]}_{target[1]}
    # step, dof, O, duration, A, tau, c, gof = int(data[0]), int(data[1]), int(data[2]), data[3], data[4], data[5], data[6], data[7]
    tau_values, A_values, c_values, gof_values, dof_values = results[0], results[1], results[2], results[3], results[4]
    t12_values = tau_values*step*np.log(2)
    bins = np.arange(5000)
    tau = 9284.04/np.log(2)
    t12 = 9284.04
    O = 2e7
    dof = 1113

    # simulation plot
    fig, ax = plt.subplots()
    exp_plot(np.column_stack([tau_values, A_values, c_values]), bins, ax, alpha=0.01)
    # exp_plot([tau, A, c], bins, ax, fmt='r-')
    plt.title('Simulation')
    plt.xlabel(f'Time [{step}s]')  # Ensure this is the case!
    plt.ylabel('Counts')

    # figure inset
    axins = zoomed_inset_axes(ax, 120, loc=1) # set zoom factor and location
    exp_plot(results, bins, axins, alpha=0.01)
    # exp_plot([tau, A, c], bins, axins, fmt='r-')
    if method == 'chi2':
        x1, x2, y1, y2 = 198, 202, 36000, 36500  # specify the limits
    elif method == 'chi2_poisson':
        x1, x2, y1, y2 = 202, 206, 36000, 36500 # they change slightly
    elif method == 'time-interval':
        x1, x2, y1, y2 = 300, 600, 30000, 40000 # adjust as necessary
    else:
        raise Exception('Invalid method statement')
    axins.set_xlim(x1, x2)  # apply the x-limits
    axins.set_ylim(y1, y2)  # apply the y-limits
    plt.yticks(visible=False) # hides ticks on the inset
    plt.xticks(visible=False)
    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5") #facecolor, edgecolor
    if save_fig:
        plt.savefig(Path(f'Figures/Simulation_plot_{factor}_{method}.pdf'), dpi=300, format='pdf')

    # tau figure
    simt_std = np.std(t12_values) # standard error
    simt_avg = np.average(t12_values)
    print(f'Tau value: {simt_avg/3600}, standard error: {simt_std/3600}, difference from actual value of {t12/3600}: {simt_avg/3600 - t12/3600}')
    print(f'The decay is {(simt_avg/t12 - 1)*100}% higher than expected')
    plt.figure()
    plt.hist(t12_values, bins=np.linspace(np.min(t12_values), np.max(t12_values), num=30))
    plt.axvline(x=simt_avg, color='k', label='Simulated average')
    plt.axvline(x=simt_avg + simt_std, color='k', linestyle='--')
    plt.axvline(x=simt_avg - simt_std, color='k', linestyle='--')
    plt.axvline(x=t12, color='r', label='Expected value')
    plt.title(r'$\tau$ distribution')
    plt.xlabel(r'$\tau$')
    plt.ylabel('Count')
    plt.legend()
    if save_fig:
        plt.savefig(Path(f'Figures/tau_plot_{factor}_{method}.pdf'), dpi=300, format='pdf')

    # c figure
    simc_std = np.std(c_values)
    simc_avg = np.average(c_values)
    print(f'c value: {simc_avg}, standard error: {simc_std}')
    plt.figure()
    plt.hist(c_values, bins=np.linspace(np.min(c_values), np.max(c_values), num=30))
    plt.axvline(x=simc_avg, color='k', label='Simulated average')
    plt.axvline(x=simc_avg + simc_std, color='k', linestyle='--')
    plt.axvline(x=simc_avg - simc_std, color='k', linestyle='--')
    plt.axvline(x=50*step, color='r', label='Expected value')
    plt.title(r'$c$ distribution')
    plt.xlabel(r'$c$')
    plt.ylabel('Count')

    # gof figure for time-interval
    if method == 'time-interval':
        plt.figure()
        gof_dof_values = np.divide(gof_values, dof_values)
        simgof_avg = np.average(gof_dof_values)
        actual_val = 2*np.euler_gamma
        actual_err = 1.606/np.sqrt(O)
        plt.hist(gof_dof_values, bins=np.linspace(np.min(gof_dof_values), np.max(gof_dof_values), num=30))
        plt.axvline(x=simgof_avg, color='k', label='Simulated average')
        plt.axvline(x=actual_val, color='r', label='Expected value')
        plt.axvline(x=actual_val + actual_err, color='r', linestyle='--')
        plt.axvline(x=actual_val - actual_err, color='r', linestyle='--')
        plt.legend()

    # chi2 figure
    else:
        plt.figure()
        avg = np.average(gof_values)
        std = np.std(gof_values)
        print(f'Chi2 value: {avg}, standard error: {std}')
        _, bins, _ = plt.hist(gof_values, bins=np.linspace(np.min(gof_values), np.max(gof_values), num=30), density=True)
        plt.plot(bins, chi2.pdf(bins, dof), 'r-', lw=2, label='chi2 pdf')
        plt.axvline(x=avg, color='k', label='Simulated average')
        plt.axvline(x=avg + std, color='k', linestyle='--')
        plt.axvline(x=avg - std, color='k', linestyle='--')
        # plt.axvline(x=gof, color='r', label='Measured value')
        plt.title(r'$\chi^2$ distribution')
        plt.xlabel('Value')
        plt.ylabel('Density')
        plt.legend()
        # plt.axis([np.min(gof_values), np.max(gof_values), 0, None])
        if save_fig:
            plt.savefig(Path(f'Figures/Chi2_plot_{factor}_{method}.pdf'), dpi=300, format='pdf')

elif plots == 'simulation distributions grouped':
    method = ['chi2', 'chi2_poisson']
    method_labels = ['$\chi^2$', '$\chi^2_\lambda$']
    model = [None, 'extending', 'non-extending']
    model_title = ['none', 'extending', 'non-extending']
    dof = 2229 # only for this specific case of simulated data - change it!
    tau = 9284.04/3600
    ### TAU FIGURE ###
    fig, ax = plt.subplots(2, 3, sharey='row', sharex='col')
    fig.add_subplot(111, frameon=False) # adds a big, single figure on top of everything else
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False) # hides all features of the big figure
    plt.xlabel(r'$t_{1/2}$') # adds labels to the big figure, so they are always nicely centered
    plt.ylabel('Count')
    plt.title(r'$t_{1/2}$ distribution')
    bins = np.arange(5000)
    for i in range(len(ax)):
        for j in range(len(ax[i])):
            folder = f'Data/{method[i]}/{model[j]}' # this changes
            # data = np.loadtxt(Path(f'{folder}/info_{factor}_{target[0]}_{target[1]}.txt'), usecols=[1])  # Loading supplementary information
            results = np.stack(np.loadtxt(Path(f'{folder}/simulated_fits_{step}_{target[0]}_{target[1]}.txt')), axis=1)
            # step, dof, O, duration, A, tau, c, gof = int(data[0]), int(data[1]), int(data[2]), data[3], data[4], data[5], data[6], data[7]
            tau_values, gof_values, dof_values = results[0]/3600, results[3], results[4]
            tau_values = tau_values*step*np.log(2)
            simtau_std = np.std(tau_values)  # standard error
            simtau_avg = np.average(tau_values)
            print(f'Average tau:  {simtau_avg}, standard deviation: {simtau_std}')
            plt.sca(ax[i, j])
            xticks = [2.577, 2.580]
            plt.xticks(xticks, ['2.577', '2.580'])
            plt.hist(tau_values, bins=np.linspace(np.min(tau_values), np.max(tau_values), num=20))
            plt.axvline(x=simtau_avg, color='k', label='Simulated average')
            plt.axvline(x=simtau_avg + simtau_std, color='k', linestyle='--')
            plt.axvline(x=simtau_avg - simtau_std, color='k', linestyle='--')
            plt.axvline(x=tau, color='r', label='Measured value')
            if i == 0 and j == 2:
                ax2 = ax[i, j].twinx()
                plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False) # hides all info of the twinned axis
                plt.ylabel(r'$\chi^2$', rotation=0)
            elif i == 1 and j == 2:
                ax2 = ax[i, j].twinx()
                plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False) # hides all info of the twinned axis
                plt.ylabel(r'$\chi^2_\lambda$', rotation=0)
            if i == 1:
                plt.title(model_title[j])

    if save_fig:
        plt.savefig(Path(f'Figures/tau_plot_all_{factor}.pdf'), dpi=300, format='pdf')

    ### CHI FIGURE ###
    fig, ax = plt.subplots(2, 3, sharey='row', sharex='col')
    fig.add_subplot(111, frameon=False) # adds a big, single figure on top of everything else
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False) # hides all features of the big figure
    plt.xlabel(r'$\chi^2$') # adds labels to the big figure, so they are always nicely centered
    plt.ylabel('Density')
    plt.title(r'$\chi^2$ distributions')
    for i in range(len(ax)):
        for j in range(len(ax[i])):
            folder = f'Data/{method[i]}/{model[j]}' # this changes
            # data = np.loadtxt(Path(f'{folder}/info_{factor}_{target[0]}_{target[1]}.txt'), usecols=[1])  # Loading supplementary information
            results = np.stack(np.loadtxt(Path(f'{folder}/simulated_fits_{step}_{target[0]}_{target[1]}.txt')),axis=1)
            # step, dof, O, duration, A, tau, c, gof = int(data[0]), int(data[1]), int(data[2]), data[3], data[4], data[5], data[6], data[7]
            tau_values, gof_values, dof_values = results[0], results[3], results[4]

            std = np.std(gof_values)
            avg = np.average(gof_values)
            print(f'Average gof:  {avg/dof}, standard deviation: {std/dof}')
            plt.sca(ax[i, j])
            xticks = [2100, 2300]
            plt.xticks(xticks, xticks)
            if i == 0 and j == 2:
                ax2 = ax[i, j].twinx()
                plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False) # hides all info of the twinned axis
                plt.ylabel(r'$\chi^2$', rotation=0)
            elif i == 1 and j == 2:
                ax2 = ax[i, j].twinx()
                plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False) # hides all info of the twinned axis
                plt.ylabel(r'$\chi^2_\lambda$', rotation=0)
            if i == 1:
                plt.title(model_title[j])
            plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
            std_range = 3 # the range of the chi2 pdf
            bins = np.arange(avg - std_range*std, avg + std_range*std, 1)
            plt.plot(bins, chi2.pdf(bins, dof), '-', lw=2, label='chi2 pdf')
            plt.axvline(x=avg, color='k', label='Simulated average')
            plt.axvline(x=avg + std, color='k', linestyle='--')
            plt.axvline(x=avg - std, color='k', linestyle='--')
            # plt.axvline(x=gof, color='r', label='Measured value')
    if save_fig:
        plt.savefig(Path(f'Figures/chi2_plot_all_{factor}.pdf'), dpi=300, format='pdf')

    ### FIGURE TIME-INTERVAL ###
    fig, ax = plt.subplots(2, 1, figsize=(5, 5))
    folder = f'Data/time-interval/extending'
    # data = np.loadtxt(Path(f'{folder}/info_{factor}_{target[0]}_{target[1]}.txt'), usecols=[1])  # Loading supplementary information
    results = np.stack(np.loadtxt(Path(f'{folder}/simulated_fits_{step}_{target[0]}_{target[1]}.txt')), axis=1)
    # step, dof, O, duration, A, tau, c, gof = int(data[0]), int(data[1]), int(data[2]), data[3], data[4], data[5], data[6], data[7]
    tau_values, gof_values, dof_values = results[0]/3600, results[3], results[4]

    tau_values = tau_values*step*np.log(2)
    ### TAU ###
    plt.sca(ax[0])
    xticks = [2.578, 2.580]
    plt.xticks(xticks, ['2.578', '2.580'])
    simtau_std = np.std(tau_values) # standard error
    simtau_avg = np.average(tau_values)
    print(f'Average tau:  {simtau_avg}, standard deviation: {simtau_std}')
    plt.hist(tau_values, bins=np.linspace(np.min(tau_values), np.max(tau_values), num=20))
    plt.axvline(x=simtau_avg, color='k', label='Simulated average')
    plt.axvline(x=simtau_avg + simtau_std, color='k', linestyle='--')
    plt.axvline(x=simtau_avg - simtau_std, color='k', linestyle='--')
    plt.axvline(x=tau, color='r', label='Measured value')
    plt.xlabel(r'$t_{1/2}$')
    plt.ylabel('Count')
    plt.title(r'$t_{1/2}$ distribution')

    ### GoF ###
    O = 33391963 # change this!!
    plt.sca(ax[1])
    xticks = [1.1535, 1.1540, 1.1545, 1.1550]
    plt.xticks(xticks, ['1.1535', '1.1540', '1.1545', '1.1550'])
    gof_dof_values = np.divide(gof_values, dof_values)
    simgof_avg = np.average(gof_dof_values)
    simgof_std = np.std(gof_dof_values)
    actual_val = 2 * np.euler_gamma
    actual_err = 1.606 / np.sqrt(O)
    print(f'Average gof:  {simgof_avg}, standard deviation: {simgof_std}')
    plt.hist(gof_dof_values, bins=np.arange(np.min(gof_dof_values), np.max(gof_dof_values), 0.0001))
    plt.axvline(x=simgof_avg, color='k', label='Simulated average')
    plt.axvline(x=actual_val, color='r', label='Actual value')
    plt.axvline(x=actual_val + actual_err, color='r', linestyle='--')
    plt.axvline(x=actual_val - actual_err, color='r', linestyle='--')
    plt.xlabel(r'$GoF/DoF$')
    plt.ylabel('Count')
    plt.title('GoF/DoF distribution')

    plt.tight_layout()
    if save_fig:
        plt.savefig(Path(f'Figures/tau_and_gof_plot_time_interval.pdf'), dpi=300, format='pdf')

elif plots == 'dead time rate dependence':
    dt = 9.116e-06

    plt.figure(figsize=figaspect(1/3))
    m_e = lambda n: n*np.exp(-n*dt)
    m_ne = lambda n: n/(1 + n*dt)
    x = np.arange(0, 5/dt, 1/(dt*100))
    plt.plot(x, m_e(x), label='extending')
    plt.plot(x, m_ne(x), label='non-extending')
    _, xmax, _, ymax = plt.axis()
    plt.plot(x, x, '--k', label=r'$m = n$')
    plt.axis([0, xmax, 0, ymax])
    plt.grid()
    plt.legend()
    xticks = np.linspace(0, 5/dt, num=5)
    yticks = np.linspace(0, 1/dt, num=5)
    # plt.xticks(xticks, ['0', r'$1/\tau$', r'$2/\tau$', r'$3/\tau$', r'$4/\tau$'])
    plt.xticks(xticks, ['0', r'$\tau^{-1}$', r'$2 \tau^{-1}$', r'$3 \tau^{-1}$', r'$4 \tau^{-1}$'])
    # plt.yticks(yticks, ['0', r'$1/(4\tau)$', r'$2/(4\tau)$', r'$3/(4\tau)$', r'$1/\tau$'])
    plt.yticks(yticks, ['0', r'$1/4 \tau^{-1}$', r'$2/4 \tau^{-1}$', r'$3/4 \tau^{-1}$', r'$\tau^{-1}$'])
    plt.xlabel(r'$n$', horizontalalignment='right', x=1)
    plt.ylabel(r'$m$', y=1, rotation=0)
    # plt.title('Observed count rate dependence on the true count rate') # is a title necessary?
    if save_fig:
        plt.savefig(Path(f'Figures/m_n_dependency.pdf'), dpi=300, format='pdf')

elif plots == 'pulser tests':
    n = [0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0] # true rates in kcps
    bg = np.array(pd.read_csv(Path(f'Data/pulser/bg_ch001.txt'), delim_whitespace=True, skiprows=4))
    bg = bg[bg[:, 1] > 0, 0] # ignore all counts with zero energy
    bg_cps = len(bg)/(int(bg[-1])*1e-8) # converting to kcps
    m = []
    for fname in n:
        curr = np.array(pd.read_csv(Path(f'Data/pulser/{fname}_ch001.txt'), delim_whitespace=True, skiprows=4))
        curr = curr[curr[:, 1] > 0, 0] # ignore all counts with zero energy
        m.append(len(curr)/(int(curr[-1])*1e-8) - bg_cps)
    n = np.array(n)*1e3
    plt.figure()
    plt.plot(n, m, 'k.')
    m_theo = lambda n, tau: n*np.exp(-n*tau)
    params, _ = curve_fit(m_theo, n, m, p0=1e-5)
    print(params)
    bins = np.linspace(0, n[-1], num=100)
    plt.plot(bins, m_theo(bins, params[0]))
    plt.title('Dead time behaviour')
    plt.xlabel(r'True rate $n$ [MHz]')
    plt.ylabel(r'Measured rate $m$ [MHz]')
    plt.xticks([0, 2000, 4000, 6000, 8000, 10000], [0, 2, 4, 6, 8, 10])
    plt.yticks([0, 500, 1000, 1500, 2000, 2500, 3000], [0, 0.5, 1, 1.5, 2, 2.5, 3])

    if save_fig:
        plt.savefig(Path(f'Figures/pulser/m_n_pulser.pdf'), dpi=300, format='pdf')

    plt.figure()
    time = np.array(pd.read_csv(Path(f'Data/pulser/2kcps.txt'), delim_whitespace=True, skiprows=4))
    ch = time[:, 1]
    time = time[time[:, 1] > 5000, 0]*1e-8  # ignore all counts with zero energy
    i = time[1:] - time[:-1]
    plt.hist(i, bins=np.linspace(np.min(i), np.max(i), num=100))

elif plots == 'spectra':
    bins = 3000
    mn = np.array(pd.read_csv(Path(f'Data/manganese_2.txt'), delim_whitespace=True, skiprows=4, usecols=[1]))
    cu = np.array(pd.read_csv(Path(f'Data/copper_1.txt'), delim_whitespace=True, skiprows=4, usecols=[1]))
    cs = np.array(pd.read_csv(Path(f'Data/two-source/measurement 2/first and second_ch001.txt'), delim_whitespace=True, skiprows=4, usecols=[1]))

    fig, ax = plt.subplots(3, 1, sharex=True)
    plt.sca(ax[0])
    plt.title('Spectra')
    plt.hist(mn, bins=np.arange(bins))
    plt.yscale('log')
    plt.ylim(1e3, None)
    ax2 = ax[0].twinx()
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)  # hides all info of the twinned axis
    plt.ylabel(r'$^{56}$Mn', rotation=0)

    plt.sca(ax[1])
    plt.hist(cu, bins=np.arange(bins))
    plt.yscale('log')
    plt.ylabel('Count')
    plt.ylim(1e3, None)
    ax2 = ax[1].twinx()
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)  # hides all info of the twinned axis
    plt.ylabel(r'$^{64}$Cu', rotation=0)

    plt.sca(ax[2])
    plt.hist(cs, bins=np.arange(bins))
    plt.yscale('log')
    plt.xlabel('Energy [keV]')
    plt.ylim(1e1, None)
    ax2 = ax[2].twinx()
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)  # hides all info of the twinned axis
    plt.ylabel(r'$^{137}$Cs', rotation=0)

    if save_fig:
        plt.savefig(Path(f'Figures/spectra.pdf'), dpi=300, format='pdf')


    # mn = np.array(pd.read_csv(Path(f'Data/manganese_1.txt'), delim_whitespace=True, skiprows=4))
    # cu = np.array(pd.read_csv(Path(f'Data/copper_1.txt'), delim_whitespace=True, skiprows=4))
    # cs = np.array(pd.read_csv(Path(f'Data/two-source/measurement 2/first and second_ch001.txt'), delim_whitespace=True, skiprows=4))
    # bg = np.array(pd.read_csv(Path(f'Data/two-source/measurement 2/background_ch001.txt'), delim_whitespace=True, skiprows=4))
    # t_mn = mn[-1, 0]
    # t_cu = cu[-1, 0]
    # t_cs = cs[-1, 0]
    # t_bg = bg[-1, 0]
    # bg = bg[:, 1]
    # mn = mn[:, 1]
    # cu = cu[:, 1]
    # cs = cs[:, 1]
    #
    # fig, ax = plt.subplots(3, 1, sharex=True)
    # bins = np.arange(0, bins, 1)
    # n_bg, _ = np.histogram(bg, bins=bins)
    #
    # plt.sca(ax[0])
    # n, x = np.histogram(mn, bins=bins)
    # n = n - n_bg/t_bg*t_mn
    # n[n[:] < 0] = 0
    # plt.hist(x[:-1], weights=n, bins=bins)
    # plt.yscale('log')
    #
    # plt.sca(ax[1])
    # n, x = np.histogram(cu, bins=bins)
    # n = n - n_bg/t_bg*t_cu
    # n[n[:] < 0] = 0
    # plt.hist(x[:-1], weights=n, bins=bins)
    # plt.yscale('log')
    #
    # plt.sca(ax[2])
    # n, x = np.histogram(cs, bins=bins)
    # n = n - n_bg/t_bg*t_cs
    # n[n[:] < 0] = 0
    # plt.hist(x[:-1], weights=n, bins=bins)
    # plt.yscale('log')

plt.show()