# Generate a raster plot of spiking

from netpyne import __gui__

if __gui__:
    import matplotlib.patches as mpatches
from ..analysis.utils import exception  # , loadData
from ..analysis.tools import loadData
from .plotter import ScatterPlotter


@exception
def plotRaster(
    rasterData=None,
    axis=None,
    timeRange=None,
    maxSpikes=1e8,
    orderBy='gid',
    popRates=True,
    popNumCells=None,
    popLabels=None,
    popColors=None,
    syncLines=False,
    colorbyPhase = None,
    legend=True,
    colorList=None,
    orderInverse=False,
    returnPlotter=False,
    **kwargs
):
    """Function to produce a raster plot of cell spiking

    NetPyNE Options
    ---------------
    include : str, int, list
        Cells and/or NetStims to return information from.

        *Default:* ``['allCells']`` includes all cells and no NetStims

        *Options:*
        (1) ``'all'`` includes all cells and all NetStims,
        (2) ``'allNetStims'`` includes all NetStims but no cells,
        (3) a *str* which matches a popLabel includes all cells in that pop,
        (4) a *str* which matches a NetStim name includes that NetStim,
        (5) an *int* includes the cell with that global identifier (GID),
        (6) a *list* of *ints* includes the cells with those GIDS,
        (7) a *list* with two items, the first of which is a *str* matching a popLabel and the second of which is an *int* (or a *list* of *ints*), includes the relative cell(s) from that population (e.g. (``['popName', [0, 1]]``) includes the first two cells in popName.

    sim : NetPyNE sim object
        The *sim object* from which to get data.

        *Default:* ``None`` uses the current NetPyNE sim object

    Parameters
    ----------
    rasterData : list, tuple, dict, str
        The data necessary to plot the raster (spike times and spike indices, at minimum).

        *Default:* ``None`` uses ``analysis.prepareRaster`` to produce ``rasterData`` using the current NetPyNE sim object.

        *Options:* if a *list* or a *tuple*, the first item must be a *list* of spike times and the second item must be a *list* the same length of spike indices (the id of the cell corresponding to that spike time).  Optionally, a third item may be a *list* of *ints* representing the number of cells in each population (in lieu of ``popNumCells``).  Optionally, a fourth item may be a *list* of *strs* representing the population names (in lieu of ``popLabels``).

        If a *dict* it must have keys ``'spkTimes'`` and ``'spkInds'`` and may optionally include ``'popNumCells'`` and ``'popLabels'``.

        If a *str* it must represent a file path to previously saved data.

    axis : matplotlib axis
        The axis to plot into, allowing overlaying of plots.

        *Default:* ``None`` produces a new figure and axis.

    timeRange : list
        Time range to include in the raster: ``[min, max]``.

        *Default:* ``None`` uses the entire simulation

    maxSpikes : int
        The maximum number of spikes to include (by reducing the max time range).

        *Default:* ``1e8``

    orderBy : str
        How to order the cells along the y-axis.

        *Default:* ``'gid'`` orders cells by their index

        *Options:* any NetPyNe cell tag, e.g. ``'pop'``, ``'x'``, ``'ynorm'`` .

    popRates : bool
        whether to include the spiking rates in the plot title and legend.

        *Default:* ``True`` includes detailed pop information on plot.

        *Options:*
        ``False`` only includes pop names.
        ``'minimal'`` includes minimal pop information.

    popNumCells : list
        A *list* of *ints* representing the number of cells in each population.

        *Default:* ``None`` puts all cells into a single population.

    popLabels : list
        A *list* of *strs* of population names.  Must be the same length as ``popNumCells``.

        *Default:* ``None`` uses generic names.

    popColors : dict
        A *dict* of ``popLabels`` and their desired color.

        *Default:* ``None`` draws from the NetPyNE default colorList.

    syncLines : bool
        Calculate synchrony measure and plot vertical lines for each spike to evidence synchrony if ``True``.

        *Default:* ``False``

    colorbyPhase : dict
        Dictionary specifying conditions to plot spikes colored by the phase of a simultaneous signal, filtered in a given range

            *Default:* ``None`` colors spikes according to other options (by populations)

            *Dictionary entries:*

            ``'signal'`` specifies the signal. Options are: ``'LFP'``, which takes the signal from the local fiel potential generated in the ongoing simulation, a numpy array of scalars (for example, an external signal used for stimulation), or an external pickle file,

            ``'fs'`` is the sampling frequency, which should be specified when the signal is obtained from external sources (pickle file or numpy array). Otherwise, it is assumed to be 1000 Hz. If the signal is specified by ``'LFP'``, then the sampling rate is obtained from the internal simulation (cfg.recordStep),

            ``'electrode'`` selects the electrode from the LFP setup. Default is electrode 1,

            ``'filtFreq'`` is a list specifying the range for filtering the signal (band-pass). For example, ``[4,8]`` to select theta rhythm. The default is a very broadband filtering (essentially, the raw signal) ``[1,500]``,

            ``'filtOrder'`` is the filter order (Butterworth) to process the signal,

            ``'pop_background'`` is a boolean option to color each population alternately with a gray background, for better visualization. The default is False,

            ``'include_signal'`` is a boolean option to plot the filtered signal below the raster plot. The default is False.

    legend : bool
        Whether or not to add a legend to the plot.

        *Default:* ``True`` adds a legend.

    colorList : list
        A *list* of colors to draw from when plotting.

        *Default:* ``None`` uses the default NetPyNE colorList.

    orderInverse : bool
        Whether or not to invert the y axis (useful if populations are defined top-down).

        *Default:* ``False`` does not invert the y-axis.

    returnPlotter : bool
        Whether to return the figure or the NetPyNE MetaFig object.

        *Default:* ``False`` returns the figure.


    Plot Options
    ------------
    showFig : bool
        Whether to show the figure.

        *Default:* ``False``

    saveFig : bool
        Whether to save the figure.

        *Default:* ``False``

    overwrite : bool
        whether to overwrite existing figure files.

        *Default:* ``True`` overwrites the figure file

        *Options:* ``False`` adds a number to the file name to prevent overwriting

    legendKwargs : dict
        a *dict* containing any or all legend kwargs.  These include ``'title'``, ``'loc'``, ``'fontsize'``, ``'bbox_to_anchor'``, ``'borderaxespad'``, and ``'handlelength'``.

    rcParams : dict
        a *dict* containing any or all matplotlib rcParams.  To see all options, execute ``import matplotlib; print(matplotlib.rcParams)`` in Python.  Any options in this *dict* will be used for this current figure and then returned to their prior settings.

    title : str
        the axis title

    xlabel : str
        label for x-axis

    ylabel : str
        label for y-axis

    s : int
        marker size

    marker : str
        marker symbol

    linewidth : int
        line width



    Returns
    -------
    rasterPlot : *matplotlib figure*
        By default, returns the *figure*.  If ``returnPlotter`` is ``True``, instead returns the NetPyNE MetaFig.


    """

    # If there is no input data, get the data from the NetPyNE sim object
    if rasterData is None:
        if 'sim' not in kwargs:
            from .. import sim
        else:
            sim = kwargs['sim']

        rasterData = sim.analysis.prepareRaster(include=['allCells'], timeRange=timeRange, maxSpikes=maxSpikes, orderBy=orderBy)

    # If the input is a dictionary, use the data inside it
    if isinstance(rasterData, dict):
        spkTimes = rasterData['spkTimes']
        spkInds = rasterData['spkInds']
        if 'popNumCells' in rasterData:
            popNumCells = rasterData['popNumCells']
        if 'popLabels' in rasterData:
            popLabels = rasterData['popLabels']
    else:  # Else it must be a list or a tuple
        spkTimes = rasterData[0]
        spkInds = rasterData[1]
        if len(rasterData) > 2:
            popNumCells = rasterData[2]
        if len(rasterData) > 3:
            popLabels = rasterData[3]

    # Ensure popNumCells and popLabels are not None
    if popNumCells is None:
        popNumCells = [len(set(spkInds))]  # assuming one population
    if popLabels is None:
        popLabels = ['Population']

    # Define the cell indices for plotting synchrony lines
    cellInds = list(set(spkInds))

    # Create the scatter plotter
    plotter = ScatterPlotter()
    fig = plotter.createFig()
    axis = plotter.getAxis(axis)
    plotter.setAttributes(title=title, xlabel=xlabel, ylabel=ylabel)

    # Ensure spkInds is ordered properly
    if orderBy == 'gid':
        spkInds = sorted(spkInds)
    else:
        # Order spkInds by some other cell tag
        spkInds = [spkInd for spkInd in sorted(spkInds, key=lambda x: sim.net.cells[x].tags[orderBy])]

    # Plot the raster
    scatter = axis.scatter(spkTimes, spkInds, s=s, marker=marker, linewidth=linewidth)

    # Add synchrony lines if needed
    if syncLines:
        for cellInd in cellInds:
            spike_times = [spkTimes[i] for i, ind in enumerate(spkInds) if ind == cellInd]
            for spike_time in spike_times:
                axis.vlines(spike_time, ymin=0, ymax=len(cellInds), colors='gray', linestyles='dotted', linewidth=linewidth)

    # Add legend
    if legend:
        legendHandles = [mpatches.Patch(color=popColors[pop], label=f'{pop} ({len([x for x in spkInds if x in range(sum(popNumCells[:i]), sum(popNumCells[:i + 1]))])})')
                         for i, pop in enumerate(popLabels)]
        axis.legend(handles=legendHandles, **kwargs.get('legendKwargs', {}))

    # Apply rcParams if provided
    if 'rcParams' in kwargs:
        for param, value in kwargs['rcParams'].items():
            plt.rcParams[param] = value

    # Show or save figure if specified
    if kwargs.get('showFig', False):
        plt.show()
    if kwargs.get('saveFig', False):
        filename = kwargs.get('filename', 'rasterPlot.png')
        plotter.saveFig(fig, filename, overwrite=kwargs.get('overwrite', True))

    return fig if not returnPlotter else plotter

