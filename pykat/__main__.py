"""Finesse command line interface

Sean Leavey
<sean.leavey@ligo.org>
"""

import sys
import os
import numpy as np
import click

from . import __version__
from .finesse import kat as katparser

@click.command(help="Python interface and tools for FINESSE")
@click.argument("file", type=click.File())
@click.option("--simulate/--no-simulate", is_flag=True, default=True,
              help="Simulate FILE in Finesse. Can be set to --no-simulate if for example "
              "you only want to display the node graph for the model specified in FILE.")
@click.option("--xstart", type=float,
              help="Set simulation start value. If specified, this overrides the xaxis start "
              "value specified in the parsed file.")
@click.option("--xstop", type=float,
              help="Set simulation stop value. If specified, this overrides the xaxis stop "
              "value specified in the parsed file.")
@click.option("--xsteps", type=int,
              help="Set number of steps to simulate between --start and --stop. If specified, "
              "this overrides the number of xaxis steps specified in the parsed file.")
@click.option("--xscale", type=click.Choice(["lin", "log"]),
              help="Set scaling for the xaxis.")
@click.option("--noxaxis", is_flag=True, default=False, help="Switch off x-axis.")
@click.option("--trace", type=click.Choice(["tem", "cavity", "mismatch", "beams", "gouy",
                                            "coupling", "modechanges", "nodes", "all"]),
              multiple=True, help="Show simulation trace results. This option can be specified "
                                  "multiple times. The following values are supported: "
                                  "'tem': list TEM modes used, "
                                  "'cavity': cavity eigenvalues and other parameters, "
                                  "'mismatch': mode mismatch parameters for the initial setup, "
                                  "'beams': beam parameters for every node, "
                                  "'gouy': Gouy phases for all spaces, "
                                  "'coupling': coupling coefficients for all components, "
                                  "'modechanges': mode matching parameter changes during calculations, "
                                  "'nodes': nodes found during cavity tracing, "
                                  "'all': all trace results.")
@click.option("--powers", type=click.Choice(["dc", "carrier", "tem00"]),
              multiple=True, help="Show powers (W) at each node in the interferometer. This "
                                  "option can be specified multiple times. The following "
                                  "values are supported: "
                                  "'dc': list dc powers, "
                                  "'carrier': list powers in the f=0 fields, "
                                  "'tem00': list powers in the TEM00 mode.")
@click.option("--maxtem", type=str, help="Set maximum transverse electric mode. Can be either "
              "an integer or 'off'.")
@click.option("--phase", type=int, help="Set Gouy phase behaviour.")
@click.option("--retrace", type=click.Choice(["force", "off"]),
              help="Set retrace behaviour: 'force' recomputes the Gaussian parameters at each "
              "node for every data point, and will trace a cavity even if it is unstable; 'off' "
              "switches off retracing even if it normally would be done.")
@click.option("--deriv-h", type=float, help="Set step size for numerical differentiation.")
@click.option("--lambda0", type=str, help="Set reference wavelength (m). Supports SI prefixes.")
@click.option("--print-frequencies", is_flag=True, default=False,
              help="Print a table listing all the frequencies used in the simulation: carriers, "
              "modulation sidebands and signal/quantum sidebands.")
@click.option("--print-noises", is_flag=True, default=False,
              help="Print a list of all the quantum noise sources being considered in the "
              "simulation and at which components and nodes it is present.")
@click.option("--ignore-block", "ignored_blocks", multiple=True,
              help="Ignore the specified block. Can be specified multiple times.")
@click.option("--plot/--no-plot", default=True, show_default=True,
              help="Display results as figure.")
@click.option("--save-figure", type=click.File("wb", lazy=False),
              help="Save image of figure to file.")
@click.option("--display-graph", is_flag=True, help="Generate and display model node graph "
              "using default system document viewer.")
@click.option("--save-input", is_flag=True, default=False,
              help="Save generated Finesse input file.")
@click.option("--save-output", is_flag=True, default=False,
              help="Save generated Finesse output file.")
@click.option("--finesse-dir", type=click.Path(exists=True, file_okay=False, dir_okay=True),
              envvar='FINESSE_DIR', help="Path to directory containing the Finesse 'kat' "
              "executable. If not specified, the environment variable FINESSE_DIR is used.")
@click.version_option(version=__version__, prog_name="Pykat")
def cli(file, simulate, xstart, xstop, xsteps, xscale, noxaxis, trace, powers, maxtem, phase,
        retrace, deriv_h, lambda0, print_frequencies, print_noises, ignored_blocks, plot,
        save_figure, display_graph, save_input, save_output, finesse_dir):
    """Base CLI command group"""
    if finesse_dir is None:
        # Use default as required by kat object.
        finesse_dir = ""
    kat = katparser(katdir=finesse_dir)
    kat.load(file.name)
    has_xaxis = hasattr(kat, "xaxis") and not noxaxis

    if ignored_blocks:
        for block in ignored_blocks:
            kat.removeBlock(block)

    if display_graph:
        from .tools.plotting.graph import NodeGraph
        nodegraph = NodeGraph(kat)
        nodegraph.view_pdf()

    if simulate:
        # Default flag for showing simulation output. Overridden below as necessary.
        output_to_show = False

        if xstart is not None or xstop is not None or xsteps is not None or xscale is not None:
            if not has_xaxis:
                click.echo("Limits can only be overridden when an xaxis is defined in FILE and "
                           "when --noxaxis is unset.", err=True)
                sys.exit(1)
            # Override xaxis.
            limits = kat.xaxis.limits
            set_limits = False
            if xstart is not None:
                limits[0] = xstart
                set_limits = True
            if xstop is not None:
                limits[1] = xstop
                set_limits = True
            if xsteps is not None:
                kat.xaxis.steps = xsteps
            if xscale is not None:
                kat.xaxis.scale = xscale
            
            if set_limits:
                kat.xaxis.limits = np.array(limits).astype(float)

        if not has_xaxis:
            kat.noxaxis = True
            if save_figure is not None:
                click.echo("Cannot plot or save figure without an xaxis defined in FILE.",
                           err=True)
                sys.exit(1)

        if maxtem is not None:
            kat.maxtem = maxtem
        if phase is not None:
            kat.phase = phase
        if retrace is not None:
            kat.retrace = retrace
        if deriv_h is not None:
            kat.deriv_h = deriv_h
        if lambda0 is not None:
            kat.lambda0 = lambda0
        if print_frequencies:
            kat.parse("frequency")
            output_to_show = True
        if print_noises:
            kat.parse("printnoises")
            output_to_show = True

        if trace:
            traceval = 0
            if "all" in trace:
                traceval = 255
            else:
                traceints = {"tem": 1, "cavity": 2, "mismatch": 4, "beams": 8,
                            "gouy": 16, "coupling": 32, "modechanges": 64, "nodes": 128}
                for tracetype in trace:
                    traceval |= traceints[tracetype]
            kat.trace = traceval
            output_to_show = True

        if powers:
            powerval = 0
            powerints = {"dc": 1, "carrier": 2, "tem00": 4}
            for powertype in powers:
                powerval |= powerints[powertype]
            kat.parse("powers %i" % powerval)
            output_to_show = True

        results = kat.run(save_output=save_output, save_kat=save_input)

        if output_to_show and results.rundata:
            click.echo(results.rundata)

        if has_xaxis and (plot or save_figure is not None):
            results.plot(show=plot, filename=save_figure)
        elif not output_to_show and not display_graph:
            click.echo("No output requested.")
