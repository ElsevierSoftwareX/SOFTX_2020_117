"""Finesse command line interface

Sean Leavey
<sean.leavey@ligo.org>
"""

from __future__ import print_function

import sys
import numpy as np
import click

from . import __version__
from .finesse import kat as katparser

@click.command(help="Python interface and tools for FINESSE")
@click.argument("file", type=click.File())
@click.option("--xstart", type=float,
              help="Simulation start value. If specified, this overrides the xaxis start "
              "value specified in the parsed file.")
@click.option("--xstop", type=float,
              help="Simulation stop value. If specified, this overrides the xaxis stop "
              "value specified in the parsed file.")
@click.option("--xsteps", type=int,
              help="Number of steps to simulate between --start and --stop. If specified, "
              "this overrides the number of xaxis steps specified in the parsed file.")
@click.option("--xscale", type=click.Choice(["lin", "log"]), help="Scaling for the xaxis.")
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
@click.option("--maxtem", type=str, help="Maximum transverse electric mode. Can be either "
              "an integer or 'off'.")
@click.option("--ignore-block", "ignored_blocks", multiple=True,
              help="Ignore the specified block. Can be specified multiple times.")
@click.option("--plot/--no-plot", default=True, show_default=True,
              help="Display results as figure.")
@click.option("--save-figure", type=click.File("wb", lazy=False),
              help="Save image of figure to file.")
@click.version_option(version=__version__, prog_name="Pykat")
def cli(file, xstart, xstop, xsteps, xscale, trace, maxtem, ignored_blocks, plot, save_figure):
    """Base CLI command group"""
    kat = katparser()
    kat.load(file.name)
    has_xaxis = hasattr(kat, "xaxis")

    if xstart is not None or xstop is not None or xsteps is not None or xscale is not None:
        if not has_xaxis:
            click.echo("Limits can only be overridden when an xaxis is defined in FILE.",
                       err=True)
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

    if maxtem:
        kat.maxtem = maxtem

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

    if ignored_blocks:
        for block in ignored_blocks:
            kat.removeBlock(block)

    results = kat.run()

    if kat.trace:
        click.echo(results.stdout)

    if has_xaxis:
        if plot or save_figure is not None:
            results.plot(show=plot, filename=save_figure)
    else:
        if save_figure is not None:
            click.echo("Cannot plot or save figure without an xaxis defined in FILE.",
                       err=True)
            sys.exit(1)
