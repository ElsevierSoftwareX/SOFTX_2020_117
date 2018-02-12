import os
import contextlib

__all__ = ["use", "context"]

# Get style directory and available styles
style_dir = os.path.split(os.path.realpath(__file__))[0]
available = [os.path.splitext(f)[0]
             for f in os.listdir(style_dir) if f.endswith(".mplstyle")]

def get_style_path(style):
    if isinstance(style, str) or hasattr(style, 'keys'):
        # If name is a single str or dict, make it a single element list.
        styles = [style]
    else:
        styles = style
    styles = list(map((lambda s: os.path.join(style_dir, s + ".mplstyle")),
                  styles))
    return styles

def use(style):
    import matplotlib.pyplot as plt
    plt.style.use(get_style_path(style))

def context(style, after_reset=False):
    import matplotlib.pyplot as plt
    return plt.style.context(get_style_path(style), after_reset)
