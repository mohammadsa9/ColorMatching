import matplotlib.pyplot as plt
from colour.plotting import *


def autolabel(rects, ax, xpos="center", p=6):
    """
    Attach a text label above each bar in *rects*, displaying its height.
    *xpos* indicates which side to place the text w.r.t. the center of
    the bar. It can be one of the following {'center', 'right', 'left'}.
    """
    ha = {"center": "center", "right": "left", "left": "right"}
    offset = {"center": 0, "right": 1, "left": -1}

    for rect in rects:
        height = rect.get_height()
        height = round(height, p)
        ax.annotate(
            "{}".format(height),
            xy=(rect.get_x() + rect.get_width() / 2, height),
            xytext=(offset[xpos] * 3, 3),  # use 3 points offset
            textcoords="offset points",  # in both directions
            ha=ha[xpos],
            va="bottom",
        )


def draw_R_style1(lines, comment=""):
    plt.legend(lines, [l.get_label() for l in lines])
    plt.gcf().canvas.set_window_title("Comparison")
    plt.xlabel("Wave Length \n\n" + comment)
    plt.ylabel("R")
    plt.gcf().set_size_inches(8, 8)
    plt.ylim(0, 1)
    plt.show()


def draw_KoverS_style1(lines):
    plt.legend(lines, [l.get_label() for l in lines])
    plt.gcf().canvas.set_window_title("Comparison")
    plt.xlabel("Wave Length")
    plt.ylabel("K/S")
    plt.gcf().set_size_inches(8, 8)
    plt.show()


def draw_CIE1931(arr=[]):
    # Plotting the *CIE 1931 Chromaticity Diagram*.
    plot_chromaticity_diagram_CIE1931(standalone=False)

    for spec in arr:
        xy_D65 = spec.getxy()
        xy = xy_D65
        x, y = xy
        plt.plot(x, y, "o-", color=spec.color)
        if spec.name != "":
            plt.annotate(
                spec.name,
                xy=xy,
                xytext=(-50, 30),
                textcoords="offset points",
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3, rad=-0.2"),
            )
    # Displaying the plot.
    render(
        standalone=True, limits=(-0.1, 0.9, -0.1, 0.9), x_tighten=True, y_tighten=True
    )


def draw_CIE1931_style2(arr=[]):
    # Plotting the *CIE 1931 Chromaticity Diagram*.
    plot_chromaticity_diagram_CIE1931(standalone=False)

    for spec in arr:
        xy_D65 = spec.getxy()
        xy = xy_D65
        x, y = xy
        plt.plot(x, y, ".", color=spec.color, markersize=5)
        if spec.name != "":
            plt.annotate(
                spec.name,
                xy=xy,
                xytext=(-50, 30),
                textcoords="offset points",
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3, rad=-0.2"),
            )
    # Displaying the plot.
    render(
        standalone=True, limits=(-0.1, 0.9, -0.1, 0.9), x_tighten=True, y_tighten=True
    )