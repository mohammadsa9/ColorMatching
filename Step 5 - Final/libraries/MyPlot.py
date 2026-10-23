import matplotlib.pyplot as plt
import matplotlib as mpl
from colour.plotting import *
from colour.plotting import override_style
from colour.colorimetry import CCS_ILLUMINANTS
import colour
import math
import random

OUTPUT = 0
Name = 1
Swatch = []

# Set the font size of plots
plt.rcParams.update(
    {"font.size": 12, "font.weight": "bold", "axes.labelweight": "bold"}
)


def array_distance(start, distance, end):
    result = []
    rg = (end - start) / distance
    for i in range(int(rg) + 1):
        result.append(int(start + i * distance))
    return result


wave_length = array_distance(400, 10, 700)


def makeout():
    global OUTPUT, Name
    OUTPUT = 1
    Name = 1


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
    plt.xlabel("Wavelength (nm)" + comment)
    plt.ylabel("Reflectance Factor")
    plt.gcf().set_size_inches(8, 8)
    # plt.ylim(0, 1)
    plt.xlim(400, 700)
    global OUTPUT, Name
    if OUTPUT:
        plt.savefig("output_img/" + str(Name) + ".jpg", bbox_inches="tight", dpi=1200)
        Name += 1
    plt.show()


def draw_R_subplot_style1(R_Samples, Method_PCC_R, array, method="PCA", mode="sub"):
    global wave_length

    M_R_minC_no = array[6]
    M_R_maxC_no = array[7]
    if M_R_minC_no == M_R_maxC_no:
        array.pop(7)
        array.pop(6)
    if mode == "sub":
        fig, axs = plt.subplots(4, 2, constrained_layout=True)
    strings = [
        "Min ΔE",
        "Max ΔE",
        "Min RMS",
        "Max RMS",
        "Max GFC",
        "Min GFC",
        "Min ΔC",
        "Max ΔC",
    ]
    row, col, index = 0, 0, 0
    for i in array:
        if mode == "sub":
            # Each Cell
            (p1,) = axs[row, col].plot(
                wave_length,
                R_Samples[i],
                color="green",
                label="R Sample " + str(i),
                linewidth=3,
            )
            (p2,) = axs[row, col].plot(
                wave_length,
                Method_PCC_R[i],
                color="black",
                label="R Interpolated " + method,
                linewidth=3,
                linestyle="dashed",
            )
            lines = [p1, p2]
            axs[row, col].legend(lines, [l.get_label() for l in lines])
            col += 1
            if col == 2:
                row += 1
                col = 0
        else:
            print(strings[index])
            (p1,) = plt.plot(
                wave_length,
                R_Samples[i],
                color="green",
                label="R Sample " + str(i),
                linewidth=3,
            )
            (p2,) = plt.plot(
                wave_length,
                Method_PCC_R[i],
                color="black",
                label="R Interpolated " + method,
                linewidth=3,
                linestyle="dashed",
            )
            lines = [p1, p2]
            draw_R_style1(lines)
            index += 1

    if mode == "sub":
        axs[0, 0].set_title("Min ΔE")
        axs[0, 1].set_title("Max ΔE")
        axs[1, 0].set_title("Min RMS")
        axs[1, 1].set_title("Max RMS")
        axs[2, 0].set_title("Max GFC")
        axs[2, 1].set_title("Min GFC")
        axs[3, 0].set_title("Min ΔC")
        axs[3, 1].set_title("Max ΔC")

        for ax in axs.flat:
            ax.set(xlabel="Wavelength (nm)", ylabel="Reflectance Factor")
            # ax.set_ylim(0, 1)
            ax.set_xlim(400, 700)
            ax.figure.set_size_inches(10, 10)

        # Hide x labels and tick labels for top plots and y ticks for right plots.
        for ax in axs.flat:
            ax.label_outer()

        plt.show()
        global OUTPUT, Name
        if OUTPUT:
            fig.savefig(
                "output_img/" + str(Name) + ".jpg", bbox_inches="tight", dpi=1200
            )
            Name += 1
    draw_R_subplot_style2(R_Samples, Method_PCC_R, array, method, mode)


def draw_R_subplot_style2(R_Samples, Method_PCC_R, array, method="PCA", mode="sub"):
    global wave_length

    M_R_minE_no = array[0]
    M_R_maxE_no = array[1]
    M_R_minRMS_no = array[2]
    M_R_maxRMS_no = array[3]
    M_R_minGFC_no = array[4]
    M_R_maxGFC_no = array[5]
    try:
        M_R_minC_no = array[6]
        M_R_maxC_no = array[7]
    except Exception:
        M_R_minC_no = array[4]
        M_R_maxC_no = array[5]

    Method_PCC_R.pop(M_R_minE_no, None)
    Method_PCC_R.pop(M_R_maxE_no, None)
    Method_PCC_R.pop(M_R_minRMS_no, None)
    Method_PCC_R.pop(M_R_maxRMS_no, None)
    Method_PCC_R.pop(M_R_minGFC_no, None)
    Method_PCC_R.pop(M_R_maxGFC_no, None)
    Method_PCC_R.pop(M_R_minC_no, None)
    Method_PCC_R.pop(M_R_maxC_no, None)
    if len(list(Method_PCC_R)) >= 8:
        plot_size = 4
    else:
        plot_size = int(len(list(Method_PCC_R)) / 2)

    random_key = random.sample(list(Method_PCC_R), plot_size * 2)
    if mode == "sub":
        fig, axs = plt.subplots(plot_size, 2, constrained_layout=True)
    row, col = 0, 0
    for i in random_key:
        if mode == "sub":
            # Each Cell
            (p1,) = axs[row, col].plot(
                wave_length,
                R_Samples[i],
                color="green",
                label="R Sample " + str(i),
                linewidth=3,
            )
            (p2,) = axs[row, col].plot(
                wave_length,
                Method_PCC_R[i],
                color="black",
                label="R Interpolated " + method,
                linewidth=3,
                linestyle="dashed",
            )
            lines = [p1, p2]
            axs[row, col].legend(lines, [l.get_label() for l in lines])

            col += 1
            if col == 2:
                row += 1
                col = 0
        else:
            print("Random Sample")
            (p1,) = plt.plot(
                wave_length,
                R_Samples[i],
                color="green",
                label="R Sample " + str(i),
                linewidth=3,
            )
            (p2,) = plt.plot(
                wave_length,
                Method_PCC_R[i],
                color="black",
                label="R Interpolated " + method,
                linewidth=3,
                linestyle="dashed",
            )
            lines = [p1, p2]
            draw_R_style1(lines)

    if mode == "sub":
        for ax in axs.flat:
            ax.set(xlabel="Wavelength (nm)", ylabel="Reflectance Factor")
            # ax.set_ylim(0, 1)
            ax.set_xlim(400, 700)
            ax.figure.set_size_inches(10, 10)
            ax.set_title("Random Sample")

        # Hide x labels and tick labels for top plots and y ticks for right plots.
        for ax in axs.flat:
            ax.label_outer()

        plt.show()
        global OUTPUT, Name
        if OUTPUT:
            fig.savefig(
                "output_img/" + str(Name) + ".jpg", bbox_inches="tight", dpi=1200
            )
            Name += 1


def plot_save(plot):
    global OUTPUT, Name
    if OUTPUT:
        plot.savefig("output_img/" + str(Name) + ".jpg", bbox_inches="tight", dpi=1200)
        Name += 1


def draw_KoverS_style1(lines):
    plt.legend(lines, [l.get_label() for l in lines])
    plt.gcf().canvas.set_window_title("Comparison")
    plt.xlabel("Wavelength (nm)")
    plt.ylabel("k/s")
    plt.xlim(400, 700)
    plt.gcf().set_size_inches(8, 8)
    global OUTPUT, Name
    if OUTPUT:
        plt.savefig("output_img/" + str(Name) + ".jpg", bbox_inches="tight", dpi=1200)
        Name += 1
    plt.show()


def draw_CIE1931(arr=[], small_point=False):
    # Plotting the *CIE 1931 Chromaticity Diagram*.
    settings = {
        "standalone": False,
        "wrap_title": False,
        # "title": "CIE 1964 Chromaticity Diagram",
        "title": "",
        "transparent_background": False,
    }
    plot_chromaticity_diagram_CIE1931(
        cmfs="CIE 1964 10 Degree Standard Observer", **settings
    )

    font_size_small = 3
    for spec in arr:
        xy_D65 = spec.getxy()
        xy = xy_D65
        x, y = xy
        text_x = -40
        text_y = 30
        if small_point:
            if spec.color.split("_")[0] == "special":
                plt.plot(
                    x,
                    y,
                    marker="v",
                    color=spec.color.split("_")[1],
                    markersize=font_size_small,
                )
            else:
                plt.plot(x, y, ".", color=spec.color, markersize=font_size_small)
        else:
            plt.plot(x, y, "o-", color=spec.color)
        if spec.name != "":
            if x < 0.3:
                text_x = 0
            elif x > 0.5:
                text_x = -80
                text_y = 20
            plt.annotate(
                spec.name,
                xy=xy,
                xytext=(text_x, text_y),
                textcoords="offset points",
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3, rad=-0.2"),
            )
    # Displaying the plot.
    fig, ax = render(
        standalone=True, limits=(-0.1, 0.9, -0.1, 0.9), x_tighten=True, y_tighten=True
    )
    global OUTPUT, Name
    if OUTPUT:
        fig.savefig("output_img/" + str(Name) + ".jpg", bbox_inches="tight", dpi=1200)
        Name += 1


def draw_CIELab(arr=[], small_point=False):

    fig = plt.figure()
    ax = plt.gca()
    # ax.grid(True)
    ax.spines["left"].set_position("zero")
    ax.spines["right"].set_color("none")
    ax.spines["bottom"].set_position("zero")
    ax.spines["top"].set_color("none")

    plt.xlim(-100, 100)
    plt.ylim(-100, 100)
    # plt.gcf().set_size_inches(8, 8)

    font_size_small = 3
    for spec in arr:
        xy_D = spec.getAB()
        xy = xy_D
        x, y = xy
        text_x = -40
        text_y = 30
        if small_point:
            if spec.color.split("_")[0] == "special":
                plt.plot(
                    x,
                    y,
                    marker="v",
                    color=spec.color.split("_")[1],
                    markersize=font_size_small,
                )
            else:
                plt.plot(x, y, ".", color=spec.color, markersize=font_size_small)
        else:
            plt.plot(x, y, "o-", color=spec.color)
        if spec.name != "":
            if y < 0:
                text_x = 40
                text_y = -30
            plt.annotate(
                spec.name,
                xy=xy,
                xytext=(text_x, text_y),
                textcoords="offset points",
                arrowprops=dict(arrowstyle="->", connectionstyle="arc3, rad=-0.2"),
            )
    # Set Label
    # Add the label as annotation. The "5" is the padding betweent the right side
    # of the axis and the label...
    ticklabelpad = mpl.rcParams["xtick.major.pad"]
    plt.annotate(
        "+a* (Red)",
        xy=(1.02, 0.52),
        ha="left",
        va="top",
        xycoords="axes fraction",
        textcoords="offset points",
    )
    plt.annotate(
        "-a* (Green)",
        xy=(-0.25, 0.52),
        ha="left",
        va="top",
        xycoords="axes fraction",
        textcoords="offset points",
    )

    plt.annotate(
        "+b* (Yellow)",
        xy=(0.48, 1.1),
        xytext=(5, -ticklabelpad),
        ha="left",
        va="top",
        xycoords="axes fraction",
        textcoords="offset points",
    )
    plt.annotate(
        "-b* (Blue)",
        xy=(0.48, -0.02),
        xytext=(5, -ticklabelpad),
        ha="left",
        va="top",
        xycoords="axes fraction",
        textcoords="offset points",
    )

    # Displaying the plot.
    # fig, ax = render()
    plt.show()
    global OUTPUT, Name
    if OUTPUT:
        fig.savefig("output_img/" + str(Name) + ".jpg", bbox_inches="tight", dpi=1200)
        Name += 1


def resetSwatch():
    global Swatch
    Swatch = []


def draw_rgb_from_XYZ(XYZ, name="", draw=True):
    global OUTPUT, Name, Swatch
    illuminant = CCS_ILLUMINANTS["CIE 1964 10 Degree Standard Observer"]["D65"]
    RGB = colour.XYZ_to_sRGB(XYZ / 100, illuminant=illuminant)
    Swatch.append(ColourSwatch("", RGB))
    if draw:
        fig, ax = plot_single_colour_swatch(
            ColourSwatch(name, RGB), text_kwargs={"size": "x-large"}
        )
        if OUTPUT:
            fig.savefig(
                "output_img/" + str(Name) + ".jpg", bbox_inches="tight", dpi=1200
            )
            Name += 1


def draw_rgb_from_all():
    global OUTPUT, Name, Swatch

    # Settings
    background_colour = "0.1"
    width = height = 1.0
    spacing = 0.25
    columns = int(math.sqrt(len(Swatch)))

    settings = {
        "width": width,
        "height": height,
        "spacing": spacing,
        "columns": columns,
        "direction": "-y",
        "text_kwargs": {"size": 8},
        "background_colour": background_colour,
        "compare_swatches": None,
    }
    settings["standalone"] = False
    fig, ax = colour.plotting.plot_multi_colour_swatches(Swatch, **settings)
    if OUTPUT:
        fig.savefig("output_img/" + str(Name) + ".jpg", bbox_inches="tight", dpi=1200)
        Name += 1


def box_plot(data):
    fig = plt.figure(figsize=(10, 7))

    # Creating plot
    plt.boxplot(data, labels=["PCC", "XYZ", "PCC + XYZ"])
    global OUTPUT, Name
    if OUTPUT:
        plt.savefig("output_img/" + str(Name) + ".jpg", bbox_inches="tight", dpi=1200)
        Name += 1
    # show plot
    plt.show()