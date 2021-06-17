def annotate_rectangles_with_values(rectangles, ax):
    # https://matplotlib.org/stable/gallery/lines_bars_and_markers/barchart.html
    for rect in rectangles:
        height = rect.get_height()
        ax.annotate(f"{height}%",
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 1),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')
