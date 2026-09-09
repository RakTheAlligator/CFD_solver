#!/usr/bin/env python3
"""Display live SIMPLE convergence data written by CFD_solver."""

import argparse
import csv
import math
from pathlib import Path


DEFAULT_SERIES = (
    ("continuity", "continuity"),
    ("x_velocity", "x-velocity"),
    ("y_velocity", "y-velocity"),
)

DIAGNOSTIC_SERIES = (
    ("velocity_change", "velocity change"),
    ("corrected_continuity", "corrected continuity"),
    ("u_linear_residual", "u linear"),
    ("v_linear_residual", "v linear"),
    ("pressure_correction_linear_residual", "p' linear"),
)


def positive_integer(text: str) -> int:
    """Parse a strictly positive command-line integer."""
    try:
        value = int(text)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be an integer") from error
    if value <= 0:
        raise argparse.ArgumentTypeError("must be strictly positive")
    return value


def parse_arguments() -> argparse.Namespace:
    """Parse the standalone plotter command line."""
    parser = argparse.ArgumentParser(
        description="Plot live SIMPLE residuals from a convergence CSV file."
    )
    parser.add_argument(
        "csv_path",
        type=Path,
        help="path to convergence.csv; the plotter waits if it does not exist yet",
    )
    parser.add_argument(
        "--refresh-ms",
        type=positive_integer,
        default=250,
        metavar="INTEGER",
        help="CSV polling interval in milliseconds (default: 250)",
    )
    parser.add_argument(
        "--diagnostics",
        action="store_true",
        help="also plot velocity change, corrected continuity, and inner linear residuals",
    )
    return parser.parse_args()


def read_rows(csv_path: Path, columns: tuple[str, ...]) -> list[tuple[int, dict[str, float]]]:
    """Read every complete finite row currently visible in the CSV file."""
    try:
        with csv_path.open("r", encoding="utf-8", newline="") as input_file:
            reader = csv.DictReader(input_file)
            required_columns = {"iteration", *columns}
            if reader.fieldnames is None or not required_columns.issubset(reader.fieldnames):
                return []

            rows: list[tuple[int, dict[str, float]]] = []
            for source_row in reader:
                try:
                    iteration = int(source_row["iteration"])
                    values = {column: float(source_row[column]) for column in columns}
                except (KeyError, TypeError, ValueError):
                    continue
                if iteration <= 0 or any(not math.isfinite(value) for value in values.values()):
                    continue
                rows.append((iteration, values))
            return rows
    except (FileNotFoundError, OSError, csv.Error):
        return []


def main() -> None:
    """Run the responsive Matplotlib live plot."""
    arguments = parse_arguments()

    try:
        import matplotlib.pyplot as plt
        from matplotlib.animation import FuncAnimation
    except ImportError as error:
        raise SystemExit(
            "matplotlib is required; install it with 'python3 -m pip install matplotlib'."
        ) from error

    series = DEFAULT_SERIES + (DIAGNOSTIC_SERIES if arguments.diagnostics else ())
    columns = tuple(column for column, _label in series)

    figure, axes = plt.subplots()
    plotted_lines = {
        column: axes.plot([], [], label=label)[0] for column, label in series
    }
    axes.set_yscale("log")
    axes.set_xlabel("Iterations")
    axes.set_ylabel("Normalized residual")
    axes.set_title("SIMPLE convergence")
    axes.set_xlim(1.0, 2.0)
    axes.set_ylim(1.0e-12, 1.0)
    axes.grid(True, which="both", linestyle=":", alpha=0.6)
    axes.legend()

    def update_plot(_frame: int):
        rows = read_rows(arguments.csv_path, columns)
        if not rows:
            for line in plotted_lines.values():
                line.set_data([], [])
            axes.set_xlim(1.0, 2.0)
            axes.set_ylim(1.0e-12, 1.0)
            return tuple(plotted_lines.values())

        iterations = [iteration for iteration, _values in rows]
        positive_values = []
        for column, _label in series:
            plot_values = []
            for _iteration, values in rows:
                value = values[column]
                if value > 0.0:
                    plot_values.append(value)
                    positive_values.append(value)
                else:
                    plot_values.append(math.nan)
            plotted_lines[column].set_data(iterations, plot_values)

        minimum_iteration = min(iterations)
        maximum_iteration = max(iterations)
        if minimum_iteration == maximum_iteration:
            axes.set_xlim(max(0.5, minimum_iteration - 0.5), maximum_iteration + 0.5)
        else:
            axes.set_xlim(minimum_iteration, maximum_iteration)

        if positive_values:
            minimum_log = math.log10(min(positive_values))
            maximum_log = math.log10(max(positive_values))
            logarithmic_span = maximum_log - minimum_log
            padding = max(0.2, 0.05 * logarithmic_span)
            axes.set_ylim(10.0 ** (minimum_log - padding), 10.0 ** (maximum_log + padding))

        figure.canvas.draw_idle()
        return tuple(plotted_lines.values())

    animation = FuncAnimation(
        figure,
        update_plot,
        interval=arguments.refresh_ms,
        blit=False,
        cache_frame_data=False,
    )
    plt.show()


if __name__ == "__main__":
    main()
